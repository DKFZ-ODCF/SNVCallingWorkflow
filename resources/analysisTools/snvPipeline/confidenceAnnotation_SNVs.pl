#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#


use strict;
use warnings;

# empirical confidence classification

if (@ARGV < 2)
{
	die "USAGE: (1)vcf formatted file for confidence classification - (2)Config file from which to get the column names - (3)optional: lower threshold for too high read depth in control (default 150, set higher for target sequencing or very deep sequencing)\n";
}

my $file = shift;
my $configfile = shift;
my $cutoff = shift;
unless (defined $cutoff)
{
	$cutoff = 150;
}

open (CF, $configfile) or die "Could not open $configfile: $!\n";

# get colunm labels from config file:
### Annotation (SNV and Indel)
#KGENOMES_COL=1K_GENOMES
#DBSNP_COL=DBSNP
#ANNOVAR_SEGDUP_COL=SEGDUP
#ANNOVAR_SELFCHAIN_COL=SELFCHAIN

my %labels = ();
my @help = ();
while (<CF>)
{
	if ($_ =~ /_COL=/)
	{
		chomp;
		@help = split ("=", $_);
		$labels{$help[0]} = $help[1];
	}
}
close CF;

foreach my $key (keys %labels)
{
	print STDERR "$key\t$labels{$key}\n";
}

#exit;
open (FH, $file) or die "Could not open $file: $!\n";

my $header = "";
while ($header = <FH>)
{
	last if ($header =~ /^\#CHROM/); # that is the line with the column names
	print "$header";
}
chomp $header;
print $header, "\tCONFIDENCE\tRECLASSIFICATION\n";
# get the columns where to look for features:
# INFO => DP4 of tumor => variant frequency and strand bias, CONTROL_INFO => DP5 of control => variant frequency , SOMATIC_GERMLINE_CLASSIFICATION => classification
my $INFO = 0;
my $CINFO = 0; 
my $CLASS = 0;
my $MAPAB = 0;
my $SEGDP = 0;
my $HSDEP = 0;
my $IN1KG = 0;
my $DBSBP = 0;
my $BLACK = 0;
my $EXCLU = 0;
my $STREP = 0;
my $REPET = 0;
my $CHAIN = 0;
@help = (split "\t", $header);
for (my $c = 0; $c < @help; $c++)
{
	# fixed column labels (from vcf_pileup_compare)
	if ($help[$c] eq "INFO")
	{
		$INFO = $c;
		print STDERR "INFO in column $c\n";
	}
	if ($help[$c] =~ "^INFO_control")
	{
		$CINFO = $c;
		print STDERR "INFO_control in column $c\n";
	}
	if ($help[$c] eq "ANNOTATION_control")
	{
		$CLASS = $c;
		print STDERR "ANNOTATION_control in column $c\n";
	}
	# from annotation
	if ($help[$c] =~ /MAP+ABILITY/)
	{
		$MAPAB = $c;
		print STDERR "MAPABILITY in column $c\n";
	}
	if ($help[$c] =~ /HISEQDEPTH/)
	{
		$HSDEP = $c;
		print STDERR "HISEQDEPTH in column $c\n";
	}
	if ($help[$c] =~ /SIMPLE_TANDEMREPEATS/)	# simple tandem repeats from Tandem Repeats Finder
	{
		$STREP = $c;
		print STDERR "SIMPLE_TANDEMREPEATS in column $c\n";
	}
	if ($help[$c] =~ /REPEAT_MASKER/)	# RepeatMasker annotation
	{
		$REPET = $c;
		print STDERR "REPEAT_MASKER in column $c\n";
	}
	if ($help[$c] =~ /DUKE_EXCLUDED/)
	{
		$EXCLU = $c;
		print STDERR "DUKE_EXCLUDED in column $c\n";
	}
	if ($help[$c] =~ /DAC_BLACKLIST/)
	{
		$BLACK = $c;
		print STDERR "DAC_BLACKLIST in column $c\n";
	}
	if ($help[$c] =~ /SELFCHAIN/)	# not used yet!
	{
		$CHAIN = $c;
		print STDERR "SELFCHAIN in column $c\n";
	}
	# column labels according to config file; if that fails, fallback hardcoded
	if ($help[$c] eq $labels{ANNOVAR_SEGDUP_COL} || $help[$c] eq "SEGDUP")
	{
		$SEGDP = $c;
		print STDERR "ANNOVAR_SEGDUP_COL in column $c\n";
	}
	if ($help[$c] eq $labels{KGENOMES_COL} || $help[$c] eq "1K_GENOMES")
	{
		$IN1KG = $c;
		print STDERR "KGENOMES_COL in column $c\n";
	}
	if ($help[$c] eq $labels{DBSNP_COL} || $help[$c] eq "DBSNP")
	{
		$DBSBP = $c;
		print STDERR "DBSNP_COL in column $c\n";
	}
}

my $line = "";	# current line

# number of reference and variant reads on forward and reverse complement strand for tumor
my ($trf, $trr, $tvf, $tvr) = (0, 0, 0, 0);
# same for control
my ($crf, $crr, $cvf, $cvr) = (0, 0, 0, 0);

# omitted
# P-values for 1) strand bias (exact test); 2) baseQ bias (t-test); 3) mapQ bias (t); 4) tail distance bias (t)
# the PV4 INFO ist not always there, but if, is it always at the end?
# but the TDbias is 1 in validated as well as in devalidated cases ...
#my ($strandbias, $baseQbias, $mapQbias, $TDbias) = (0, 0, 0, 0);

# absolute coverage at position
my $depthT = 0;
my $depthC = 0;

# variant supporting reads in control incl. low quality
my $varsuppC = 0;

# fraction of variant reads in tumor
my $fr_var_tum = 0;
# fraction of variant reads in control: the classification called "somatic" is < 1/30 of control reads support variant
my $fr_var_ctrl = 0;	#=round(float((DP5[2]+DP5[3]))/sum(DP5),2)

my $confidence=10;	# start with maximum value and subtract something for each "bad" feature
			# "high" = 9-10, "medium" = 6-8, low <= 5
my $mapp = 0;	# Mapability
my $is_repeat = 0;	# true if SNP conicides with any of the suspicious repeat classes simple, low and satellite, which are partially redundant to other annotations
my $is_STR = 0;	# also short tandem repeats from Tandem Repeats Finder are suspicious and can conincide with other classes
my $is_weird = 0;	# coindicende with known artefact-producing regions
my $in1KG = 0;
my $indbSNP = 0;
my $precious = 0;
my $class = "";	# for potential re-classification (e.g. low coverage in control and in dbSNP => probably germline)

while (<FH>)
{
	$confidence=10;	# start with maximum value
	# reset global variables
	$in1KG = 0;
	$indbSNP = 0;
	$precious = 0;
	$is_repeat = 0;
	$is_STR = 0;
	$is_weird = 0;
	$line = $_;
	chomp $line;
	@help = split ("\t", $line);
	$class = $help[$CLASS];	# start with original classification

# 1) external information of if these SNPs have already been found (incl. false positives from 1000 genomes!)
	# 1000 genomes
	#if ($help[$IN1KG] ne "." && $help[$IN1KG] ne "0")
	if ($help[$IN1KG] =~ /MATCH=exact/)	# verified somatic SNVs have MATCH=position instead of exact, so it would be bad to reclassify them as SNP_support_germline!
	{
		$in1KG = 1;	# good for germline, bad for "somatic"
	}
	# dbSNP
	#if ($help[$DBSBP] ne "." && $help[$DBSBP] ne "0")	# verified somatic SNVs have MATCH=position instead of exact, so it would be bad to reclassify them as SNP_support_germline!
	if ($help[$DBSBP] =~ /MATCH=exact/)
	{
		$indbSNP = 1;	# good for germline, bad for "somatic"
		# precious!
		#INFO=<ID=CLN,Number=0,Type=Flag,Description="SNP is Clinical(LSDB,OMIM,TPA,Diagnostic)">
		#INFO=<ID=PM,Number=0,Type=Flag,Description="SNP is Precious(Clinical,Pubmed Cited)">
		if ($help[$DBSBP] =~ /;CLN;/ || $help[$DBSBP] =~ /;PM;/)	# more likely somatic
		{
			$precious = 1;
		}
	}
# 2) annotations of regions that cause problems: some classes of repeats from RepeatMasker track, segmental duplications,
# (cf. Reumers et al. 2012, Nature Biotech 30:61), external blacklists, mapability 
	# simple repeats and low complexity (not the same as homopolymer, but similar enough); some satellites are not annotated in blacklist ...
	if ($help[$REPET] =~ /Simple_repeat/ || $help[$REPET] =~ /Low_/ || $help[$REPET] =~ /Satellite/)
	{
		$is_repeat = 1;
		$confidence-=2;
	}
	# other repeat elements to penalize at least a bit
	elsif ($help[$REPET] ne "." && $help[$REPET] ne "0")
	{
		$confidence--;
	}
	# simple tandem repeats most often coincide with other bad features - do not penalize twice
	if ($help[$STREP] ne "." && $help[$STREP] ne "0")
	{
		$is_STR = 1;
		unless ($is_repeat)
		{
			$confidence-=2;
		}
	}
	# Segmental Duplications are less effective than homopolymers, short tandem repeats and microsatellites, do not penality twice
	if ($help[$SEGDP] ne "." && $help[$SEGDP] ne "0" && ! $is_repeat && ! $is_STR)
	{
		$confidence-=2;	# bad region
		$is_weird = 1;
	}
	# Duke excluded and ENCODE DAC blacklist, only consider if not already annotated as suspicious repeat
	if (($help[$EXCLU] ne "." && $help[$EXCLU] ne "0") || ($help[$BLACK] ne "." && $help[$BLACK] ne "0"))
	{
		$confidence-=3;	# really bad region, usually centromeric repeats
		$is_weird = 1;
	}
	# HiSeqDepth: regions "attracting" reads; often coincide with tandem repeats and CEN/TEL, not always with low mapability
	if ($help[$HSDEP] ne "." && $help[$HSDEP] ne "0")
	{
		$confidence-=3;	# really really bad region!
		$is_weird = 1;
	}
	# Mapability is 1 for unique regions, 0.5 for regions appearing twice, 0.33... 3times, ...
	# Everything with really high number of occurences is artefacts
	# does not always correlate with the above regions
	# is overestimating badness bc. of _single_ end read simulations
	$mapp = $help[$MAPAB];
	if ($mapp eq ".")	# in very rare cases (CEN), there is no mapability => ".", which is not numeric but interpreted as 0
	{
		$confidence-=5;
	}
	else
	{
		if ($mapp < 0.5)	# 0.5 does not seem to be that bad: region appears another time in the genome and we have paired end data!
		{
			$confidence--;
			$is_weird = 1;	# something _is_ weird already there and known SNPs might be artefacts
			if ($mapp < 0.4)	#  3-4 times appearing region is worse but still not too bad
			{
				$confidence--;
				if ($mapp < 0.25)	# > 4 times appearing region
				{
					$confidence--;
				}
				if ($mapp < 0.1)	# > 5 times is bad
				{
					$confidence-=2;
					if ($mapp < 0.05)	# these regions are clearly very bad (Lego stacks)
					{
						$confidence-=3;
					}
				}
			}
		}
	}
	# if others have found the SNP already, it may be interesting despite low score - but only if it's not a weird region.
	# if a position gets up from score 4 to 6 by giving +2 for presence in dbSNP, it's very likely an artefact reported to dbSNP
	if (($in1KG || $indbSNP) && ! $is_weird)
	{
		$confidence++;
	}

# 3) information from the calls and germline comparisons: coverage, strand bias, variant support, ... => can lead to reclassification
# what we do not know: "custering", misplaced indels, in read end region (!overlapping reads) (cf. Reumers et al.)
	# parse INFO fields
	($trf, $trr, $tvf, $tvr) = $help[$INFO] =~/DP4=(\d+),(\d+),(\d+),(\d+);/;
	# the PV4 INFO may be missing!
	#($strandbias, $baseQbias, $mapQbias, $TDbias) = $help[$INFO] =~ /PV4=(.+),(.+),(.+),(.+)/;
	#if (defined $TDbias)	# SNPs near ends may be from low quality or even overlapping reads
	#{
		#if ($TDbias > 0.2)
		#{
			#$confidence--;
			#if ($TDbias > 0.5)	# a real bad one will have 1
			#{
			#	$confidence-=2;
			#}
		#}
		# strand bias is empirically estimated below, not clear how bad it is
		#if ($strandbias > 0.5)
		#{
		#	$confidence--;
		#}
	#}
	($crf, $crr, $cvf, $cvr) = $help[$CINFO] =~/DP5=(\d+),(\d+),(\d+),(\d+),\d+;/;
	($varsuppC) = $help[$CINFO] =~/TSR=(\d+)/;	# incl. low quality reads
	$depthT = ($trf + $trr + $tvf + $tvr);	# sum in DP4 - can become zero after filter_PEoverlap.py! thus add pseudocount
	if ($depthT < 1)
	{
		$depthT = 1;
	}
	$fr_var_tum = ($tvf + $tvr)/$depthT;
	$depthC = ($crf + $crr + $cvf + $cvr);	# this does not include low base quality reads - might be 0!
	if ($depthC)
	{
		$fr_var_ctrl = ($cvf + $cvr)/$depthC;
	}
	else
	{
		$fr_var_ctrl = 0;
	}
	# regions have an excessive number of reads in control, may not be recognized by the mapability and High Seq Depth tracks
	# (unclear if these tracks are from "normal" hg19 or 1KG reference with decoy sequences!)
	# tumor could have deletion here
	if ($depthC > $cutoff)
	{
		$confidence-=2;
		if ($depthC > $cutoff*2)	# >>300 reads at lego stack regions
		{
			$confidence-=2;
		}
	}
        # total coverage too low
	#if ($depthT <= 5 || $depthC <= 5)	# low coverage is bad in general
	#{
	#	$confidence-=3;	# "low", very bad!
	#}
	# split up in 2, to subtract 3 for tumor low coverage bc. that's so bad and additional penalty for low control depth
	if ($depthT <= 5)
	{
		$confidence-=3;
	}
	if ($depthC <= 5)	# not _that_ bad as depthT low - will very likely become reclassified as "lowCov_SNP_support_germline"
	{
		$confidence-=2;
	}
	if ($tvf < 1 || $tvr < 1)	# strand bias
	{
		$confidence-=2;	# "moderate"
		if ($fr_var_tum <= 0.1)	# and too low variant coverage
		{
			$confidence-=3;	# "low"
		}
	}
	if ($tvf + $tvr <= 5)	# low variant support
	{
		$confidence-=2;	# "moderate"
	}
	if ($tvf == 1 && $tvr == 1)   # only 1 read from each strand (might be PCR error of overlapping reads!) ### TODO: But overlaps are removed???
	{
		$confidence-=3;	# "low"
	}
	if ($tvf + $tvr < 6 && $fr_var_tum < 0.1)  # too low variant coverage
	{
		$confidence-=3;	# "low"
	}
	if ($help[$CLASS] =~/somatic/)	# this includes multi_somatic
	{
		# the overlapping reads filter may have reduced the variant read to 1 - before, such SNVs were never considered
		if ($tvf + $tvr < 2)
		{
			$confidence-=3;	# very low!
		}
		if ($fr_var_ctrl > 0)	# might be germline or artefact! 1/30 reads required for calling germline
		{
			$class = "LQVSIG";
			$confidence-=2;	# "moderate"
			if ($fr_var_tum < 0.3)
			{
				$confidence-=3;	# "low"
			}
			#if ($tvf < 1 || $tvr < 1)	# already penalized in general above with -2
			#{
			#	$confidence-=3;	# "low"
			#}
		}
		if ($depthC && $varsuppC / $depthC > 1 / 30)	# present in germline but with low qual, so likely seq error ### WAS: $varsuppC > 2. THIS WAS TOO STRINGENT FOR HIGH COVERAGE DATA
		{
			$class = "LQVSIG";
			$confidence-=3;	# "low"
		}
		# an SNV that is in dbSNP but not "precious" or/and in 1 KG with high frequency is probably germline,
		# especially if control coverage is low (may be due to library problems)
		if ($in1KG || ($indbSNP && ! $precious))
		{
			$class = "SNP_support_germline";	# missed germline?
			if ($depthC < 10)	# Matthias' suggestion that below 10 (and not <= 5 as before) the probability to miss the variant is >95%
			{
				$class = "lowCov_SNP_support_germline";
			}
		}
	}
	# in the classification, 1/30 reads required for calling germline
        if ($help[$CLASS] =~ /germline/)
	{
		if ($fr_var_ctrl < 0.3)	# too far from 50%, but that depends on coverage. I'm not convinced that qualifies for "low" / -3
		{
			$confidence-=2;	# "low"
		}
		if ($in1KG || ($indbSNP && ! $precious))	# but this supports again - number of reads may be low!
		{
			$class .= "_SNP_support";
		}
		if ($depthC <= 10)	# very probably germline due to skewed distribution at low coverage
		{
			$class .= "_lowCov";	# => can end up as "germline_SNP_support_lowCov"
		}
	}
	if ($help[$CLASS] =~ /LQVSIG/)	# low quality variant support in germline: usually a weird region; a rare case anyways
	{
		if ($in1KG || ($indbSNP && ! $precious))
		{
			$class = "SNP_support_germline";
		}
		else	# weird region
		{
			$confidence-=2;
		}
	}
	if ($help[$CLASS] =~ /unclear/)
	# < 1/30 of the high quality supporting the variant AND other bases present
	# In rare cases, unclear could be LOH in tumor at a multiallelic SNP position, but we
	# cannot check that without SNV calling in control (done in LOH-BAF pipeline)
	{
		if ($varsuppC > 1)
		# < 1/30 of the high quality supporting the variant, or only low quality ones
		# these positions are mostly in regions of extremely high coverage (HiSeq depth track)
		# bc. of targetted sequencing, a coverage threshold cannot be used here.
		{
			$class = "LQVSIG";
		}
		if ($in1KG || ($indbSNP && ! $precious))
		{
			$class = "SNP_support_germline";
		}
		else	# just a weird region
		{
			$confidence-=2;
		}
	}
	if ($help[0] =~ /Y/ && $help[5] < 100)	# SNV on chrY would have to be homozygous, too
						# het SNV on male chrX is another weird thing (unless a subpopulation had amplified chrX and mutated)
						# any SNV on chrY in female is an artefact
	{
		$confidence-=3;	# "low"
	}
	# to make sure that the confidence score is positive:
	if ($confidence < 1)
	{
		$confidence = 1;
	}
	# and not > 10
	if ($confidence > 10)
	{
		$confidence = 10;
	}
	if ($precious)
	{
		$class.= "_precious";
	}
	print $line, "\t$confidence\t$class\n";
}
close FH;
exit;
