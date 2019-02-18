#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

use strict;
use warnings;
use List::Util qw(min max);

# empirical confidence classification

if (@ARGV < 2)
{
	die "USAGE: (1)vcf formatted file for confidence classification - (2)Config file from which to get the column names\n";
}

my $file = shift;
my $configfile = shift;

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
print $header, "\tCONFIDENCE\tCLASSIFICATION\n";
# get the columns where to look for features:
# INFO => DP4 of tumor => variant frequency and strand bias
my $INFO = 0;
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
my $ExAC  = 0;
my $EVS   = 0;
my $CILC  = 0;
@help = (split "\t", $header);
for (my $c = 0; $c < @help; $c++)
{
	# fixed column labels (from vcf_pileup_compare)
	if ($help[$c] eq "INFO")
	{
		$INFO = $c;
		print STDERR "INFO in column $c\n";
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
	if ($help[$c] eq "ExAC")
   	{
   		$ExAC = $c;
   		print STDERR "ExAC in column $c\n";
   	}
 	if ($help[$c] eq "EVS")
   	{
       	$EVS = $c;
       	print STDERR "EVS in column $c\n";
    }
    if ($help[$c] eq "CountInLocalControl")
    {
    	$CILC = $c;
    	print STDERR "CountInLocalControl in column $c\n";
    }
}

my $line = "";	# current line

# number of reference and variant reads on forward and reverse complement strand for sample (tumor)
my ($trf, $trr, $tvf, $tvr) = (0, 0, 0, 0);

# absolute coverage at position
my $depthT = 0;

# fraction of variant reads in tumor
my $fr_var_tum = 0;
my $confidence=10;	# start with maximum value and subtract something for each "bad" feature
			# "high" = 9-10, "medium" = 6-8, low <= 5
my $mapp = 0;	# Mapability
my $is_repeat = 0;	# true if SNP conicides with any of the suspicious repeat classes simple, low and satellite, which are partially redundant to other annotations
my $is_STR = 0;	# also short tandem repeats from Tandem Repeats Finder are suspicious and can conincide with other classes
my $is_weird = 0;	# coindicende with known artefact-producing regions
my $in1KG = 0;
my $indbSNP = 0;
my $precious = 0;
my $class = "";	# for germline/somatic classification (e.g. in dbSNP => probably germline)

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
	$class = "somatic";	# start with default somatic

# 1) external information of if these SNPs have already been found (incl. false positives from 1000 genomes!)
	# 1000 genomes
	#if ($help[$IN1KG] ne "." && $help[$IN1KG] ne "0")
	if ($help[$IN1KG] =~ /MATCH=exact/)	# verified somatic SNVs have MATCH=position instead of exact, so it would be bad to classify them as SNP_support_germline!
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
	### Nov 15, 2013: mouse /icgc/ngs_share/assemblies/mm10/databases/UCSC/wgEncodeCrgMapabilityAlign100mer_lifted_from_mm9.bedGraph.gz
	### has overlapping regions => use lowest score of the concatenated ones, as done in confidenceAnnotation_indels.pl
		if ($mapp =~ /&/)
		{
			my @mappab = split ("&", $mapp);
			$mapp = min @mappab;	# just chose one - here: min
		}
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
	# an SNV that is in dbSNP but not "precious" or/and in 1 KG with high frequency is probably germline
	if ($in1KG || ($indbSNP && ! $precious))
	{
		$class = "SNP_support_germline";
	}

# 3) information from the calls: coverage, strand bias, variant support, ...
# what we do not know: "custering", misplaced indels, in read end region (!overlapping reads) (cf. Reumers et al.)
	# parse INFO fields
	($trf, $trr, $tvf, $tvr) = $help[$INFO] =~/DP4=(\d+),(\d+),(\d+),(\d+);/;
	$depthT = ($trf + $trr + $tvf + $tvr);	# sum in DP4 - better take DP?
	$fr_var_tum = ($tvf + $tvr)/$depthT;
	# for whole genomes, the absolute values for low coverage/allele frequency make sense, but for
	# exomes etc. with high coverage, fractions would make much more sense
	if ($depthT <= 5)
	{
		$confidence-=3;
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
	if ($tvf == 1 && $tvr == 1)   # only 1 read from each strand (might be PCR error of overlapping reads!)
	{
		$confidence-=3;	# "low"
	}
	if ($tvf + $tvr < 6 && $fr_var_tum < 0.1)  # too low variant coverage
	{
		$confidence-=3;	# "low"
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
