#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

# "somatic call": tumor vcf from mpileup+bcftools+filter <-> control from normal pileup
# call somatic if variant is in tumor and control
# call germline if the variant is only there in tumor
# report numbers of high-quality bases for reference and variant allele, separated by strand (as mpileup for tumor)
# correct status: if only low quality bases support the variant, only 1 low-quality base for a different allele, ...
# also include a general nucleotide counter

use strict;
use warnings;
use feature "switch";
use Math::CDF qw(pbinom);

if (@ARGV < 2)
{
	die "$0 (1)vcf file for tumor - (2)normalpileup for control - optional(3)header (default no)\n";
}

my $tumor = shift;
my $control = shift;
my $makeheader = shift;

open (T, $tumor) or die "Could not open $tumor: $!\n";
open (C, $control) or die "Could not open $control: $!\n";

# hardcoded:
# at least 1/30th of reads must support the variant in control to call it germline
my $fraction = 30;
# mpileup counts only the bases with score >= 13, so we want to use that cutoff, too
my $cutoff = 13;
# for RNA, cutoff for calling monoallelic (% of high quality reads). Not 100% bc. of sqcing errors and normal contamination; at least 10 reads
my $limit = 0.85;
my $rcutoff = 10;
# allow 1 high quality base that is neiter reference nor variant
my $allowedothers = 1;

# change this for RNA or 4fold comparison
my $somatic = "somatic";
my $germline = "germline";
my $unclear = "unclear";
my $notcovered = "not_covered";
my $ctrlname = "control";

# p-value when assuming binomial distribution: gives probability to observe maximally that many variant supporting reads in the germline under the condition that germline is heterozygous variant
my $pval;


# current lines of tumor and control
my $lineT = "";
my $lineC = "";
# splitting the lines by tabs:
my @tum;
my @ctrl;
# count all kinds of things
# variant allele taken from tumor must be present in same way in control
my $varT = "";
# for indels, to get exact string of inserted/deleted bases
my $difference = 0;
# different allele than reference or variant present in control
my $ctrT = 0;	# lines in tumor
my $ctrC = 0;	# lines in control
my $match = 0;	# the pairs found
my $trueSNP = 0;	# clear SNP positions
my $multi = 0;	# tri- or quadallelic SNPs, most are probably sequencing errors
# these have several alleles:
my @alleles = ();
my $is_multi = 0;
# number of insertions and deletions
my $insert = 0;
my $delet = 0;
# for all variants, assign status:
my $somatic_ctr = 0;
my $germline_ctr = 0;
my $unclear_ctr = 0;	# suspicious cases in which there is another allele than reference or variant present in control
my $artefact_ctr = 0;	# likely artefacts: "LQVCIG"
my $expr_mono = 0;
my $missingC = 0;	# no report for the tumor position in control pileup
# get the variant allele in upper case letters
my $varuc = "";
# also in lower case to count occurences on reverse complement in control
my $varlc = "";
# get the bases for control pileup
my $basestr = "";
# and the corresponding quality scores
my $qualstr = "";
my @quals = ();
# a single base or an indel to compare to the variant
my $indel = "";
# check if a . or , is followed by + or -, then it's part of the indel
my $plusmin = "";
#the ASCII is re-encoded into Sanger score by (ord($_)-33)
my $qstr = 0;
# variant allele frequency (fraction)
my $allelefrq = 0;
# numbers of reads with base quality >= cutoff supporting: reference forward, r. revcomp, variant forward, variant revcomp - and others
my ($rfw, $rrc, $vfw, $vrc, $other) = (0, 0, 0, 0, 0);
# total numbers of reads in above categories _without_ looking at quality
my ($rfwNQ, $rrcNQ, $vfwNQ, $vrcNQ, $otherNQ) = (0, 0, 0, 0, 0);
# explicit counting of AGCTNacgt bases w.r.t. quality
my @bases = qw (A C G T N a c g t n);
# the reference never appears literally in the pileup, only as . and ,
my $refuc = "";
my $reflc = "";
# assigned status (germline, somatic, unclear, no_coverage)
my $newanno = "";
# has status changed?
my $changed = 0;
# for indels:
my $is_indel = 0;
# number of ins/del in control pileup
my $num = 0;
# to advance correctly in basestring
my $offset = 0;
# if the indel is > 9:
my $plus = 0;

# finally start parsing!
# first get header line
my $header = "";
if (defined $makeheader && $makeheader ne "no")
{
	while ($header = <T>)
	{
		last if ($header =~ /^\#CHROM/); # that is the line with the column names
		print "$header";
	}
	chomp $header;
	print "$header\tINFO_$ctrlname(VAF=variant_allele_fraction;TSR=total_variant_supporting_reads_incl_lowqual)\tANNOTATION_$ctrlname\n";
}
my $ccoord = -1; # make sure that ccoord in first iteration is smaller than tcoord
my $tcoord = 0;

# one file is longer than the other, nevermind which
#(if we did pileup for only the tumor SNV positions, control will be shorter; otherwise longer)
TUMORFILE_LOOP: while ($lineT=<T>)
{
	# after reading over the comments, see above, we start with a line containing chromosome and coordinate etc.
	# in tumor:
	#chr1	16566	.	A	C	10.7	.	DP=11;AF1=0.6642;CI95=0.5,0.75;DP4=0,4,0,3;MQ=30;PV4=1,0.14,0.096,0.012	PL:GT:GQ	29,3,0:1/1:5	16,0,82:0/1:22	
	# in control:
	#chr1	16566	A	42	,.,,.,..,,,,,.,,,,,,,,,,,c,,,,.,c,,cc,^1,^!,^!,^!,	IIII6I-)CFIII+9I,IIIFII)IHB/IGIII9IIIBI52;
	# remove newline
	$ctrT++;
	chomp $lineT;
	@tum = split (/\s+/, $lineT);
	$tcoord = $tum[1];	# start coordinate
	if ($tcoord < $ccoord)	# no matching control line => tumor position not covered in control (never true in 1st iteration)
	{
		$missingC++;
		# print the tumor line with according "flag"
		# to have same number of columns, put "." in place of chr, pos, ref, number of reads, variants, qual
		print $lineT, "\tDP=0;DP5=0,0,0,0,0;DP5all=0,0,0,0,0;ACGTNacgtnHQ=0,0,0,0,0,0,0,0,0,0;ACGTNacgtn=0,0,0,0,0,0,0,0,0,0;VAF=0;TSR=0\t$notcovered\n";
		next;	# read next line from tumor
	}
	if ($tcoord == $ccoord)	# matching pair found!
	{
		$match++;
		chomp $lineC;
		# split lines
		@ctrl = split (/\s+/, $lineC);
		$ccoord = $ctrl[1];
		# do the evaluation in a subroutine
		evaluate_pos();
		next;	# read next line from tumor
	}
	# else missing in tumor - can never be the case!
	while ($lineC = <C>)	# when we are here tcoord is higher than ccoord; read new control line
	{
		chomp $lineC;
		$ctrC++;
		# split lines
		@ctrl = split (/\s+/, $lineC);
		$ccoord = $ctrl[1];
		if ($tcoord == $ccoord)	# matching pair found!
		{
			$match++;
			# do the evaluation in a subroutine
			evaluate_pos();
			next TUMORFILE_LOOP; # and read next line from tumor
		}
		if ($tcoord < $ccoord)	# new ccoord is higher than tcoord =>  tumor position not covered in control (should never be the case!)
		{
			$missingC++;
			print $lineT, "\tDP=0;DP5=0,0,0,0,0;DP5all=0,0,0,0,0;ACGTNacgtnHQ=0,0,0,0,0,0,0,0,0,0;ACGTNacgtn=0,0,0,0,0,0,0,0,0,0;VAF=0;TSR=0\t$notcovered\n";
			next TUMORFILE_LOOP; # read next line from tumor
		}
		# else missing in tumor - can never be the case!
	}
	# when both if conditions in the inner loop fail: next iteration (read control line)
	# when this is reached the control file has ended
	# everything remaining in the tumor file is not covered in control
	# print the tumor line with according "flag"
	print $lineT, "\tDP=0;DP5=0,0,0,0,0;DP5all=0,0,0,0,0;ACGTNacgtnHQ=0,0,0,0,0,0,0,0,0,0;ACGTNacgtn=0,0,0,0,0,0,0,0,0,0;VAF=0;TSR=0\t$notcovered\n";
	$missingC++;
}
close C;
close T;
# if there are lines left in control pileup, these are not interesting anyways
print STDERR "$ctrT lines of tumor and $ctrC of control processed, $match are paired. After filtering, $trueSNP normal SNPs, $multi multiallelic SNPs, $insert insertions, $delet deletions. $germline_ctr $germline, $somatic_ctr $somatic, $unclear_ctr with different allele, $artefact_ctr 'LQVCIG' possible artefacts, $missingC missing in control.\n";

# check if the variant is also present in the control and assign corresponding status
sub evaluate_pos
{
	# tumor: the variant to search in control, reads forward and reverse complement, read depth
	#chr1	16566	.	A	C	10.7	.	DP=11;AF1=0.6642;CI95=0.5,0.75;DP4=0,4,0,3;MQ=30;PV4=1,0.14,0.096,0.012	PL:GT:GQ	29,3,0:1/1:5	16,0,82:0/1:22
	# simple procedure for SNPs does not work for insertions or deletions
	# 1) the "DP" is at a different position
	# 2) the variant can't be detected by simple string search in control pileup
	# => special treatment below
	$varT = $tum[4];
	# control pileup
	# . and , stand for the reference, which also needs to be counted explicitely as a base
	$refuc = uc ($tum[3]);
	$reflc = lc ($tum[3]);
	#print STDERR "reference: $refuc $reflc\n";
	# quality scores are in the last column; $#ctrl gives index of last column
	$qualstr = $ctrl[$#ctrl];
	# split up:
	@quals = split ("", $qualstr);
	# bases are now in the second-to-last column:
	$basestr = $ctrl[$#ctrl-1];
	# remove the mapping qualities for starting reads (^1 etc.), and the $ signs for read ends
	$basestr =~ s/\^.{1}//g;
	$basestr =~ s/\$//g;
	#print STDERR "after ^ removal - $basestr\n";
	# explicit counting of AGCTNacgt bases w.r.t. quality
	my %bases_highq;
	my %bases_all;
	if ($lineT !~ /INDEL/)
	{
		# finally: look at control, variant count (variant base is given by $varT=$tum[4])
		#chr1	16566	A	42	,.,,.,..,,,,,.,,,,,,,,,,,c,,,,.,c,,cc,^1,^!,^!,^!,	IIII6I-)CFIII+9I,IIIFII)IHB/IGIII9IIIBI52;
		$is_indel = 0;
		if ($varT =~ /,/)	# A complicated case: heterozygous with 2 (or even 3!) possibilites where none of the alleles is reference
					# The most frequent base, which is probably the only true variant (reference genome is wrong or
					# the alternative bases are sequencing errors), is listed first
		{
			$is_multi = 1;
			$multi++;
			@alleles = split (",", $varT);
			# get the most frequent variant allele in both upper and lower case to search
			$varuc = uc ($alleles[0]);
			$varlc = lc ($alleles[0]);
		}
		else	# the simple case, just counting
		{
			$trueSNP++;
			$is_multi = 0;
			# get the variant allele in both upper and lower case to search
			$varuc = uc ($varT);
			$varlc = lc ($varT);
		}
	}
	else	# indel: special treatment
	# example insertion:
	#chr22   17333902        .       caaaaaaaaaaa    cAaaaaaaaaaaa   73.6    .    INDEL;DP=21;AF1=1;CI95=0.8,1;DP4=0,1,8,8;MQ=58;PV4=1,1,0.4,0.18 ...
	# example deletion:
	#chr22   24031480        .       TGGGGGGGG       TGGGGGGG        69.2    .    INDEL;DP=12;AF1=1;CI95=0.8,1;DP4=0,0,5,5;MQ=50 
	{
		# might also be multiallelic, in this case just take the most frequent one:
		$tum[4] =~ /(\w+)/;
		$tum[4] = $1;
		$is_indel = 1;
		if (length ($tum[3]) < length ($tum[4]))	# insertion.
		# The inserted nucleotide(s) start at 2nd position of variant: caaaaaaaaaaa <-> cAaaaaaaaaaaa
		# in the control pileup:
		#chr22   17333902        c       36      ,,,...+1A,+1a,+1a.+1A,+1a,+1a,+1a.+1A.+1A,+1a,+1a,+1a,+1a.+1A,+1a,+1a,+1a,+1a,+1a,+1a.+1A,+1a.+1A.+1A.+1A.+1A.+1A,+1a.+1A,+1aa  EIII#F<6GA9FIG80'#II:4#-#IGEIIGIAH##
		# if there were >1 nucl. inserted: longer variant
		#e.g.
		#chr22	16052167	.	aaaacaaacaaacaaacaaacaaacaaacaaac	aAAACaaacaaacaaacaaacaaacaaacaaacaaac	99	.	INDEL;DP=12;AF1=0.9999;CI95=0.7,1;DP4=0,0,3,5;MQ=42	PL:GT:GQ	84,6,0:1/1:45	0,0,0:1/1:39	66,6,0:1/1:45	76,9,0:1/1:48	29,3,0:1/1:42
		# in control:
		#chr22	16052167	a	23	,,,..,,.,..,.,.,....+4AAAC,+4aaac.+3AACC	DIGGDEIGH=5IGFDGEEGH2GI
		# the actual variant is there twice, the +3AACC won't be detected in any way
		# here, TAT was inserted, but control has a:
		#chr22	23103420	.	TTAt	TTATAt	99	.	INDEL;DP=19;AF1=0.9999;CI95=0.7,1;DP4=0,0,6,4;MQ=58
		#chr22	23103420	T	17	,,..,..,,,...a.AA	GGIIF?ACIACEH?DIG
		{
			# get the length of the insertion for correct string to search:
			$difference = length ($tum[4]) - length ($tum[3]);
			$varuc = uc (substr($tum[4], 1, $difference));
			$varlc = lc ($varuc);
			# to search for exactly that string in the pileup - otherwise would
			# consider "A", "+1A" and "-1A" as support for insertion (or deletion) of A!
			$varuc = "+".$difference.$varuc;
			$varlc = "+".$difference.$varlc;
			#print STDERR "Insertion of $varuc\n";
			$insert++;
		}
		else	# deletion. The deleted nucleotide(s) start at 2nd position of reference: TGGGGGGGG <->  TGGGGGGG
		# in the control pileup:
		#chr22   24031480        T       10      *,-1g.-1G,-1g.-1G,-1g.-1G,-1g,-1g.-1G   #A1>8I&6(D
		# if there were >1 nucl. deleted: longer variant:
		#chr22	16080452	.	atatattatat	atatat	99	.	INDEL;DP=21;AF1=0.518;CI95=0.4,0.6;DP4=3,10,2,4;MQ=39;PV4=1,0.16,1,1	PL:GT:GQ	0,3,29:0/1:3	16,0,142:0/1:19	17,0,116:0/1:20	97,6,0:1/1:5	48,0,82:0/1:51
		#control:
		#chr22	16080452	a	19	.$,,.,,,,,,,.,,,-5tatat,,-5tatat,-5tatat.	#@DDIGIGIIG2IIIHI7D
		# looks complicated:
		#chr22	34027139	.	ACTGTCT	ACT	57.6	.	INDEL;DP=17;AF1=0.2164;CI95=0.2,0.3;DP4=10,2,1,1;MQ=60;PV4=0.4,1,1,0.21	PL:GT:GQ
		# control here unclear:
		#chr22	34027139	A	12	.....,,.,,.g	C=CD@E#G##E#
		# another problem:
		#chr22	36274088	.	caca	cca	97.1	.	INDEL;DP=20;AF1=0.4782;CI95=0.4,0.6;DP4=5,1,8,1;MQ=57;PV4=1,0.28,1,1	PL:GT:GQ	22,3,0:0/1:3	59,0,23:0/1:26	17,0,40:0/1:20	53,0,22:0/1:25	0,6,57:0/0:5
		# the control pilup does _not_ report "-1a." but * - does that mean that the reference c itself is deleted?
		#chr22	36274088	c	23	,...****.*.*.*...,..,*.
		{
			# get the length of the deletion for correct string to search:
			$difference = length ($tum[3]) - length ($tum[4]);
			$varuc = uc (substr($tum[3], 1, $difference));
			$varlc = lc ($varuc);
			# to search for exactly that string in the pileup:
			$varuc = "-".$difference.$varuc;
			$varlc = "-".$difference.$varlc;
			#print STDERR "Deletion of $varuc\n";
			$delet++;
		}
	}
	# look for the variant(s) in control pileup
	$offset = 0;
	$plus = 0;
	for (my $p=0; $p < @quals; $p++)
	{
		#print STDERR "Pos. $p - ";
		# the ASCII is re-encoded into Sanger score by (ord($_)-33)
		$qstr = ord($quals[$p])-33;
		#print STDERR "ASCII: $quals[$p] => $qstr for ";
		# the base or symbol at the current position
		# always check if there is a + or - after a . or ,
		# bc. that signals beginning of indel (and the . or , belongs to the indel!)
		$indel = substr ($basestr, ($p+$offset), 1);
		$plusmin = substr ($basestr, ($p+$offset+1), 1);
		if ($plusmin eq "+" || $plusmin eq "-") 
		{
			# get next pos: the number of inserted/deleted. May be > 9! But at least not > 99 for gapped aligners
			# that don't do split read mapping (N operations)
		 #chrX    3696728 .       tgt     tGTCCGTAATGgt   70.4    .       INDEL;DP=16;AF1=1;CI95=0.5,1;DP4=0,0,1,2;MQ=43  PL:GT:GQ        110,9,0:1/1:63  chrX    3696728 t       15      ..,$.,.,,,,,.+10GTCCGTAATG.+10GTCCGTAATG.,      ```````````Q```     DP=15;DP5=5,8,2,0,0     germline        2

			$num = int(substr ($basestr, ($p+$offset+2), 1));
			if (substr ($basestr, ($p+$offset+3), 1) =~ /\d+/)
			{
				$plus = 1;
				# 1 and 0 must be 10 (not 1!)
				$num .= (substr ($basestr, ($p+$offset+3), 1));
			}
			#print STDERR "Indel: $indel $num ";
			# to search for exactly that string in the pileup:
			$indel = substr ($basestr, ($p+$offset+3+$plus), $num);
			$indel = $plusmin.$num.$indel;
			# correct offset: the + or -, the number, and the inserted/deleted base
			# the . or , following an indel is _not_ the next base but part of the indel
			$offset+=($num+2+$plus);
			$plus = 0;
		}
		#print STDERR "base $indel\n";
		# now consider qualites:
		if ($qstr >= $cutoff)
		{
			#print STDERR "$indel - highQ\n";
			given ($indel)
			{
				# the reference never appears literally, only as . and ,
				when (".") { $rfw++; $bases_highq{$refuc}++; }
				when (",") { $rrc++; $bases_highq{$reflc}++; }
				when ($varuc) { $vfw++; $bases_highq{$indel}++; }
				when ($varlc) { $vrc++; $bases_highq{$indel}++; }
				default    { $other++; $bases_highq{$indel}++; }
			}
		}
		# also count low quality bases
		given ($indel)
		{
			when (".") { $rfwNQ++; $bases_all{$refuc}++; }
			when (",") { $rrcNQ++; $bases_all{$reflc}++; }
			when ($varuc) { $vfwNQ++; $bases_all{$indel}++;}
			when ($varlc) { $vrcNQ++; $bases_all{$indel}++;}
			default    { $otherNQ++; $bases_all{$indel}++;}
		}
	}

	# depth for only high quality reads
	my $hqdepth = $vfw+$vrc+$rfw+$rrc+$other;
	if ($hqdepth > 0)	# yes it can be zero!
	{
		$allelefrq = sprintf ("%.2f", (($vfw+$vrc)/$hqdepth));
	}
	else
	{
		$allelefrq = 0;
	}
	# total reads supporting the variant in control, incl. low quality
	my $supp_var = $vfwNQ + $vrcNQ;
	# how many reads for which base with base quality / in general
	my $ACGTNacgtnHQ = "";
	my $ACGTNacgtn = "";
	for (my $b = 0; $b < 10; $b++)
	{
		if (defined $bases_highq{$bases[$b]})
		{
			$ACGTNacgtnHQ.=$bases_highq{$bases[$b]}.",";
		}
		else
		{
			$ACGTNacgtnHQ.= "0,";
		}
		if (defined $bases_all{$bases[$b]})
		{
			$ACGTNacgtn.=$bases_all{$bases[$b]}.",";
		}
		else
		{
			$ACGTNacgtn.= "0,";
		}
	}
	# remove last , from the count strings
	chop $ACGTNacgtnHQ;
	chop $ACGTNacgtn;

# 	if ($is_indel)	# whenever anything, even low quality, supports the indel, it's germline
# 	{
# 		if ($vfwNQ > 0 || $vrcNQ > 0)
# 		{
# 			$newanno = $germline;
# 			$germline_ctr++;
# 		}
# 		elsif ($other)	# not support for variant but a different high-qual base that is neither ref nor var
# 		{
# 			$newanno = $unclear;
# 			$unclear_ctr++;
# 		}
# 		else	# neither variant nor 3rd option found
# 		{
# 			$newanno = $somatic;
# 			$somatic_ctr++;
# 		}
# 	}
# 	elsif (! $is_multi)	# only if a high qual base supports the variant, it's germline
# new 10.12.2013: no special treatment for indels anymore!

	if (! $is_multi)	# only if a high qual base supports the variant, it's germline
	{
		#if ($vfw > 0 || $vrc > 0) Barbara's original Sept 2011     $rfw,$rrc
		#if ((($rfw+$rrc)>50) && ((($vfw + $vrc)/($rfw+$rrc+$vfw + $vrc)) < 0.05))
		#if ((($rfw+$rrc)>100) && (($vfw + $vrc) <= 5))
		my $allowedGermlineVariants = int(($rfw+$rrc)/$fraction);
		# this is problematic in case of e.g. CLL when tumor blood cells are in control blood => false negative
		# => only call germline if the VAF is ~ 0.5 (Binomial distribution?)
		if (($vfw + $vrc) > $allowedGermlineVariants)
		{
			$newanno = $germline;
			$germline_ctr++;
		}
		elsif ($other > $allowedothers)	# no support for variant but a different high-qual base that is neither ref nor var
				# to reduce number of these, filter in a way such as "if < 10% of all reads are 'other', still call it somatic"
		{
			$newanno = $unclear;
			$unclear_ctr++;
			# the reference allele was always there in tested cases, even if not the most frequent allele!
		}
		else	# neither variant nor 3rd option found, or fewer than $allowedGermlineVariants
		{
			# likely artefacts: neither somatic nor germline, but per chance the same sequencing error in both samples
			# LQVCIG: "low quality variant call in germline"
			# DP=11;AF1=0.6642;CI95=0.5,0.75;DP4=0,4,0,3;MQ...
			my ($dp4f, $dp4r) = $tum[7] =~ /;DP4=\d+,\d+,(\d+),(\d+);/;
			# 1. the quality score assigned by bcftools is < 20 and the sum of support in control DP5 is > 3
			if ($tum[5] < 20 && $supp_var > 3)
			{
				$newanno = "LQVCIG";
				$artefact_ctr++;
			}
			# 2. the quality score assigned by bcftools is < 3, only one variant on one of the strands, i.e. DP4[2] == 1 or DP4[3] == 1, and the sum of support in control DP5 is > 1
			elsif ($tum[5] < 3 && $supp_var > 1 && ($dp4f == 1 || $dp4r == 1))
			{
				$newanno = "LQVCIG";
				$artefact_ctr++;
			}
			# 3. variant is only on one strand, but control has > 1 support on the _opposite_ strand in DP5: DP4=5,2,0,4; tumor and DP5all=21,6,0,8,0 for control, so seems like seq error since only low qual reads in control show the variant
			elsif (($dp4f == 0 && $vrcNQ > 1) || ($dp4r == 0 && $vfwNQ > 1))
			{
				$newanno = "LQVCIG";
				$artefact_ctr++;
			}
			else	# it passes as somatic
			{
				$newanno = $somatic;
				$somatic_ctr++;
			}
		}
	}
	else	# multiallelic; if the most frequent base is _not_ found, call it "somatic". The other base(s) probably sequencing error
		# best test again later with pileup of tumor; other bases in control will be reported in DP5 anyways, just not separated
	{
		if ($vfw > 0 || $vrc > 0)
		{
			$newanno = "multi_".$germline;
			$germline_ctr++;
		}
		elsif ($other > $allowedothers)	# to be consistent
		{
			$newanno = "multi_".$unclear;
			$unclear_ctr++;
		}
		else
		{
			$newanno = "multi_".$somatic;
			$somatic_ctr++;
		}
	}
	$pval = pbinom($vfw+$vrc, $rfw+$rrc+$vfw+$vrc+$other, 0.5);
	print $lineT, "\tDP=$ctrl[$#ctrl-2];DP5=$rfw,$rrc,$vfw,$vrc,$other;DP5all=$rfwNQ,$rrcNQ,$vfwNQ,$vrcNQ,$otherNQ;ACGTNacgtnHQ=${ACGTNacgtnHQ};ACGTNacgtn=$ACGTNacgtn;VAF=$allelefrq;TSR=$supp_var;PBINOM=$pval\t$newanno\n";
	# reset global variables
	($rfw, $rrc, $vfw, $vrc, $other) = (0, 0, 0, 0, 0);
	($rfwNQ, $rrcNQ, $vfwNQ, $vrcNQ, $otherNQ) = (0, 0, 0, 0, 0);
	return 0;
}
