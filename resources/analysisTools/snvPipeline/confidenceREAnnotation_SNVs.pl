#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

use strict;
use warnings;
# use Getopt::Std;
use Getopt::Long;
use List::Util qw(min max);
use POSIX qw(strftime);

# empirical confidence classification, version with option to have no penalties for "bad control events" and output reasons
# modifications inserted to get a vcf file that is conform to the ICGC pancancer regulations

### Options including default values (also for header)

my $makehead = 1;
my $panCanOut;
my $file;
my $cutoff = 150;
my $controlbad = 0;
my $print_info = 0;
my $newpun = 5;
my $refgenome = "hs37d5,ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz";
my $center = "DKFZ";
my @additionalHeader = "";
my $help;
my $pid = "NA";
my $somout;
my $round = 0;
my $runlowmaf = 0;
my $punPCR = 1;
my $punSEQ = 1;
my $runExome = 0;



GetOptions (	'makehead|m=i'		=> \$makehead,			# Set to 1 if you want to produce a pancancer conform head
				'panCanOut|o=s'		=> \$panCanOut,			# Outfile name including path for PanCan, will only be used if set
				'infile|i=s'		=> \$file,				# Infile, vcf file after annotations (will be piped in when using the standard pipeline)
				'cutoff|t=i'		=> \$cutoff,			# lower threshold for too high read depth in control (default 150; set to 500 for exomes etc.), only effective with -c 1
				'punishControl|c=i'	=> \$controlbad,		# lower confidence score if there are bad control events (default 0)
				'printpenalty|p=i'	=> \$print_info,		# print info on penalties into additional column 'PENALTIES' at the end (default 0)
				'refgenome|r=s'		=> \$refgenome,			# reference genome used for calling ID, path (default hs37d5,ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz)
				'center|z=s'		=> \$center,			# Center (unclear if the center where the data was produced or the center of the calling pipeline; default DKFZ)
				'pid|d=s'			=> \$pid,				# Patient ID (default NA)
				'newPunish|n=f'		=> \$newpun,		# Punish if alternative is supported with low allele frequency and the alternative allel is present in the control (default = 4 times)
				'addHead|s=s{,}'	=> \@additionalHeader,	# String with additional header line infer multiple times for multiple additional lines
				'somout|f=s'		=> \$somout,			# File with the high confidence somatic SNVs
				'round|a=i'			=> \$round,				# Round of iteration for filtering (default = 0)
				'punPCR|r=i'		=> \$punPCR,			# Run the PCR bias filter (should always be 1)
				'punSEQ|e=i'		=> \$punSEQ,			# Run the sequencing bias filer (should be 0 for exomes)
				'runlowmaf|l=i'		=> \$runlowmaf,			# Set to 1 if you want to run the low maf punishment
				'runExome|x=i'		=> \$runExome,			# Run on exome, will turn off the high control coverage punishment and the PCR bias filter
				'help|h|?'			=> \$help				# print help/usage
) or die "Could not get the options!\n";

if (defined $help || ! defined $file)
{
	die "USAGE: $0 [options]
	-m/--makehead [0/1]		Set to 1 if you want to produce a pancancer conform head (default 1)
	-o/--panCanOut		Outfile name including path for PanCan, will only be used if set will be zipped when ending with .gz and if bgzip is installed
	-i/--infile FILE		SNV file in vcf format for confidence (re)classification, required (set to - for pipe from STDIN)
	-t/--cutoff INT			lower threshold for too high read depth in control (default 150; set to 500 for exomes etc.), only effective with -c 1
	-c/--punishControl [0/1]	lower confidence score if there are bad control events (default 0)
	-p/--printpenalty [0/1]		print info on penalties into additional column 'PENALTIES' at the end (default 0)
	-r/--refgenome PATH		reference genome used for calling \"ID,path\" (default hs37d5,ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz)
	-z/--center NAME		Center (unclear if the center where the data was produced or the center of the calling pipeline; default DKFZ)
	-d/--pid NAME			Patient ID (default NA)
	-n/--newpun			Punish if alternative is supported with low allele frequency and the alternative allel is present in the control (cutoff freqContr/freqTumor, default = 5, set to 0 to turn off)
	-s/--addHead		tring with additional header line infer multiple times for multiple additional lines
	-f/--somout			File with the high confidence somatic SNVs
	-a/--round			Round of iteration for filtering (default = 0)
	-r/--runPCR			Run the PCR bias filter (should always be 1)
	-e/--runSEQ			Run the sequencing bias filer (should be 0 for exomes)
	-l/--runlowmaf		Set to 1 if you want to run the low maf punishment
	-x/--runExome		Run on exome, will turn off the high control coverage punishment and the PCR bias filter
	-h/--help			help\n";
}
if($file =~ /vcf\.gz$/){
	open(FH, "zcat $file |") or die "Could not open $file: $!\n";
}else{
	open (FH, $file) or die "Could not open $file: $!\n";
}
if(defined $panCanOut){
	if($panCanOut =~ /\.gz$/){
		open(PANOUT, "| bgzip > $panCanOut") or die "The outfile for panCan was set: $panCanOut but could not be opened\n";
	}else{
		open(PANOUT, ">$panCanOut") or die "The outfile for panCan was set: $panCanOut but could not be opened\n";
	}
}

my $additionalHeader;
if(defined $additionalHeader[0]){
	$additionalHeader = join("\n", @additionalHeader);
	$additionalHeader .= "\n";
}else{
	$additionalHeader = "";
}

my $date = strftime "%Y%m%d", localtime;
my @refgenome = split(",", $refgenome);

my $pancanhead = "##fileformat=VCFv4.1
##fileDate=$date
##pancancerversion=1.0
##reference=<ID=$refgenome[0],Source=$refgenome[1]>;
##center=\"$center\"
##workflowName=DKFZ_SNV_workflow
##workflowVersion=1.0.0";
$pancanhead .= $additionalHeader;
$pancanhead .= "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Indicates if record is a somatic mutation\">
##INFO=<ID=GERMLINE,Number=0,Type=Flag,Description=\"Indicates if record is a germline mutation\">
##INFO=<ID=UNCLEAR,Number=0,Type=Flag,Description=\"Indicates if the somatic status of a mutation is unclear\">
##INFO=<ID=VT,Number=1,Type=String,Description=\"Variant type, can be SNP, INS or DEL\">
##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency in primary data, for each ALT allele, in the same order as listed\">
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership\">
##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"RMS Mapping Quality\">
##INFO=<ID=1000G,Number=0,Type=Flag,Description=\"Indicates membership in 1000Genomes\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth at this position in the sample\">
##FORMAT=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-forward, ref-reverse, alt-forward and alt-reverse bases\">
##FILTER=<ID=RE,Description=\"variant in UCSC_27Sept2013_RepeatMasker.bed.gz region and/or SimpleTandemRepeats_chr.bed.gz region, downloaded from UCSC genome browser and/or variant in segmental duplication region, annotated by annovar\">
##FILTER=<ID=BL,Description=\"variant in DAC-Blacklist from ENCODE or in DUKE_EXCLUDED list, both downloaded from UCSC genome browser\">
##FILTER=<ID=DP,Description=\"<= 5 reads total at position in tumor\">
##FILTER=<ID=SB,Description=\"Strand bias of reads with mutant allele = zero reads on one strand\">
##FILTER=<ID=TAC,Description=\"less than 6 reads in Tumor at position\">
##FILTER=<ID=dbSNP,Description=\"variant in dbSNP135\">
##FILTER=<ID=DB,Description=\"variant in 1000Genomes, ALL.wgs.phase1_integrated_calls.20101123.snps_chr.vcf.gz or dbSNP\">
##FILTER=<ID=HSDEPTH,Description=\"variant in HiSeqDepthTop10Pct_chr.bed.gz region, downloaded from UCSC genome browser\">
##FILTER=<ID=MAP,Description=\"variant overlaps a region from wgEncodeCrgMapabilityAlign100mer.bedGraph.gz:::--breakPointMode --aEndOffset=1 with a value below 0.5, punishment increases with a decreasing mapability\">
##FILTER=<ID=SBAF,Description=\"Strand bias of reads with mutant allele = zero reads on one strand and variant allele frequency below 0.1\">
##FILTER=<ID=FRQ,Description=\"variant allele frequency below 0.05\">
##FILTER=<ID=TAR,Description=\"Only one alternative read in Tumor at position\">
##FILTER=<ID=UNCLEAR,Description=\"Classification is unclear\">
##FILTER=<ID=DPHIGH,Description=\"Too many reads mapped in control at this region\">
##FILTER=<ID=DPLOWC,Description=\"Only 5 or less reads in control\">
##FILTER=<ID=1PS,Description=\"Only two alternative reads, one on each strand\">
##FILTER=<ID=ALTC,Description=\"Alternative reads in control\">
##FILTER=<ID=ALTCFR,Description=\"Alternative reads in control and tumor allele frequency below 0.3\">
##FILTER=<ID=FRC,Description=\"Variant allele frequency below 0.3 in germline call\">
##FILTER=<ID=YALT,Description=\"Variant on Y chromosome with low allele frequency\">
##FILTER=<ID=VAF,Description=\"Variant allele frequency in tumor < $newpun times allele frequency in control\">
##FILTER=<ID=BI,Description=\"Bias towards a PCR strand or sequencing strand\">
##SAMPLE=<ID=CONTROL,SampleName=control_$pid,Individual=$pid,Description=\"Control\">
##SAMPLE=<ID=TUMOR,SampleName=tumor_$pid,Individual=$pid,Description=\"Tumor\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCONTROL\tTUMOR\n";

if ($makehead == 1 && defined $panCanOut)
{
print PANOUT $pancanhead;
}

my $header = "";
while ($header = <FH>)
{
	last if ($header =~ /^\#CHROM/); # that is the line with the column names
	print "$header";
}
chomp $header;

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
my $SEQERR;
my $SEQERR1;
my $PCRERR1;
my $PCRERR;
my $JUSTIF = 0;
my ($CONFIDENCE, $RECLASSIFICATION) = (0, 0);
my @help = (split "\t", $header);

# find out if there were already columns called CONFIDENCE and RECLASSIFICATION
my ($cex, $rcex, $rjust);
foreach (@help)
{
	if ($_ eq 'CONFIDENCE')
	{
		$cex = 1;
	}
	if ($_ eq 'RECLASSIFICATION')
	{
		$rcex = 1;
	}
	if ($_ eq 'PENALTIES')
	{
		$rjust = 1;
	}
}

# if there were none, append
if (! $cex)
{
	push(@help, 'CONFIDENCE');
}
if (! $rcex)
{
	push(@help, 'RECLASSIFICATION');
}
if (! $rjust && $print_info)
{
	push(@help, 'PENALTIES');
}

for (my $c = 0; $c < @help; $c++)
{
	# fixed column labels (from vcf_pileup_compare)
	if ($help[$c] eq "CONFIDENCE")
	{
		$CONFIDENCE = $c;
		print STDERR "CONFIDENCE in column $c\n";
	}
	if ($help[$c] eq "RECLASSIFICATION")
	{
		$RECLASSIFICATION = $c;
		print STDERR "RECLASSIFICATION in column $c\n";
	}
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
	if ($help[$c] =~ /MAPABILITY/)
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
	if ($help[$c] eq "SEGDUP")
	{
		$SEGDP = $c;
		print STDERR "ANNOVAR_SEGDUP_COL in column $c\n";
	}
	if ($help[$c] eq "1K_GENOMES")
	{
		$IN1KG = $c;
		print STDERR "KGENOMES_COL in column $c\n";
	}
	if ($help[$c] eq "DBSNP")
	{
		$DBSBP = $c;
		print STDERR "DBSNP_COL in column $c\n";
	}
	if ($help[$c] eq "seqBiasPresent")
	{
		$PCRERR = $c;
		print STDERR "PCR bias in column $c\n";
	}
	if ($help[$c] eq "seqingBiasPresent")
	{
		$SEQERR = $c;
		print STDERR "Sequencing bias in column $c\n";
	}
	if ($help[$c] eq "PENALTIES")
	{
		$JUSTIF = $c;
		print STDERR "Penalty justification in column $c\n";
	}
	if ($help[$c] eq "seqBiasPresent_1")
	{
		$PCRERR1 = $c;
		print STDERR "PCR bias in column $c\n";
	}
	if ($help[$c] eq "seqingBiasPresent_1")
	{
		$SEQERR1 = $c;
		print STDERR "Sequencing bias in column $c\n";
	}
}

if(defined $PCRERR && defined $SEQERR && $round == 2 && (!defined $PCRERR1 || !defined $SEQERR1))
{
	die "There was no column with the first bias filereing present\n";
}

if(defined $PCRERR && defined $SEQERR && $round == 1)
{
	$help[$PCRERR] .= "_1";
	$help[$SEQERR] .= "_1";
}
elsif(defined $PCRERR && defined $SEQERR && $round == 2)
{
	$help[$PCRERR] .= "_2";
	$help[$SEQERR] .= "_2";
}

### Check if running on an exome and set the parameters accordingly (turn of PCR bias filter)
if($runExome == 1){
	$punPCR = 0;
}

### Print header
print join("\t", @help);

print "\n";

if(defined $somout)
{
	open(SOM, ">$somout") or die "Could not open the outfile for somatic SNVs $somout\n";
	print SOM join("\t", @help), "\n";
}

####################-END-HEADER-######################

my $line = "";	# current line

# number of reference and variant reads on forward and reverse complement strand for tumor
my ($trf, $trr, $tvf, $tvr) = (0, 0, 0, 0);
# same for control
my ($crf, $crr, $cvf, $cvr) = (0, 0, 0, 0);

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
my $is_STR = 0; # also short tandem repeats from Tandem Repeats Finder are suspicious and can conincide with other classes
my $is_weird = 0;	# coindicende with known artefact-producing regions
my $in1KG = 0;
my $indbSNP = 0;
my $precious = 0;
my $class = ""; # for potential re-classification (e.g. low coverage in control and in dbSNP => probably germline)

while (<FH>)
{
	$confidence=10; # start with maximum value
	my $reasons = "";	# collect info on which penalties came into effect
	my $dbsnpPos;
	my $dbsnpId;
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
	$class = $help[$CLASS]; # start with original classification

### For pancancer
    #	next if($help[4] =~ /,/);	# For now throw out all multiple alternative bases
	my ($gttum) = $help[9] =~ m/(^\d\/\d)/; # genotype tumor as originally from mpileup
	my ($gtcontr) = $help[$CINFO] =~ /VAF\=(\d+\.*\d*)/;	# VAF in control will be genotype 1/1 if > 0.7 else will be 0/1
	my %infofield;
	my %filterfield;
	my %DP4tumor;
	my @DP4control;

	$infofield{"VT"} = "SNP";	# All variants here are SNVs, no indels in this file

# 1) external information of if these SNPs have already been found (incl. false positives from 1000 genomes!)
	# 1000 genomes
	#if ($help[$IN1KG] ne "." && $help[$IN1KG] ne "0")
	if ($help[$IN1KG] =~ /MATCH=exact/)	# verified somatic SNVs have MATCH=position instead of exact, so it would be bad to reclassify them as SNP_support_germline!
	{
		$in1KG = 1;	# good for germline, bad for "somatic"
		$infofield{"1000G"} = "1000G";
	}
	# dbSNP
	#if ($help[$DBSBP] ne "." && $help[$DBSBP] ne "0")	# verified somatic SNVs have MATCH=position instead of exact, so it would be bad to reclassify them as SNP_support_germline!
	if ($help[$DBSBP] =~ /MATCH=exact/)
	{
		$indbSNP = 1;	# good for germline, bad for "somatic"
		$infofield{"DB"} = "DB";
		$dbsnpPos = $help[1];
		($dbsnpId) = $help[$DBSBP] =~ /ID\=(rs\d+)/;
		# precious!
		#INFO=<ID=CLN,Number=0,Type=Flag,Description="SNP is Clinical(LSDB,OMIM,TPA,Diagnostic)">
		#INFO=<ID=PM,Number=0,Type=Flag,Description="SNP is Precious(Clinical,Pubmed Cited)">
		if ($help[$DBSBP] =~ /;CLN;/ || $help[$DBSBP] =~ /;PM;/)	# more likely somatic
		{
			$precious = 1;
		}
	}
# Punish for biases round 1
	if(defined $PCRERR && defined $SEQERR && $round == 1)
	{
		if($help[$PCRERR] !~ /^[\.0]$/ && $help[$SEQERR] !~ /^[\.0]$/ && $punPCR == 1 && $punSEQ == 1){
			$confidence-=3;
			$reasons.="bias_filter_round1_PCR_and_SEQ(-3)";
			$filterfield{"BI"} = 1;
		}elsif($help[$PCRERR] !~ /^[.0]$/ && $punPCR == 1)
		{
			$confidence-=3;
			$reasons.="bias_filter_round1_PCR(-3)";
			$filterfield{"BI"} = 1;
		}elsif($help[$SEQERR] !~ /^[.0]$/ && $punSEQ == 1)
		{
			$confidence-=3;
			$reasons.="bias_filter_round1_SEQ(-3)";
			$filterfield{"BI"} = 1;
		}
	}
# Punish for biases round 1 in round 2
	if(defined $PCRERR1 && defined $SEQERR1 && $round == 2)
	{
		if($help[$PCRERR1] !~ /^[\.0]$/ && $help[$SEQERR1] !~ /^[\.0]$/ && $punPCR == 1 && $punSEQ == 1){
			$confidence-=3;
			$reasons.="bias_filter_round1_PCR_and_SEQ(-3)";
			$filterfield{"BI"} = 1;
		}elsif($help[$PCRERR1] !~ /^[.0]$/ && $punPCR == 1)
		{
			$confidence-=3;
			$reasons.="bias_filter_round1_PCR(-3)";
			$filterfield{"BI"} = 1;
		}elsif($help[$SEQERR1] !~ /^[.0]$/ && $punSEQ == 1)
		{
			$confidence-=3;
			$reasons.="bias_filter_round1_SEQ(-3)";
			$filterfield{"BI"} = 1;
		}
	}
# Punish for biases round 2
	if(defined $PCRERR && defined $SEQERR && $round == 2)
	{
		if($help[$PCRERR] !~ /^[\.0]$/ && $help[$SEQERR] !~ /^[\.0]$/ && $punPCR == 1 && $punSEQ == 1){
			$confidence-=3;
			$reasons.="bias_filter_round2_PCR_and_SEQ(-3)";
			$filterfield{"BI"} = 1;
		}elsif($help[$PCRERR] !~ /^[.0]$/ && $punPCR == 1)
		{
			$confidence-=3;
			$reasons.="bias_filter_round2_PCR(-3)";
			$filterfield{"BI"} = 1;
		}elsif($help[$SEQERR] !~ /^[.0]$/ && $punSEQ == 1)
		{
			$confidence-=3;
			$reasons.="bias_filter_round2_SEQ(-3)";
			$filterfield{"BI"} = 1;
		}
	}
# 2) annotations of regions that cause problems: some classes of repeats from RepeatMasker track, segmental duplications,
# (cf. Reumers et al. 2012, Nature Biotech 30:61), external blacklists, mapability
	# simple repeats and low complexity (not the same as homopolymer, but similar enough); some satellites are not annotated in blacklist ...
	if ($help[$REPET] =~ /Simple_repeat/ || $help[$REPET] =~ /Low_/ || $help[$REPET] =~ /Satellite/)
	{
		$is_repeat = 1;
		$confidence-=2;
		$reasons.="Simple_repeat(-2)";
		$filterfield{"RE"} = 1;
	}
	# other repeat elements to penalize at least a bit
	elsif ($help[$REPET] ne "." && $help[$REPET] ne "0")
	{
		$confidence--;
		$reasons.="Other_repeat(-1)";
		$filterfield{"RE"} = 1;
	}
	# simple tandem repeats most often coincide with other bad features - do not penalize twice
	if ($help[$STREP] ne "." && $help[$STREP] ne "0")
	{
		$is_STR = 1;
		unless ($is_repeat)
		{
			$confidence-=2;
			$reasons.="Tandem_repeat(-2)";
			$filterfield{"RE"} = 1;
		}
	}
	# Segmental Duplications are less effective than homopolymers, short tandem repeats and microsatellites, do not penality twice
	if ($help[$SEGDP] ne "." && $help[$SEGDP] ne "0" && ! $is_repeat && ! $is_STR)
	{
		$confidence-=2; # bad region
		$is_weird = 1;
		$reasons.="Segmental_dup(-2)";
		$filterfield{"RE"} = 1;
	}
	# Duke excluded and ENCODE DAC blacklist, only consider if not already annotated as suspicious repeat
	if (($help[$EXCLU] ne "." && $help[$EXCLU] ne "0") || ($help[$BLACK] ne "." && $help[$BLACK] ne "0"))
	{
		$confidence-=3; # really bad region, usually centromeric repeats
		$is_weird = 1;
		$reasons.="Blacklist(-3)";
		$filterfield{"BL"} = 1;
	}
	# HiSeqDepth: regions "attracting" reads; often coincide with tandem repeats and CEN/TEL, not always with low mapability
	if ($help[$HSDEP] ne "." && $help[$HSDEP] ne "0")
	{
		$confidence-=3; # really really bad region!
		$is_weird = 1;
		$reasons.="Hiseqdepth(-3)";
		$filterfield{"HSDEPTH"} = 1;
	}
	# Mapability is 1 for unique regions, 0.5 for regions appearing twice, 0.33... 3times, ...
	# Everything with really high number of occurences is artefacts
	# does not always correlate with the above regions
	# is overestimating badness bc. of _single_ end read simulations

	$mapp = $help[$MAPAB];
	if ($mapp eq ".")	# in very rare cases (CEN), there is no mapability => ".", which is not numeric but interpreted as 0
	{
		$confidence-=5;
		$reasons.="Not_mappable(-5)";
		$filterfield{"MAP"} = 1;
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
		my $reduce = 0;
		if ($mapp < 0.5)	# 0.5 does not seem to be that bad: region appears another time in the genome and we have paired end data!
		{
			$confidence--;
			$reduce++;
			$reasons.="Low_mappability($mapp=>";
			$is_weird = 1;	# something _is_ weird already there and known SNPs might be artefacts
			if ($mapp < 0.4)	#  3-4 times appearing region is worse but still not too bad
			{
				$confidence--;
				$reduce++;
				if ($mapp < 0.25)	# > 4 times appearing region
				{
					$confidence--;
					$reduce++;
				}
				if ($mapp < 0.1)	# > 5 times is bad
				{
					$confidence-=2;
					$reduce+=2;
					if ($mapp < 0.05)	# these regions are clearly very bad (Lego stacks)
					{
						$confidence-=3;
						$reduce+=3;
					}
				}
			}
			$reasons.="-$reduce)";
			$filterfield{"MAP"} = 1;
		}
	}
	# if others have found the SNP already, it may be interesting despite low score - but only if it's not a weird region.
	# if a position gets up from score 4 to 6 by giving +2 for presence in dbSNP, it's very likely an artefact reported to dbSNP
	#~ if (($in1KG || $indbSNP) && ! $is_weird)
	#~ {
		#~ $confidence++;
		#~ $reasons.="SNP_known_noweirdregion(+1)";
	#~ }

# 3) information from the calls and germline comparisons: coverage, strand bias, variant support, ... => can lead to reclassification
# what we do not know: "custering", misplaced indels, in read end region (!overlapping reads) (cf. Reumers et al.)
	# parse INFO fields
	($trf, $trr, $tvf, $tvr) = $help[$INFO] =~/DP4=(\d+),(\d+),(\d+),(\d+);/;
	my ($DP4tumor) = $help[$INFO] =~/DP4=(\d+,\d+,\d+,\d+);/;
	my ($MAPqual) = $help[$INFO] =~ /MQ=(\d+);/;
	$infofield{"MQ"} = "MQ=$MAPqual";
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
	my ($DP4control) = $help[$CINFO] =~/DP5=(\d+,\d+,\d+,\d+),\d+;/;
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
	### DP4 based allele frequency for tumor and normal for AF info fields
	my $alfrtum = sprintf('%.2f', $fr_var_tum);
	my $altfctrl = sprintf('%.2f', $fr_var_ctrl);
	$infofield{"AF"} = "AF=$altfctrl,$alfrtum";

    # Punish if alternative allele frequency is low and variant is also seen in control (try first with factor of tumor/control)
    if ($newpun != 0 && $altfctrl > 0 && $class eq "somatic")
    {
	if (($alfrtum/$altfctrl) < $newpun)
	{
	    $confidence-=3;
	    $reasons.="Alternative_allel_freq_not_".$newpun."_times_bigger_in_tum_than_in_control(-3)";
	    $filterfield{"VAF"} = 1;
	}
    }

	# regions have an excessive number of reads in control, may not be recognized by the mapability and High Seq Depth tracks
	# (unclear if these tracks are from "normal" hg19 or 1KG reference with decoy sequences!)
	# tumor could have deletion here
	if ($depthC > $cutoff && $runExome != 1)
	{
		$confidence-=2;
		if ($depthC > $cutoff*2)	# >>300 reads at lego stack regions
		{
			$confidence-=2;
			$reasons.="Controlcoverage>2*$cutoff(-4)";
			$filterfield{"DPHIGH"} = 1;
		}
		else
		{
			$reasons.="Controlcoverage>$cutoff(-2)";
			$filterfield{"DPHIGH"} = 1;
		}
	}
	# total coverage too low
	#if ($depthT <= 5 || $depthC <= 5)	# low coverage is bad in general
	#{
	#	$confidence-=3; # "low", very bad!
	#}
	# split up in 2, to subtract 3 for tumor low coverage bc. that's so bad and additional penalty for low control depth
	if ($depthT <= 5)
	{
		$confidence-=3;
		$reasons.="Tumor<5reads(-3)";
		$filterfield{"DP"} = 1;
	}
	if ($depthC <= 5)	# not _that_ bad as depthT low - will very likely become reclassified as "lowCov_SNP_support_germline"
	{
		$confidence-=3;
		$reasons.="Controlcoverage<=5(-3)";
		$filterfield{"DPLOWC"} = 1;
	}
	if ($tvf < 1 || $tvr < 1)	# strand bias
	{
		$confidence-=2; # "moderate"
		if ($fr_var_tum <= 0.1) # and too low variant coverage
		{
			if($runlowmaf != 1){
				$confidence-=3; # "low"
				$reasons.="Strandbias_and_MAF<=0.1(-3)";
			}else{
				$confidence-=1; # "low"
				$reasons.="Strandbias_and_MAF<=0.1_and_low_maf_filter(-1)";
			}
			$filterfield{"SBAF"} = 1;
		}
		else
		{
			$reasons.="Strandbias(-2)";
			$filterfield{"SB"} = 1;
		}
	}
	if ($tvf + $tvr <= 5)	# low variant support
	{
		if($runlowmaf != 1){
			$confidence-=2; # "moderate"
			$reasons.="ALT<=5(-2)";
		}else{
			$confidence-=0; # "moderate"
			$reasons.="ALT<=5_and_low_maf_filter(-0)";
		}
		$filterfield{"TAC"} = 1;
	}
	if ($depthT >= 150 && ($tvf + $tvr <=2) && $runlowmaf == 1){
		$confidence -= 3;
		$reasons.="Coverage_above_200_and_less_than_3_variant_reads(-3)";
	}
	if ($tvf == 1 && $tvr == 1)   # only 1 read from each strand (might be PCR error of overlapping reads!) ### TODO: But overlaps are removed??? - Still unreliable!
	{
		if($runlowmaf != 1){
			$confidence-=3; # "low"
			$reasons.="Only_one_ALT_per_strand(-3)";
		}else{
			$confidence-=1; # "low"
			$reasons.="Only_one_ALT_per_strand_and_low_maf_filter(-1)";
		}
		$filterfield{"1PS"} = 1;
	}
	if ($tvf + $tvr < 5 && $fr_var_tum < 0.05 && $runlowmaf != 1)  # too low variant coverage
	{
		$confidence-=3; # "low"
		$reasons.="ALT<6_and_MAF<=0.1(-3)";
		$filterfield{"FRQ"} = 1;
	}
	elsif ($tvf + $tvr < 5 && $fr_var_tum < 0.05 && $runlowmaf == 1)  # too low variant coverage
	{
		$confidence-=1; # "low"
		$reasons.="ALT<6_and_MAF<=0.1_and_low_maf_filter(-1)";
		$filterfield{"FRQ"} = 1;
	}
	if ($help[$CLASS] =~/somatic/)	# this includes multi_somatic
	{
		# the overlapping reads filter may have reduced the variant read to 1 - before, such SNVs were never considered
		if ($tvf + $tvr < 2)
		{
			$confidence-=3; # very low!
			$reasons.="Somatic_ALT<2(-3)";
			$filterfield{"TAR"} = 1;
		}
		if ($fr_var_ctrl > 0)	# might be germline or artefact! 1/30 reads required for calling germline
		{
			$class = "LQVSIG";
			if ($controlbad)
			{
				$confidence-=2; # "moderate"
				$reasons.="Somatic_Control_${fr_var_ctrl}_ALT(-2)";
				$filterfield{"ALTC"} = 1;
			}
			if ($fr_var_tum < 0.3 && $controlbad)
			{
				$confidence-=3; # "low"
				$reasons.="Somatic_MAF<03_and_Control_${fr_var_ctrl}_ALT(-3)";
				$filterfield{"ALTCFR"} = 1;
			}
			#if ($tvf < 1 || $tvr < 1)	# already penalized in general above with -2
			#{
			#	$confidence-=3; # "low"
			#}
		}
		if ($depthC && $varsuppC / $depthC > 1 / 30)	# present in germline but with low qual, so likely seq error ### WAS longer before: $varsuppC > 2. THIS WAS TOO STRINGENT FOR HIGH COVERAGE DATA ### has already have been punished above!
		{
			unless ($class eq "LQVSIG")
			{
				$class = "LQVSIG";
				if ($controlbad)
				{
					$confidence-=3;
					$reasons.="Somatic_Control>1/30_ALT(-3)";
					$filterfield{"ALTC"} = 1;
				}
			}
			# TODO: subtract only 2 because Rosario's exome recurrent example with LQVSIG score 7 is true! (contamination!)
		}
		# an SNV that is in dbSNP but not "precious" or/and in 1 KG with high frequency is probably germline,
		# especially if control coverage is low (may be due to library problems)
		if ($in1KG || ($indbSNP && ! $precious))
		{
			$class = "SNP_support_germline";	# missed germline?
			if ($depthC < 10)	# Matthias' suggestion that below 10 (and not <= 5 as before) the probability to miss the variant is >95%
			{
				$class = "lowCov_SNP_support_germline";
				$filterfield{"DB"} = 1;
			}
		}
	}
	# in the classification, 1/30 reads required for calling germline
	if ($help[$CLASS] =~ /germline/)
	{
		if ($fr_var_ctrl < 0.3) # too far from 50%, but that depends on coverage. I'm not convinced that qualifies for "low" / -3
		{
			if ($controlbad)
			{
				$confidence-=2; # "low"
				$reasons.="Germline_ALT<0.3(-2)";
				$filterfield{"FRC"} = 1;
			}
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
	if ($help[$CLASS] =~ /LQVSIG/)	# low quality variant support in germline: usually a weird region; a rare case anyways. BUG in original script ("LQVSIG"): it's called LQVCIG with C (for count) instead of S (support) in original classification!
	{
		if ($in1KG || ($indbSNP && ! $precious))
		{
			$class = "SNP_support_germline";
			$filterfield{"DP"} = 1;
		}
		else	# weird region
		{
			if ($controlbad)
			{
				$confidence-=2;
				$reasons.="Original_LQVSIG(-2)";
				$filterfield{"ALTC"} = 1;
			}
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
			$filterfield{"DP"} = 1;
		}
		else	# just a weird region
		{
			if ($controlbad)
			{
				$confidence-=2;
				$reasons.="Unclear_no_known_SNP(-2)";
				$filterfield{"ALTC"} = 1;
			}
		}
	}
	if ($help[0] =~ /Y/ && $help[5] < 100)	# SNV on chrY would have to be homozygous, too
	# het SNV on male chrX is another weird thing (unless a subpopulation had amplified chrX and mutated)
	# any SNV on chrY in female is an artefact
	{
		$confidence-=3; # "low"
		$reasons.="Y_mpileup_score<100(-3)";
		$filterfield{"YALT"} = 1;
	}
	#To make sure that changes in the raw filter do not influence the final result we punish them with -3
	if($runlowmaf != 1 && (($help[5] < 3 && ($alfrtum < 0.05 || ($tvf + $tvr < 5))) || ($help[5] < 20 && (($tvf == 0 || $tvr == 0) && ($tvf + $tvr <= 3)))))
	{
		$confidence-=3;
		$reasons.="raw_filter_punishment(-3)";
		$filterfield{"FRQ"} = 1;
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

	$help[$CONFIDENCE] = $confidence;
	$help[$RECLASSIFICATION] = $class;

### Make here different outputs, first the one for the pancancer, then if wanted add the others
	my @panout = @help[0 .. 4];
	$help[6] =~ s/^\.//;
	$panout[5] = ".";
	if($confidence > 7)
	{
		$panout[6] = "PASS";	# FILTER column
		$help[6] = "PASS";
	}
	else
	{
		foreach(my @filteroptions = ("RE","BL","DP","SB","TAC","dbSNP","DB","HSDEPTH","MAP","SBAF","FRQ","TAR","UNCLEAR","DPHIGH","DPLOWC","1PS","ALTC","ALTCFR","FRC","YALT","VAF","BI"))
		{
			if(exists $filterfield{$_} && $filterfield{$_} == 1){
				$panout[6] .= "$_;";
				if($help[6] !~ /^$_$/ && $help[6] !~ /^$_;/ && $help[6] !~ /;$_;/ && $help[6] !~ /;$_$/){
					$help[6] .= "$_;";
				}
			}
		}
		$panout[6] =~ s/;$//;
		$help[6] =~ s/;$//;
	}
	### Definition of genotypes still crude, could be more precise

	if($help[$CLASS] =~ /germline/ && $class ne "lowCov_SNP_support_germline" )
	{
		$infofield{"GERMLINE"} = "GERMLINE";
		if($gtcontr >= 0.75){$gtcontr = "1/1";}
		else{$gtcontr = "0/1";}
	}
	elsif($class eq "lowCov_SNP_support_germline" || $help[$CLASS] !~ /somatic/)
	{
		$infofield{"UNCLEAR"} = "UNCLEAR";
		if($gtcontr >= 0.75){$gtcontr = "1/1";}
		elsif($gtcontr >= 0.3){$gtcontr = "0/1";}
		else{$gtcontr = "0/0";}
	}
	else
	{
		$infofield{"SOMATIC"} = "SOMATIC";
		$gtcontr = "0/0";
#		if($gttum eq "0/0"){$gttum = "0/1";} # To change that the genotype ot the alternative is at least 0/1
	}

	foreach(my @infooptions = ("SOMATIC","GERMLINE","UNCLEAR","VT","AF","MQ","DB","1000G"))
	{
		if(exists $infofield{$_}){$panout[7] .= "$infofield{$_};";}
	}
	if(defined $dbsnpId && defined $dbsnpPos)
	{
		$panout[2] = $dbsnpId."_".$dbsnpPos;
		$help[2] = $dbsnpId;
	}
	$panout[7] =~ s/;$//;
	$panout[8] = "GT:DP:DP4";
	$panout[9] = "$gtcontr:$depthC:$DP4control";
	$panout[10] = "$gttum:$depthT:$DP4tumor";


### Print
	if(defined $panCanOut){
		print PANOUT join("\t", @panout), "\n";
	}
	if ($print_info)
	{
		if(!defined $reasons || $reasons eq "")
		{
			$reasons = ".";
		}
		$help[$JUSTIF] = "$reasons";
	}
	print join ("\t", @help);
	print "\n";

	if(defined $somout && $confidence > 7 && $help[$CLASS] eq "somatic" && $help[$RECLASSIFICATION] !~ /lowCov_SNP_support_germline/){print SOM join("\t", @help), "\n";}
}
close FH;
exit;
