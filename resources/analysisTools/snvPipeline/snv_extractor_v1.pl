#!usr/bin/perl

# By Ivo Buchhalter
# Hard coded script to extract snvs in the snv pipeline
# Extracts: somatic snvs -> snvs_PID_somtatic_snvs_conf_?_to_10.vcf
#           somatic coding snvs -> snvs_PID_somtatic_coding_snvs_conf_?_to_10.vcf
#           germline coding snvs -> snvs_PID_germline_coding_snvs_conf_?_to_10.vcf

use strict;
use warnings;
use Getopt::Long;

my $infile;
my $minconf=8;
my $region = 0;
my $pid;
my $extractsyn = 0;
my $extractNcRNA = 0;
my $bgzip = "bgzip";
my $tabix = "tabix";
GetOptions (	"infile=s"	 		=> \$infile,		# vcf file, can be bgzipped
				"minconf=i"			=> \$minconf,		# minimum confidence score
				"region=s"			=> \$region,		# file with the region when only a region should be extracted
				"ncRNA=i"			=> \$extractNcRNA,	# If non coding RNAs (functional) should be extracted
				"pid=s"				=> \$pid,
				"synonymous=i"		=> \$extractsyn,
				"bgzip=s"			=> \$bgzip,
				"tabix=s"			=> \$tabix
) or die "Could not get the options!\n";

if($region ne "0" && (!-f $region || $infile !~ /\.gz$/ || !-f $infile.".tbi")){die "region-file: $region is not a valid file or, infile: $infile is not zipped or there is no index for the infile $infile.tbi\n";}
if(!-f $infile){die "The provided infile: $infile is not a valid file\n";}
if(-f $infile && $infile =~ /\.gz$/){open(IN, "zcat $infile |");}
else{open(IN, "<$infile") or die "Could not open the infile: $infile\n";}

if($region ne "0"){print "Region file provided: $region\n";}

my $outsom = $pid."_somatic_snvs_conf_".$minconf."_to_10.vcf";
my $outsomcod = $pid."_somatic_functional_snvs_conf_".$minconf."_to_10.vcf";
my $outgermcod = $pid."_germline_functional_snvs_conf_".$minconf."_to_10.vcf";
my $outsyn = $pid."_somatic_functional_and_synonymous_snvs_conf_".$minconf."_to_10.vcf";
my $outNcRNA = $pid."_somatic_functional_ncRNA_snvs_conf_".$minconf."_to_10.vcf";

open(SOM, ">$outsom") or die "Could not open the file $outsom\n";
open(COD, ">$outsomcod") or die "Could not open the file $outsomcod\n";
open(GER, ">$outgermcod") or die "Could not open the file $outgermcod\n";
if($extractsyn == 1){open(SYN, ">$outsyn") or die "Could not open the file $outsyn\n";}
if($extractNcRNA == 1){open(NCRNA, ">$outNcRNA") or die "Could not open the file $outNcRNA\n";}

my $head;
while(<IN>)
{
	chomp;
	$head=$_;
	last if($_ =~ /^#CHR/);
}
close IN;

my @head=split("\t", $head);
my %col;


my $i = 0;
while($i < @head)
{
	if($head[$i] eq "EXONIC_CLASSIFICATION"){$col{"EXONIC_CLASSIFICATION"} = $i; print "EXONIC_CLASSIFICATION in column ", $i+1,"\n";}
	if($head[$i] eq "ANNOVAR_FUNCTION"){$col{"ANNOVAR_FUNCTION"} = $i; print "ANNOVAR_FUNCTION in column ", $i+1,"\n";}
	if($head[$i] eq "RECLASSIFICATION"){$col{"RECLASSIFICATION"} = $i; print "RECLASSIFICATION in column ", $i+1,"\n";}
	if($head[$i] eq "ANNOTATION_control"){$col{"ANNOTATION_control"} = $i; print "ANNOTATION_control in column ", $i+1,"\n";}
	if($head[$i] eq "CONFIDENCE"){$col{"CONFIDENCE"} = $i; print "CONFIDENCE in column ", $i+1,"\n";}
	$i++;
}

print SOM $head, "\n";
print COD $head, "\n";
print GER $head, "\n";
if($extractsyn == 1){print SYN $head, "\n";}
if($extractNcRNA == 1){print NCRNA $head, "\n";}

if($region ne "0"){open(IN, "$tabix $infile -B $region |") or die "Could not open the file with tabix and regions\n";}
elsif($infile =~ /\.gz$/){open(IN, "zcat $infile |");}
else{open(IN, "<$infile");}

while(<IN>)
{
	chomp;
	next if($_ =~ /^#/);
	my @line = split("\t", $_);
	next if($line[$col{"CONFIDENCE"}] < $minconf);
	if($line[$col{"ANNOTATION_control"}] eq "somatic" && $line[$col{"RECLASSIFICATION"}] !~ /lowCov_SNP_support_germline/){print SOM $_, "\n";}
	
	if($line[$col{"ANNOTATION_control"}] eq "somatic" && $line[$col{"RECLASSIFICATION"}] !~ /lowCov_SNP_support_germline/ && $line[$col{"ANNOVAR_FUNCTION"}] !~ /ncRNA/ && ($line[$col{"EXONIC_CLASSIFICATION"}] =~ /nonsynonymous/ || $line[$col{"EXONIC_CLASSIFICATION"}] =~ /stopgain/ ||$line[$col{"EXONIC_CLASSIFICATION"}] =~ /stoploss/ || $line[$col{"ANNOVAR_FUNCTION"}] =~ /splicing/)){
		print COD $_, "\n";
	}
	if($extractNcRNA == 1 && $line[$col{"ANNOTATION_control"}] eq "somatic" && $line[$col{"RECLASSIFICATION"}] !~ /lowCov_SNP_support_germline/ && ($line[$col{"ANNOVAR_FUNCTION"}] =~ /ncRNA_exonic/ || $line[$col{"ANNOVAR_FUNCTION"}] =~ /ncRNA_splicing/)){
		print NCRNA $_, "\n";
	}
	if($extractsyn == 1 && $line[$col{"ANNOTATION_control"}] eq "somatic" && $line[$col{"RECLASSIFICATION"}] !~ /lowCov_SNP_support_germline/ && ($line[$col{"EXONIC_CLASSIFICATION"}] =~ /^synonymous/ )){
		print SYN $_, "\n";
	}
	
	if($line[$col{"ANNOTATION_control"}] eq "germline" && $line[$col{"ANNOVAR_FUNCTION"}] !~ /ncRNA/ && ($line[$col{"EXONIC_CLASSIFICATION"}] =~ /nonsynonymous/ || $line[$col{"EXONIC_CLASSIFICATION"}] =~ /stopgain/ ||$line[$col{"EXONIC_CLASSIFICATION"}]  =~ /stoploss/ || $line[$col{"ANNOVAR_FUNCTION"}] =~ /splicing/)){
		print GER $_, "\n";
	}
}

close IN;
close SOM;
close COD;
close GER;
