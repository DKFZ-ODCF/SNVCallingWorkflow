#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

# count how many / what percentage of somatic SNVs of certain confidence score (default >= 8) are in dbSNP and 1KG

use strict;
use warnings;

my $usage = "USAGE: 1.snv file (can be the pre-filtered one from rainfall plots) - optional 2.min score (default 8)\n";
if (@ARGV < 1)
{
	die $usage;
}

my $file = shift;
my $minscore = shift;

unless (defined $minscore)
{
	$minscore = 8;
}
open (my $FH, $file) or die "Could not open $file: $!\n";

# count all, somatic, s. with minconfidence, in dbSNP, in 1KG
my $all = 0;
my $somatic = 0;
my $scored = 0;
my $dbSNP = 0;
my $oneKG = 0;
my $perc_dbSNP = 0;
my $perc_1KG = 0;
my @help = ();


my $header = "";
my ($DBSBP, $KG, $CONFIDENCE, $CLASSIFICATION);

while(!eof($FH)){
	my $line = readline($FH) || die "Could not read line: $?";
	chomp $line;
	if ($line =~ /^#/)
	{
		#print STDOUT "HELLO\n";
		if ($line =~/^#CHROM/){
			chomp $line;

			my @head=split("\t", $line);
			my %col;

			my $i = 0;
			while($i < @head)
			{
				if($head[$i] eq "DBSNP"){
					$col{"DBSNP"} = $i;
					print STDERR "DBSNP in column ", $i+1,"\n";
				}
				if($head[$i] eq "1K_GENOMES"){
					$col{"1K_GENOMES"} = $i;
					print STDERR "1K_GENOMES in column ", $i+1,"\n";
				}
				if($head[$i] eq "RECLASSIFICATION"){
					$col{"RECLASSIFICATION"} = $i;
					print STDERR "RECLASSIFICATION in column ", $i+1,"\n";
				}
				if($head[$i] eq "ANNOTATION_control"){
					$col{"ANNOTATION_control"} = $i;
					print STDERR "ANNOTATION_control in column ", $i+1,"\n";
				}
				if($head[$i] eq "CONFIDENCE"){
					$col{"CONFIDENCE"} = $i;
					print STDERR "CONFIDENCE in column ", $i+1,"\n";
				}
				$i++;
			}

			$DBSBP = $col{"DBSNP"};
			$KG = $col{"1K_GENOMES"};
			$CONFIDENCE = $col{"CONFIDENCE"};
			$CLASSIFICATION = $col{"RECLASSIFICATION"};
		}
	}
	else{
		$all++;
		@help = split ("\t", $line);
		if ($help[$CLASSIFICATION] =~ /somatic/)
		{
			$somatic++;
			if ($help[$CONFIDENCE] >= $minscore)
			{
				$scored++;
				if ($help[$DBSBP] =~ /MATCH=exact/)
				{
					$dbSNP++;
				}
				if ($help[$KG] =~ /MATCH=exact/)
				{
					$oneKG++;
				}
			}
		}
	}
}
close $FH;

if ($scored > 0)
{
	$perc_dbSNP = sprintf ("%.2f", ($dbSNP/$scored*100));
	$perc_1KG = sprintf ("%.2f", ($oneKG/$scored*100));
	print "file\tall\tSomatic\tS_with_minscore$minscore\tS_in_dbSNP\tS_in_1KG\tperc_dbSNP\tperc_1KG\n";
	print "$file\t$all\t$somatic\t$scored\t$dbSNP\t$oneKG\t$perc_dbSNP\t$perc_1KG\n";
}
else	# there might be no ...
{
	die "!!! There are no SNVs with confidence >= $minscore in the input, it makes no sense to do further analyses with the current settings on these data!!!\n";
}
exit;
