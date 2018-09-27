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
open (FH, $file) or die "Could not open $file: $!\n";

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
while ($header = <FH>)
{
	last if ($header =~ /^\#CHROM/); # that is the line with the column names
}
chomp $header;

my $DBSBP = 13;
my $KG = 14;
my $CONFIDENCE = 28;
my $CANNOTATION = -1;
my $RECLASSIFICATION = -1;
my $CLASSIFICATION;

@help = split(/\t/, $header);
for (my $c = 0; $c < @help; $c++)
{
	# fixed column labels (from vcf_pileup_compare)
	if ($help[$c] eq "CONFIDENCE")
	{
		$CONFIDENCE = $c;
		print STDERR "CONFIDENCE in column $c\n";
	}
	if ($help[$c] eq "1K_GENOMES")
	{
		$KG = $c;
		print STDERR "INFO in column $c\n";
	}
	if ($help[$c] eq "DBSNP")
	{
		$DBSBP = $c;
		print STDERR "DBSNP_COL in column $c\n";
	}
	if ($help[$c] eq "ANNOTATION_control")
    {
    	$CANNOTATION = $c;
    	print STDERR "ANNOTATION_control in column $c\n";
    }
    if ($help[$c] eq "RECLASSIFICATION")
    {
        $RECLASSIFICATION = $c;
        print STDERR "RECLASSIFICATION in column $c\n";
    }
}

if ($CANNOTATION > 0) {
    $CLASSIFICATION = $CANNOTATION;
} else {
    $CLASSIFICATION = $RECLASSIFICATION;
}

# I know that it's evilly hardcoded!
while (<FH>)
{
	if ($_ =~ /^#/)	# skip header lines
	{
		next;
	}
	$all++;
	@help = split ("\t", $_);
	if ($help[$CLASSIFICATION] =~ /somatic/)
	{
		$somatic++;
		if ($help[$CONFIDENCE] >= $minscore)
		{
			$scored++;
			#print STDERR $_;
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

close FH;
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
