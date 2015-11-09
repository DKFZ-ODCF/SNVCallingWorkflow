#!/usr/bin/env perl

# tumor variant frequency for high confidence somatic SNVs as input for MAF plots
# output: PID (number!)	chrom	tumor variant frequency

use strict;
use warnings;

if (@ARGV < 2)
{
	die "1. vcf file to calculate tumor variant frequency for MAF plot - 2.minimal coverage cutoff (default 0) 3.minimal confidence thereshold (default 8)\n";
}
my $file = shift;
open (FH, $file) or die "Could not open $file: $!\n";

################################################################################
my $header = "";
while ($header = <FH>)
{
	last if ($header =~ /^\#CHROM/); # that is the line with the column names
}
chomp $header;

my $CONFIDENCE = 28;
my $INFO = 7;
my $DBSBP = 0;

my @help = split(/\t/, $header);
for (my $c = 0; $c < @help; $c++)
{
	# fixed column labels (from vcf_pileup_compare)
	if ($help[$c] eq "CONFIDENCE")
	{
		$CONFIDENCE = $c;
		print STDERR "CONFIDENCE in column $c\n";
	}
	if ($help[$c] eq "INFO")
	{
		$INFO = $c;
		print STDERR "INFO in column $c\n";
	}
	if ($help[$c] eq "DBSNP")
	{
		$DBSBP = $c;
		print STDERR "DBSNP_COL in column $c\n";
	}
}

###############################################################################


my $mincov = shift;
unless (defined $mincov)
{
	$mincov = 0;
}
my $minconfidence = shift;
unless (defined $minconfidence)
{
	$minconfidence = 8;
}

#INFO => DP4 of tumor => depth and variant frequency
my ($trf, $trr, $tvf, $tvr) = (0, 0, 0, 0);
@help = ();
my $depth = 0;
# tumor variant frequency from DP4 fields (pos. 7 of vcf file)
my $tumvarfrq = 0;

while (<FH>)
{
	if ($_ =~ /^#/)
	{
		next;
	}
	chomp;
	@help = split ("\t", $_);
	if ($help[$CONFIDENCE] >= $minconfidence)
	{
		($trf, $trr, $tvf, $tvr) = $help[$INFO] =~/DP4=(\d+),(\d+),(\d+),(\d+);/;
		$depth = $tvf+$tvr+$trf+$trr;
		if ($depth >= $mincov)
		{
			$tumvarfrq = sprintf ("%.2f", (($tvf+$tvr)/($tvf+$tvr+$trf+$trr)));
			if($help[$DBSBP] =~ /MATCH=exact/)
			{
				print "1\t$help[0]\t$tumvarfrq\t$depth\n";
			}
			else
			{
				print "0\t$help[0]\t$tumvarfrq\t$depth\n";
			}
		}
	}
}
close FH;
exit;
