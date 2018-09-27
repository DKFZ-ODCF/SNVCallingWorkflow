#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

use strict;
use warnings;

my $infile = $ARGV[0];

if($infile =~ /\.gz/){open(IN, "zcat $infile |") or die "Could not open the $infile to detect the confidence column\n";}
else{open(IN, "<$infile") or die "Could not open the $infile to detect the confidence column\n";}

my $line;
while(<IN>)
{
	chomp;
	$line=$_;
	last if($_ =~ /^#CHROM\s/);
}
close IN;
my $i = 0;
my @line = split("\t", $line);
while($i <= @line)
{
	last if($line[$i] =~ /^CONFIDENCE$/);
	$i++;
}
print $i;