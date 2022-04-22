#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

use strict;
use warnings;

my @bins;
my $max = 50000;
my $median;
my $outfile = $ARGV[1];

my $infile = $ARGV[0];
open(my $IN, "<$infile");

my $count=0;
my @head;
while (!eof($IN)) {
	defined(my $line = readline($IN))
		|| die "Could not read from '$infile': $!";
	print $line;
	if ($line =~ /^#CHR/) {
		@head = split("\t", $line);
		last;
	}
}

my $j = 0;
my $dpcol;
foreach (@head) {
	if ($_ =~ /^INFO_control/) {
		$dpcol = $j;
	}
	$j++;
}

while (!eof($IN)) {
	defined(my $line = readline($IN))
		|| die "Could not read from '$infile': $!";
	print $line;
	my @l = split("\t", $line);
	my ($dp) = $l[$dpcol] =~ /^DP=(\d+);/;
	next if ($dp == 0);
	if ($dp <= $max) {
		$bins[$dp]++;
	} else {
		$bins[$max]++;
	}
	$count++;
}

my $mid = $count/2;

my $medcount = 0;
my $i = 0;
while ($i < @bins) {
	if(!defined $bins[$i]){
		$bins[$i] = 0;
	}
	$medcount += $bins[$i];
	last if($medcount >= $mid);
	$i++;
}

$median = 5*$i;
close $IN;

open(my $OUT, ">$outfile");
print $OUT $median, "\n";
close $OUT;
