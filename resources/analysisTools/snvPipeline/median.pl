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

open(IN, "<$ARGV[0]");

my $count=0;
my @head;
while(<IN>){
	print $_;
	if($_ =~ /^#CHR/){
		@head = split("\t", $_);
		last;
	}
}

my $j = 0;
my $dpcol;
foreach(@head){
	if($_ =~ /^INFO_control/){
		$dpcol = $j;
	}
	$j++;
}

while(<IN>){
	print $_;
	my @l = split("\t", $_);
	my ($dp) = $l[$dpcol] =~ /^DP=(\d+);/;
	next if($dp == 0);
	if($dp <= $max){
		$bins[$dp]++;
	}else{
		$bins[$max]++;
	}
	$count++;
}

my $mid = $count/2;

my $medcount = 0;
my $i = 0;
while($i < @bins){
	if(!defined $bins[$i]){
		$bins[$i] = 0;
	}
	$medcount += $bins[$i];
	last if($medcount >= $mid);
	$i++;
}

$median = 5*$i;
close IN;
open(OUT, ">$outfile");
print OUT $median, "\n";
close OUT;
