#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#
use strict;
use warnings;
use v5.10;

my $snv_of = shift @ARGV;
my $indel_of = shift @ARGV;

open(SNV, ">$snv_of") || die "Could not open $snv_of for writing ($!)";
open(INDEL, ">$indel_of")  || die "Could not open $indel_of for writing ($!)";

while(<>) {
  if (/^\#/) {
    print SNV;
    print INDEL;
    next;
  }
  if ((split(/\t/))[7] =~ /^INDEL/) {
    print INDEL;
  } else {
    print SNV;
  }
}
close SNV;
close INDEL;
