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

open(my $SNV, ">$snv_of") || die "Could not open $snv_of for writing ($!)";
open(my $INDEL, ">$indel_of")  || die "Could not open $indel_of for writing ($!)";

while(!eof(\*STDIN)) {
  defined (my $line = readline(\*STDIN))
      || die "Could not read from STDIN: $!";
  if ($line =~/^\#/) {
    print $SNV $line;
    print $INDEL $line;
    next;
  }
  if ((split(/\t/, $line))[7] =~ /^INDEL/) {
    print INDEL $line;
  } else {
    print SNV $line;
  }
}
close $SNV;
close $INDEL;
