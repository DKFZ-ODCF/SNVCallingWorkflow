#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

use strict;
use warnings;
use v5.10;
use List::Util qw(min);
use constant DEBUG => 0;

my $ref;
my $alt;
my $info;
my $homlen;
my $homseq;
my ($ri, $rl, $ai);
my @fields;
my @alts;
my @homlens;

while (<>) {
  if (/^#/) {
    print $_;
    next;
  }
  chomp;
  DEBUG && say $_;
  @fields = split(/\t/);
  $ref = $fields[3];
  @alts = split(/,/,$fields[4]);
  $ri = $rl = length($ref) -1;
  @homlens = ();
  foreach $alt (@alts) {
    $ri = $rl;
    $ai = length($alt) -1;
    
    while ($ri > 0 && $ai > 0) {
      last if uc(substr($ref, $ri, 1)) ne uc(substr($alt, $ai, 1));
      $ri--;
      $ai--;
    }
    #say $rl-$ri;
    push(@homlens, $rl - $ri);
    #say "@homlens";
  }

  $homlen = min(@homlens);
  $fields[7] .= ";HOMLEN=$homlen";
  if ($homlen) {
    $homseq = uc(substr($ref, $rl-$homlen+1));
    $fields[7] .= ";HOMSEQ=$homseq";
    $ref = substr($ref, 0, -$homlen);
    $fields[3] = $ref;
    foreach $alt (@alts) {
      $alt = substr($alt, 0, -$homlen);
    }
    $fields[4] = join ',',@alts;
  }
  say join "\t", @fields;
}

