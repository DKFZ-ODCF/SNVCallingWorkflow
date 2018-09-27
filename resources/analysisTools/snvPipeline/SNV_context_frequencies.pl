#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

use strict;
use warnings;
use v5.10;

if (@ARGV < 2)
{
	die "vcf file to generate transition transversion rates from - minimal confidence score (recommended: 9)\n";
}

my $file = shift;
my $minconfidence = shift;

open (FH, $file) or die "Could not open $file: $!\n";
my $header;
while ($header = <FH>)
{
	last if ($header =~ /^\#CHR/); # that is the line with the column names
}

chomp($header);
my @columns = split(/\t/, $header);

my $CLASSIFICATION;
my %header_hash = map { $_ => 1 } @columns;
if (exists $header_hash{ANNOTATION_control}) {
   $CLASSIFICATION = "ANNOTATION_control";
} else {
   $CLASSIFICATION = "RECLASSIFICATION";
}

my %fields;
my %counts;
say join "\t", (qw(REF ALT preceeding following type exonic value));
my ($ref, $alt, $pre, $fo, $ty, $ex);

while (<FH>)
{
  chomp;
  # make hash of array with header lines => header names are keys
  @fields{@columns} = split(/\t/);
  # skip low confidence calls
  next if ($fields{CONFIDENCE} < $minconfidence);
  # skip unclear, LQVSIG, and multi_ (bc. of eq instead of =~ //)
  next unless ($fields{$CLASSIFICATION} eq 'germline' or $fields{$CLASSIFICATION} eq 'somatic');
  $fields{type} = ($fields{$CLASSIFICATION} eq 'somatic') ? 'somatic' : 'germline';
  $fields{exonic} = ($fields{ANNOVAR_FUNCTION} eq 'exonic') ? 'exonic' : 'non-exonic';
  $fields{REF} = uc($fields{REF});
  $fields{ALT} = uc($fields{ALT});
  $fields{ALT} =~ s/,.+//;	# get rid of multiallelic - but that was already skipped above ?!
  #$fields{SEQUENCE_CONTEXT} = uc($fields{SEQUENCE_CONTEXT});	# is upper case anyways
  # the entry looks lie this: TCCGCCTCTC,ATTCCCTTCC => preceeding base is C, following is A
  $fields{SEQUENCE_CONTEXT} =~ /.{9}(.)(.)(.).{9}/ or die "SHit with seq context";
  if ($fields{REF} eq 'G' || $fields{REF} eq 'A')
  {
    $fields{REF} =~ tr/ACGT/TGCA/;
    $fields{ALT} =~ tr/ACGT/TGCA/;
    ($fields{preceeding}, $fields{following}) = ($3,$1);
    $fields{preceeding}=~ tr/ACGT/TGCA/;
    $fields{following} =~ tr/ACGT/TGCA/;
  }
  else
  {
    ($fields{preceeding}, $fields{following}) = ($1,$3);
  }
  $counts{join ':', @fields{(qw(REF ALT preceeding following type exonic))}}++;
}
close FH;
for $ref ('T','C') {
  for $alt ('A','C','G','T') {
    next if ($alt eq $ref);
    for $pre ('A','C','G','T') {
      for $fo ('A','C','G','T') {
        for $ty ('somatic', 'germline') {
          for $ex ('exonic', 'non-exonic') {
            if (!defined($counts{join ':', $ref,$alt,$pre,$fo,$ty,$ex})) {
              $counts{join ':', $ref,$alt,$pre,$fo,$ty,$ex} = 0;
            }
          }
        }
      }
    }
  }
}
foreach my $key (keys %counts) {
  say join "\t", split(':', $key), $counts{$key};
}
