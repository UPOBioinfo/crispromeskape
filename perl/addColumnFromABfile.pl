#!/usr/bin/perl
# AJPerez, 2022
# Create an output to add to the metadata table, from a file with 2 columns: metadata ID.
use strict;

my ($F, $f) = @ARGV; # File with all IDs in the first colum (only IDs that we want; ./prokka/mapping_ids.tsv), and file with 2 columns: ID & Metadata
my %g;

# Abaumannii strains
my @AB;
open in, $F;
while (<in>) {
  chomp;
  my ($ab) = split/\t/;
  push @AB, $ab;
  $g{$ab} = "";
}
close in;

open in, $f || die "Error opening file\n";
while (<in>) {
  chomp;

  my ($ab, $c2) = split/\t/;
  $g{$ab} = $c2;
}
close in;

foreach my $ab (@AB) {
  print "$ab\t$g{$ab}\n";
}

exit;
