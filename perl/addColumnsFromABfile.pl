#!/usr/bin/perl
# AJPerez, 2022
# Creates an output with several columns to add to the metadata table, from a list of files, each with a list of IDs.

use strict;

my ($F, @f) = @ARGV; # File with all IDs in the first colum (only IDs that we want; ./prokka/mapping_ids.tsv), and files with lists of IDs
my %g;

# Abaumannii strains
my @AB;
open in, $F;
while (<in>) {
  chomp;
  my ($ab) = split/\t/;
  push @AB, $ab;
  @{$g{$ab}} = ();
}
close in;

foreach my $f (@f) {
  chomp $f;
  open in, $f || die "Error opening file\n";
  chomp(my @l = <in>);
  push @{$g{$_}}, $f for @l;
  close in;
}

foreach my $ab (@AB) {
  print "$ab\t";
  print join ",", @{$g{$ab}};
  print "\n";
}

exit;
