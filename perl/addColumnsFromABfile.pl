#!/usr/bin/perl
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
