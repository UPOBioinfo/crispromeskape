#!/usr/bin/perl
# AJPerez, 2022
# Create gene presence/absence matrix

use strict;

my $F = $ARGV[0] || "roary/ab2/clustered_proteins";
my $A = $ARGV[1] || "strains.ab";
my $BIN = $ARGV[2] || 0;
my $CAS = $ARGV[3];
my @ab;
my %cas;

# Gather strains IDs
open in, $A;
while (<in>) {
  chomp;
  push @ab, $_;
}
close in;

# Gather cas genes to remove
open in, $CAS;
while (<in>) {
  chomp;
  $cas{$_} = 1;
}
close in;

# Print header
print "#genes/strains\t";
print join "\t", sort @ab;
print "\n";

# Run through clustered_protein
open in, $F;
while (<in>) {
  chomp;

  my ($c1, @c) = split/\t/;
  my ($g, $r) = split/: /, $c1;
  push @c, $r;

  next if $cas{$r};
  my %g;
  for (@c) {
    $_ =~ s/_[0-9]+//;
    $g{$_}++;
  }

  print $r; # reference or gene_name
  foreach my $a (sort @ab) {
    $g{$a} = 0 if !$g{$a};
    $g{$a} = 1 if $g{$a} > 1 && $BIN == 1;
    print "\t$g{$a}";
  }
  print "\n";

}
close in;

exit;

