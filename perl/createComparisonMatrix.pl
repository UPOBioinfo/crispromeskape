#!/usr/bin/perl
# AJPerez, 2022
# Create matrix for heatmaps comparing genes vs genomes

use strict;

my $G1 = $ARGV[0] || "pan_omp_9090.id";
my $G2 = $ARGV[1] || "vir_phages.id";
my $C  = $ARGV[2] || "roary/ab2/clustered_proteins";

my %g1;
my %g2;
my %ref;
my %good;

open in, $G1;
chomp (my @G1 = <in>);
foreach (@G1) { $g1{$_} = 1 }
close in;

open in, $G2;
chomp (my @G2 = <in>);
foreach (@G2) { $g2{$_} = 1 }
close in;

# clustered_proteins
open in, $C;
while (<in>) {
  chomp;

  my ($c1, @c) = split/\t/;
  my ($c1) = (split/: /, $c1)[1];

  next unless $g1{$c1} || $g2{$c1};
  push @c, $c1;
  @c = map { s/_[0-9]+//; $_ } @c; # gather strain ab
  my %c = map {$_ => 1} @c; @{$ref{$c1}} = sort keys %c; # remove redundancy
  $good{$c1} = 1;
#  $good{$c1} = 1 if scalar(@c) > 100 && scalar(@c) < 5000; # filter by number of strains
}
close in;

# pares
print "#g1/g2\t";
print join "\t", @G2;
print "\n";
foreach my $g1 (@G1) {
  next unless $good{$g1};
  print $g1;
  foreach my $g2 (@G2) {
    next unless $good{$g2};
    my %seen = map { $_ => 1 } @{$ref{$g1}};
    my @common = grep { $seen{$_} } @{$ref{$g2}};
    print "\t";
    
  # Number of strains with the pair
  # my $p = scalar(@common);

  # Percent regart to G2
  my $p = (scalar(@common) / scalar(@{$ref{$g2}})) * 100;

  # Print out
  printf "%.2f", $p;
  }
  print "\n";
}

exit;
