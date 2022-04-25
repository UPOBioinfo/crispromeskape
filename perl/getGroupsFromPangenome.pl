#!/usr/bin/perl
# AJPerez, 2022
# Count appearances of each gene in the pangenome

use strict;

my $C  = $ARGV[0] || "roary/ab2/clustered_proteins";
my %ref;

# clustered_proteins
open in, $C;
while (<in>) {
  chomp;

  my ($c1, @c) = split/\t/;
  my ($c1) = (split/: /, $c1)[1];
  
  push @c, $c1;
  @c = map { s/_[0-9]+//; $_ } @c; # gather strain ab
  my %c = map {$_ => 1} @c; @c = sort keys %c; # remove redundancy

  print "$c1\t";
  print scalar(@c);
  print "\n";
}
close in;

exit;
