#!/usr/bin/perl
# countGenesAndShareGenesFromPangenome.pl
# AJPerez, 29/06/2021
# Inputs: clusters from Roary, and I (number of times the IQR)
# Count genes and average number of shared genes for each genome in the pangenome

use strict;
use Statistics::Descriptive;

my $F = $ARGV[0] || "clustered_proteins";
my $I = $ARGV[1] || 1.5;
my %ng; # number of genes/orthologs/clusters
my %sh; # number of shared genes
my $n = 1; # row (genes analyzed)

# Read clusters
open in, $F;
while (<in>) {
  chomp;

  my ($c1, @c) = split/\t/;
  $c1 =~ s/[^:]+: //;
  unshift @c, $c1;

  @c = map { s/_[0-9]+//; $_ } @c; # only strains
  my %c = map {$_ => 1} @c; my @c = sort keys %c; # remove redundancy

  # Count
  #print "gene $n\n"; $n++;
  for (my $x = 0; $x <= $#c; $x++) {
    $ng{$c[$x]}++;
    next if $x == $#c; # next last strain
    for (my $y = $x + 1; $y <= $#c; $y++) {
      $sh{$c[$x]}{$c[$y]}++;
    }
  }
}
close in;

# Print output
my @ng;
my @sh;
my %out; my %outm;
foreach my $str (sort keys %ng) {
  my $o;
  my @o;
  my $oo;
  $out{$str} .= "$ng{$str}\t";
  push @ng, $ng{$str};
  foreach my $str2 (sort keys %ng) {
    if ($sh{$str}{$str2}) {
      $oo = $sh{$str}{$str2};
    } elsif ($sh{$str2}{$str}) {
      $oo = $sh{$str2}{$str};
    } elsif ($str eq $str2) {
      $oo = $ng{$str};
    } else {
      $oo = 0;
    }
    $o .= "\t$oo";
    push @o, $oo;
  }

  my $stat = Statistics::Descriptive::Full->new();
  $stat->add_data(@o);
  my $sh = sprintf "%.2f", $stat->mean();
  push @sh, $sh;
  $out{$str} .= $sh;
  $outm{$str} .= $o;
}

# Filtered
my $ng_threshold = &threshold(@ng);
my $sh_threshold = &threshold(@sh);

# Check if filtered
open FILE, ">$F\_ngenes.tsv";
open MATRIX, ">$F\_matrix.tsv";
open REMOVE, ">removables.id";
print MATRIX "#\t";
print MATRIX join "\t", sort keys %ng;
print MATRIX "\n";
foreach my $str (sort keys %ng) {
  my $remove;
  print FILE "$str\t$out{$str}\t";
  print MATRIX "$str";
  my ($c1, $c2) = split/\t/, $out{$str};
  if ($c1 < $ng_threshold) {
    print FILE "1";
    $remove = 1;
  } 
  if ($c2 < $sh_threshold) {
    print FILE "2";
    $remove = 1;
  }
  print FILE "0" if !$remove;
  print REMOVE "$str\n" if $remove;
  print MATRIX "$outm{$str}";
  print MATRIX "\n";
  print FILE "\n";
}
close FILE;
close MATRIX;
close REMOVE;

print "Threshold (number of genes): $ng_threshold\n";
print "Threshold (shared genes): $sh_threshold\n";
print "Files $F\_matrix.tsv, $F\_ngenes.tsv and removables.id created\n";

exit;

##############
# Subrutines #
##############
sub threshold () {
  my @arr = @_;

  my $stat = Statistics::Descriptive::Full->new();
  $stat->add_data(@arr);
  my $q1 = $stat->quantile(1);
  my $q3 = $stat->quantile(3);
  return $q1 - $I * ($q3 - $q1);
}
