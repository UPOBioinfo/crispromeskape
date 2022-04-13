#!/usr/bin/perl
use strict;

my $G1 = $ARGV[0]; # File 1
my $G2 = $ARGV[1]; # File 2;
my $C  = $ARGV[2] || "roary/ab2/clustered_proteins";

my (@G1_) = split/\//, $G1;
my (@G2_) = split/\//, $G2;
$G1_[$#G1_] =~ s/\.[^\.]+//;
$G2_[$#G2_] =~ s/\.[^\.]+//;
my $outfile = "./diffs/$G1_[$#G1_]\_$G2_[$#G2_].tsv";

my %g1;
my %g2;

open in, $G1;
chomp (my @G1 = <in>);
my $NG1 = scalar(@G1);
close in;

open in, $G2;
chomp (my @G2 = <in>);
my $NG2 = scalar(@G2);
close in;

# clustered_proteins
open FILE, ">$outfile";
print FILE "#ID\tN1\tN2\tFreq1\tFreq2\tDiff\tDiff_abs\n";
open in, $C;
while (<in>) {
  chomp;

  my ($c1, @c) = split/\t/;
  my ($c1) = (split/: /, $c1)[1];

  push @c, $c1;
  @c = map { s/_[0-9]+//; $_ } @c; # gather strain ab
  my %c = map {$_ => 1} @c; @c = sort keys %c; # remove redundancy

  my %seen = map { $_ => 1 } @c;
  my $nc1 = scalar(grep { $seen{$_} } @G1);
  my $nc2 = scalar(grep { $seen{$_} } @G2);
  
  printf FILE "%s\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\n", $c1, $nc1, $nc2, $nc1/$NG1, $nc2/$NG2, ($nc1/$NG1)-($nc2/$NG2), abs(($nc1/$NG1)-($nc2/$NG2));
}
close in;
close FILE;

exit;
