#!/usr/bin/perl
use strict;

my $AB = $ARGV[0] || "ab";
my %n;
my %gns;
my $nunk = 1;

# Get Sma3s genenames
`mv pangenome_references_$AB\_uniprot_bacteria_go.tsv pangenome_references_$AB\_uniprot_bacteria_go.tsv.bak`;
open in, "pangenome_references_$AB\_uniprot_bacteria_go.tsv.bak";
while (<in>) {
  chomp;

  my ($c1, $c2, @c) = split/\t/;
  my $gn = $c2;
  if ($gn) {
    $n{$c2}++;
    $gn .= "__$n{$c2}" if $n{$c2} > 1;
    $gns{$c1} = $gn;
  } else {
    $gns{$c1} = "unknown$nunk";
    $nunk++;
  }
}
close in;

# Run through clustered_proteins
`mv clustered_proteins clustered_proteins.bak`;
open FILE, ">clustered_proteins";
open in, "clustered_proteins.bak";
while (<in>) {
  chomp;
  
  my ($c1, @c) = split/\t/;
  my ($g, $id) = split/: /, $c1;
  if ($gns{$id} =~ /^unknown[0-9]+$/ && $g !~ /^group_[0-9]+$/) {
    print FILE "_$g: $id\t";
    $gns{$id} = "_$g";
  } elsif ($gns{$id}) {
    print FILE "$gns{$id}: $id\t";
  } else {
    next;
  }
  print FILE join "\t", @c;
  print FILE "\n";
}
close in;
close FILE;

open FILE, ">pangenome_references_$AB\_uniprot_bacteria_go.tsv";
open in, "pangenome_references_$AB\_uniprot_bacteria_go.tsv.bak";
while (<in>) {
  chomp;

  my ($c1, $c2, @c) = split/\t/;
  print FILE "$c1\t$gns{$c1}\t";
  print FILE join "\t", @c;
  print FILE "\n";
}
close in;
close FILE;

`gzip pangenome_references_$AB\_uniprot_bacteria_go.tsv.bak clustered_proteins.bak`;

exit;
