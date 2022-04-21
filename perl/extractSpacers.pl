#!/usr/bin/perl
use strict;

my $S = $ARGV[0] || "sp"; # or dr or flank
my $P = $ARGV[1] || 4; # score
my $F = $ARGV[2] || "."; # path

$S = "CRISPRspacer" if $S eq "sp";
$S = "CRISPRdr" if $S eq "dr";
$S = "LeftFLANK|RightFLANK" if $S eq "flank";

# Run through GFF files
my @ab = `ls -d $F/prokka/[a-z][a-z][0-9][0-9][0-9][0-9][0-9]`;
foreach my $ab (@ab) {
  chomp $ab;
  $ab =~ m/\/([^\/]+)$/;
  my $ab2 = $1;
  
  # CRISPR
  # TSV files (evidence codes)
  my %ec;
  open in, "$F/ccfinder/$ab2/TSV/Crisprs_REPORT.tsv";
  while (<in>) {
    chomp;
    my ($id, $p1, $p2, $ec) = (split/\t/)[1,5,6,26];
    next unless $ec >= $P; # Evidence code
    $id = "$id\_$p1\_$p2";
    $ec{$id} = 1;
  }
  close in;

  # GFF files
  my @gff = `ls -d $F/ccfinder/$ab2/GFF/*.gff`;
  foreach my $gff (@gff) {
    open in, $gff || die "Error! Problem reading $gff\n";
    while (<in>) {
      chomp;
      
      next unless /(\t.+){8}/;
      my ($element,$annot) = (split/\t/)[2,8];
      $annot =~ /Parent=([^;]+)/;
      my $idc = $1;

      next unless $element =~ /^$S$/ && $ec{$idc} == 1;
      #$annot =~ /Name=([^;]+)/;
      $annot =~ /ID=([^;]+)/;
      my $name = $1;
      $name =~ s/spacer_//;
      $name = "$ab2\_$idc\_$name";
      $annot =~ /sequence=([^;]+)/;
      my $seq = $1;
      print ">$name\n$seq\n";
    }
    close in;
  }
}

exit;
