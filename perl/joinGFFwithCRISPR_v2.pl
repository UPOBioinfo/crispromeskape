#!/usr/bin/perl
use strict;

my $F = $ARGV[0] || "/mnt/data/research/aba";
my %gn;
my %dr;
my $G = "ssrA_dr";
$gn{"CRISPR"} = "CRISPR";
$gn{$G} = $G;

# Read pangenome annotation
my $roary = `ls $F/roary/??2/clustered_proteins`;
open in, $roary;
while (<in>) {
  chomp;
  my ($id, @ids) = split/\t/;
  $id =~ s/^(.+): (.+)//;
  push @ids, $2; # reference id
  $gn{$_} = $1 for @ids;
}
close in;

# Read direct repeat (ssrA)
my $id_old;
my $p1_old;
my $p2_old;
open in, "$F/ssra_3p_filtered.blast";
while (<in>) {
  chomp;
  my ($id, $p1, $p2) = split/\t/;
  my $s = "+";
  if ($p1 > $p2) {
    ($p1, $p2) = ($p2, $p1);
    $s = "-";
  }

  next if ($id eq $id_old && ($p1 == $p1_old || $p2 == $p2_old)); # pass repeated shorter sequences
  push @{$dr{$id}}, "$G\t$p1\t$p2\t$s";
  $id_old = $id;
  $p1_old = $p1;
  $p2_old = $p2;
}
close in;

# Run through GFF files
mkdir("$F/genes") unless(-d "$F/genes");
my @ab = `ls -d $F/prokka/[a-z][a-z][0-9][0-9][0-9][0-9][0-9]`;
open FILE, ">$F/genes/aba_genes_$G.tsv" || die "Error! Problem creating gff files\n";
foreach my $ab (@ab) {
  chomp $ab;
  $ab =~ m/\/([^\/]+)$/;
  my $ab2 = $1;
  my %l;
  
  # Prokka
  open in, "$ab/$ab2.gff" || die "Error! Problem reading $ab/$ab2.gff\n";
  while (<in>) {
    chomp;
    
    next unless /(\t.+){8}/;
    my ($id, $source, $type, $p1, $p2, $sign, $annot) = (split/\t/)[0,1,2,3,4,6,8];
    next if $source =~ /^minced:/;

    $annot =~ /ID=([^;]+)/;
    my $idg = $1;
    push @{$l{$id}}, "$idg\t$p1\t$p2\t$sign";
    if ($type eq "tRNA") { $annot =~ /product=([^;]+)/; $gn{$idg} = $1; }
    if ($type eq "tmRNA") { $annot =~ /gene=([^;]+)/; $gn{$idg} = $1; }
  }
  close in;

  # CRISPR
  # TSV files (evidence codes)
  my %ec;
  open in, "$ab/../../ccfinder/$ab2/TSV/Crisprs_REPORT.tsv" || die "Error! Problem reading ccfinder files\n";
  while (<in>) {
    chomp;
    my ($id, $ec) = (split/\t/)[4,26];
    next unless $ec >= 4;
    $ec{$id} = $ec;
  }
  close in;

  # GFF files
  my @gff = `ls -d $F/ccfinder/$ab2/GFF/*.gff`;
  foreach my $gff (@gff) {
    open in, $gff || die "Error! Problem reading $gff\n";
    while (<in>) {
      chomp;
      
      next unless /(\t.+){8}/;
      my ($id,$element,$p1,$p2,$sign,$annot) = (split/\t/)[0,2,3,4,6,8];
      $annot =~ /ID=([^;]+)/;
      my $idc = $1; $idc =~ s/Crispr_//;
      next unless $element eq "CRISPR" && $ec{$idc};
      push @{$l{$id}}, "CRISPR\t$p1\t$p2\t$ec{$idc}";
    }
    close in;
  }

  # Create output
  foreach my $id (keys %l) {
    @{$l{$id}} = (@{$l{$id}}, @{$dr{$id}}) if ($dr{$id}); # ADD direct repeats from ssrA
    print FILE "$ab2\t$id\t";
    my @idg = ();
#    my @pos = ();
    foreach my $l (sort { (split(/\t/, $a))[1] <=> (split(/\t/, $b))[1] } @{$l{$id}}) {
      my ($idg, $p1, $p2, $sign) = split/\t/, $l;
#      push @pos, "$p1-$p2;$sign";
      push @idg, "$gn{$idg}\{$idg $p1 $p2 $sign\}";
    }
#    print FILE join ",", @pos;
#    print FILE "\t";
    print FILE join ",", @idg;
    print FILE "\n";
  }
}
close FILE;

exit;
