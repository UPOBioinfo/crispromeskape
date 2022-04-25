#!/usr/bin/perl
use strict;
# collectGFFMetadata.pl
# AJPerez, 30/06/2021
# Inputs: metadata from NCBI (dataformat tsv genome --package ac_dh.zip > assembly_isolates_ac.tsv),
#   GFF file from GenBank (datasets download genome taxon 'acinetobacter calcoaceticus' --filename ac.zip)
#   and Pathogen detection table (including BioSample metadata)
# GCA_/GCF_ is removed from the IDs, to get more data from Pathogen Detection

my $MD = $ARGV[0] || "assembly_isolates_ec.tsv";
my $PD = $ARGV[1] || "pathogen_isolates_ec.tsv";
my $IS = $ARGV[2] || "isolation_sources.tsv";
my $ID = $ARGV[3] || "prokka/mapping_ids.tsv";
my $RM = $ARGV[4] || "roary/ec/clustered_proteins_ngenes.tsv";
$MD =~ /(..\.tsv)$/; my $OUT = "metadata_$1";
my %is; my %di; my %id; my %id2; my %rm; my %pd; my %prokka; my %assembly; my %md;
my %is2; my %di2;
my @names;
my $pd_header;
my %L;
my %cr;

# Collect strain identifiers
print "Collecting strain identifiers...";
open in, $ID;
while (<in>) {
  chomp;

  my ($c1, $c2) = split/\t/;
  $id2{$c1} = $c2; # iD complete, con prefix GCA_/GCF_
  $c2 =~ s/[A-Z]{3}_//;
  $id{$c1} = $c2;
}
close in;

# Count number of protein-coding genes by Prokka (AB)
print "ok\nCollecting faa from prokka...";
my @prokka = `ls /mnt/beegfs/ajperez/pangenomes/ab/prokka/*/*.faa /mnt/beegfs/ajperez/pangenomes/ab/prokka/removables/*/*.faa`;
foreach my $p (@prokka) {
  chomp $p;

  $p =~ /\/([a-z]{2}[0-9]{5})\./;
  my ($ab) = $1;
  my $n = `grep ">" $p | wc -l`; chomp $n;
  $prokka{$ab} = $n;

  # inicializate variables
  $rm{$ab} = "";
  $pd{$id{$ab}} = "";
  @{$cr{$ab}} = ();
  $assembly{$id{$ab}} = ""; 
}

# Collect isolation sources (IS)
print "ok\nCollecting isolation sources...";
open in, $IS;
while (<in>) {
  chomp;

  my ($c1, $c2, $c3, $c4, $c5, $c6) = (split/\t/)[1,2,3,4,5];
  my $c0 = "$c1\t$c2";
  $is{$c0} = $c3;
  $di{$c0} = $c4;
  $is2{$c0} = $c5;

print "$c0\t$c1\t$c2\n";
}
close in;

# Collect ngenes (AB)
print "ok\nCollecting ngenes...";
open in, $RM;
while (<in>) {
  chomp;

  my ($c1, $c2, $c3, $c4) = split/\t/;
  $rm{$c1} = "$c2\t$c3\t$c4";
}
close in;

# Run thorough the GFF files (ID)
print "ok\nCollecting GFFs...";
my @GFF = `ls /mnt/beegfs/ajperez/pangenomes/ab/ncbi_dataset/data/*/*gff`;
foreach my $file (@GFF) {
  chomp $file;

  my ($as) = (split/\//, $file)[2];
  $as =~ s/[A-Z]{3}_//;
  open in, $file;
  while (<in>) {
    chomp;

    my ($r, $md) = (split/\t/)[2,8];
    next unless $r eq "region";
    my (@md) = split/;/, $md;
    foreach my $m (@md) {
      my ($n, $v) = split/=/, $m;
      $md{$as}{$n} = $v;
      $L{$as} = lc $md{$as}{$n} if $n eq "isolation-source"; # useful when assembly isn't in PD
      push @names, $n unless grep /^$n$/, @names;
    }

    last;
  }
  close in;
}

# Collect pathogen detection metatable
print "ok\nCollecting pathogen detection metadata...";
open in, $PD;
while (<in>) {
  chomp;

  if (/^#/) {
    $pd_header = $_;
  } else {
    my ($is,$id,$hd) = (split/\t/)[7,13,16];
    $id =~ s/[A-Z]{3}_//;
    #next unless $id{$id};
    $pd{$id} = $_;
    $is = lc $is;
    
    #$L{$id} = lc $md{$id}{"isolation-source"}; # it is more better above
    if (exists $md{$id}{"isolation-source"}) {
      if ($L{$id} && $is && $L{$id} ne $is) {
        $L{$id} = "$L{$id} ## $is";
      } elsif ($is) {
        $L{$id} = $is;
      } # else, the original from GFF file
    }
    $L{$id} .= "\t$hd";
  }
}
close in;

# Collect CRISPR-Cas systems from cctyper
print "ok\nCollecting cctyper...";
#my @cct = `ls /mnt/beegfs/ajperez/pangenomes/ab/cctyper/*/CRISPR_Cas.tab`;
my @cct = `ls /mnt/beegfs/ajperez/pangenomes/ab/cctyper/*/cas_operons.tab`;
foreach my $file (@cct) {
  chomp $file;
  
  my ($id) = (split/\//, $file)[1];
  open in, $file;
  while (<in>) {
    chomp;

    #next if /^Contig\t/;
    next if /^Contig\t/;
    #my ($cr) = (split/\t/)[6];
    my ($cr, $int, $ada, $gs) = (split/\t/)[4,5,6,9];
    $int =~ s/\%//; $ada =~ s/\%//;
    $gs =~ s/[\'\[\]\s]//g;
    my (@gs) = split/,/, $gs;
    @gs = join ";", sort @gs;
    push @{$cr{$id}}, "$cr $int $ada $gs";
  }
  close in;
}

# Run through the metadata table
print "ok\nCollecting Assembly metadata...";
open in, $MD;
my $header = <in>;
chomp $header;
while (<in>) {
  chomp;

  my ($id) = (split/\t/)[17];
  #next unless $id2{$id};
  $id =~ s/[A-Z]{3}_//;
  #next unless $rm{$id{$id}};
  $assembly{$id} = $_;
}
close in;

# Final output
print "ok\nFinal calculations...";
open FILE, ">$OUT";
print FILE "ID\t$header\t";
print FILE join "\t", @names;
print FILE "\t$pd_header\tProkka_proteins\tNgenes\tNshared\tRemovable\tCollection_year\tIsolation_source\tDisease\tHOST\tBSI?\tCCType\tInterference\tAdaptation\tCCGenes\n";
my ($n_pd) = $pd_header =~ s/\t/\t/g; 
my ($n_as) = $header =~ s/\t/\t/g; 
foreach my $ab (sort keys %prokka) {
  my $id = $id{$ab};
 
  $assembly{$id} = "\t"x$n_as if !$assembly{$id}; 
  print FILE "$ab\t$assembly{$id}";
  foreach my $n (@names) {
    print FILE "\t$md{$id}{$n}";
  }

  $pd{$id} = "\t"x$n_pd if !$pd{$id};
  my ($year) = (split/\t/, $pd{$id})[29];
  $year =~ /([0-9]{4})/;
  $year = $1;
  $L{$id} .= "\t" if ($L{$id} && $L{$id} !~ /\t/); # is without pathogen detection line
  $is{$L{$id}} = "other" if !$is{$L{$id}}; # class 'other' for isolation source
  $di{$L{$id}} = "na" if !$di{$L{$id}}; # class 'na' for disease
  $is2{$L{$id}} = "na" if !$is2{$L{$id}}; # class 'other' for isolation source
  $di2{$L{$id}} = "na" if !$di2{$L{$id}}; # class 'na' for 
  print FILE "\t$pd{$id}\t$prokka{$ab}\t$rm{$ab}\t$year\t$is{$L{$id}}\t$di{$L{$id}}\t$is2{$L{$id}}";

  if (@{$cr{$ab}}) {
    my %x;
    for (my $x = 0; $x <= $#{$cr{$ab}}; $x++) {
      my (@c) = split/ /, ${$cr{$ab}}[$x];
      for (my $y = 0; $y <= 3; $y++) {
        push @{$x{$y}}, $c[$y];
        #@{$x{$y}} = sort @{$x{$y}} if $y == 0;
      }
    }
    for (my $y = 0; $y <= 3; $y++) {
      print FILE "\t";
      print FILE join ",", @{$x{$y}};
    }
  } else {
      print FILE "\t"x4;
  }
  print FILE "\n";
}
close in;
close FILE;

print "ok\nFile $OUT created\n\n";

exit;
