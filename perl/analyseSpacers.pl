#!/usr/bin/perl
# AJPerez, 2022
# Calculate proportions of genomes with unique spacer/phages in two different clusters
use strict;

my $SP = $ARGV[0] || "ab";
my $TP = $ARGV[1] || "ifb";
my %cp;  # clustered_proteins
my %c;   # clusters
my %c_c; # clusters together
my %v;   # viral genes (all)
my %vd;  # viral genes (diff)
my %v2;  # final viral genes
my @cc = ("c1", "c2");
my %r;   # result

# Gather clusters
foreach my $cc (@cc) {
  open in, "$SP/$cc\_$TP";
  chomp (my @C = <in>);
  foreach (@C) { $c{$cc}{$_} = 1 }
  close in;
}
%c_c = (%{$c{$cc[0]}}, %{$c{$cc[1]}});
my $ncc1 = scalar keys %{$c{$cc[0]}};
my $ncc2 = scalar keys %{$c{$cc[1]}};

# Gather viral genes
open in, "$SP/virus_sma3s_db_btw1_0.id";
chomp (my @C = <in>);
foreach (@C) { $v{$_} = 1 }
close in;

# Gather clustered_proteins
my $CP = "$SP/roary/$SP" . "2/clustered_proteins";
open in, $CP;
while (<in>) {
  chomp;

  my ($c1, @c) = split/\t/;
  my ($c1) = (split/: /, $c1)[1];

  push @c, $c1;
  @c = map { s/_[0-9]+//; $_ } @c; # gather strain ab
  my %i = map {if ($c_c{$_}) { $_ => 1 } else { () } } @c; @{$cp{$c1}} = sort keys %i; # remove redundancy
}
close in;

# Gather specific viral genes from diffs folder (I)
open in, "$SP/diffs/c1_$TP\_c2_$TP.tsv";
while (<in>) {
  chomp;

  my ($id, $f1, $f2) = split/\t/;
  next unless $v{$id};

  # cluster 1
  if ($f2 == 0) { # if it is not found in c2
    $vd{"c1"}{$id} = 1;
  }
  # cluster 2
  if ($f1 == 0) { # if it is not found in c1
    $vd{"c2"}{$id} = 1;
  }
}
close in;

# Read spacers
my %spacers; # all spacers
my %sp; # spacers by cluster
my %spref; # references by cluster
open in, "spacers_vsall_$SP\_95.tsv"; # if we remove _95, 90% is taken
while (<in>) {
  chomp;

  my ($id, $ref, $sp) = (split/\t/)[0,8,10];
  next unless $ref;
  my (@ref) = split/ /, $ref;
  my (@sp) = split/ /, $sp;

  # Check strains belonging only to one cluster
  my @s1 = grep { $c{"c1"}{$_} } @sp;
  my @s2 = grep { $c{"c2"}{$_} } @sp;

  next unless ((@s1 && !@s2) || (@s2 && !@s1));
  my $cc = "c1";
  my $s1 = join " ", @s1;
  if (@s2 && !@s1) {
    $cc = "c2";
    $s1 = join " ", @s2;
  }

  # Only viral genes
  my @v1 = grep { $vd{$cc}{$_} } @ref; # only frequent viral genes in cluster
  next unless (@v1);

  # Gather columns from spacers, but now only viral genes (v1) and strains from cluster (s1)
  my $genes = join " ", @v1;
  $spacers{$id} = "$genes\t$s1";
  
  # Gather spacer IDs and viral genes for later analysis
  push @{$sp{$cc}}, $id;
  foreach my $v1 (@v1) { $spref{$cc}{$v1} = 1 }
}
close in;

# Compare spacers from clusters
open FILE, ">$SP\_$TP\_0_elements.tsv";
foreach my $cc1 (@cc) {
  my $cc2 = "c2";
  if ($cc1 eq "c2") { $cc2 = "c1" }
  foreach my $id (@{$sp{$cc1}}) {
    my ($v1, $s1) = split/\t/, $spacers{$id};
    my (@v1) = split/ /, $v1;
    my (@s1) = split/ /, $s1;
    next if (grep { $spref{$cc2}{$_} } @v1);

    # Assign as spacer (cluster 1 or 2)
    foreach (@s1) { $r{$cc1}{$_}->{spacer} = 1 }
    
    # Gather strains from viral genes
    my @str;
    foreach (@v1) {
      @str = (@str, @{$cp{$_}});
      print FILE "$_\t$cc1\tspacer\t", $#s1 + 1, "\t", join " ", @s1, "\n";
      next if (!@str);
      print FILE "$_\t$cc1\tphage\t", $#str + 1, "\t", join " ", @str, "\n";
    }
    foreach (@str) { $r{$cc1}{$_}->{virus} = 1 } 

    # Assign as together
    my %s1 = map { $_ => 1 } @s1;
    my @tog = grep { $s1{$_} } @str;
    foreach (@tog) {
      $r{$cc1}{$_}->{together} = 1;
    }

    # Finally print togethers
    next if (!@tog);
    foreach (@v1) {
      print FILE "$_\t$cc1\ttogether\t", $#tog + 1, "\t", join " ", @tog, "\n";
    }
    #print "$cc1\t$id\n"; # spacers by cluster
  }
}
close FILE;

# Final table
open FILE, ">$SP\_$TP\_0_suppl.tsv";
my %o = ("c1" => 0, "c2" => 0);
my %s = ("c1" => 0, "c2" => 0); 
my %p = ("c1" => 0, "c2" => 0);
my %b = ("c1" => 0, "c2" => 0);
my %t = ("c1" => 0, "c2" => 0);
foreach my $cc (@cc) {
  foreach my $c1 (keys %{$c{$cc}}) {
    if ($r{$cc}{$c1}->{together}) {
      $t{$cc}++;
    } elsif ($r{$cc}{$c1}->{spacer} && $r{$cc}{$c1}->{virus}) {
      $b{$cc}++;
    } elsif ($r{$cc}{$c1}->{spacer}) {
      $s{$cc}++;
    } elsif ($r{$cc}{$c1}->{virus}) {
      $p{$cc}++;
    } else {
      $o{$cc}++;
    }
  }
}

print FILE <<OUT;
Strains\tcluster_1\tcluster_2
no_hits\t$o{"c1"}\t$o{"c2"}
together\t$t{"c1"}\t$t{"c2"}
both\t$b{"c1"}\t$b{"c2"}
spacer\t$s{"c1"}\t$s{"c2"}
phage\t$p{"c1"}\t$p{"c2"}
total\t$ncc1\t$ncc2
OUT

close FILE;

exit;
