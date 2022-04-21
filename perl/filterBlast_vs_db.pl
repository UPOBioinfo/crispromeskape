#!/usr/bin/perl
use strict;

# Filter Blast output:
# blast -query pan_gn.faa -db omp.fasta -evalue 1e-05 
#   -outfmt '6 qseqid sseqid pident qcovs qcovhsp length qlen slen evalue qstart qend sstart send' -max_target_seqs 1'
# AJPerez, 2019/10/19
# Updated, 2022/01/19

# Thresholds
my $ID = $ARGV[1] || 90; # identity
my $SC = $ARGV[2] || 90; # subject (DB) coverage

# Check blast file
die "Please. give me a blast output as argument" if !$ARGV[0];
open INFILE, $ARGV[0];
while (<INFILE>) {
  chomp;

  my ($qseqid, $sseqid, $pident, $qcovs, $qcovhsp, $length, $qlen, $slen, $evalue, $qstart, $qend, $sstart, $send) = split/\t/;
  my $scov = ($length / $slen) * 100;

  next unless ($scov >= $SC && $pident >= $ID);
  print "$qseqid\n";
}
close INFILE;

exit;
