#!/usr/bin/perl
# removeAssemblies.pl
# AJPerez, 15/07/2021
# Check assembly files and remove those from GenBank when RefSeq is available

$FOLDER = "ncbi_dataset/data/"; # Folder with genome files

my @F = `ls $FOLDER`;
chomp (@F);
for my $f (@F) {
  my $f1 = $f;
  $f =~ s/\.[0-9]+//; # remove version number

  my $f2 = $f;
  $f2 =~ s/GCA_/GCF_/;

  if ($f1 =~ /^GCA_/ && grep/^$f2\.[0-9]+$/, @F) {
    `rm -rf $FOLDER$f1`;
    print "removing $f1...\n";
  }
}

exit;
