#!/usr/bin/perl
use strict;

my $file_sp = $ARGV[1]; ###spacers

my $file_dr = $ARGV[0]; ####repeats

open (SP, "$file_sp");

open (DR, "$file_dr");

my @lista_def;

while (my $linea=<DR>){
	chomp $linea;
	my @lista=split (/\t/, $linea);
	my @lista2=split (/,/, $lista[2]);
	push @lista_def, @lista2;
	#print "$lista[2]\n";
}

my %hash;

while (my $line=<SP>){
	chomp $line;
	my @list=split (/\t/, $line);
	foreach my $var (@lista_def){
		chomp $var;
		#print "$var\n";
		if ($line=~/$var/){
			$hash{$list[0]}{"TYPE"}="REP";
			#print "$list[0]\n";
			next;
		}
	}
	$hash{$list[0]}{"LIST"}=$list[2];
	#print "$list[0]\n";
}

foreach my $aa (keys %hash){
	chomp $aa;
	if ($hash{$aa}{TYPE}=~/REP/){
		next;
	}
	print "$aa\t$hash{$aa}{LIST}\n";
}