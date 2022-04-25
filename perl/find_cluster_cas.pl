#!/usr/bin/perl
use strict;

my $file_annot = $ARGV[0];
my $matrix= $ARGV[1];

open (IN, "$file_annot");
open (OUT_LISTA, ">lista_id_cas.txt");

my %hash;

while (my $linea=<IN>){
	chomp $linea;
	next if($linea=~/#ID/);
	my @lista=split (/\t/, $linea);
	$hash{$lista[0]}{"SM"}="";
	if ($linea=~/\tcas[1-9](.+)\t|\tcas[A-E](.+)\t/ or $linea=~/CRISPR/){
		#print "$linea\n";
		$hash{$lista[0]}{"SM"}="sma3s";
		#print "$lista[0]\tCRISPR_SMA3S\n";
		print OUT_LISTA "$lista[0]\n";
	}
}

my @strains_ccfinfer=`ls /mnt/beegfs/ajperez/pangenomes/kp/ccfinder/`; ##evidence level column 27

foreach my $var (@strains_ccfinfer){
	chomp $var;
	print "$var\t";
	############################################################
	open (IN2, "/mnt/beegfs/ajperez/pangenomes/kp/ccfinder/$var/TSV/Crisprs_REPORT.tsv");
	open (OUT, ">/mnt/beegfs/ajperez/pangenomes/kp/ccfinder/$var/TSV/Crisprs_REPORT_def.gff");
	while (my $line=<IN2>){
		chomp $line;
		my @valores=split (/\t/, $line);
		if ($valores[26]==4){
			my $ori=".";
			my $ini=$valores[5];
			my $fin=$valores[6];
			if ($valores[8]=~/Forward/){
				$ori="+";
				$ini=$valores[5]-10000;
			} elsif ($valores[8]=~/Backward/){
				$ori=".";
				$fin=$valores[6]+10000;
			}
			print OUT "$valores[1]\tCCFINDER_EVL4\tCAS\t$ini\t$fin\t\.\t$ori\t\.\t$valores[11]\n";
			print "$ini\t$fin\t$valores[11]\t";
		}
	}
	close IN2;
	close OUT2;
	my $gff=$var.".gff";
	my $gff_def=$var."_def.gff";
	open (IN3, "/mnt/beegfs/ajperez/pangenomes/kp/prokka/$var/$gff");
	open (OUT2, ">/mnt/beegfs/ajperez/pangenomes/kp/prokka/$var/$gff_def");
	while (my $line=<IN3>){
		chomp $line;
		last if ($line=~/##FASTA/);
		print OUT2 "$line\n";
	}
	close IN3;
	close OUT2;
	my $out=$var."_result";
	`bedtools intersect -b /mnt/beegfs/ajperez/pangenomes/kp/ccfinder/$var/TSV/Crisprs_REPORT_def.gff -a /mnt/beegfs/ajperez/pangenomes/kp/prokka/$var/$gff_def > /mnt/beegfs/ajperez/pangenomes/kp/ccfinder/$var/TSV/$out`;
	open (IN4, "/mnt/beegfs/ajperez/pangenomes/kp/ccfinder/$var/TSV/$out");
	while (my $linea=<IN4>){
		chomp $linea;
		if ($linea=~/note=CRISPR/){
			$linea=~s/(.+)note=//g; $linea=~s/CRISPR(.+)/CRISPR/g;
			print "$linea ";
		} else {
			$linea=~s/(.+)ID=//g; $linea=~s/;(.+)//g;
			print OUT_LISTA "$linea\n";
			print "$linea ";
		}
	}
	print "\t";
	close IN4;
	############################################
	open (IN5, "/mnt/beegfs/ajperez/pangenomes/kp/cctyper/$var/cas_operons.tab");
	open (OUT3, ">/mnt/beegfs/ajperez/pangenomes/kp/cctyper/$var/cas_operons_def.tab");
	while (my $line=<IN5>){
		chomp $line;
		next if ($line=~/^Contig(.+)/);
		my @valors=split (/\t/, $line);
		print "$valors[2]\t$valors[3]\t$valors[4]\t";
		print OUT3 "$valors[0]\tCCTYPER\tCAS\t$valors[2]\t$valors[3]\t\.\t\.\t\.\t$valors[4]\n";
	}
	close IN5;
	close OUT3;
	my $out_cctyper=$var."_result_cctyper";
	`bedtools intersect -b /mnt/beegfs/ajperez/pangenomes/kp/cctyper/$var/cas_operons_def.tab -a /mnt/beegfs/ajperez/pangenomes/kp/prokka/$var/$gff_def > /mnt/beegfs/ajperez/pangenomes/kp/ccfinder/$var/TSV/$out_cctyper`;
	open (IN6, "/mnt/beegfs/ajperez/pangenomes/kp/ccfinder/$var/TSV/$out_cctyper");
	while (my $linea=<IN6>){
		chomp $linea;
		if ($linea=~/note=CRISPR/){
			$linea=~s/(.+)note=//g; $linea=~s/CRISPR(.+)/CRISPR/g;
			print "$linea ";
		} else {
			$linea=~s/(.+)ID=//g; $linea=~s/;(.+)//g;
			print OUT_LISTA "$linea\n";
			print "$linea ";
		}
	}
	print "\n";
	close IN6;

	###########################################
}

close OUT_LISTA;

`sort lista_id_cas.txt | uniq >lista_id_cas_def.txt`;

open (IN7, "lista_id_cas_def.txt");

open (OUT_CAS, ">list_casgenes_KP.txt");

while (my $line=<IN7>){
	chomp $line;
	#print "$line\n";
	my $busq=`grep -w $line $matrix | cut -f1 | cut -d ":" -f1`;
	print OUT_CAS "$busq";
}

`sort list_casgenes_KP.txt | uniq >list_casgenes_KP_def.txt`;
