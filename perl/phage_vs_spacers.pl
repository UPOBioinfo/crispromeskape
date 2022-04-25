#!/usr/bin/perl
use strict;

my $MD = $ARGV[0] || "phages_ncbi.fasta"

open (IN, "$MD");

my $id; 
my $seq;
my %hash; 
my $file;

while (my $linea=<IN>){
	chomp $linea;
	if ($linea=~/^>(.+)/){
		$linea=~s/>//g;
		$linea=~s/\s(.+)//g;
		$id=$linea;
		$seq="";
		$file=$id.".fasta";
		open (OUT, ">$file");
		print OUT ">$id\n";
	} else {
		print OUT "$linea\n";
		$seq.=$linea;
		$hash{$id}=$seq;
	}
}

close IN;

foreach my $var (keys %hash){
	$file=$var.".fasta";
	my $out=$var.".out";
	my $outfil=$var."_filtrado.out";
	#open (OUT, ">/home/arubval/Proyecto_Ministerio/spacers_vs_db/phage_vs_spab/phage/$file");
	#print OUT ">$var\n$hash{$var}\n";
	#close OUT;
	`makeblastdb -in $file -dbtype nucl`;
	`blastn -task blastn-short -query spacers_ab_filtrate.fasta -db $file -outfmt '6 qseqid sseqid pident qcovs qcovhsp length qlen slen evalue qstart qend sstart send' -max_target_seqs 100 -num_threads 40 | awk '\$3 >= 95 && \$4 >= 95' > $out`;
	open (FILE, "/home/arubval/Proyecto_Ministerio/spacers_vs_db/phage_vs_spab/Resultados/$out");
	open (FILTRADO, ">/home/arubval/Proyecto_Ministerio/spacers_vs_db/phage_vs_spab/Resultados/$outfil");
	while (my $line=<FILE>){
		chomp $line;
		#print "$line\n";
		my @list=split (/\t/, $line);
		if ($list[3]>=95 and $list[2]>=95 and $list[4]>=95){
			print FILTRADO "$list[0]\t$list[1]\t$list[2]\t$list[4]\n";
		}
	}
	close FILE;
}
