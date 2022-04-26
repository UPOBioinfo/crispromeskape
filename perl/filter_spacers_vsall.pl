#!/usr/bin/perl
use strict;

open (IN, "spacers_vs_virus_sa.txt");

my $ID = $ARGV[0] || "95";
my $COV = $ARGV[1] || "95";

my %hash;

while (my $line=<IN>){
	chomp $line;
	#print "$line\n";
	my @list=split (/\t/, $line);
	my $length=$list[5];
	my $slen=$list[7];
	if ($list[3]>=$ID and $list[2]>=$COV){
		$hash{$list[0]}{"VIRUS"}=$list[1];
		#print "$list[0]\t$list[2]\t$cov\n";
	} else {
		$hash{$list[0]}{"VIRUS"}="";
		next;
	}
}

close IN;

open (IN2, "spacers_vs_plasmid_sa.txt");

while (my $line=<IN2>){
	chomp $line;
	my @list=split (/\t/, $line);
	my $length=$list[5];
	my $slen=$list[6];
	if ($list[3]>=$ID and $list[2]>=$COV){
		$hash{$list[0]}{"PLASMID"}=$list[1];
		#print "$list[0]\t$list[2]\t$cov\n";
	} else {
		$hash{$list[0]}{"PLASMID"}="";
		next;
	}
}

close IN2;

open (IN3, "spacers_vs_ICE_sa.txt");

while (my $line=<IN3>){
	chomp $line;
	my @list=split (/\t/, $line);
	my $length=$list[5];
	my $slen=$list[7];
	if ($list[3]>=$ID and $list[2]>=$COV){
		$hash{$list[0]}{"ICE"}=$list[1];
		#print "$list[0]\t$list[2]\t$cov\n";
	} else {
		$hash{$list[0]}{"ICE"}="";
		next;
	}
}

close IN3;

open (IN4, "spacers_vs_CIME_sa.txt");

while (my $line=<IN4>){
	chomp $line;
	my @list=split (/\t/, $line);
	my $length=$list[5];
	my $slen=$list[7];
	if ($list[3]>=$ID and $list[2]>=$COV){
		$hash{$list[0]}{"CIME"}=$list[1];
		#print "$list[0]\t$list[2]\t$cov\n";
	} else {
		$hash{$list[0]}{"CIME"}="";
		next;
	}
}

close IN4;

open (IN5, "spacers_vs_IME_sa.txt");

while (my $line=<IN5>){
	chomp $line;
	my @list=split (/\t/, $line);
	my $length=$list[5];
	my $slen=$list[7];
	if ($list[3]>=$ID and $list[2]>=$COV){
		$hash{$list[0]}{"IME"}=$list[1];
		#print "$list[0]\t$list[2]\t$cov\n";
	} else {
		$hash{$list[0]}{"IME"}="";
		next;
	}
}

close IN5;

open (IN6, "spacers_vs_AICE_sa.txt");

while (my $line=<IN6>){
	chomp $line;
	my @list=split (/\t/, $line);
	my $length=$list[5];
	my $slen=$list[7];
	if ($list[3]>=$ID and $list[2]>=$COV){
		$hash{$list[0]}{"AICE"}=$list[1];
		#print "$list[0]\t$list[2]\t$cov\n";
	} else {
		$hash{$list[0]}{"AICE"}="";
		next;
	}
}

close IN6;

open (IN7, "spacers_vs_T4SS_sa.txt");

while (my $line=<IN7>){
	chomp $line;
	my @list=split (/\t/, $line);
	my $length=$list[5];
	my $slen=$list[7];
	if ($list[3]>=$ID and $list[2]>=$COV){
		$hash{$list[0]}{"T4SS"}=$list[1];
		#print "$list[0]\t$list[2]\t$cov\n";
	} else {
		$hash{$list[0]}{"T4SS"}="";
		next;
	}
}

close IN7;

open (ID, "id_no_duplicados_sa.txt");

my @id;

while (my $line2=<ID>){
	chomp $line2;
	push @id, $line2;
}

open (CL, "clustered_proteins");

my %cluster;
my %cluster_strain;

while (my $line6=<CL>){
	chomp $line6;
	my $referencia=$line6;
	$referencia=~s/\t(.+)//g;
	$referencia=~s/(.+) //g;
	$line6=~s/(.+): //g;
	$cluster_strain{$referencia}=$line6;
	#print "$referencia\n";
	my @secuencia=split (/\t/, $line6);
	foreach my $var5 (@secuencia){
		$cluster{$var5}=$referencia;
	}
}

open (PAN, "sa_sp_filtrado.blastshort"); #ab_sp_filtrado.blastshort_v2_def

my %pan;

while (my $line3=<PAN>){
	chomp $line3;
	my @lista=split (/\t/, $line3);
	$pan{$lista[0]}=$lista[2];
}

my %pangedef;
my %straindef;

foreach my $var2 (keys %pan){
	#print "$pan{$var2}\n";
	my @pangenome=split (/,/, $pan{$var2});
	#my @ref;
	my %pangev2;
	foreach my $var3 (@pangenome){
		my $busq= $cluster{$var3};
		#if ($busq ~~ @ref) {
		#	next;
		#} else {
		#	push @ref, $busq;
		#}
		$pangev2{$busq}=1;
	}
	my @ref=sort keys %pangev2;
	#print "$var2\t@ref\n";
	foreach my $var4 (@ref){
		chomp $var4;
		$pangedef{$var2}.=" ".$var4;
		$pangedef{$var2}=~s/\t//g;
		#my @rsa_strain;
		my @lista4=split (/\t/, $cluster_strain{$var4});
		my %pangev1;
		foreach my $var6 (@lista4){
			$var6=~s/_(.+)//g;
			$pangev1{$var6}=1;
			#if ($var6 ~~ @rsa_strain) {
			#	next;
			#} else {
			#	push @rsa_strain, $var6;
			#}
		}
		my @rsa_strain=sort keys %pangev1;
		foreach my $var7 (@rsa_strain){
			$straindef{$var2}.=" ".$var7;
			$straindef{$var2}=~s/\t//g;
		}
	}
}

open (ASOC, "aso_duplicados_sa.txt");

my %aso;

while (my $line4=<ASOC>){
	chomp $line4;
	my @lista2=split (/\t/, $line4);
	$aso{$lista2[0]}=$lista2[1];
	#print "$lista2[1]\n";
}

open (ST, "index_sa_sp_ccfinder_uniq.tsv");

my %st;

while (my $line5=<ST>){
	chomp $line5;
	my @lista3=split (/\t/, $line5);
	$st{$lista3[0]}=$lista3[1];
	#print "$lista2[1]\n";
}

print "Strain\tPhage\tPlasmid\tICE\tCIME\tIME\tAICE\tT4SS\tReferences\tStrain_proto\tStrain_spacer\n";

foreach my $var (@id){
	chomp $var;
	my @strain= split (/\s/, $st{$var});
	$pangedef{$var}=~s/\s//;
	$straindef{$var}=~s/\s//;
	print "$var\t$hash{$var}{VIRUS}\t$hash{$var}{PLASMID}\t$hash{$var}{ICE}\t$hash{$var}{CIME}\t$hash{$var}{IME}\t$hash{$var}{AICE}\t$hash{$var}{T4SS}\t$pangedef{$var}\t$straindef{$var}\t";
	my @listado= split (/\s/, $aso{$var});
	foreach my $ident (@listado){
		my @id_2=split (/\s/, $st{$ident});
		foreach my $ident2 (@id_2){
			if ($ident2 ~~ @strain) {
				next;
			} else {
				push @strain, $ident2;
			}
		}
	}
	print "@strain\n";
}
