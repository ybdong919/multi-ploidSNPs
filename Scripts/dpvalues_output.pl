#!/usr/bin/perl
my $input1="All_SNP_Genotypes_DP.txt";
my $input2="part_Clean_SNP_Genotypes.txt";


open (INPUT1,'<',$input1) or die;
my @input1=<INPUT1>;
my $firstline= $input1[0];
chomp $firstline;
my @firstlinecells= split /\t/, $firstline;
my $num= @firstlinecells;
for (my $i=3; $i< $num; $i++){
	chomp $firstlinecells[$i];
	$firstlinecells[$i] = $firstlinecells[$i]."_GT\t".$firstlinecells[$i]."_DP";
}
my $newfirstline= join "\t", @firstlinecells;
shift @input1;

open (INPUT2,'<',$input2) or die;
my @input2=<INPUT2>;
shift @input2;
my @items;
my $u=0;
foreach (@input2){
	chomp $_;
	if (/(\d+__len__\d+\t\d+)\t.+/){
	$items[$u]=$1;
	$u +=1;	
	}
}

my %hash;
foreach (@input1){
	chomp $_;
	if (/(\d+__len__\d+\t\d+)\t.+/) {
		$hash{$1}=$_;
	}
}

my $w=0;
my @filteredlines;
foreach (@items) {
    my $value = $hash{$_};
	if (defined ($value)){
		$filteredlines[$w]=$value;
		$w +=1;	
	}
}

my @newotherlines;
my $x=0;
foreach my $line(@filteredlines) {
    chomp $line;
	$line =~ s/\|/\t/g;
	$line =~ s/\tNA/\t0\t0/g;
	$newotherlines[$x]=$line;
	$x +=1;
}
my $outputfile= "Output_results/Nu_genotypes_DP.txt";
open (OUTPUT, '>', "$outputfile") or die;
print OUTPUT "$newfirstline\n";
foreach (@newotherlines){
	print OUTPUT "$_\n";
}
close OUTPUT;
close INPUT1;
close INPUT2;