#!/usr/bin/perl

open (ALLCONTIGS,'<',"allcontigs.fa") or die;
my @allcontigs=<ALLCONTIGS>;
my $num= @allcontigs;
my @contigtitles;
my $x=0;
for (my $i=0; $i<$num; $i +=2){
	chomp $allcontigs[$i];
	$contigtitles[$x]= substr($allcontigs[$i], 1);
	$contigtitles[$x]=~s/\s+$//;
	$x +=1;
}

open (INPUT,'<',"readsnum_namefiles") or die;
my @filenames=<INPUT>;
my %hash;
foreach (@contigtitles) {
	chomp $_;
	#print "$_\n";
	$hash{$_}=$_;
}
my @conkey= keys %hash;

foreach (@filenames){
	chomp $_;
	open (FILE,'<',$_) or die;
	my @lines=<FILE>;
	my %contigs;
	foreach (@lines){
		chomp $_;	
		my @cells = split /\s+/, $_;
		$contigs{$cells[2]}=$cells[1];
	}
	
    foreach (@conkey){
		chomp $_;
		my $var= $contigs{$_};
		if (defined($var)){
			$hash{$_} .= "\t".$contigs{$_};		
		}else {
			$hash{$_} .="\t0";		
		}
	}
}

my @titleline;
$titleline[0]="Contigs";
my $j=1;
foreach (@filenames){
    my @parts= split /_L001_/, $_;
	$titleline[$j]= $parts[0];
	$j +=1;
}

my $output= "Output_results/Readsnum_of_contig.txt";
open (OUTPUT,'>',$output) or die;
foreach (@titleline){
	print OUTPUT "$_\t";
}
print OUTPUT "\n";
foreach (@contigtitles){
	chomp $_;
	
    my $val= $hash{$_};
	print OUTPUT "$val\n";
}

close OUTPUT;
close ALLCONTIGS;
close FILE;




