#!/usr/bin/perl -w
#$ -S /usr/bin/perl
#$ -N perl_findGT
#$ -j y
#$ -cwd
#$ -R y
use strict;

##### remove SNP sites with missing according to missing-threshold-setting value
open(THRESHOLD,'<',"Threshold_set/Missing_threshold.txt") or die;
my @miss_thres =<THRESHOLD>;
chomp $miss_thres[1];
my $miss_threshold = $miss_thres[1];
#print $miss_threshold;

my $p=0;
my @inputfiles= glob 'Output_results/N*SNP*.txt';
my @outputfiles;
for (@inputfiles){
    my $textname = $_;
	$textname =~ s/^Output_results\//Output_results\/Clean_/;
    $outputfiles[$p]= $textname;
	$p +=1;   
}

for (my $a=0; $a<$p; $a++){
    #print $inputfiles[$a];
	open(INPUTGT,'<', $inputfiles[$a]) or die;
	my @inputgt1 = <INPUTGT>;
	my $title_gt1 = $inputgt1[0];
	chomp $title_gt1;
	
	shift @inputgt1;
    my @inputgt2;
	my $i=0;
	foreach (@inputgt1){
		chomp $_;
		my @splits = split /\t/,$_;
		my $num=@splits;
		#print $num;
		my $mis_num =0;
		for(my $x=3; $x<$num; $x++){      
			if (($splits[$x] eq "NA")|| ($splits[$x] eq "0")) {
				$mis_num +=1;
			}             		
		}
		if ($mis_num <= $miss_threshold) {
			$inputgt2[$i]= $_; 
            $i +=1;			
			}
	}
	open(OUTPUTGT,'>', $outputfiles[$a]) or die;
	print OUTPUTGT $title_gt1;
	print OUTPUTGT "\n";
	foreach (@inputgt2){
	    print OUTPUTGT "$_\n";
	}	
	close INPUTGT;
	close OUTPUTGT;	
}
close THRESHOLD;


