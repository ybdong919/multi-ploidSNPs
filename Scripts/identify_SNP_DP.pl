#!/usr/bin/perl -w
use strict;

############################# step1 find GT in all examples ####################################################################

open(CONTIGS,'<','align_contigs') or die;
my @contigs = <CONTIGS>;
close CONTIGS;

my %list;
foreach (@contigs){
    chomp $_;
    if (/(.+\d+)\t(\w+)/) {
        $list{$1} = $2;         ### using hash key's unique trait to delete duplications 
    }    
}


##############################################################

open(STEP1,'<', 'names_of_step4') or die;
my @namefile_step1 = <STEP1>;
#print @namefile_step1;
close STEP1;

my @list =keys %list;

foreach (@namefile_step1){
    chomp $_;
    open(FILE,'<', $_) or die;
    my @input =<FILE>;
    close FILE;
    
    my %data;
    foreach (@input){
      chomp $_;
      if (/^(.+\t\d+)\t(.+)/) {
         my $keys = $1;
        $data{$keys}=$2; 
      }  
    }

    foreach (@list){
        chomp $_;
        my $var = $data{$_};
        if (defined($var)) {
            $list{$_} .= "\t".$data{$_};
        }else {
            $list{$_} .="\tNA\t0";
        }
    } 
}


###################################### delete Quality columns ####################
######################################  ################
my @k2 = keys %list;
foreach (@k2){
    my @line = split /\t/, $list{$_};
    #print $line[2]."\n";
    my $stri;
    foreach (@line){
        if (/^[^0-9]/){
          $stri .= "$_\t"; 
        }    
    }
    $list{$_} = $stri;   
    #print $stri."\n";
    #print $av_qt."\n";
}

##############################################################################################
#########################################   delete consensus sites###################

my $j = @namefile_step1;
my @k = keys %list;
   
foreach (@k){
     my @line = split /\t/, $list{$_};
     my $con;
     
     shift @line;
     foreach (@line){
		my @gtdpcell=split /\|/,$_;        ### seperate genotype value and DP value. 
		my $gtvalue=$gtdpcell[0];
        if ($gtvalue ne 'NA') {
            $con = $gtvalue;
            last;
        }        
     }
     #print @line;
     #print "\n";
	 
	 ####################### add at 10.8.2015
	 chomp $con;
	 my @con_base=split//,$con;
	 my $con_base1=$con_base[0];
	 my $con_base_homo=1;
	 foreach (@con_base){
	    if ($con_base1 ne $_){
		    $con_base_homo=0;
			last;
		}
	 }
	 
	 if ($con_base_homo==1){
	    my $num =0;
		foreach (@line){
			my @gtdpcell=split /\|/,$_;       ### seperate genotype value and DP value.
			my $gtvalue=$gtdpcell[0];
			my @gtvalue_base=split//,$gtvalue;
			my $gtvalue_base1=$gtvalue_base[0];
			my $gtvalue_homo=1;
			foreach (@gtvalue_base){
			    if ($gtvalue_base1 ne $_){
				   $gtvalue_homo=0;
				   last;
				}
			}
			if (($gtvalue_homo==1 && $gtvalue_base1 eq $con_base1 )|| $gtvalue eq "NA") {
				$num += 1;
			}   
		}
		if ($num == $j) {
			delete $list{$_};
		}   	 
	 } 
	 
	 if ($con_base_homo==0){
	    my $num =0;
		foreach (@line){
			my @gtdpcell=split /\|/,$_;       ### seperate genotype value and DP value.
			my $gtvalue=$gtdpcell[0];
			if ($gtvalue eq $con || $gtvalue eq "NA") {
				$num += 1;
			}   
		}
		if ($num == $j) {
			delete $list{$_};
		}   
	 }
     
}

###############################################################
my @titles;
my $x = 0;
foreach (@namefile_step1){
    if (/step4_(.+S\d+)/) {
        $titles[$x] = $1;
        $x += 1;
    }
}

open(OUTPUT,'>', 'All_GTs.txt') or die;

print OUTPUT "CHROM\tPOS\tREF\t";
foreach (@titles){
    print OUTPUT "$_\t";    
}
print OUTPUT "\n";
my $k;
my $v;
while (($k,$v)= each %list) {      
    print OUTPUT "$k\t$v\n";
}
close OUTPUT;

############################# step2 sort GT ####################################################################

my $input_gt = "All_GTs.txt";
open(INPUTGT,'<', $input_gt) or die;
my @inputgt = <INPUTGT>;
my $title_gt = $inputgt[0];
chomp $title_gt;
shift @inputgt;

my %sorthash;
foreach (@inputgt){
    if (/(\d+)__len__\d+\t(\d+)\t.+/) {       
        $sorthash{$1}{$2}=$_;             #definate a multi-dimension hash
    }       
}

my $i=0;
my @value;
foreach my $key1 (sort {$a<=>$b} keys %sorthash){                         # using Multi-dimension hash for sorting 
    foreach my $key2 (sort {$a <=> $b} keys %{$sorthash{$key1}}){
        $value[$i]=$sorthash{$key1}{$key2};
        $i +=1;
    }
}

open(OUTPUTGT,'>',"All_SNP_Genotypes_DP.txt") or die;
print OUTPUTGT $title_gt;
print OUTPUTGT "\n";
foreach (@value){
   chomp $_;
   print OUTPUTGT "$_\n";
}

close INPUTGT;
close OUTPUTGT;

############################# Plus:  Seperate genotype value and DP value into two files ###################
my @inputcomble = ("All_SNP_Genotypes_DP.txt");
my @outputunique = ("All_SNP_Genotypes.txt");
my $yy=0;
foreach (@inputcomble){
	
    open(HAHAINPUT,'<', $_) or die;         
    my @samples =<HAHAINPUT>;
	my $firstline = $samples[0];
    chomp $firstline;
	shift @samples;
	my @genotypes;
	#my @dpvalues;
	my $xx=0;
	foreach (@samples){
	    chomp $_;
		my @cells= split /\t/, $_;
		my @firstthree= splice (@cells, 0,3);
		my $numcells = @cells;
		my @gtval;
		my @dpval;
		for (my $i=0; $i< $numcells; $i++) {
			chomp $cells[$i];
			my @mint= split /\|/, $cells[$i];
			$gtval[$i]=$mint[0];
			#$dpval[$i]=$mint[1];
		
		}
		$genotypes[$xx]=join "\t", @firstthree, @gtval;
		#$dpvalues[$xx]=join "\t", @firstthree, @dpval; 
		$xx +=1;
	}
	
	open (OUTHAHA,'>', "Only_SNP_Genotypes.txt") or die;
	print OUTHAHA "$firstline\n";
	foreach (@genotypes){
		chomp $_;
		print OUTHAHA "$_\n";
	}
	close HAHAINPUT;
	close OUTHAHA;
		
	my $oldnam = "Only_SNP_Genotypes.txt";
	rename $oldnam => $outputunique[$yy];
    $yy +=1;
}


############################# step3 remove duplicated GT-data ####################################################################

#open(THRESHOLD,'<',"Threshold_set/Missing_threshold.txt") or die;
#my @miss_thres =<THRESHOLD>;
#chomp $miss_thres[1];
#my $miss_threshold = $miss_thres[1];
#print $miss_threshold;

my $input_gt1 = "All_SNP_Genotypes.txt";
open(INPUTGT,'<', $input_gt1) or die;
my @inputgt1 = <INPUTGT>;
my $title_gt1 = $inputgt1[0];
chomp $title_gt1;
shift @inputgt1;

my %del_dup;
foreach (@inputgt1){
    if (/(.+\d+\t\d+)\t(.+)/) {
        $del_dup{$2}=$1;
    }   
}

#my @keys = keys %del_dup;
#foreach (@keys){
#    chomp $_;
#    my @splits_keys = split /\t/,$_;
#    my $mis_num =0;
#    foreach my $unit(@splits_keys){      
#        if ($unit eq "NA") {
#            $mis_num +=1;
#        }          
#    }
#    if ($mis_num > $miss_threshold) {
#        delete $del_dup{$_};    
#        }
#}

open(OUTPUTGT,'>',"GT_data_without_duplication.txt") or die;
print OUTPUTGT $title_gt1;
print OUTPUTGT "\n";
my $k1;
my $v1;
while (($k1,$v1) = each %del_dup) {
    print OUTPUTGT "$v1\t$k1\n";
}
close INPUTGT;
close OUTPUTGT;
close THRESHOLD;

############################# step4 remove GT-data if SNP position <=0 ############################################

#open(THRESHOLD,'<',"Threshold_set/SNP_position_threshold.txt") or die;
#my @SNP_thres =<THRESHOLD>;
#chomp $SNP_thres[1];
#my $SNP_threshold = $SNP_thres[1];
#print $SNP_threshold;

my $input_gt2 = "GT_data_without_duplication.txt";
open(INPUTGT,'<', $input_gt2) or die;
my @inputgt2 = <INPUTGT>;
my $title_gt2 = $inputgt2[0];
chomp $title_gt2;
shift @inputgt2;

my %mulhash;
my %outhash;
foreach (@inputgt2){
    if (/(\d+)__len__(\d+)\t(\d+)\t.+/) {       
        $mulhash{$1}{$3}=$2;             #definate a multi-dimension hash
        $outhash{$1}{$3}=$_;
    }       
}

my $ib=0;
my @value1;
foreach my $key1 (sort {$a<=>$b} keys %mulhash){               # using Multi-dimension hash for sorting and delete SNPs located at both ends of each contig
    
    foreach my $key2 (sort {$a <=> $b} keys %{$mulhash{$key1}}){
        
            if (($key2 > 0)&&($key2 < ($mulhash{$key1}{$key2} - 0))) {
                 $value1[$ib]=$outhash{$key1}{$key2};
                 $ib +=1;
            }            
    }
}

open(OUTPUTGT,'>',"part_Clean_SNP_Genotypes.txt") or die;
print OUTPUTGT $title_gt2;
print OUTPUTGT "\n";
foreach (@value1){
   chomp $_;
   print OUTPUTGT "$_\n";
}

close INPUTGT;
close OUTPUTGT;
close THRESHOLD;

############################ step5 homologous-ployplody(>=3) genotype ##################################
#################################   diplody genotype may be heterozygous,while ployplody must be homologous, in order to solve homoeology problem. 
open (FIVEINPUT, '<', "part_Clean_SNP_Genotypes.txt") or die;
open (FIVEOUTPUT,'>',"part_Clean_SNP_Homo_Genotypes.txt")or die;
my @fiveinput=<FIVEINPUT>;
my $five_titleline=$fiveinput[0];
chomp $five_titleline;
shift @fiveinput;
my $fiveseries=0;
my @homoline;
foreach my $eachline(@fiveinput){
     chomp $eachline;
     my @line=split/\t/,$eachline;
	 my $num=@line;
	 my $hetegenotype=0;
	 LINE:for (my $i=3;$i<$num;$i++){
			chomp $line[$i];
			if ($line[$i]=~/^NA/){       # if 'NA', jump to next loop
				next;
			}
			my @bases=split//,$line[$i];
			my $basenumb=@bases;
			if ($basenumb > 2){           # if not 'NA' and ploid >2
			    my $base1=$bases[0];		 
				foreach (@bases){
					if ($base1 ne $_){
						$hetegenotype=1;
						last LINE;			    
					}		        
				}	  
			}
			   
	 }
     if ($hetegenotype==0){
	    $homoline[$fiveseries]=$eachline;
        $fiveseries += 1;		
	 }
}

print FIVEOUTPUT "$five_titleline\n";
foreach (@homoline){
    chomp $_;
	print FIVEOUTPUT "$_\n";
}
close FIVEINPUT;
close FIVEOUTPUT;

############ step6 polyploidy-homologous haplotype ########################################

open (SIXINPUT, '<', "part_Clean_SNP_Homo_Genotypes.txt") or die;
open (SIXOUTPUT,'>',"part_Clean_SNP_Homo_haplotype.txt")or die;
my @sixinput=<SIXINPUT>;
my $six_titleline = $sixinput[0];
chomp $six_titleline;
shift @sixinput;
foreach my $eachline (@sixinput){
    chomp $eachline;
	my @cells=split/\t/, $eachline;
	my $num = @cells;
    for (my $i=2;$i< $num; $i++){
	    chomp $cells[$i];
	    my @bases=split//, $cells[$i];
		if ($bases[0] eq "N"){
		   $cells[$i]=0;
		} else {
		$cells[$i]=$bases[0];
		}
	}	
    $eachline = join "\t", @cells;
}

print SIXOUTPUT "$six_titleline\n";
foreach (@sixinput){
    chomp $_;
	print SIXOUTPUT "$_\n";
}
close SIXINPUT;
close SIXOUTPUT;







