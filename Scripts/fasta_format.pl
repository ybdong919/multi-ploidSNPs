
############################# Transform into FASTA format ################################
######   1) add "Loci number" column, delete "CONTIG", "POS" and "REF" columns       #####
######   2) 0 change into -                                                          #####
######   3) transpose between row and column                                         #####
######   4) transform into FASTA format                                              #####
##########################################################################################
#!/usr/bin/perl -w
#$ -S /usr/bin/perl
#$ -N perl_findGT
#$ -j y
#$ -cwd
#$ -R y
use strict;
my $q=0;
my @textfiles= glob 'Output_results/*SNP_hap.txt';
my @fastafiles;
for (@textfiles){
    my $textname = $_;
	$textname =~ s/txt$/fasta/;
    $fastafiles[$q]= $textname;
	$q +=1;   
}

for (my $round=0; $round< $q;$round +=1){

	##### 1) add "Loci number" column, delete "CONTIG", "POS" and "REF" columns
	
	my $in = $textfiles[$round];
	open(NUM,'<', $in) or die;
	open(NUMOUT,'>', "FASTA_step1.txt") or die;

	my @num =<NUM>;                          
	$num[0] = "LOCI_NO\t".$num[0];         ### add "Loci number" column
	my $j = @num;
	for (my $i=1; $i < $j; $i++){
	    $num[$i]= $i."\t".$num[$i];
	}

	my @out1;                            ### delete "CONTIG", "POS" and "REF" columns
	my $m=0;
	foreach (@num){
	    my @row = split /\t/, $_;
	    splice @row, 1, 3;
	    $out1[$m] = join "\t", @row;
	    $m += 1;   
	}
	foreach (@out1){                  ### 0 change into -
	    $_ =~s/0/-/g;
	}

	foreach (@out1){
	    chomp $_;
	    print NUMOUT "$_\n";

	}
	close NUM;
	close NUMOUT;

	##### 2) transposes the data between row and column
	
	open(TRANSPOSE,'>', "FASTA_step2_transpose.txt") or die;
	open(DATA,'<', "FASTA_step1.txt") or die;
	my @data =<DATA>;
	foreach (@data){
	    $_ =~s/\s+$//g;
	}
	my $old_rows = @data;

	my @for_columns = split /\t/, $data[0];
	my $old_columns = @for_columns;
	my $wholestring;
	foreach (@data){
	   chomp $_; 
	   $wholestring .= $_."\t"; 
	}

	my @all_cells = split /\t/, $wholestring;
	my $a=@all_cells;
	my @sheet;
	for (my $il=0; $il<$old_columns; $il++){
	    for (my $nl=$il; $nl<$a; $nl += $old_columns){
		$sheet[$il].= $all_cells[$nl]."\t";     
	    }   
	}

	foreach (@sheet){
	    print TRANSPOSE "$_\n";
	}
	close TRANSPOSE;
	close DATA;

	###### 3) output FASTA format

	open(FASTAFORM,'>', $fastafiles[$round]) or die;
	open(OUME,'<', "FASTA_step2_transpose.txt") or die;
	my @oume =<OUME>;
	shift @oume;

	my @fastaform;
	my $u=0;
	foreach (@oume){
	    chomp $_;
	    my @ouline = split /\t/,$_;
	    $fastaform[$u] = ">".$ouline[0];
	    shift @ouline;
	    my $seq = join "",@ouline;
	    chomp $seq;
	    my $leng = length $seq;
	    my @part;
	    my $p = 0;
	    for (my $i=0; $i<$leng; $i += 60){
		$part[$p]=substr($seq, $i, 60);
		$p += 1;
	    }
	    foreach my $part (@part){
	      $fastaform[$u+1] .= $part."\n"; 
	    }
	    chomp $fastaform[$u+1];
	    $u += 2;
	}

	foreach (@fastaform){
	  chomp $_;
	  print FASTAFORM "$_\n";
	}
	
	close FASTAFORM;
	close OUME;
	unlink qw (FASTA_step1.txt FASTA_step2_transpose.txt);
}
