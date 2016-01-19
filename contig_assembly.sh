#!/bin/bash
# De novo assemble contigs. 
# If no commandline parameter, all fastq sequencing files will be used;
# If using the file name including list of some fastq sequencing-file 
# names as commandline parameter, only these sequencing files will be 
# used for denovo assembly contigs. 
sta=$(date)
ls ./Input_data/*fastq > datafilenames.txt
selectedfile=${1:-'datafilenames.txt'}
#echo $selectedfile

samnum=$(cat $selectedfile|wc -l)
samnum=$(printf "%.0f" $(echo "scale=2;$samnum*0.6"|bc))
#echo $samnum

while read -r line
do  
   mv $line ./    
done < "$selectedfile"

for file1 in *fastq
do
  fastx_collapser -Q33 -i $file1 -o $file1.fx 
done

ls -L *.fx > fx-list.txt
minia fx-list.txt 100 $samnum 300000000 fx                  
mv fx.contigs.fa allcontigs.fa
mv *fastq ./Input_data/
cp allcontigs.fa ./Output_results/Nu_contigs.fasta
rm *.fx fx* 
rm datafilenames.txt
echo "Generating contigs is finished"
fin=$(date)
echo -e "Start time is $sta.\nFinish time is $fin."