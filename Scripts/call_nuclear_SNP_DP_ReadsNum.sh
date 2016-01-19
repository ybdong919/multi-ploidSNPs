#!/bin/bash

echo "Calling nuclear SNPs is beginning"

### generate hash array: $key is the name of "*R1_001.fastq" (equal to below $file1); $value is the ploid of sequencing sample ### 
declare -A ploids
filename="$1"
while read -r line
do
   name=$line
   key=${name%%/*}
   echo "$key"
   value=${name##*/}
   echo "$value"
   ploids+=(["$key"]="$value")      
done < "$filename"

bowtie2-build allcontigs.fa bt2ref
ref=allcontigs.fa
bwa index -a bwtsw $ref
samtools faidx $ref
java -jar ./Scripts/picard-tools-1.137/picard.jar CreateSequenceDictionary REFERENCE=allcontigs.fa OUTPUT=allcontigs.dict

for file1 in *R1_001.fastq
do
   file2=$(echo ${file1} | sed 's/R1_001/R2_001/') 
   bowtie2 -x bt2ref -1 $file1 -2 $file2 -S $file1.sam --no-unal
   samtools view -Sbt ${ref}.fai ${file1}.sam > ${file1}.bam
   samtools sort ${file1}.bam ${file1}.sorted
   
   #samtools index ${file1}.sorted.bam
   java -jar ./Scripts/picard-tools-1.137/picard.jar AddOrReplaceReadGroups I=${file1}.sorted.bam O=${file1}.sorted.grp.bam LB=whatever PL=illumina PU=whatever SM=${file1}
   java -jar ./Scripts/picard-tools-1.137/picard.jar MarkDuplicates INPUT=${file1}.sorted.grp.bam OUTPUT=${file1}.sorted.grp.md.bam METRICS_FILE=${file1}.duplicates ASSUME_SORTED=TRUE CREATE_INDEX=TRUE
   #echo "${ploids["$file1"]}"
   java -jar ./Scripts/GenomeAnalysisTK.jar -T UnifiedGenotyper -R allcontigs.fa -I ${file1}.sorted.grp.md.bam -ploidy ${ploids["$file1"]} -stand_call_conf 10 -stand_emit_conf 10 -out_mode EMIT_ALL_CONFIDENT_SITES -o ${file1}.vcf
      
   samtools view -F 4 ${file1}.sorted.bam | cut -f3 | uniq -c > ${file1}.readsnum   ## get reads count of every samples matching each contig
   echo -e "${file1}.readsnum" >> readsnum_namefiles                                ## get reads count of every samples matching each contig
   
   echo -e "${file1}.vcf" > dirfile
   perl ./Scripts/screen_sampleSNP_DP.pl
     
done
perl ./Scripts/contig_reads_num.pl                   ## get reads count of every samples matching each contig
perl ./Scripts/identify_SNP_DP.pl
perl ./Scripts/dpvalues_output.pl                  ## get DP value of each SNP, each sample
mv part_Clean_SNP_Genotypes.txt ./Output_results/Nu_SNP_genotypes.txt
mv part_Clean_SNP_Homo_Genotypes.txt ./Output_results/Nu_homo_SNP_genotypes.txt
mv part_Clean_SNP_Homo_haplotype.txt ./Output_results/Nu_homo_SNP_hap.txt

rm *.bai *.bam *.sam
rm *step1*.vcf *step2*.vcf 
rm All_SNP_Genotypes.txt All_SNP_Genotypes_DP.txt 
rm GT_data_without_duplication.txt
rm *.vcf *.bt2 *.fai *name* dirfile *align_contigs All_GTs.txt
#rm pre_Clean_SNP_hap.txt readsnum_namefiles
rm *.readsnum 
rm *.duplicates *.idx allcontigs.dict
