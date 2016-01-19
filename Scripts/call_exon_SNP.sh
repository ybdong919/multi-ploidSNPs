#!/bin/bash
#$ -S /bin/bash
#$ -N getting_genotype
#$ -j y
#$ -cwd
#$ -R y
#export PATH=$PATH:/home/AAFC-AAC/dongy/mpp/bowtie2-2.2.3/:/home/AAFC-AAC/dongy/mpp/samtools:/home/AAFC-AAC/dongy/mpp/samtools/bcftools
### generate hash array: $key is the name of "*R1_001.fastq" (equal to below $file1); $value is the ploid of sequencing sample ### 
echo "Calling exon SNPs is beginning"

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



bowtie2-build exoncontigs.fa bt2ref
ref=exoncontigs.fa
bwa index -a bwtsw $ref
samtools faidx $ref
java -jar ./Scripts/picard-tools-1.137/picard.jar CreateSequenceDictionary REFERENCE=allcontigs.fa OUTPUT=allcontigs.dict

for file1 in *R1_001.fastq
do

   #echo $file1 
   file2=$(echo ${file1} | sed 's/R1_001/R2_001/') 
   #echo $file2 
 
   bowtie2 -x bt2ref -1 $file1 -2 $file2 -S $file1.sam --no-unal 
   samtools view -Sbt ${ref}.fai ${file1}.sam > ${file1}.bam
   samtools sort ${file1}.bam ${file1}.sorted
   java -jar ./Scripts/picard-tools-1.137/picard.jar AddOrReplaceReadGroups I=${file1}.sorted.bam O=${file1}.sorted.grp.bam LB=whatever PL=illumina PU=whatever SM=${file1}
   java -jar ./Scripts/picard-tools-1.137/picard.jar MarkDuplicates INPUT=${file1}.sorted.grp.bam OUTPUT=${file1}.sorted.grp.md.bam METRICS_FILE=${file1}.duplicates ASSUME_SORTED=TRUE CREATE_INDEX=TRUE
   #echo "${ploids["$file1"]}"
   java -jar ./Scripts/GenomeAnalysisTK.jar -T UnifiedGenotyper -R allcontigs.fa -I ${file1}.sorted.grp.md.bam -ploidy ${ploids["$file1"]} -out_mode EMIT_ALL_CONFIDENT_SITES -o ${file1}.vcf
   
   rm *.bai *.bam *.sam
   echo -e "${file1}.vcf" > dirfile
   perl ./Scripts/screen_sampleSNP_DP.pl
   rm *step1*.vcf *step2*.vcf
    
done
perl ./Scripts/identify_SNP_DP.pl
mv part_Clean_SNP_Genotypes.txt ./Output_results/NE_SNP_genotypes.txt
mv part_Clean_SNP_Homo_Genotypes.txt ./Output_results/NE_homo_SNP_genotypes.txt
mv part_Clean_SNP_Homo_haplotype.txt ./Output_results/NE_homo_SNP_hap.txt
rm All_SNP_Genotypes.txt All_SNP_hap.txt
rm GT_data_without_duplication.txt 
rm *.vcf *.bt2 *.fai *name* dirfile *align_contigs All_GTs.txt
rm pre_Clean_SNP_hap.txt


