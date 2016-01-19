#!/bin/bash
#$ -S /bin/bash
#$ -N getting_genotype
#$ -j y
#$ -cwd
#$ -R y
#export PATH=$PATH:/home/AAFC-AAC/dongy/mpp/bowtie2-2.2.3/:/home/AAFC-AAC/dongy/mpp/samtools:/home/AAFC-AAC/dongy/mpp/samtools/bcftools

sta=$(date)
mv ./Input_data/*.fastq ./
#./Scripts/generate_nuclear_contigs.sh
./Scripts/call_nuclear_SNP_DP_ReadsNum.sh ./Threshold_set/Ploid.txt

perl ./Scripts/fasta_format.pl
#rm allcontigs.fa
echo "Calling SNPs is finished"
mv ./*.fastq ./Input_data/

perl ./Scripts/missing_level.pl
fin=$(date)
echo -e "miniAdaptive running is over.\nStart time is $sta.\nFinish time is $fin."
