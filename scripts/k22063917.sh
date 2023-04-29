#!/usr/bin/env bash
mkdir AdvancedBioinformatics
cd AdvancedBioinformatics/
mkdir data meta results logs
ls -lf
ls -lls

mkdir untrimmed_fastq trimmed_fastq
cd untrimmed_fastq/
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
mv *fastq.qz ~/AdvancedBioinformatics/data/untrimmed_fastq

bgzip -d NGS0001.R1.fastq.qz > NGS0001.R1.fastq
bgzip -d NGS0001.R2.fastq.qz > NGS0001.R2.fastq



cd ..
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

fastqc -t 4 *.fastq

mkdir ~/AdvancedBioinformatics/results/fastqc_untrimmed_reads

mv *fastqc* ~/AdvancedBioinformatics/results/fastqc_untrimmed_reads/

cd ../../data/untrimmed_fastq/



trimmomatic PE -threads 4 -phred33 ~/AdvancedBioinformatics/data/untrimmed_fastq/NGS0001.R1.fastq ~/AdvancedBioinformatics/data/untrimmed_fastq/NGS0001.R2.fastq -baseout ~/AdvancedBioinformatics/data/trimmed_fastq/NGS0001_trimmed_R ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10   TRAILING:25 MINLEN:50

cd ../trimmed_fastq/

fastqc -t 4 *P

cd ../../results/

mkdir fastqc_trimmed_reads

cd ../data/trimmed_fastq/
mv *fastqc* ~/AdvancedBioinformatics/results/fastqc_trimmed_reads/

cd ../../data/
mkdir reference
mv ~/AdvancedBioinformatics/data/hg19.fa.gz ~/AdvancedBioinformatics/data/reference/

bwa index ~/AdvancedBioinformatics/data/reference/hg19.fa.gz

mkdir ~/AdvancedBioinformatics/data/aligned_data


bwa mem -t 4 -v 1 -R '@RG\tID:11V6WR.111.D1375ACXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-NGS0001-blood\tDT:2017-02-23\tPU:11V6WR1' -I 250,50  ~/AdvancedBioinformatics/data/reference/hg19.fa.gz ~/AdvancedBioinformatics/data/trimmed_fastq/NGS0001_trimmed_R_1P ~/AdvancedBioinformatics/data/trimmed_fastq/NGS0001_trimmed_R_2P > ~/AdvancedBioinformatics/data/aligned_data/NGS0001.sam

samtools view -h -b NGS0001.sam > NGS0001.bam
samtools sort NGS0001.bam > NGS0001_sorted.bam
samtools index NGS0001_sorted.bam
ls

picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt

samtools index NGS0001_sorted_marked.bam
samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam

samtools index NGS0001_sorted_filtered.bam

samtools flagstat NGS0001_sorted_filtered.bam
samtools idxstats NGS0001_sorted_filtered.bam
samtools stats NGS0001_sorted_filtered.bam

picard CollectInsertSizeMetrics I=NGS0001_sorted_filtered.bam O=insert_size_metrics.txt H=insert_size_histogram.pdf M=0.5

ls

cp annotation.bed ~/AdvancedBioinformatics/data/aligned_data/annotation.bed2
bedtools coverage -a  NGS0001_sorted_filtered.bam -b annotation.bed2


zcat ~/AdvancedBioinformatics/data/reference/hg19.fa.gz > ~/AdvancedBioinformatics/data/reference/hg19.fa
samtools faidx ~/AdvancedBioinformatics/data/reference/hg19.fa

freebayes --bam ~/AdvancedBioinformatics/data/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference ~/AdvancedBioinformatics/data/reference/hg19.fa --vcf ~/AdvancedBioinformatics/results/NGS0001.vcf

bgzip ~/AdvancedBioinformatics/results/NGS0001.vcf
tabix -p vcf ~/AdvancedBioinformatics/results/NGS0001.vcf.gz

vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ~/AdvancedBioinformatics/results/NGS0001.vcf.gz > ~/AdvancedBioinformatics/results/NGS0001_filtered.vcf

bedtools intersect -header -wa -a ~/AdvancedBioinformatics/results/NGS0001_filtered.vcf -b ../data/annotation.bed > ~/AdvancedBioinformatics/results/NGS0001_filtered_hg19.vcf

bgzip ~/AdvancedBioinformatics/results/NGS0001_filtered_hg19.vcf

tabix -p vcf ~/AdvancedBioinformatics/results/NGS0001_filtered_hg19.vcf.gz

#compresses some file to make space in the virtual machine
bgzip ~/AdvancedBioinformatics/data/untrimmed_fastq/NGS0001.R1.fastq
bgzip ~/AdvancedBioinformatics/data/untrimmed_fastq/NGS0001.R2.fastq
bgzip ~/AdvancedBioinformatics/data/trimmed_fastq/NGS0001_trimmed_R_1P
bgzip ~/AdvancedBioinformatics/data/trimmed_fastq/NGS0001_trimmed_R_1U
bgzip ~/AdvancedBioinformatics/data/trimmed_fastq/NGS0001_trimmed_R_2U
bgzip ~/AdvancedBioinformatics/data/trimmed_fastq/NGS0001_trimmed_R_2P
bgzip ~/AdvancedBioinformatics/data/aligned_data/NGS0001.bam
bgzip ~/AdvancedBioinformatics/data/aligned_data/NGS0001.sam
bgzip ~/AdvancedBioinformatics/data/aligned_data/NGS0001_sorted.bam
bgzip ~/AdvancedBioinformatics/data/aligned_data/NGS0001_sorted.bam.bai

#unpack annovar
tar -zxvf annovar.latest.tar.gz

#download Annovar db
cd annovar
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/

./convert2annovar.pl -format vcf4 ~/AdvancedBioinformatics/results/NGS0001_filtered_hg19.vcf.gz > ~/AdvancedBioinformatics/results/NGS0001_filtered_hg19.avinput

./table_annovar.pl ~/AdvancedBioinformatics/results/NGS0001_filtered_hg19.avinput humandb/ -buildver hg19 -out ~/AdvancedBioinformatics/results/NGS0001_filtered_hg19 -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operati
