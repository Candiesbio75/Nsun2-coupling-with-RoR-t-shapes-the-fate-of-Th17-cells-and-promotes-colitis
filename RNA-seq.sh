#!/bin/sh

data=./RNAseq/ywl-Th17-CKO-T-1
rawdata=ywl-Th17-CKO-T-1
hisat2_index=./Mouse/ENSEMBLE_68/index/genome.fa
gtf=./Mouse/ENSEMBLE_68/GTF/Mus_musculus.GRCm38.68.chr.gtf
single=./Mouse/ENSEMBLE_68/GeneList/single-transcript-chr

thread=6

####################################################
cd $data
mkdir clean_data
cd clean_data
mkdir fastqc

##
fastq1=${rawdata}/*1.fq.gz
fastq2=${rawdata}/*2.fq.gz

fastqc ${fastq1} -t $thread -o ./fastqc
fastqc ${fastq2} -t $thread -o ./fastqc

cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG  -o cutadapt_fastq1.fastq.gz -p cutadapt_fastq2.fastq.gz ${fastq1} ${fastq2}
fastqc cutadapt_fastq1.fastq.gz -t $thread -o ./fastqc
fastqc cutadapt_fastq2.fastq.gz -t $thread -o ./fastqc

java -Xmx4g -jar /software/biosoft/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 cutadapt_fastq1.fastq.gz cutadapt_fastq2.fastq.gz ./trim_qua_1.fastq.gz ./unpaired_1.fastq.gz ./trim_qua_2.fastq.gz ./unpaired_2.fastq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:35
rm  cutadapt_fastq1.fastq.gz cutadapt_fastq2.fastq.gz
fastqc ./trim_qua_1.fastq.gz -t $thread -o ./fastqc
fastqc ./trim_qua_2.fastq.gz -t $thread -o ./fastqc

#######################mapping
cd ${data}
#mkdir hisat2
#hisat2 -p 12 -N 1 --dta -x $hisat2_index -1 ${data}/clean_data/trim_qua_1.fastq.gz -2 ${data}/clean_data/trim_qua_2.fastq.gz -S ./hisat2/hisat2.sam --un-conc ./hisat2

cd hisat2
gzip un-conc-mate.1 un-conc-mate.2
echo "raw read" > readsCount
zcat ${fastq1}|awk 'NR%4==2'|wc -l >> readsCount
echo "trimmed read" >> readsCount
zcat ${data}/clean_data/trim_qua_1.fastq.gz|awk 'NR%4==2'|wc -l >> readsCount
echo "reads with no -q:" > readsCount
cat hisat2.sam|grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> readsCount

samtools view -hSb hisat2.sam > hisat2.bam
samtools view hisat2.bam -q 20 -h > hisat2_20.sam
cat hisat2_20.sam|grep -E "^@|NH:i:1$|NH:i:1[^0-9]"|awk '$1~/^@/ || $1~/^[A-Z]00/'> uniqmap.sam
echo "reads with -q:" >> readsCount
cat hisat2_20.sam|grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> readsCount
rm hisat2_20.sam
rm hisat2.sam

echo "reads with uniqmap:" >> readsCount
cat uniqmap.sam|cut -f1|awk '!x[$0]++'|wc -l >> readsCount
echo "pair reads with uniqmap:" >> readsCount
cat uniqmap.sam|awk '$7~/=/'|cut -f1|awk '!x[$0]++'|wc -l >> readsCount

samtools view -S uniqmap.sam -b -o uniqmap.bam
bedtools bamtobed -i uniqmap.bam -split > uniqmap.bed
rm uniqmap.sam


intersectBed -a uniqmap.bed -b $single -wa -wb -f 0.5 > tmp.SG2
echo "" >> readsCount
echo "" >> readsCount
echo "" >> readsCount

echo "mapping reads:" >> readsCount
cat uniqmap.bed|awk '!x[$4]++'|wc >> readsCount
echo "ncRNA+mRNA reads:" >> readsCount
cat tmp.SG2|awk '!a[$4]++'|wc >> readsCount
echo "mRNA region reads:" >> readsCount

cat tmp.SG2|awk '$13=="protein_coding"||$13=="mRNA"'|cut -f14|sed -e 's/_/\t/'|cut -f1|awk -v OFS="\t" '{num[$1]++}END{for(i in num){print i,num[i]}}' >> readsCount
echo "" >> readsCount
echo "" >> readsCount
echo "" >> readsCount


##htSeq
cd ${data}/hisat2
samtools sort -n uniqmap.bam -@ 6 -o sort.bam
samtools index sort.bam

mkdir ./htSeq
htseq-count -m union -f bam -s no sort.bam $gtf > ./htSeq/union_no.out
rm sort.bam  hisat2.bam

samtools sort uniqmap.bam -@ 6 -o sort.bam
samtools index sort.bam
bamCoverage -b sort.bam --normalizeUsing RPKM -o sort.bw

































































