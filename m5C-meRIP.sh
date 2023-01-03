#!/bin/sh
single=./Mouse/ENSEMBLE_68/GeneList/single-transcript-chr

genome="./Mouse/ENSEMBLE_68/index/genome.fa"
rRNA="./Mouse/ENSEMBLE_68/45S/45S_Mus.fasta"
gtf="./Mouse/ENSEMBLE_68/GTF/Mus_musculus.GRCm38.68.chr.gtf"


path=ywl-Th17-m5C-IP
thread=5

R1=ywl-Th17-m5C-IP_1.fq.gz
R2=ywl-Th17-m5C-IP_2.fq.gz


#################pre-processing
cd $path
mkdir -p clean_data/fastqc
cd clean_data
fastqc $R1 -t $thread -o ./fastqc
fastqc $R2 -t $thread -o ./fastqc
cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -U 3 -o cutadapt_1.fastq.gz -p cutadapt_2.fastq.gz  $R1  $R2

fastqc cutadapt_1.fastq.gz -t $thread -o ./fastqc
fastqc cutadapt_2.fastq.gz -t $thread -o ./fastqc

java -Xmx4g -jar trimmomatic-0.36.jar PE -phred33 cutadapt_1.fastq.gz cutadapt_2.fastq.gz ./trim_qua_1.fastq.gz ./unpaired_1.fastq.gz ./trim_qua_2.fastq.gz ./unpaired_2.fastq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20
rm  cutadapt_1.fastq.gz cutadapt_2.fastq.gz
fastqc ./trim_qua_1.fastq.gz -t $thread -o ./fastqc
fastqc ./trim_qua_2.fastq.gz -t $thread -o ./fastqc

###############rRNA filtering
cd $path
mkdir rRNA_mapping
hisat2 -p $thread -N 1 --dta -x $rRNA  -1 ./clean_data/trim_qua_1.fastq.gz -2 ./clean_data/trim_qua_2.fastq.gz -S rRNA_mapping/hisat2.sam --un-conc ./rRNA_mapping/2_delrRNA 2> mapping_stat

gzip ./rRNA_mapping/2_delrRNA.1
gzip ./rRNA_mapping/2_delrRNA.2

fastqc ./rRNA_mapping/2_delrRNA.1.gz -t $thread -o ./clean_data/fastqc
fastqc ./rRNA_mapping/2_delrRNA.2.gz -t $thread -o ./clean_data/fastqc
rm ./rRNA_mapping/hisat2.sam

mv ./rRNA_mapping/2_delrRNA.1.gz ./clean_data
mv ./rRNA_mapping/2_delrRNA.2.gz ./clean_data


#################mapping
cd $path
mkdir hisat2
hisat2 -p $thread -N 1 --dta -x $genome -1 ./clean_data/2_delrRNA.1.gz -2 ./clean_data/2_delrRNA.2.gz --un-conc ./ -S hisat2/hisat2_output.sam

##bam2bed
cd hisat2
echo "reads with no -q:" >readsCount
cat hisat2_output.sam|grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >>readsCount
samtools view -S hisat2_output.sam -b -o hisat2_output.bam
rm hisat2_output.sam
echo "reads mapped without -q:" >>readsCount
samtools view  hisat2_output.bam |awk '$3!="*"'|awk '!x[$1]++'|wc -l >>readsCount

samtools view hisat2_output.bam -q 20 -h|grep -E "^@|NH:i:1$|NH:i:1[^0-9]"|awk '$1~/^@/ || $1~/^[A-Z]/' > uniqmap.sam
echo "reads with uniqmap:" >>readsCount
cat uniqmap.sam|cut -f1|awk '!x[$0]++'|wc -l >>readsCount
echo "pair reads with uniqmap:" >>readsCount
cat uniqmap.sam|awk '$7~/=/'|cut -f1|awk '!x[$0]++'|wc -l >>readsCount

samtools view -S uniqmap.sam -b -o uniqmap.bam
rm uniqmap.sam
bedtools bamtobed -i uniqmap.bam -split > uniqmap.bed

intersectBed -a uniqmap.bed -b $single -wa -wb -f 0.5 > tmp.bed
echo "" >> readsCount
echo "" >> readsCount
echo "" >> readsCount

echo "mapping reads:" >> readsCount
cat uniqmap.bed|awk '!x[$4]++'|wc >> readsCount
echo "ncRNA+mRNA reads:" >> readsCount
cat tmp.bed |awk '!a[$4]++'|wc >> readsCount
echo "ncRNA reads:" >> readsCount
cat tmp.bed |awk '$13!="mRNA"'|awk '!a[$4]++'|wc >> readsCount


echo "mRNA region reads:" >> readsCount
cat tmp.bed|awk '$13=="mRNA"||$13=="protein_coding" '|cut -f14|sed -e 's/_/\t/'|cut -f1|awk -v OFS="\t" '{num[$1]++}END{for(i in num){print i,num[i]}}' >> readsCount
echo "" >> readsCount
echo "" >> readsCount
echo "" >> readsCount

rm tmp.bed



