#!/bin/sh

thread=5
####################需要准备的文件
genome_index=./Mouse/ENSEMBLE_68/index/genome.fa
gtf=./Mouse/ENSEMBLE_68/GTF/Mus_musculus.GRCm38.68.chr.gtf
single=./Mouse/ENSEMBLE_68/GeneList/single-transcript-chr
size=./Mouse/ENSEMBLE_68/GeneList/mm-size

#######################################################
#
#ATAC
#
#######################################################
sample=ywl-ATAC-CKO-1
path=./ATAC/${sample}

cd $path
mkdir clean_data
mkdir ./clean_data/fastqc
cd clean_data
R1=./$sample/${sample}_R1.fq.gz
R2=./$sample/${sample}_R2.fq.gz 

#Step1 =============================
##fastqc
# ==================
cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o 1_cutadapt_1.fastq.gz -p 1_cutadapt_2.fastq.gz $R1 $R2

cutadapt -a CTGTCTCTTATA -A CTGTCTCTTATA -o 2_cutadapt_1.fastq.gz -p 2_cutadapt_2.fastq.gz  1_cutadapt_1.fastq.gz  1_cutadapt_2.fastq.gz

java -Xmx4g -jar /software/biosoft/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 2_cutadapt_1.fastq.gz 2_cutadapt_2.fastq.gz ./3_trim_qua_1.fastq.gz ./unpaired_1.fastq.gz ./3_trim_qua_2.fastq.gz ./unpaired_2.fastq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:35

fastqc $R1 $R2 -t $thread -o ./fastqc
fastqc $R1 $R2 -t $thread -o ./fastqc
fastqc 2_cutadapt_1.fastq.gzgz -t $thread -o ./fastqc
fastqc 2_cutadapt_2.fastq.gzgz -t $thread -o ./fastqc
fastqc 3_trim_qua_1.fastq.gz -t $thread -o ./fastqc
fastqc 3_trim_qua_2.fastq.gz -t $thread -o ./fastqc
rm  1_cutadapt_1.fastq.gz 1_cutadapt_2.fastq.gz  2_cutadapt_1.fastq.gz 2_cutadapt_2.fastq.gz


##Step2 =============================
#mapping using bowtie
# ==================
cd $path
mkdir bowtie2
cd bowtie2

bowtie2 -p $thread -x $genome_index -1 ../clean_data/3_trim_qua_1.fastq.gz -2 ../clean_data/3_trim_qua_2.fastq.gz -t -q -N 1 -L 25 -X 2000 --no-mixed --no-discordant | samtools view -bS -1 -h - -o 0_align.bam 
echo "reads with no -q:" > readsCount
samtools view 0_align.bam|grep -v "^@"|awk '$3~/^chr/'|cut -f1|awk '!x[$0]++ '|wc -l >> readsCount

#Step3 =============================
# Remove unmapped, mate unmapped
# not primary alignment, reads failing platform
# ==================
samtools view -h -F 1804 -b 0_align.bam | samtools sort - -o 3_filter.bam
if [ -e 3_filter.bam.bam ]; then mv 3_filter.bam.bam 3_filter.bam; fi
echo "reads after Remove unmapped, mate unmapped, not primary alignment, reads failing platform :" >> readsCount
samtools view 3_filter.bam|grep -v "^@"|awk '$3~/^chr/'|cut -f1|awk '!x[$0]++'|wc -l >> readsCount

#Step4 ======================== 
# Mark duplicates 
# ======================
java -Xmx6G -jar picard.jar  MarkDuplicates  INPUT=3_filter.bam   OUTPUT=4_tmp.bam  METRICS_FILE=log_dup_qc_matrix VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false

mv 4_tmp.bam 3_filter.bam 

#step5 ============================
# Remove duplicates
# Index final position sorted BAM
# ============================
samtools view -F 1804 -b 3_filter.bam  > 4_final.bam
echo "reads after deduplication:" >> readsCount
samtools view 4_final.bam|grep -v "^@"|awk '$3~/^chr/'|cut -f1|awk '!x[$0]++'|wc -l >> readsCount

#step6 ============================= 
# Compute library complexity
# ============================= 
# sort by position and strand
# Obtain unique count statistics
# PBC File output
echo 'TotalReadPairs  DistinctReadPairs   OneReadPair  TwoReadPairs   NRF=Distinct/Total  PBC1=OnePair/Distinct PBC2=OnePair/TwoPair' >> readsCount
bedtools bamtobed -i 3_filter.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | grep -v 'chrM' | sort | uniq -c | awk -v OFS='\t' 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{print "TotalReadPairs="mt,"DistinctReadPairs="m0,"OneReadPair="m1,"TwoReadPairs="m2,"NRF="m0/mt,"PBC1="m1/m0,"PBC2="m1/m2}' >> readsCount
rm 3_filter.bam


#step7 ============================= 
# dechrMT
# ============================= 
samtools view 4_final.bam -h|grep  -v 'chrM'|samtools view -bS - |samtools sort - -o - > 5_deMT.bam
rm  4_final.bam
echo "Delete chrMT infor:" >> readsCount
samtools view 5_deMT.bam|grep -v "^@"|awk '$3~/^chr/'|cut -f1|awk '!x[$0]++'|wc -l >> readsCount


#step8 ============================= 
# unique mapping
# ============================= 
samtools view -h -q 20 -b 5_deMT.bam > 6_last.bam
samtools index 6_last.bam
rm 5_deMT.bam

echo "read with quality > Q20" >> readsCount
samtools view 6_last.bam|grep -v "^@"|awk '$3~/^chr/'|cut -f1|awk '!x[$0]++'|wc -l >> readsCount

#step9 ============================= 
# ============================= 
outpath=$path/bowtie2
bamfile=$path/bowtie2/6_last.bam
bamindex=$path/bowtie2/6_last.bam.bai

Rscript ./software/ATAC/ATACseq/R/shiftbam.R $outpath $bamfile $bamindex

#step10 ============================= 
# bw
# ============================= 
bamCoverage -p 12 --normalizeUsing RPKM -b 7_shifted.bam -o  8_norm.bw 

#step11 ============================= 
# call peak
# ============================= 
################bam2bed call peak
cd $path
mkdir macs
cd macs
bedtools bamtobed -i ../bowtie2/7_shifted.bam  -split > ../bowtie2/7_shifted.bed
macs2 callpeak -t ../bowtie2/7_shifted.bed  -n $sample -f BED --nomodel -g mm -B -q 0.05 --keep-dup all

#####step10 QC
cd $path
mkdir quality
Rscript ./ATACseq/R/QC_mm10.R $path/quality $path/bowtie2/6_last.bam  $path/bowtie2/6_last.bam.bai $sample 






