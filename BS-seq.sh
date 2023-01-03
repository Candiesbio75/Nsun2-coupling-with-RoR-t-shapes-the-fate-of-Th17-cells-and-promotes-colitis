#!/bin/sh
thread=12

call=./meRanTK-1.2.0/meRanCall
align=./meRanTK-1.2.0/meRanGh

#######################################################
path=./BSseq/ywl-Th17-CKO-BS-1
sample=ywl-Th17-CKO-BS-1

R1=./ywl-Th17-CKO-BS-1/ywl-Th17-CKO-BS-1_R1.fq.gz
R2=./ywl-Th17-CKO-BS-1/ywl-Th17-CKO-BS-1_R2.fq.gz

#####################trim
cd ${path}
mkdir -p clean_data/fastqc
cd clean_data

fastqc ${R1} -t ${thread} -o ./fastqc
fastqc ${R2} -t ${thread} -o ./fastqc
cutadapt -a AGATCGGAAG -A AGATCGGAAGA -o cutadapt_${sample}_1.fastq.gz -p cutadapt_${sample}_2.fastq.gz ${R1} ${R2}

fastqc cutadapt_${sample}_1.fastq.gz -t ${thread} -o ./fastqc
fastqc cutadapt_${sample}_2.fastq.gz -t ${thread} -o ./fastqc

java -Xmx4g -jar trimmomatic-0.36.jar PE -phred33 cutadapt_${sample}_1.fastq.gz cutadapt_${sample}_2.fastq.gz ./${sample}_trim_qua_1.fastq.gz ./unpaired_1.fastq.gz ./${sample}_trim_qua_2.fastq.gz ./unpaired_2.fastq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:35

rm  cutadapt_${sample}_1.fastq.gz cutadapt_${sample}_2.fastq.gz
fastqc ./${sample}_trim_qua_1.fastq.gz ./${sample}_trim_qua_2.fastq.gz -t ${thread} -o ./fastqc

gunzip ./${sample}_trim_qua_1.fastq.gz
gunzip ./${sample}_trim_qua_2.fastq.gz

#####################################################
##DHFR Bisulfite mapping#####################################################
htseq_index="./Dhfr_hx/Dhrf.fasta"
meRanT_Dhfr_index="./Dhfr_hx/meRanT_index/Dhrf_C2T"
map=./Dhfr_hx/Dhrf.map

cd ${path}
mkdir DHFR_meRanTK-T
cd DHFR_meRanTK-T
d1=../clean_data/${sample}_trim_qua_1.fastq
d2=../clean_data/${sample}_trim_qua_2.fastq

meRanT align -o ./meRanGhResult -f $d2 -r $d1 -t $thread -S tmp.sam -un -ds -i2g $map -x $meRanT_Dhfr_index  -mbp -fmo -mmr 0.01

cd meRanGhResult
meRanCall -f $htseq_index -p $thread -bam tmp_sorted.bam -gref -o out.txt -mBQ 20 -mr 0
rm *.sam *unmapped.fastq  *_sorted*.bam *_sorted.bam.bai

##Bisulfite mapping to mouse
########
htseq_index="./mouse/index/genome.fa"
meRanTK_index="./mouse/ENSEMBLE_68/meRanTK-1.2.0/hisat2_index/"
gtf="./mouse/ENSEMBLE_68/GTF/Mus_musculus.GRCm38.68.chr.gtf"


cd ${path}
mkdir meRanTK
cd meRanTK
d1=../clean_data/${sample}_trim_qua_1.fastq
d2=../clean_data/${sample}_trim_qua_2.fastq

$align align -o ./meRanGhResult -f $d2 -r $d1 -t ${thread} -S tmp_${sample}.sam -un -ds -id $meRanTK_index -GTF $gtf -mbp -fmo -mmr 0.01


mv ./meRanGhResult/*1_unmapped.fastq ${path}/clean_data/3_deHuman_1.fastq
mv ./meRanGhResult/*2_unmapped.fastq ${path}/clean_data/3_deHuman_2.fastq

gzip ${path}/clean_data/3_deHuman_1.fastq
gzip ${path}/clean_data/3_deHuman_2.fastq

fastqc  ${path}/clean_data/3_deHuman_1.fastq.gz ${path}/clean_data/3_deHuman_2.fastq.gz -t $thread -o ${path}/clean_data/fastqc

##m5C calling
cd  ${path}/meRanTK/meRanGhResult
##########################################
if [[  ! -e ${sample}_sorted.bam.bai  &&   -e tmp_${sample}.sam   ]]; then
head -n29 tmp_${sample}.sam > ${sample}.sam
awk -v F="\t" 'NF==18 && length($10)==length($11) && $1~/^[A-Z]00/ && $0!~/::/' tmp_${sample}.sam|perl validate_CIGAR_zt.pl - >> ${sample}.sam
samtools view ${sample}.sam -Sb -o ${sample}.bam
samtools sort ${sample}.bam ${sample}_sorted
samtools index ${sample}_sorted.bam
rm ${sample}.sam
$call -f $htseq_index -p ${thread} -bam ${sample}_sorted.bam -gref -o ${sample}.txt -mBQ 20 -mr 0 -fdr 1

elif [ -e tmp_${sample}_sorted.bam.bai ]; then
$call -f $htseq_index -p ${thread} -bam tmp_${sample}_sorted.bam -gref -o ${sample}.txt -mBQ 20 -mr 0 -fdr 1

elif [ -e ${sample}_sorted.bam.bai ]; then
$call -f $htseq_index -p ${thread} -bam ${sample}_sorted.bam -gref -o ${sample}.txt -mBQ 20 -mr 0 -fdr 1

fi


