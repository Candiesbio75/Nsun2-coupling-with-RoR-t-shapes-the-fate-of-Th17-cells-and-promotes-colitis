#!/bin/sh
thread=12

####################
genome_index=./Mouse/ENSEMBLE_68/index/genome.fa
gtf=./Mouse/ENSEMBLE_68/GTF/Mus_musculus.GRCm38.68.chr.gtf
single=./Mouse/ENSEMBLE_68/GeneList/single-transcript-chr
size=./Mouse/ENSEMBLE_68/GeneList/mm-size

##################total trans
genebody=./Mouse/ENSEMBLE_68/GeneList/Mus_musculus.GRCm38.68_genebody_totalTrans.bed
############################################
# awk -v OFS='\t' '{if($6=="-"){print $1,$3-1000,$3+1000,$6,"proximal_promoter"}else{print $1,$2-1000,$2+1000,$6,"proximal_promoter"}}' $genebody >  ./Mouse/ENSEMBLE_68/GeneList/Mus_musculus.GRCm38.68_proximal_promoter_totalTrans.bed
#############################################
TSS=./Mouse/ENSEMBLE_68/GeneList/Mus_musculus.GRCm38.68_proximal_promoter_totalTrans.bed



#######################################################
#
#CHIP
#
#######################################################
sample=ywl-TH17-CHIP-10
R1=ywl-TH17-CHIP-10_R1.fq.gz
R2=ywl-TH17-CHIP-10_R2.fq.gz

path=./CHIP/ywl-TH17-CHIP-10


cd $path
mkdir clean_data
mkdir ./clean_data/fastqc
cd clean_data
cutadapt -a AGATCGGAAGAGCACACGTCT -A AGATCGGAAGAGCGTCGTGTAG -o 1_cutadapt_1.fastq.gz -p 1_cutadapt_2.fastq.gz $R1 $R2
cutadapt -u 10 -u -10 -U 10 -U -10 -o 2_cutadapt_1.fastq.gz -p 2_cutadapt_2.fastq.gz 1_cutadapt_1.fastq.gz 1_cutadapt_2.fastq.gz
java -Xmx4g -jar trimmomatic-0.36.jar PE -phred33 2_cutadapt_1.fastq.gz 2_cutadapt_2.fastq.gz ./3_trim_qua_1.fastq.gz ./unpaired_1.fastq.gz ./3_trim_qua_2.fastq.gz ./unpaired_2.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35

########################quality control
fastqc $R1 $R2 -t $thread -o ./fastqc
fastqc 2_cutadapt_1.fastq.gz -t $thread -o ./fastqc
fastqc 2_cutadapt_2.fastq.gz -t $thread -o ./fastqc
fastqc 3_trim_qua_1.fastq.gz -t $thread -o ./fastqc
fastqc 3_trim_qua_2.fastq.gz -t $thread -o ./fastqc
rm  1_cutadapt_1.fastq.gz 1_cutadapt_2.fastq.gz  2_cutadapt_1.fastq.gz 2_cutadapt_2.fastq.gz

###################rRNA filtering
cd $path
mkdir rRNA_mapping
cd  rRNA_mapping

rRNA=./Mouse/ENSEMBLE_68/45S/45S_Mus.fasta
bowtie2 -p $thread -x $rRNA -1 ../clean_data/3_trim_qua_1.fastq.gz -2 ../clean_data/3_trim_qua_2.fastq.gz  --un-conc ./ -S bowtie.sam 2>rRNA_mapping

gzip un-conc-mate.1 
gzip un-conc-mate.2

mv un-conc-mate.1.gz $path/clean_data/derRNA.1.gz
mv un-conc-mate.2.gz $path/clean_data/derRNA.2.gz

##################mapping 
cd $path
mkdir bowtie2
cd bowtie2
bowtie2 -p $thread -x $genome_index -1 ../clean_data/derRNA.1.gz -2 ../clean_data/derRNA.2.gz -t -q -N 1 -L 25 | samtools view -bS -1 -q 20 -h - -o unique.bam  

# ####################unique mapping read 
echo "reads with uniqmap:" > readsCount
samtools view  unique.bam |grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> readsCount
echo "pair reads with uniqmap:" >> readsCount
samtools view  unique.bam |grep -v "^@"|awk '$7~/=/'|cut -f1|awk '!x[$0]++'|wc -l >> readsCount

# #################sort the unique mapping bam and then index
samtools sort unique.bam -o sort.bam 
samtools index sort.bam 

# ####################deduplication
java -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true METRICS_FILE=sort.matrix INPUT=sort.bam OUTPUT=sorted.nodup.bam ASSUME_SORTED=true
samtools index sorted.nodup.bam

echo "deduplication infor:" >> readsCount
samtools flagstat sorted.nodup.bam >> readsCount
rm sort.bam   sort.bam.bai  

######################生成big wig文件用于IGV
genomeCoverageBed -bga -ibam ./sorted.nodup.bam -g $size -split |sort - -k1,1 -k2,2n > sorted.nodup.bedgraph

awk 'NR==FNR{sum+=($3-$2)*$4;len+=$3-$2}NR>FNR{print $1,$2,$3,$4*len/sum}' sorted.nodup.bedgraph sorted.nodup.bedgraph > norm_sorted.nodup.bedgraph
bedGraphToBigWig norm_sorted.nodup.bedgraph $size norm_sorted.nodup.bw

################bam2bed call peak
cd $path
mkdir macs
cd macs
bedtools bamtobed -i ../bowtie2/sorted.nodup.bam  -split > ../bowtie2/sorted.nodup.bed 
macs2 callpeak -t ../bowtie2/sorted.nodup.bed  -n $sample -f BED --nomodel -g mm -B -q 0.05


########Quality
cd $path
mkdir quality
cd quality

####TSS
computeMatrix reference-point  --referencePoint TSS  -p $thread  \
-b 10000 -a 10000    \
-R $gtf  \
-S ../bowtie2/*.bw  \
--skipZeros  -o matrix_TSS.gz  \
--outFileSortedRegions genes.bed

plotHeatmap -m matrix_TSS.gz  -out TSS_Heatmap.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m matrix_TSS.gz  -out TSS_Profile.pdf --plotFileFormat pdf --perGroup --dpi 720


#genebody
computeMatrix scale-regions  -p $thread  \
-R $gtf  \
-S ../bowtie2/*.bw \
-b 1000 -a 1000  \
--regionBodyLength \
--skipZeros -o matrix_body.gz

plotHeatmap -m matrix_body.gz  -out Body_Heatmap.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m matrix_body.gz -out Body_Profile.pdf --plotFileFormat pdf --perGroup --dpi 720


##########peak annotation
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ../macs/*peaks.narrowPeak > homer_peaks.tmp
./homer/bin/annotatePeaks.pl  homer_peaks.tmp mm10 -gtf $gtf 1>peakAnn.xls  2> annLog.txt

#########motif annotation
./homer/bin/findMotifsGenome.pl homer_peaks.tmp  mm10 motifDir -len 8,10,12 












































