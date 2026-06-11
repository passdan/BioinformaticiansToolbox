#!/bin/bash

## Useful shortcuts
export workingdir=~/RNAseq-Processing
export sing_folder=~/Share/Singularity_images
threads=2

## The commands you want to run
mkdir -p $workingdir/star

#list=("sample1" "sample2" "sample3")
list=("SRR5222797_10pc" "SRR5222798_10pc" "SRR5222799_10pc")

for i in ${list[@]}
do
# map forward and reverse reads to genome
# If input data is gzipped (.fastq.gz) inculde the additional parameter:   --readFilesCommand zcat
singularity exec --bind /home/ubuntu/Share/Day5/REFS:/data/REFS $sing_folder/rnaseq-mini_latest.sif \
       STAR   --outMultimapperOrder Random \
       --outSAMmultNmax 1 \
       --runThreadN $threads  \
       --runMode alignReads \
       --outSAMtype BAM Unsorted \
       --quantMode GeneCounts \
       --outFileNamePrefix $workingdir/star/${i}-unsort. \
       --genomeDir /data/REFS \
       --readFilesIn $workingdir/fastq/${i}-trim_1.fastq $workingdir/fastq/${i}-trim_2.fastq
done       

