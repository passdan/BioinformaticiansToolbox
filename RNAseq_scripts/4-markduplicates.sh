#!/bin/bash

## Useful shortcuts
export workingdir=~/RNAseq-Processing
export sing_folder=~/Share/Singularity_images
threads=2

mkdir -p $workingdir/markdup

#list=("sample1" "sample2" "sample3")
list=("SRR5222797_10pc" "SRR5222798_10pc" "SRR5222799_10pc")


for i in ${list[@]}
do

singularity exec $sing_folder/rnaseq-mini_latest.sif \
	samtools sort -@ $threads -o $workingdir/star/${i}.sorted.bam $workingdir/star/${i}-unsort.Aligned.out.bam

## REMOVE DUPLICATES ##
singularity exec $sing_folder/rnaseq-mini_latest.sif \
	picard MarkDuplicates \
		I=$workingdir/star/${i}.sorted.bam \
		O=$workingdir/markdup/${i}.rmdup.bam \
		M=$workingdir/markdup/${i}.metrics.rmdup.txt \
		REMOVE_DUPLICATES=true \
		VALIDATION_STRINGENCY=SILENT
done
