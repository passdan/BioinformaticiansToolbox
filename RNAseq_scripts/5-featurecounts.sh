#!/bin/bash

## Useful shortcuts
export workingdir=~/RNAseq-Processing
export sing_folder=~/Share/Singularity_images
threads=2

mkdir -p $workingdir/featureCounts

#list=("sample1" "sample2" "sample3")
list=("SRR5222797_10pc" "SRR5222798_10pc" "SRR5222799_10pc")

for i in ${list[@]}
do

singularity exec --bind /home/ubuntu/Share/Day5/REFS:/data/REFS $sing_folder/rnaseq-mini_latest.sif \
	featureCounts \
	-T $threads -p -F GTF -t exon -g gene_id \
	-a /data/REFS/Arabidopsis_thaliana.TAIR10.53.gtf \
	-o $workingdir/featureCounts/${i}.markdup.featurecount \
	$workingdir/star/${i}.sorted.bam
done
