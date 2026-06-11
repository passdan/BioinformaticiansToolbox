#!/bin/bash

## Useful shortcuts
export workingdir=~/RNAseq-Processing
export sing_folder=~/Share/Singularity_images
export refdir=/home/ubuntu/Share/Day5/REFS
threads=2

## Salmon image. If it isn't in your Singularity_images folder, pull it once with:
##   singularity pull $sing_folder/salmon_1.10.3.sif docker://combinelab/salmon:1.10.3
salmon_img=$sing_folder/salmon_1.10.3.sif

mkdir -p $workingdir/salmon

#list=("sample1" "sample2" "sample3")
list=("SRR5222797_10pc" "SRR5222798_10pc" "SRR5222799_10pc")


## Quantify each sample (this is the bit you run).
for i in ${list[@]}
do
	singularity exec --bind $refdir:/data/REFS $salmon_img \
		salmon quant -i /data/REFS/salmon_index -l A \
			-1 $workingdir/fastq/${i}-trim_1.fastq \
			-2 $workingdir/fastq/${i}-trim_2.fastq \
			-p $threads --gcBias \
			-o $workingdir/salmon/${i}
done

## Each sample now has $workingdir/salmon/<sample>/quant.sf
## In R: tximport() these to gene level (using a tx2gene table from the GTF),
## then DESeqDataSetFromTximport() -> standard DESeq2 / SARTools as in the analysis session.
