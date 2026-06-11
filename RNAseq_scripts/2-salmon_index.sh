#!/bin/bash

## Useful shortcuts
export workingdir=~/RNAseq-Processing
export sing_folder=~/Share/Singularity_images
export refdir=/home/ubuntu/Share/Day5/REFS
threads=2

## Salmon image. If it isn't in your Singularity_images folder, pull it once with:
##   singularity pull $sing_folder/salmon_1.10.3.sif docker://combinelab/salmon:1.10.3
salmon_img=$sing_folder/salmon_1.10.3.sif

## Build the index (-k 31 suits reads >= ~75bp
singularity exec --bind $refdir:/data/REFS $salmon_img \
      salmon index -t /data/REFS/Arabidopsis_thaliana.TAIR10.cdna.all.fa \
          -i /data/REFS/salmon_index -k 31 -p $threads
