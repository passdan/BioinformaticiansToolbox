#!/bin/bash

## Useful shortcuts
export workingdir=~/RNAseq-Processing
export sing_folder=~/Share/Singularity_images
threads=2


#list=("sample1" "sample2" "sample3")
list=("SRR5222797_10pc" "SRR5222798_10pc" "SRR5222799_10pc")

for i in ${list[@]}
do
## The commands you want to run
	# fastqc the raw data 
	singularity exec $sing_folder/fastqc_v0.11.9_cv8.sif fastqc -t $threads $workingdir/fastq/${i}_1.fastq
	singularity exec $sing_folder/fastqc_v0.11.9_cv8.sif fastqc -t $threads $workingdir/fastq/${i}_2.fastq

	# Run qc with fastp
	singularity exec $sing_folder/fastp_0.23.1.sif fastp \
		--in1 $workingdir/fastq/${i}_1.fastq \
	    	--in2 $workingdir/fastq/${i}_2.fastq \
	        --out1 $workingdir/fastq/${i}-trim_1.fastq \
	        --out2 $workingdir/fastq/${i}-trim_2.fastq \
		--cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 30 \
	        --qualified_quality_phred 30 --unqualified_percent_limit 30 \
        	--n_base_limit 5 --length_required 60 --detect_adapter_for_pe \
		-w $threads

	# fastqc the outputs
	singularity exec $sing_folder/fastqc_v0.11.9_cv8.sif fastqc -t $threads $workingdir/fastq/${i}-trim_1.fastq
	singularity exec $sing_folder/fastqc_v0.11.9_cv8.sif fastqc -t $threads $workingdir/fastq/${i}-trim_2.fastq

done
