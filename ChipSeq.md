# A Brief Introduction to ChIPSeq

Here we'll process through the main stages of a ChIPSeq analysis. I have only included simple default parameters for QC and mapping in this tutorial so we can focus on ChIP, but you should read the full QC and genome alignment tutorials to parameterise these correctly.

### Source Data
```
# Download C. elegans genome
wget https://ftp.ensembl.org/pub/release-113/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna_rm.chromosome.I.fa.gz
```

### Samples
Distribution of histone modifications across the genome in C. elegans sperm vs. oocytes vs. early embryos

BioProject Accession Number: [PRJNA475794](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115704)
```
# filename              GEO             description
SRR7297994.fastq.gz     GSM3187950      H3K27me3OocyteChIP1
SRR7297995.fastq.gz     GSM3187951      H3K27me3OocyteChIP2
SRR7297996.fastq.gz     GSM3187952      H3K27me3SpermChIP1
SRR7297997.fastq.gz     GSM3187953      H3K27me3SpermChIP2

SRR7297998.fastq.gz     GSM3187954      H3K36me3OocyteChIP1
SRR7297999.fastq.gz     GSM3187955      H3K36me3OocyteChIP2
SRR7298000.fastq.gz     GSM3187956      H3K36me3SpermChIP1
SRR7298001.fastq.gz     GSM3187957      H3K36me3SpermChIP2

SRR7298003.fastq.gz     GSM3187959      H3K4me3OocyteChIP1
SRR7298004.fastq.gz     GSM3187960      H3K4me3OocyteChIP2
SRR7298005.fastq.gz     GSM3187961      H3K4me3SpermChIP2
SRR7298006.fastq.gz     GSM3187962      H3K4me3SpermChIP2

SRR7298009.fastq.gz     GSM3187965      OocyteInput1
SRR7298010.fastq.gz     GSM3187966      OocyteInput2
SRR7298011.fastq.gz     GSM3187967      SpermInput1
SRR7298012.fastq.gz     GSM3187968      SpermInput2


```


## ChIPSeq Data Processing
### Pre-Processing
QC
```
fastp \
	--in fastq/Sample1.fastq.gz \
	--out fastq/Sample1-trim.fastq.gz \

fastqc -t $threads fastq/Sample*.fastq.gz
```
Check that your data looks good quality and reasonable.


## Genome Alignment (one sample)
Typically either BWA or Bowtie2 is used for standard mapping of reads to the reference genome. BWA (or BWA-mem) often for short reads (i.e. 50-100bp) and bowtie2 for longer illumina reads, however these are not fixed rules!

#### Index genome
First we need to convert the reference genome into and indexed format. Only need to do this once ever! 
```
bowtie2-build c_elegans.genomic.fa.gz c_elegans_index
```
#### Align reads against the reference genome 
```
# One sample
bowtie2 -x c_elegans_index \
    -p 4 \
    -q Sample1-trim.fastq.gz \
    -S Sample1.sam
```
#### Format the outputs for downstream processing
```
samtools view -bS Sample1-align.sam | samtools sort -o Sample1.bam
samtools index Sample1.bam
```
#### Make a visualisation for genome browsers
```
# Generate bedGraph (intermediate format)
genomecov -ibam Sample1.bam \
    -g c_elegans.genomic.fa.gz \
    -bg > Sample1.bedGraph

# Convert bedGraph to bigWig
bedGraphToBigWig Sample1.bedGraph \
    c_elegans.genomic.fa.gz \
    Sample1.bw
```


## Genome Alignment Multiple samples loop
Of course, we'll usually have many samples to process together. Here is a simple loop of the above steps. Make sure to have the genome indexed first, and modify the code to reflect the location on your filesystem.
```
#!/bin/bash

# List of sample names (not the .fastq.gz part)
samples=("Sample1" "Sample2" "Sample3")

REFERENCE_INDEX="c_elegans_index"
REFERENCE_GENOME="c_elegans.genomic.fa.gz"

# Number of threads to use
THREADS=4

# Loop through each sample
for sample in "${samples[@]}"; do
    echo "Processing ${sample}"
    
    singularity exec docker://staphb/fastp \
        fastp \
    	    --in fastq/${sample}.fastq.gz \
	        --out fastq/${sample}-trim.fastq.gz \

    # Alignment step
    singularity exec docker://staphb/bowtie2 \
        bowtie2 -x "${REFERENCE_INDEX}" \
            -p "${THREADS}" \
            -q "${sample}-trim.fastq.gz" \
            -S "${sample}.sam"
    
    # Convert SAM to sorted BAM
    singularity exec docker://staphb/samtools \
        samtools view -bS "${sample}.sam" | samtools sort -o "${sample}.bam"
    singularity exec docker://staphb/samtools \    
        samtools index "${sample}.bam"
    
    # Generate genome coverage files
    singularity exec docker://staphb/bedtools \
        genomecov -ibam "${sample}.bam" \
            -g "${REFERENCE_GENOME}" \
            -bg > "${sample}.bedGraph"
        
    # Convert bedGraph to bigWig
    singularity exec docker://staphb/bedtools \
        bedGraphToBigWig "${sample}.bedGraph" \
            "${REFERENCE_GENOME}" \
            "${sample}.bw"
    
    # Optional: Clean up intermediate files
    rm "${sample}.sam" "${sample}.bedGraph"
done
```

### Review our outputs
Lets use multiqc to summarise all of out samples and review our alignment statistics
```
singularity exec --bind `pwd`:`pwd` --pwd `pwd` docker://ewels/multiqc:latest multiqc .
```

### Genome visualisation
Lets use our wig files (bw - bigwig) to look at the mapped reads.

In the browser go to https://igv.org/app/, select _C elegans_ and import the `.bw` files.

## Peak Calling
Now that we have seen that we have some peaks in our data, lets technically identify them

### MACS2
```
GENOME_SIZE="100286401"
mkdir -p macs2

singularity exec docker://dceoy/macs2 \
    macs2 callpeak \
        -t Sample1_chip.bam \
        -c Sample1_input.bam \
        -f BAM \
        -g ${GENOME_SIZE} \
        -n Sample1 \
        --outdir macs2 \
        -q 0.01 \
        --broad \
        2> ${OUTPUT_DIR}/logs/macs2_log.txt
```
Now lets take our annotated peaks and add them to the IGV browser to see agreement (or not)

A good practice is to run additional peak calling algorythms and use those in agreement as your default set.

### HOMER
```
# 1. Make Tag Directory from BAM
makeTagDirectory Sample1_tag_directory \
    Sample1.bam \
    -format BAM \
    -illuminaReads  # Use this for Illumina sequencing data

makeTagDirectory Sample1_input_tag_directory \
    Sample1_Input.bam \
    -format BAM \
    -illuminaReads

# 2. Call Peaks 
# For Histone Modifications (Broad Peaks)
findPeaks Sample1_tag_directory \
    -style histone \S
    -o Sample1_histone_peaks.txt \
    -i input_tag_directory  

# 3. Convert HOMER peaks to BED format
pos2bed.pl Sample1_histone_peaks.txt > Sample1_peaks.bed
```
We can check this in the genome browser too. How well do the two peak calling methods agree?

## Annotation & Interpretation
We now have processed our raw data through the essential steps and it is now ready for further downstream analysis. 