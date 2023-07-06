# Processing RNAseq with a reference genome

Three main steps are taken to analyse a series of RNAseq fastq sample files (following quality filtering): 
- Alignment: STAR – (Spliced Transcript Alignments to a Reference) is an alignment package which functions similarly to standard genome alignments but is designed for short regions of RNA that could span intron-exon junctions and with low compute requirements. STAR outputs a bam format file which contains the locations where all the reads in your dataset have aligned and the genes they cover.
- Counting: FeatureCounts is a simple package that takes the positions of mapped reads and outputs a file quantifying the expression of each gene or exon (based on parameter choices). At this point raw read counts are hard to interpret due to likely different levels of sequencing achieved per sample and methodological biases. 
    - One step prior to counting is marking duplicates that arise from data generation for further information, or so that they can be removed. This used to be very common and investigated, but these days it is not too common as it is found to introduce more errors. We can run this step just for our information using the picard tool MarkDuplicates.
- Differential Gene Analysis: Contrasting the expression profile of the samples is typically done with one of two R packages: Deseq2 or EdgeR (the mac vs windows of the RNAseq fight), however a multitude of alternatives exist. These packages perform the normalization and statistical steps of contrasting samples as defined in a metadata file stating your experimental design (replicates, tissue type, treatment etc). The output here is a range of significant genes, ordinance and cluster analysis of sample similarity, and various quality control figures.

Following these three steps, there are an almost infinite number of tools and packages to look deeper into your data, find experimentally specific insights, and prior published data to contrast against.

## Data

This data comes from a paper looking at the chromatin organisation within the Arabidopsis genome (Genome-wide chromatin mapping with size resolution reveals a dynamic sub-nucleosomal landscape in Arabidopsis - [https://doi.org/10.1371/journal.pgen.1006988](https://doi.org/10.1371/journal.pgen.1006988))

Full data is available here: https://www.ebi.ac.uk/ena/browser/view/PRJNA369530

We will be using  scripts to run these steps. In the Share/Day4 folder you will find the following that you can use to base your analysis, however make sure you’re tuning it to your own file structure and file names. 

So far we have used only a small dataset to quickly practice the steps but now we’ll be using full sized RNAseq samples. This is because otherwise it causes the programs to think it’s bad data and causes errors. 

In the Share/Day4 folder there are three pairs of RNAseq files from an Arabidopsis RNAseq study. In the folder Share/REFS there is a reference genome, and a gtf file. The step 2 “star index genome” has already been run for you (you don’t need to do this!)
```
$ ls Share/Day4/Processing:
1-QC.sh  
2-star_index_genome.sh  (already done, don’t repeat!)
3-star.sh  
4-markduplicates.sh  
5-featurecounts.sh 
```
```
$ ls Share/Day4/fastqs
SRR5222797_1.fastq    SRR5222797_2.fastq
SRR5222798_1.fastq    SRR5222798_2.fastq
SRR5222799_1.fastq    SRR5222799_2.fastq

$ ls Share/REFS
Arabidopsis_thaliana.TAIR10.47.gtf
Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa
……… <Lots of other index files for star to function that you don’t need to touch!>
```

### EXCERCISES

We will be using a docker container with all of the tools in. All steps are done with this one container, for example for the first QC processing :
```
$ docker run \
	-u $(id -u):$(id -g) \
	-v $(pwd):/data \
	-v /home/ubuntu/Share/REFS/RNAseqREFS:/data/REFS \
	-it passdan/rnaseq-mini ./1-QC.sh
```

Using the pre-made scripts perform the steps on three pairs of fastq files. There are examples of all of these files in the Share/Day4 directory which you should copy into your own folder. You may need to edit them to represent your own working folder and filenames:

0. Copy the folder ```~/Share/Day4/RNAseq-Processing``` to your local directory and enter it.
1. Open and read the 5 scripts to understand what their functions are
2. QC and trim your sample data (script 1)
3. Use your trimmed data as inputs to run star (script 3)
4. [optional] Use the outputs from star to run mark duplicates to mark duplicates and optionally remove them (script 4)
5. Use featureCounts to count abundance of mapped reads to each gene (script 5)

<details>
    <summary>

### Extension: Create a MultiQC report for  the outputs

</summary>

6. Run multiQC on the processed directory using this full command (you don’t need to give any additional parameters):

```
$  docker run --rm -v $(pwd):/in -w /in -it ewels/multiqc:v1.12 .
```
</details>
	
These outputs are now ready to put into R and perform Differential Gene Expression Analysis!
