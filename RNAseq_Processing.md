# Processing RNAseq with a reference genome

Three main steps are taken to analyse a series of RNAseq fastq sample files (following quality filtering): 
- Alignment: 
  - Whole genome: e.g. STAR An alignment package which functions similarly to standard genome alignments but is designed for short regions of RNA that could span intron-exon junctions and with low compute requirements. STAR outputs a bam format file which contains the locations where all the reads in your dataset have aligned and the genes they cover.
  - Psuedo-alignment: e.g. Salmon
- Counting: FeatureCounts takes the positions of mapped reads in a genome file and outputs a file quantifying the expression of each gene or exon (based on parameter choices). N.b. not required with psuedo-alignment as it is completed in the process
    - One ***potential*** step prior to counting is removing duplicates that arise from data generation. This used to be very common and investigated, **but now is not advised** as it is found to introduce more errors. We can run this step just for our information using the picard tool MarkDuplicates.
- Differential Gene Analysis: Contrasting the expression profile of the samples is typically done with one of two R packages: Deseq2 or EdgeR however a multitude of alternatives exist. These packages perform the normalization and statistical steps of contrasting samples as defined in a metadata file (replicates, tissue type, treatment etc). 

Following these three steps, there are an unending number of tools and packages to look deeper into your data, find experimentally specific insights, and prior published data to contrast against.

## Data

This data comes from a paper looking at the chromatin organisation within the Arabidopsis genome (Genome-wide chromatin mapping with size resolution reveals a dynamic sub-nucleosomal landscape in Arabidopsis - [https://doi.org/10.1371/journal.pgen.1006988](https://doi.org/10.1371/journal.pgen.1006988))

Full data is available here: https://www.ebi.ac.uk/ena/browser/view/PRJNA369530

### Basic genome alignment process
#### QC
```
fastp \
	--in1 fastq/Sample_1.fastq \
	--in2 fastq/Sample_2.fastq \
	--out1 fastq/Sample-trim_1.fastq \
	--out2 fastq/Sample-trim_2.fastq

fastqc -t $threads fastq/Sample-trim_*.fastq
```

#### Index the reference genome
```
STAR 	\
	--runThreadN 8 \
	--runMode genomeGenerate \
	--genomeDir  data/REFS \
	--genomeFastaFiles ReferenceGenome.fa \
	--sjdbGTFfile ReferenceGeneLocations.gtf \
	--sjdbOverhang 99
```
Set `--sjdbOverhang` to your read length minus 1 (here 99 for 100 bp reads); the default of 100 is a reasonable fallback if read lengths vary.

#### Map reads against reference genome
```
STAR \
    --outMultimapperOrder Random \
    --runThreadN $threads  \
    --runMode alignReads \
    --quantMode GeneCounts \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix star/Sample. \
    --genomeDir data/REFS \
    --readFilesIn fastq/Sample-trim_1.fastq fastq/Sample-trim_2.fastq
```
With `--outSAMtype BAM SortedByCoordinate`, STAR writes a sorted BAM named `star/Sample.Aligned.sortedByCoord.out.bam`. Note that `--quantMode GeneCounts` *also* writes a per-gene count table `star/Sample.ReadsPerGene.out.tab`, so STAR alone already gives you counts — featureCounts below is the second of two possible counting routes.

#### Optional: mark duplicates (for QC information only)
```
picard MarkDuplicates \
		I=star/Sample.Aligned.sortedByCoord.out.bam \
		O=markdup/Sample.markdup.bam \
		M=markdup/Sample.metrics.txt \
		REMOVE_DUPLICATES=false
```
For RNAseq we set `REMOVE_DUPLICATES=false` — we *mark* duplicates so MultiQC can report the rate, but we don't remove them, because in RNAseq high-expression genes legitimately produce many identical fragments and removing them biases the counts. So the counting step below reads the STAR BAM, not the dedup one.

Count features/genes
```
featureCounts \
	-T $threads -p --countReadPairs -F GTF -t exon -g gene_id \
	-a data/REFS/ReferenceGeneLocations.gtf \
	-o featureCounts/Sample.featurecount \
	star/Sample.Aligned.sortedByCoord.out.bam
```
Note: In Subread ≥ 2.0.2 `-p` only declares the data is paired; add `--countReadPairs` to count fragments rather than individual reads.


### Pseudo-alignment - Alignment to transcriptome (Salmon / kallisto)

The STAR → featureCounts route above aligns every read to genomic coordinates, which is accurate but slow and disk-heavy. A faster mainstream alternative is **pseudo-alignment** (or "lightweight mapping") with common tools **Salmon** or **kallisto**. Rather than a full base-by-base alignment, these match each read directly to the *transcriptome* and quantify in a fraction of the time and memory, while handling multi-mapping reads and isoforms more gracefully.

```
# build an index once from the transcriptome (cDNA) fasta
salmon index -t Arabidopsis_thaliana.TAIR10.cdna.fa -i salmon_index

# quantify each sample
salmon quant -i salmon_index -l A \
    -1 fastq/Sample-trim_1.fastq -2 fastq/Sample-trim_2.fastq \
    -p $threads --validateMappings -o salmon/Sample
```

This produces per-transcript estimates (`quant.sf`), which you summarise to gene level with the **tximport** R package and feed straight into DESeq2 — no BAM, sorting, MarkDuplicates or featureCounts step needed.

### Which to choose
- pseudo-alignment when you only need gene/transcript *quantification* (most differential-expression studies); 
- full alignment with STAR when you also need the BAM itself — for novel-junction discovery, variant calling, or coverage visualisation.

### Data

In the `~/Share/Day5/RNAseq-Processing` folder there are three pairs of RNAseq files from an Arabidopsis RNAseq study. In the folder `Share/Day5/REFS` there is a reference genome, and a gtf file. The step 2 “star index genome” has already been run for you (you don’t need to do this!)
```
$ ls Share/Day5/RNAseq-Processing:
1-QC.sh
2-salmon_index.sh
2-star_index_genome.sh
3-salmon_quant.sh
3-star_quant.sh
4-markduplicates.sh
5-featurecounts.sh
```
```
$ ls Share/Day5/RNAseq-Processing/fastqs
SRR5222797_1.fastq    SRR5222797_2.fastq
SRR5222798_1.fastq    SRR5222798_2.fastq
SRR5222799_1.fastq    SRR5222799_2.fastq

$ ls Share/Day5/REFS
Arabidopsis_thaliana.TAIR10.53.gtf
Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa
……… <Lots of other index files for star to function that you don’t need to touch!>
```

### EXERCISES

Using the pre-made scripts perform the steps on three pairs of fastq files. There are examples of all of these files in the ```~/Share/Day5/RNAseq-Processing``` directory which you should copy into your own folder. You may need to edit them to represent your own working folder and filenames

Note: Here, we will use a local version of the singularity image files to avoid downloading uniquely each (alternatively replace with ```docker://passdan/rnaseq-mini```).

0. Copy the folder ```~/Share/Day5/RNAseq-Processing``` to your local directory (```cp -r```) and enter it with cd.
1. Choose if you want to run Genome alignment or psuedo-alignment
2. Perform the processing steps:
   1. QC and trim your sample data (script 1)
   2. Use your trimmed data as inputs to run star or salmon (script 3)
   3. If using star, use featureCounts to count abundance of mapped reads to each gene (script 5)
3. Review the outputs to see what files you have created!

4. Create a MultiQC report for  the outputs

You can run multiQC on the processed directory using this full command (you don’t need to give any additional parameters):

```
singularity exec --bind `pwd`:`pwd` --pwd `pwd` docker://multiqc/multiqc:latest multiqc .
```
	
These outputs are now ready to put into R and perform Differential Gene Expression Analysis!
