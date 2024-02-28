# Genome Assembly integrating Long and Short read data 

Here we are going to work with a 4.6Mb E. coli genome using a pipeline of a de novo assembly of middle-sized genomes at a low cost and with high quality using Pacific Biosciences and Illumina reads. Neither high quality short read, or the lower quality but long-read data can be suitably assembled for publication or downstream comparative genomics.

However, by incorporating the short reads with the long reads and using the strengths of each, we are able to obtain a significantly more complete genome.

The process and steps below are mostly identical for Nanopore or Pac Bio data and are used as an example. However, there are constant improvements in this field and new software is being developed and released regularly! Keep an eye out for developments!

## Contents
- Data
- Assembly
  - Short/long hybrid assembly
    - Spades Assembler
  - Long Read assembly
    - Wtdbg2
    - Flye
  - Long read correction and cleaning
    - FMLRC
    - Medaka
- Assessing genome completion
  - QUAST
  - BUSCO
- Exercises

## Data

A set of paired genomic illumina data and corresponding PacBio data for the E. Coli is in the ```~/Shared_folder/Day3/longRead``` folder.

Illumina Data
```
ERR022075_1.fastq.gz (~35m)
ERR022075_2.fastq.gz (~35m)
```
Pac Bio Reads 
```
ERR022075_PacBio.fastq -  Raw Pac Bio data (~50k)
```

We will be using a docker container from dockerhub named:
```
chrishah/fmlrc-wtdbg2-plus:v062022
chrishah/flye:2.9-b1778
```

## Assembly
Using only high quality short read data often results in smaller islands of non-overlapping data i.e. many shorter contigs, especially in complex or repetitive genome regions.

Long read data spans greater lengths but is largely of a lower accuracy and typically will fail to assemble very successfully as variation and errors are unable to agree.

There are a number of processes available to combine the data together and use the strengths of each to nullify the flaws. Here we will use the high quality Illumina data to “correct” or “polish” the long read scaffolds and then use those new high quality long reads as a trusted input into a standard assembly. 

## Short/long Hybrid Assembly
### Spades Assembler

Having already used Spades to perform a short-read assembly, it is straight forward to integrate long read data too. Spades has a number of options for working with long read data, and below is a small selection of the many useful parameters.
```
Usage: spades.py [options] -o <output_dir>
Input data:
-1 <filename>		file with forward paired-end reads
-2 <filename>		file with reverse paired-end reads
--sanger  <filename>	file with Sanger reads
--pacbio  <filename>	file with PacBio reads
--nanopore  <filename>	file with Nanopore reads
--trusted-contigs  <filename>	file with trusted contigs
--untrusted-contigs  <filename>	file with untrusted contigs

Pipeline options:
--only-error-correction	runs only read error correction (without assembling)
--only-assembler	runs only assembling (without read error correction)
--careful	tries to reduce number of mismatches and short 

Advanced options:
-t/--threads <int> 	number of threads [default: 16]
-m/--memory <int> 	RAM limit for SPAdes in Gb (terminates if exceeded) [default: 250]
--tmp-dir <dirname> 	directory for temporary files [default: <output_dir>/tmp]
-k <int,int,...>	comma-separated list of k-mer sizes (must be odd and less than 128) [default: 'auto']
--cov-cutoff <float>	coverage cutoff value (a positive float number, or 'auto', or 'off') [default: 'off']
```

A basic spades command looks like:
```
$ spades.py -o spades-default \
    -1 Sample1_1.fastq.gz -2 Sample1_2.fastq.gz \
	-t 4 -m 8 --only-assembler
```

## Long Read assembly
There are many different long-read assemblers. Two we will look at are wtdbg2 and Flye.
### Wtdbg2 
Source: https://github.com/ruanjue/wtdbg2

Wtdbg2 is a de novo sequence assembler for long noisy reads produced by PacBio or Oxford Nanopore Technologies (ONT). It assembles raw reads without error correction and then builds the consensus from intermediate assembly output and can assemble the human and even the 32Gb Axolotl genome at a speed tens of times faster than CANU and FALCON while producing contigs of comparable base accuracy.

During assembly, wtdbg2 chops reads into 1024bp segments, merges similar segments into a vertex and connects vertices based on the segment adjacency on reads. The resulting graph is called fuzzy Bruijn graph (FBG). It is similar to De Bruijn graph but permits mismatches/gaps and keeps read paths when collapsing k-mers. The use of FBG distinguishes wtdbg2 from the majority of long-read assemblers. 

Wtdbg2 has two key components: an assembler wtdbg2 and a consenser wtpoa-cns. Executable wtdbg2 assembles raw reads and generates the contig layout and edge sequences in a file "ctg.lay.gz". Executable wtpoa-cns takes this file as input and produces the final consensus in FASTA. A typical workflow looks like this:

The basic command is as such, assigning parameters for number of threads (-t), technology (-x), genome size (-g), and output format (-o), and a second step

Note for the -x parameter for choosing the technology
- ```ont``` for Oxford Nanopore
- ```rs``` for PacBio RSII
- ```sq``` for PacBio Sequel
- ```ccs``` for PacBio CCS reads

It is a two step process:
```
$ wtdbg2 -x rs -g 4.6m -t 4 -i ERR022075_PacBio.fastq -fo wtdgb2_raw
$ wtpoa-cns -t 4 -i wtdgb2_raw.ctg.lay.gz -fo wtdgb2_raw_assembly.fasta
```


### Flye 
Source: https://github.com/fenderglass/Flye

Flye is a de novo assembler for single molecule sequencing reads, such as those produced by PacBio and Oxford Nanopore Technologies. Flye is using repeat graph as a core data structure. In difference to de Bruijn graphs (which require exact k-mer matches), repeat graphs are built using approximate sequence matches, and can tolerate higher noise of SMS reads.

The edges of repeat graph represent genomic sequence, and nodes define the junctions. Each edge is classified into unique or repetitive. The genome traverses the graph (in an unknown way), so as each unique edge appears exactly once in this traversal. Repeat graphs reveal the repeat structure of the genome, which helps to reconstruct an optimal assembly.
 
It can accept either raw long read data or corrected. By default, expected error rates are <30% for raw, < 3% for corrected, and <1% for HiFi.

A simple default command for running flye with raw PacBio reads are:
```
$ flye --pacbio-raw ERR022075_PacBio.fastq --out-dir raw_flye_assemb --threads 4
or
$ flye --pacbio-corr ERR022075_PacBio_corrected.fasta --out-dir corr_flye_assemb --threads 4
```

There are a number of parameters which can improve assembly

## Long read correction and cleaning 
### FMLRC
Source: https://github.com/holtjma/fmlrc

FMLRC or FM-index Long Read Corrector, is a generic tool for performing correction of long read sequencing using short-read sequencing data. 

Given a BWT of the short-read sequencing data, FMLRC will build an FM-index and use that as an implicit de Bruijn graph. Each long read is then corrected independently by identifying low frequency k-mers in the long read and replacing them with the closest matching high frequency k-mers in the implicit de Bruijn graph. In contrast to other de Bruijn graph based implementations, FMLRC is not restricted to a particular k-mer size and instead uses a two pass method with both a short "k-mer" and a longer "K-mer". This allows FMLRC to correct through low complexity regions that are computational difficult for short k-mers. 

Note: there is a FMLRC2 version now available that works exactly the same but runs in half the time. However this is written in the rust language which is not universally available.

#### Prepare the data
We need to create an index of the short reads prior to correcting. It’s a complex command with some conversions in there, but the fundamentals are running ropebwt2 and fmlrc-convert. I don’t understand why this isn’t an automated step! 

__I’ve done this step as it is very memory intensive and to save time__

```
$ mkdir temp
$ awk "NR % 4 == 2" ERR022075_*.fastq | sort -T temp | tr NT TN | ropebwt2 -LR | tr NT TN | fmlrc-convert -f ecoli_illumina_msbwt.npy
```
ecoli_illumina_msbwt.npy is a fast-access index of your short read data.

#### Polish the reads
Now we can use fmlrc to correct the long-read data using the short-read index:
```
$ fmlrc2 -t 4 ecoli_illumina_msbwt.npy \
ERR022075_PacBio.fastq \
ERR022075_PacBio_fmlrc2_corr.fasta
```

__You now have some nice corrected long reads! Note that FMLRC2 outputs a fasta file (note quality data!)__

### Medaka
Source: https://github.com/nanoporetech/medaka 

Medaka is different in that it doesn’t combine short read data with long read, but is __a Nanopore specific data cleaning tool__ to create consensus sequences and from nanopore sequencing data. This is performed using neural networks applied to a pileup of individual sequencing reads against a draft assembly (to remove stochastic errors).

Their quote: "It outperforms graph-based methods operating on basecalled data, and can be competitive with state-of-the-art signal-based methods whilst being much faster."

You feed it a draft genome such as the output of wtdbg2, alongside your basecalls.
```
$ medaka_consensus -i raw_nano_data \
-d draft_genome \
-o medaka_output \
-t 4 \
-m r941_min_high_g303
```
Note that -m is the flowcell and sequencing chemistry that you used

## Assessing genome completion
Two tools that can be used to judge the assembly are Quast and BUSCO. Quast evaluates the number of contigs and size of assembly, whereas BUSCO looks for expected genes and whether they are complete or duplicated.

### QUAST
To run quast is very simple, with just the command and the genome of interest:
```
quast.py my_genome.fasta
```

However you can also run it on multiple together and do direct comparisons:
```
quast.py *fasta
```

To run Quast on one sample in a docker container:
```
singularity exec reslp/quast:5.0.2 \
      quast.py my_genome.fasta
```

Quast makes a simple text output of the results, but my favourite view is the html, so I recommend downloading the whole result folder and exploring it.

### BUSCO
We will test our genome assembly against a gammaproteobacteria database. This can be downloaded straight from busco if you know the code, or list them using:
```
busco --list-datasets
```

To run BUSCO on my_genome.fasta in a singularity container:
```
singularity exec docker://ezlabgva/busco:v5.3.2_cv1 \
    busco --in my_genome.fasta  \
          --out busco-my_genome \
          -l gammaproteobacteria_odb10 \
          --mode genome \
          -c 4 -f 
```
However if running busco on multiple genomes then it's easier to put it in a loop! (Here named ```busco_loop.sh```)
```
#!/bin/bash
for i in *fa*
do
    name=${i%.*}

    busco --in $i  \
          --out busco-$name \
          -l gammaproteobacteria_odb10 \
          --mode genome -c 4 -f
done
```
and run the script once:
```
singularity exec docker://ezlabgva/busco:v5.3.2_cv1 \
        ./busco_loop.sh
```


Inside the output folder there is a file named ```short_summary_busco-………txt``` with the relevant information i.e.:

        C:94.2%[S:94.0%,D:0.2%],F:2.2%,M:3.6%,n:452

        426     Complete BUSCOs (C)
        425     Complete and single-copy BUSCOs (S)
        1       Complete and duplicated BUSCOs (D)
        10      Fragmented BUSCOs (F)
        16      Missing BUSCOs (M)
        452     Total BUSCO groups searched


### MultiQC

MultiQC is a package for collecting lots of Quality Control outputs together to create a simple visual overview and is great for comparing multiple samples together. It can combine raw data, assemblies, and outputs across all sequence data. 

Source: https://multiqc.info

It can be run on the current directory without any further parameters (don't miss the ```.``` on the end of the line!)
```
singularity exec docker://multiqc/multiqc:latest multiqc .
```

Download the html that is generated to view the outputs


## Exercises
The processing of this data takes quite a long time! After all, it is Whole Genome _E. coli_ long and short read data! For time we can jump straight to the outputs of these steps, however you may want to generate these outputs yourself with more time
1. Copy the long read data to your home folder and use the docker before each command

Singularity generic command:
```
singularity exec docker://maintainer/containerName:version
```
Optional docker containers for different program running:
```
reslp/spades:3.15.3
passdan/fmlrc2:latest
staphb/wtdbg2:2.5
staphb/flye:latest
reslp/quast:5.0.2
ezlabgva/busco:v5.3.2_cv1
```
Note: All steps have been completed and are available in ```~/Shared_folder/Day3/longReads/Assemblies``` if you would prefer to look at the outputs rather than process the data. Step one (short read only spades assembly) specifically has been completed for you as it require huge resources.

### Assembly Exercises:
Estimated times are using 4 CPUs: 
1. [__Skip this for now, and do the long read assemblies!__] Using spades, assemble the 1% illumina short reads on their own and evaluate the resulting assembly (~15 minutes) 
2. Using wtdbg2 (~5 minutes) or flye  (~20 minutes), assemble the raw PacBio data alone
3. Use FMLRC2 (~2 minutes when using the pre-generated index) to correct the PacBio dataset with the Illumina short reads 
4. Use wtdbg2  (~5 minutes) or flye (~20 minutes) again to assemble the now high-quality long reads
5. [Optional] Use spades in hybrid mode to include illumina and corrected PacBio reads together (~20 minutes)

### Evaluating Assembly Exercises

In the folder ```longRead/Assemblies``` you have the outputs from the Assemblies above:
```
spades-illumina.fasta               Illumina alone Assembly
wtdbg_assembly.cns.fa               WTDBG2 alone Assembly
flye_raw_assembly.fasta		    Flye raw data Assembly
polished_wtdbg_assembly.cns.fa      WTDBG2 Assembly after FMLRC polishing
polished_flye_assembly.fasta	    Flye data Assembly after FMLRC polishing
spades-ill-and-pb.fasta             Spades default hybrid mode (Illumina & pacBio comibined)
spades-ill-and-pb-careful.fasta     Spades Hybrid mode with 'careful' parameter
```

1. Use quast.py with all assembled genome fasta files together to judge the assemblies
2. Use busco to determine the completeness of the assemblies
3. Determine which assembly is best!


<details>
  <summary>

  ### Extension: Annotate the genome

  </summary>

### Prokka
Source: https://github.com/tseemann/prokka

In this session we have been doing an assembly of a bacterial genome which is pretty simple to annotate (none of those annoying introns or big repeats).

We'll learn about annotating complex genomes later, but for bacteria then a simple prokka command does 99% of the work automatically! This is all down to the incredible work of Torsten Seemann and creating Prokka which is a great piece of software!


 1. Lets annotate the polished flye assembly
 ```
singularity exec docker://staphb/prokka:latest \
        prokka polished_flye_assembly.fasta
```
2. Read the gtf file that was generated. Try using grep to find the 16S gene, or another favourite gene of yours!
3. Later, we will look at visualising genome annotations!

</details>
