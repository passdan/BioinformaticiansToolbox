# Physalia Winter School in Bioinformatics - 2024
## Instructors: Dr Daniel Pass & Dr Christoph Hahn

---

## OVERVIEW
This course introduces participants to Next Generation Sequencing biology, both understanding the concepts and handling the data. We will cover a broad range of software and analysis from quality assessment of sequencing runs, through assembling and annotating small genomes, RNAseq and differential gene expression, and phylogenetics with NGS data. 

Primarily focussing on the most popular Illumina data, we will also look at the different requirements and opportunities utilising long read data (Nanopore/PacBio). This course will also include use of the linux command line and docker for bioinformatic analysis as modern standard approaches.
 
## FORMAT
The course is structured in modules over five days with each session including an introductory lecture, class discussion of key concepts, and practical hands-on sessions.

We also provide short self-study preparation materials on basic linux usage if you have not worked in linux before, to have a foundation to aid you in the course.

Live sessions involve both mirroring exercises with the instructor to demonstrate a skill as well as applying these skills on your own to complete individual exercises. After and during each exercise, interpretation of results will be discussed as a group.

## Schedule
 
**Note**: These teaching materials are designed to be understood alongside the files found on the teaching server presentations shared during the course.

If you're not on a course and some steps don't make sense you may want to go sign up for a course!

### Day 1
- Accessing the bioinformatics cloud Image 
- [Review of Linux basics](Introduction_to_Linux.md)
- [NGS data and Quality Control - How well did my sequencing run work?](NGS_QualityControl.md)
- [Linux methods for multiple sample handling](Looping_in_Linux.md)

### Day 2
- [Using Docker & Singularity for reproducible bioinformatics](Using_Containers.md)
- [Short-read Genome Assembly]()
- [Assessing Assembly Quality]()

### Day 3
- [Assembly with long read data - Nanopore/PacBio & hybrid assemblies](LongReadAssembly.md)
- [Genome Annotation]()
- [Genome Visualisation]()

### Day 4
- [Phylogenetics and Phylogenomics]()
- Bioinformatic analysis with pipelines
  - [Snakemake]()
  - [Nextflow](A_breif_view_on_nextflow.md)
- [Job Queueing & SLURM](Queueing_with_SLURM.md)

### Day 5
- [RNAseq Processing and Transcriptomics](RNAseq_Processing.md)
- [Differential Gene Analysis](RNAseq_DifferentialGeneAnalysis.md)
- [Downstream RNAseq Interpretation]()
