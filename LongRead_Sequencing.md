# Long-Read Sequencing: Data, QC and Assembly

Long-read sequencing has changed dramatically in the last few years. Where it was once a low-accuracy technology used only to *scaffold* short-read assemblies, modern long reads are accurate enough to assemble genomes on their own, frequently producing single-contig (or even telomere-to-telomere) chromosomes.

This page covers the two dominant platforms (Oxford Nanopore and PacBio), how to assess and filter the data, and a modern toolkit for aligning, assembling, polishing and evaluating long-read genomes. The older workflow of correcting noisy long reads with short reads (see [LongReadAssembly.md](LongReadAssembly_errorProneData.md)) is still useful for legacy or low-quality data, but is no longer the default for high-quality reads.

## Contents
- [Contents](#contents)
- [The technologies: ONT vs PacBio](#the-technologies-ont-vs-pacbio)
- [Read QC and filtering](#read-qc-and-filtering)
  - [NanoPlot](#nanoplot)
  - [seqkit stats](#seqkit-stats)
  - [Filtlong and chopper](#filtlong-and-chopper)
- [Reference available: Aligning long reads with minimap2](#reference-available-aligning-long-reads-with-minimap2)
- [De novo assembly](#de-novo-assembly)
  - [Flye](#flye)
  - [hifiasm](#hifiasm)
  - [Legacy assemblers](#legacy-assemblers)
- [Polishing and correction](#polishing-and-correction)
  - [Medaka](#medaka)
  - [Short-read polishing](#short-read-polishing)
- [Assessing the assembly](#assessing-the-assembly)
  - [QUAST](#quast)
  - [BUSCO and compleasm](#busco-and-compleasm)
  - [Mapping reads back](#mapping-reads-back)
- [Choosing a workflow](#choosing-a-workflow)
- [Data for these exercises](#data-for-these-exercises)
  - [Primary dataset: *E. coli* (HiFi + ONT, same DNA)](#primary-dataset-e-coli-hifi--ont-same-dna)
  - [Optional extension: latest ONT (R10.4.1)](#optional-extension-latest-ont-r1041)
- [Exercises](#exercises)
- [Container quick reference](#container-quick-reference)

## The technologies: ONT vs PacBio

Both platforms read single, native DNA molecules and produce reads thousands to millions of bases long, but they get there very differently.

**Oxford Nanopore (ONT)** passes DNA through a protein pore and measures changes in electrical current, which a basecaller (now **Dorado**) translates into sequence. Key points:
- Current chemistry is **R10.4.1**, basecalled in `fast`, `hac` (high accuracy) or `sup` (super-high accuracy) modes. Modern `sup` simplex reads reach Q20+ (≥99% accuracy), and **duplex** reads higher still.
- ONT excels at **ultra-long** reads (100kb–1Mb+), which span repeats that nothing else can, and runs on portable devices (MinION) through to high-throughput (PromethION).
- Accuracy depends heavily on the basecalling model, so always note which model produced your data — downstream tools need it.

**PacBio** reads a circularised molecule many times and builds a consensus:
- **HiFi** (also called CCS) reads are the modern standard: ~10–25 kb at ~99.9% accuracy (Q30+), produced on the Revio and Sequel II/IIe instruments. HiFi reads are accurate enough that assemblies usually need **no polishing at all**.
- **CLR** (Continuous Long Read) is the older, noisy (~85–90%) PacBio mode. You will still see it in legacy datasets but it is essentially retired.

| | ONT (R10.4.1, sup) | PacBio HiFi | Legacy (CLR / old ONT) |
| -- | -- | -- | -- |
| Typical read length | 10–100kb+ (ultra-long possible) | 10–25kb | 10–50kb |
| Raw accuracy | Q20+ (~99%) | Q30+ (~99.9%) | ~85–90% |
| Polishing needed? | Light (Medaka) | Usually none | Short-read correction/polishing |
| Best for | Repeats, structural variants, field work | Reference-grade assemblies, phasing | — (superseded) |

The practical implication: **match your tools and flags to your data type and quality**. Feeding HiFi reads into a "raw noisy read" mode wastes their accuracy, and feeding noisy reads into a HiFi mode will fail.

## Read QC and filtering

Standard short-read QC tools (FastQC) are not designed for long reads. Use long-read-aware tools to inspect read-length and quality distributions before assembling.

### NanoPlot
Source: https://github.com/wdecoster/NanoPlot

NanoPlot produces read-length histograms, quality plots and summary statistics (N50, total yield, median quality). `NanoComp` from the same author compares several runs side by side.

```
singularity exec docker://staphb/nanoplot:latest \
    NanoPlot --fastq reads.fastq.gz --outdir nanoplot_out --N50 --loglength
```
Download the HTML report to explore the read-length and quality distributions.

### seqkit stats
Source: https://github.com/shenwei356/seqkit

A fast way to get key numbers (count, total bases, N50, min/max/mean length) for any FASTA/FASTQ:
```
singularity exec docker://staphb/seqkit:latest \
    seqkit stats -a reads.fastq.gz
```

### Filtlong and chopper
Sources: https://github.com/rrwick/Filtlong  &  https://github.com/wdecoster/chopper

Long-read assemblies improve when you remove the shortest, lowest-quality reads. **Filtlong** filters by length and quality, optionally keeping only the best fraction of data:
```
# Keep reads >1kb and the best 90% by quality-weighted length
singularity exec docker://staphb/filtlong:latest \
    filtlong --min_length 1000 --keep_percent 90 reads.fastq.gz \
    | gzip > reads.filt.fastq.gz
```

**chopper** (Replacing NanoFilt) is a fast alternative that filters by quality and length in a stream:
```
zcat reads.fastq.gz \
    | chopper -q 10 -l 1000 \
    | gzip > reads.chopped.fastq.gz
```
Tip: don't over-filter. With high-yield ONT/HiFi runs, discarding the shortest reads helps; with low coverage, keep as much as you can.

## Reference available: Aligning long reads with minimap2
Source: https://github.com/lh3/minimap2

**minimap2** is the universal long-read aligner — used for mapping reads to a reference, mapping reads to an assembly (for polishing or QC), and all-vs-all overlap. The right preset (`-x`) matters:

| Preset | Use |
| -- | -- |
| `map-ont` | ONT reads vs reference |
| `lr:hq` | High-quality long reads (modern ONT sup / HiFi) vs reference |
| `map-hifi` | PacBio HiFi reads vs reference |
| `map-pb` | Legacy PacBio CLR reads vs reference |
| `asm5` / `asm10` / `asm20` | Assembly-to-assembly (≈5/10/20% divergence) |

minimap2 outputs SAM, which you pipe straight into **samtools** to make a sorted, indexed BAM:
```
# Map ONT reads to a reference and produce a sorted, indexed BAM
singularity exec docker://staphb/minimap2:latest \
    minimap2 -ax map-ont -t 4 reference.fasta reads.fastq.gz \
  | singularity exec docker://staphb/samtools:latest \
    samtools sort -@ 4 -o reads.sorted.bam -

singularity exec docker://staphb/samtools:latest samtools index reads.sorted.bam
```
For HiFi reads, swap the preset to `-ax map-hifi`. You can then load the BAM and reference into IGV to inspect coverage and structural variation.

## De novo assembly

### Flye
Source: https://github.com/mikolmogorov/Flye

Flye is a versatile long-read assembler that can work with all types of long read data (ONT and PacBio, raw or high-quality) and one of the most used. It builds a *repeat graph*, tolerating noise, and works well from bacteria to large eukaryotes.

| Flag | Data |
| -- | -- |
| `--nano-raw` | Older/`fast`-basecalled ONT (>5% error) |
| `--nano-hq` | Modern ONT `sup`/Q20 (R10.4.1) — use this for most ONT today |
| `--pacbio-hifi` | PacBio HiFi reads |
| `--pacbio-raw` | Legacy PacBio CLR |
| `--nano-corr` / `--pacbio-corr` | Reads already error-corrected with short reads |

```
# Modern ONT (sup-basecalled R10.4.1) E. coli assembly
singularity exec docker://staphb/flye:latest \
    flye --nano-hq reads.filt.fastq.gz --out-dir flye_ont --threads 4

# PacBio HiFi assembly
singularity exec docker://staphb/flye:latest \
    flye --pacbio-hifi hifi_reads.fastq.gz --out-dir flye_hifi --threads 4
```
The assembly is `flye_ont/assembly.fasta`, with a graph (`assembly_graph.gfa`) you can view in **Bandage**.

### hifiasm
Source: https://github.com/chhylp123/hifiasm

**hifiasm** is the de facto assembler for PacBio HiFi and is **haplotype-resolved** — for a diploid genome it can output both haplotypes rather than a collapsed consensus. Since release 0.21+ it also assembles modern ONT R10 reads with the `--ont` flag.

```
# PacBio HiFi
singularity exec docker://staphb/hifiasm:latest \
    hifiasm -o sample.asm -t 4 hifi_reads.fastq.gz

# Modern ONT R10 reads
singularity exec docker://staphb/hifiasm:latest \
    hifiasm --ont -o sample.asm -t 4 ont_reads.fastq.gz
```
hifiasm outputs assembly graphs in **GFA** format. The primary contigs are `sample.asm.bp.p_ctg.gfa`; convert to FASTA with:
```
awk '/^S/{print ">"$2; print $3}' sample.asm.bp.p_ctg.gfa > sample.p_ctg.fasta
```
For phased output, `*.bp.hap1.p_ctg.gfa` and `*.bp.hap2.p_ctg.gfa` hold the two haplotypes.

### Legacy assemblers
You may still meet these in older pipelines:
- **Canu** — accurate but slow; historically the standard, now usually replaced by Flye/hifiasm.
- **wtdbg2 (redbean)** — very fast on noisy reads but superseded; see [LongReadAssembly_errorProneData.md](LongReadAssembly_errorProneData.md).
- **Raven**, **NextDenovo** — fast ONT assemblers used in some pipelines.

## Polishing and correction

Polishing improves per-base accuracy of a draft assembly. **How much you need depends entirely on the input data**: HiFi assemblies are usually left as-is, modern ONT gets a light Medaka pass, and only noisy/legacy data needs short-read correction.

### Medaka
Source: https://github.com/nanoporetech/medaka

Medaka is ONT's neural-network consensus tool. It maps the reads back to the draft and corrects systematic basecalling errors. Medaka **auto-selects the model** from the basecaller information in the reads, so you usually no longer need to specify it by hand (you may see that in older tutorials):
```
singularity exec docker://staphb/medaka:latest \
    medaka_consensus -i reads.fastq.gz -d flye_ont/assembly.fasta \
                     -o medaka_out -t 4
```
If auto-selection can't determine the model, set it explicitly with `-m` using the basecaller-matched name, e.g. `-m r1041_e82_400bps_sup_v5.2.0` for R10.4.1 `sup` data. (The old `r941_*` models are for the retired R9.4.1 chemistry.) The polished assembly is `medaka_out/consensus.fasta`.

### Short-read polishing
For lower-accuracy assemblies, high-quality Illumina reads can correct residual errors:
- **NextPolish** / **POLCA** (from MaSuRCA) — current short-read polishers.
- **Pilon** — older and slower, still common in teaching material.
- **FMLRC2** — corrects the *reads* with a short-read index before assembly (the approach used in [LongReadAssembly.md](LongReadAssembly_errorProneData.md)).

Note: polishing HiFi assemblies with short reads can *introduce* errors and is generally discouraged. Only look to short-read polishing when the long-read data is noisy.

## Assessing the assembly

Two complementary questions: is the assembly *contiguous* (QUAST) and is it *complete* (BUSCO/compleasm)?

### QUAST
Source: https://github.com/ablab/quast

QUAST reports contig counts, total length, N50 and (with a reference) misassemblies. Run it on several assemblies at once to compare:
```
singularity exec docker://staphb/quast:latest \
    quast.py flye_ont/assembly.fasta flye_hifi/assembly.fasta sample.p_ctg.fasta \
             -o quast_compare
```
The HTML report (`quast_compare/report.html`) is the nicest way to compare assemblies side by side.

### BUSCO and compleasm
Sources: https://busco.ezlab.org  &  https://github.com/huangnengCSU/compleasm

**BUSCO** checks for the presence of near-universal single-copy orthologues expected in a lineage, reporting Complete / Duplicated / Fragmented / Missing. Datasets are now based on **OrthoDB v12** (`*_odb12`); list what's available with `busco --list-datasets`.
```
singularity exec docker://ezlabgva/busco:v5.8.2_cv1 \
    busco --in assembly.fasta \
          --out busco_assembly \
          -l bacteria_odb12 \
          --mode genome -c 4 -f
```
**compleasm** is a faster (3–14×) reimplementation that uses the same OrthoDB lineages via miniprot. It's a good drop-in when BUSCO is slow on large genomes:
```
compleasm run -a assembly.fasta -o compleasm_out -l bacteria_odb12 -t 4
```
Both report a completeness summary like `C:99.2%[S:98.9%,D:0.3%],F:0.2%,M:0.6%` — aim for high Complete and low Missing/Duplicated.

### Mapping reads back
A quick sanity check is to map the reads back to the assembly with minimap2 (above) and look at coverage:
```
singularity exec docker://staphb/samtools:latest samtools flagstat reads.sorted.bam
```
Even coverage with few unmapped reads suggests a well-represented assembly; large zero-coverage regions or pile-ups can flag misassemblies.

## Choosing a workflow

A reasonable default by data type:

- **PacBio HiFi** → `hifiasm` (or Flye `--pacbio-hifi`) → QUAST + BUSCO. No polishing.
- **Modern ONT (R10.4.1 sup)** → Filtlong → Flye `--nano-hq` → Medaka → QUAST + BUSCO.
- **ONT + HiFi together / large genomes** → hifiasm with HiFi plus ultra-long ONT for scaffolding.
- **Legacy noisy data (CLR / old ONT)** → short-read correction (FMLRC2) → Flye → short-read polish, as in [LongReadAssembly.md](LongReadAssembly_errorProneData.md).

## Data for these exercises

The data in `~/Share/Day3/longRead` is drawn from public repositories so you can reproduce or extend the exercises with your own machine after the course. There is the raw data, and also some pre-computed assemblies to allow you to jump straight to analysis of the outputs to save time.

### Primary dataset: *E. coli* (HiFi + ONT, same DNA)

The main set comes from Tvedte *et al.* 2021, [*Comparison of long-read sequencing technologies in interrogating bacteria and fly genomes*](https://doi.org/10.1093/g3journal/jkab083) (BioProject [PRJNA602597](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA602597)). The **same *E. coli* DNA** was sequenced with PacBio HiFi, PacBio CLR, ONT and Illumina, and there is a ground-truth reference assembly (`GCF_014117345.2`) — ideal for comparing HiFi vs ONT vs hybrid against a known answer.

> Note on chemistry: the **HiFi** reads are representative of current PacBio output, but the **ONT** reads are **R9.4.1** (one generation behind today's R10.4.1). They are fine for learning the assembly/polishing mechanics; for the *latest* ONT accuracy see the optional extension below.


### Optional extension: latest ONT (R10.4.1)

To show off **current** ONT accuracy, use Ryan Wick's [*Perfect bacterial genome* tutorial](https://github.com/rrwick/Perfect-bacterial-genome-tutorial/wiki/Sample-data) data — *S. aureus* JKD6159, R10.4 reads basecalled with a modern model, plus a reference. (Different organism, and no matched HiFi, so it's a standalone "modern ONT" comparison rather than part of the HiFi-vs-ONT contrast.)

## Exercises

Using the E. coli long-read data in `~/Share/Day3/longRead` choose **either ONT or PB** to analyse. You can go back and run the other sample(s) after the session.

For time, you may want to skip running the assembly in-session and use the pre-completed output from step 3 to save time.

1. [3 minutes] Run **NanoPlot** (or `seqkit stats`) on the long reads. Note: It does give some warnings about googe-chrome not being present but ignore it, and just download the output folder to look at on your own computer. Open the report first. What is the read N50 and median quality? Is this HiFi-grade or noisier data?
2. [3 minutes] Filter the reads with **Filtlong** (`--min_length 1000 --keep_percent 90`) and compare the stats before and after.
3. [15 minutes] Assemble the reads with **Flye**, choosing the flag that matches the data type. How many contigs do you get?
4. [7 minutes] [If ONT] Polish the Flye assembly with **Medaka** and note whether the model was auto-selected.
5. [Extension] Map the reads back to your assembly with **minimap2 + samtools**, and run `samtools flagstat`. What fraction of reads mapped?
6. [3 minutes each] Evaluate the assembly with **QUAST** and **BUSCO** (`bacteria_odb12`). Compare against the short-read-only or hybrid assemblies provided — what improved?

<details>
  <summary>Extension</summary>

7. Assemble the same data with a second tool (**hifiasm** for HiFi, or a second Flye mode) and use QUAST to compare contiguity and BUSCO to compare completeness. Which assembly would you take forward, and why?
8. Load the assembly graph (`assembly_graph.gfa` or the hifiasm `.gfa`) into **Bandage** [(**download here to your local computer**)](https://rrwick.github.io/Bandage/) and look at its structure — is it a clean circle, or tangled by repeats?
</details>

## Container quick reference
All tools above run from containers on dockerhub:
```
staphb/nanoplot:1.42.6
staphb/seqkit:latest
staphb/filtlong:latest
staphb/minimap2:latest
staphb/samtools:latest
staphb/flye:latest
staphb/hifiasm:latest
staphb/medaka:latest
staphb/quast:latest
ezlabgva/busco:v5.8.2_cv1
```
Generic pattern:
```
singularity exec docker://maintainer/containerName:version  command [options]
```
