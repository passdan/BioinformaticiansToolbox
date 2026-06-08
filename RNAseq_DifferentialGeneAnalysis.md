## Differential Gene Expression Analysis with DESeq2 

Following alignment to a reference genome OR de novo transcriptome annotation, we can perform differential gene analysis.

Differential Gene Analysis is contrasting the expression profile of the samples is typically done with one of two R packages: Deseq2 or EdgeR (the mac vs windows of the RNAseq fight), however a multitude of alternatives exist. These packages perform the normalization and statistical steps of contrasting samples as defined in a metadata file stating your experimental design (replicates, tissue type, treatment etc). The output here is a range of significant genes, ordinance and cluster analysis of sample similarity, and various quality control figures.

Following these three steps, there are an almost infinite number of tools and packages to look deeper into your data, find experimentally specific insights, and prior published data to contrast against.

### Contents


### Differential Gene Expression analysis using Deseq2 & Sartools

We will be using DESeq2 to perform our analysis but it is a massive package with a huge number of options, parameters and steps ([See the full manual here](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). 

Today for a rapid overview we will use a R wrapper package named SARTools (Statistical Analysis of RNA-Seq Tools) to automate the core steps and produce a graphical output. It also creates the deseq object for further analysis using standard deseq2 commands too. The full guide to SARTools can be found on the [github page](https://github.com/PF2-pasteur-fr/SARTools).  

I have also provided a ”mapping file” to annotate genes with common names for some easy investigation! This is just a reduced version of the GTF file you originally mapped with. You have been provided with it for today’s data but in the future you can generate it yourself with the following command (changing the “Homo_sapiens.GRCh38.100” part for your species):
```
$ echo "Id      Gene" > Homo_sapiens.GRCh38.100.map.txt
$ grep -P '\tgene\t' Homo_sapiens.GRCh38.100.gtf | cut -f 9 |cut -f1,3 -d';' | sed -r 's/^.*"(EN.*)".*"(.*)"$/\1\t\2/' >>  Homo_sapiens.GRCh38.100.map.txt
```

## Data

A set of nine Human neuronal differential RNAseq samples have been sequenced, consisting of 6x control samples (3x two different individuals) and 3x a deletion mutant of the 1q21.1 cytogenetic region of the human genome, and 3x of a duplication of this region. 

Deletions and alterations to this region has a range of impacts on neuropsychiatric disorders and is under active study. See the [omim page]( https://www.omim.org/entry/612474) and the [Wikipedia page](https://en.wikipedia.org/wiki/1q21.1_deletion_syndrome) is also remarkably informative! One potential relationship is to Autism Spectrum Disorders (ASD) and looking at the effect deletion/duplication of this region has on related genes can be of interest. A list of genes identified as being symptomatically-related have been put in today’s folder named ASD_list.txt from the [SFARI database](https://gene.sfari.org/database/gene-scoring/)

# Exercises
The data for this session is in: ```$ ~/Share/Day5/RNAseq-Analysis``` and we will be using scripts found in the R-scripts folder.

When doing the analysis we will generate comparisons between Control, Duplication, and Deletion of the gene region. We also will select from removed or non-removed duplicates. Pick one comparison to focus on for these questions:

## Performing a Differential Gene Expression analysis

You can either open the DESeq2/Sartools script in the RStudio server at today's IP address and ```:8787``` (i.e 123.4.5.67:8787) or run R on the command line once you have edited the script.

You can also open it directly on your own computer (not recommended during a workshop as it will take time to install some packages). 
1. Copy the folder RNAseq-Analysis to your local directory, and then enter it.
```
$ cp -r ~/Share/Day5/RNAseq-Analysis .
```
2. View and edit the Rscript to confirm the file locations and choose your parameters and testing options.

In Rstudio run line by line until you reach the ```check.parameters``` function.

<details>
  <summary>
  
  ### Using R on the command line

  </summary>

Lets first test the parameters by processing the first lines of the script:
```
$ head -n64 Sartools-template-deseq2.r > parameter_test.r
```
```
$ docker run --rm -u $(id -u):$(id -g) -v $(pwd):/in \
    -w /in chrishah/r-sartools-plus:2b95eaa \
    Rscript parameter_test.r
```
</details>

3. If everything comes back successful (no errors!), then we can run the full script. On Rstudio continue until you reach the end of the "tables" function (before starting heatmaps).

<details>
  <summary>
  
  ### Using R on the command line

  </summary>

```
$ docker run --rm -u $(id -u):$(id -g) -v $(pwd):/in \
    -w /in chrishah/r-sartools-plus:2b95eaa \
    Rscript Sartools-template-deseq2.r
```

</details>


4. This will have created a html output and a folder of tables. Download and inspect the outputs.

## Generating Heatmaps
Heatmaps are a common way of interrogating and displaying your results. And can be pretty!

We will continue the script in Rstudio, but you could also directly run the heatmap generating scripts with standard R.


<details>
  <summary>
  
  ### Using R on the command line

  </summary>

These each generate a pdf named ```Rplots.pdf```.
```
$ docker run --rm -u $(id -u):$(id -g) -v $(pwd):/in -w /in \
    chrishah/r-sartools-plus:2b95eaa \
    Rscript R-scripts/RNAseq-SimpleHeatmap.r
```
```
$ docker run --rm -u $(id -u):$(id -g) -v $(pwd):/in -w /in \
    chrishah/r-sartools-plus:2b95eaa \
    Rscript R-scripts/RNAseq-ComplexHeatmap.r
```

</details>

Firstly, pick a comparison you are most interested in (Deletion vs Duplication vs Control)
1. Open and look at the most differentially expressed genes in the outputted tables (on your own computer with excel/Google Sheets etc)
2. Use the provided script to generate a basic heatmap of the top 100 most DEGs
3. Use the complexHeatmap script to generate a more advanced heatmap which includes annotating your samples with the metadata.

<details>
  <summary>
  
  ### EXTENSION: Filter genes to known ASD relevant selection
  
  </summary>

4. Use the ASD list file to extract genes that are of known importance, and use that data as input to the heatmap creation. We can use grep to extract the gene counts of interest from our tables (Note using head to extract the header first, and >> to append the counts to the file)

Example code, where A & B is your choice of Deletion/Duplication/Control:
```
$ head -n1 AvsB.complete.txt > AvsB.ASD.txt
$ grep -f gene_list.txt AvsB.complete.txt >> AvsB.ASD.txt
```

You can now use that file for your heatmap generation using the same method as above


</details>

