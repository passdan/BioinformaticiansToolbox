## Differential Gene Expression Analysis with DESeq2 

Following alignment to a reference genome OR de novo transcriptome annotation, we can perform differential gene analysis.

Differential Gene Analysis is contrasting the expression profile of the samples is typically done with one of two R packages: Deseq2 or EdgeR (the mac vs windows of the RNAseq fight), however a multitude of alternatives exist. These packages perform the normalization and statistical steps of contrasting samples as defined in a metadata file stating your experimental design (replicates, tissue type, treatment etc). The output here is a range of significant genes, ordinance and cluster analysis of sample similarity, and various quality control figures.

Following these three steps, there are an almost infinite number of tools and packages to look deeper into your data, find experimentally specific insights, and prior published data to contrast against.

### Contents


### Differential Gene Expression analysis using Deseq2 & Sartools

We will be using DESeq2 to perform our analysis but for a rapid overview we will use SARTools to automate the core steps and produce a graphical output, but it also creates the deseq objects for further analysis. The full guide to SARTools can be found on the github page PF2-pasteur-fr/SARTools: Statistical Analysis of RNA-Seq Tools.

I have also provided a ”mapping file” to annotate genes with common names for some easy investigation! This is just a reduced version of the GTF file you originally mapped with. You have been provided with it for today’s data but in the future you can generate it yourself with the following command (changing the “Homo_sapiens.GRCh38.100” part for your species):
```
$ echo "Id      Gene" > Homo_sapiens.GRCh38.100.map.txt
$ grep -P '\tgene\t' Homo_sapiens.GRCh38.100.gtf | cut -f 9 |cut -f1,3 -d';' | sed -r 's/^.*"(EN.*)".*"(.*)"$/\1\t\2/' >>  Homo_sapiens.GRCh38.100.map.txt
```

## Data

A set of nine Human neuronal differential RNAseq samples have been sequenced, consisting of 6x control samples (3x two different individuals) and 3x a deletion mutant of the 1q21.1 cytogenetic region of the human genome, and 3x of a duplication of this region. 

Deletions and alterations to this region has a range of impacts on neuropsychiatric disorders and is under active study. See the [omim page]( https://www.omim.org/entry/612474) and the [Wikipedia page](https://en.wikipedia.org/wiki/1q21.1_deletion_syndrome) is also remarkably informative! One potential relationship is to Autism Spectrum Disorders (ASD) and looking at the effect deletion/duplication of this region has on related genes can be of interest. A list of genes identified as being symptomatically-related have been put in today’s folder named ASD_list.txt from the [SFARI database](https://gene.sfari.org/database/gene-scoring/)

# Exercises
The data for this session is in: ```$ ~/Share/Day4``` and we will be using scripts found in the R-scripts folder

When doing the analysis we will generate comparisons between Control, Duplication, and Deletion of the gene region. We also will select from removed or non-removed duplicates. Pick one comparison to focus on for these questions:

## Performing a Differential Gene Expression analysis

We can use R on the command line to run the DESeq2/Sartools script (you could also run this in RStudio directly on your own computer). 
1. Copy the folder RNAseq-Analysis to your local directory, and then enter it.
```
$ cp -r ~/Share/Day4/RNAseq-Analysis .
```
2. View and edit the Rscript to confirm the file locations and choose your parameters and testing i.e. Which condition you are testing or whether to remove duplicates or not (we can run with default for now). 

Lets first test the parameters by processing the first lines of the script:
```
$ head -n64 Sartools-template-deseq2.r > parameter_test.r
```
```
$ docker run --rm -u $(id -u):$(id -g) -v $(pwd):/in \
    -w /in chrishah/r-sartools-plus:2b95eaa \
    Rscript parameter_test.r
```
If everything comes back successful (no errors!), then we can run the full script
```
$ docker run --rm -u $(id -u):$(id -g) -v $(pwd):/in \
    -w /in chrishah/r-sartools-plus:2b95eaa \
    Rscript Sartools-template-deseq2.r
```
3. This will have created a html output and a folder of tables. Download and inspect the outputs.

## Generating Heatmaps
We are going to use RStudio to investigate the data we have generated. If you already use RStudio on your own computer you can do that, or you can [connect to it on the server with these instructions](https://docs.google.com/document/d/1SlwJ1okSSg0TuIT8M9nIooK9bHIi60gEMljR4TJ5-rY/edit?usp=sharing).

You could also directly run the heatmap generating scripts with standard R (this generates a pdf named ```Rplots.pdf```) like so:
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

Firstly, pick a comparison you are most interested in (Deletion vs Duplication vs Control)
1. Open and look at the most differentially expressed genes in the outputted tables (on your own computer with excel/Google Sheets etc)
2. Use the provided script ```SimpleHeatmap.r``` to generate a basic heatmap of the top 100 most DEGs
3. Use the script ```ComplexHeatmap.r``` to generate a more advanced heatmap which includes annotating your samples with the metadata.

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

# Funtional and Network Analysis

There are a huge number of tools for performing further analysis based on your RNAseq Differentially Expressed Gene results. Here are just a few tools that can help explore your data further and gain novel insights.

#### Converting between common gene IDs
- [bioDBnet - Biological Database Network](https://biodbnet-abcc.ncifcrf.gov/db/db2db.php)
- [g:Convert Gene ID conversion](https://biit.cs.ut.ee/gprofiler/convert)

#### Whole gene set ontology annotation and enrichment
- [gProfiler g:Profiler – functional enrichment analysis](https://biit.cs.ut.ee/gprofiler/)
- [KEGG KEGG: Kyoto Encyclopedia of Genes and Genomes](https://www.genome.jp/kegg/)
  - Note that kegg uses some different geneIDs and can take some conversion to get it working properly!
- [NeVOmics GitHub - bioinfproject/bioinfo: NeVOmics](https://github.com/bioinfproject/bioinfo)
- [WebGestalt](http://www.webgestalt.org/)

#### Interaction Networks
- [StringDB](https://string-db.org/)
- [GOnet GOnet](https://tools.dice-database.org/GOnet/)
- Cytoscape 
  - [Bingo plugin  Cytoscape App Store - BiNGO](http://apps.cytoscape.org/apps/bingo)
  - [EnrichmentMap plugin Cytoscape App Store - EnrichmentMap](http://apps.cytoscape.org/apps/enrichmentmap)

#### Other visualisation tools
- [Morpheus Morpheus](https://software.broadinstitute.org/morpheus/)
- [diVenn DiVenn 2.0](https://divenn.tch.harvard.edu/)
- [Complex venn Interactive Venn Diagrams](http://www.interactivenn.net/)

#### Whole packages
- [IPA Ingenuity Pathway Analysis | QIAGEN Digital Insights](https://digitalinsights.qiagen.com/products-overview/discovery-insights-portfolio/analysis-and-visualization/qiagen-ipa/)
- [Shiny-Seq https://szenitha.shinyapps.io/shiny-seq3/](https://szenitha.shinyapps.io/shiny-seq3/)
- [BeavR BEAVR: A Browser-based tool for the Exploration And Visualization of RNAseq data](https://github.com/developerpiru/BEAVR#loading-your-data-into-beavr)

## Exercises

When doing the analysis we generated comparisons between Control, Duplication, and Deletion of the gene region. Pick one comparison to focus on for these questions. (Control vs Duplication, Control vs Deletion, or Deletion vs Duplication).

For these downstream analysis you can use the “annotated” csv file (containing just GeneID, Fold Change & pAdjusted) which will be the most easy to use (This was generated by just selecting the most useful columns from the larger data files)

1. Open your chosen comparison in excel/Sheets (or similar tool). You can choose from upregulated genes, downregulated genes, or the complete file (up and down regulated)
2. Filter the genes to an appropriate number. This could be by P value or fold change or a combination:
    - gProfiler:  ~ 300 genes
    - goNet:       ~300 genes
    - String:       ~100 genes
    - diVenn:       any really!

3. Explore your data using one (or more) of these 4 tools. Try to generate a figure that supports you writing a statement:
        This shows that [Deletion/Duplication] results in [increased/decreased]...

4. Put a figure in the shared document. It doesn’t have to be attractive, a basic screenshot is fine! We just want to see the amount of variation that we can generate within a short time! 

Shared Document: [RNAseq Downstream Analysis - Shared Results](https://docs.google.com/presentation/d/1ZJhtYOjzVINXjKvA-Kcbf0ib-SrIYVs1yq7_iuw1DYc/edit?usp=sharing)
