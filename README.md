# 16S Metagenomic Analysis Pipeline

## Overview

This pipeline provides a streamlined workflow for analyzing 16S rRNA gene sequencing data using DADA2 and Snakemake. It is designed to process data from Illumina NovaSeq platforms, automating the journey from raw sequencing reads to taxonomic classification and community analysis.

## Features
- This pipeline is inspired by a wonderful [ErnakovichLab tutorial][1] for Dada2 analysis of NovaSeq datasets

- Utilizes DADA2 for high-resolution sample inference
- Implements a Snakemake workflow for reproducibility and scalability
- Optimized for Illumina NovaSeq data
- Handles paired-end sequencing data
- Performs quality control using FastQC
- Generates comprehensive QC reports with MultiQC
- Performs quality filtering, denoising, and chimera removal
- Generates diagnostic plots along the way for troubleshooting and parameter tweaking
- Generates Amplicon Sequence Variants (ASVs)
- Assigns taxonomy using a [SILVA reference database][2] 

## Requirements

- Snakemake 8.16
- Conda
- R (version 4.0 or higher)
- DADA2 R package
- Other R dependencies (listed in the envs/*.yml files)

## Installation

1. clone this repository
2. Install all the environments in ```/workflow/envs``` directory
## Before running the pipeline
- In config file ("config/config.yml") specify relative path from Snakefile.py to your metadata.csv and the name of your experiment.
- In config file under Dada2_global specify how your forward and reverse reads are distinguished (default is _1.fq.gz and _2.fq.gz for forward and reverse read respectively).
- Place your raw reads into "raw_reads/<name_of_your_experiment>/" folder
- use snakemake profiles to adjust execution of HPC

## Running the pipeline

### if running on a local machine
```bash
cd 16S_metagenomics_pipeline/workflow
conda create -f envs/snakemake.yml
conda activate snakemake
snakemake -s Snakefile.py --profile profiles/default --cores <number_of_cores>
```

### if running on a HPC using slurm workload manager
```bash
cd 16S_metagenomics_pipeline/workflow
conda create -f envs/snakemake.yml
conda activate snakemake
snakemake -s Snakefile.py --workflow-profile profiles/16S_pipeline/
```

## Rules and configuration explanation:
- For more information about each step see 
- For detailed description of Dada2 arguments see [DADA2 documentation][3]

### Dada2_global:
- Specify how your forward and reverse reads are distinguished
- Specify whether to run the pipeline in multithread or not 

### Dada2_removeNs:
- Using Dada2 will remove sequences with ambiguous bases (Ns) in order to help cutadapt to remove primers in the next step
- Number of allowed Ns per sequence can be adjusted in config file

### Dada2_remove_primers:
- Using Cutadapt will remove primer sequences from the reads
- Can remove reads with less than minimal length after primer removal
- Deals with incorrect G calls at the end of the reads caused by two-color chemistry of NovaSeq
- Settings and primer sequences can be changed in config file

### Dada2_filter_and_trim:
- Trims the reads where their quality drops (truncLen and truncQ) and removes low-quality reads (maxEE). The parameters can be changed in the config file
- Specify at which base you want to trim forward and reverse read in the config file
- Plots of random 20 sequences and their quality before and after trimming will be generated. Inspect and adjust settings as necessary
- Be aware that excessive trimming will prevent read merging in the later steps. (If quality of the forward read is good enough do not trim it at all (trunc_length_FWD: 0), this will allow more trimming on the reverse read)

### Dada2_learn_error_rates:
- At this step Dada2 tries to distinguish sequencing errors from biological variability
- Unfortunately at the time of writing this, Dada2 has not been adjusted for the quality scores binning performed by NovaSeq
- [ErnakovichLab tutorial][1] compares multiple methods to overcome this limitation
- In this workflow I have implemented option 1, originally proposed by [JacobRPrice][4]
- Inspect error rate plots after running this script. These will be in ../results/Dada2/experiment/02_filter/error_rates. Make sure no sudden drops are present

### Dada2_sequence_table:
- Here Dada2 will merge and denoise sequences and generates sequence table
- In config file you can change pooling option. Default = FALSE
- From ErnakovichLab description:
    - FALSE = Sequence information is not shared between samples. Fast processing time, less sensitivity to rare taxa
    - TRUE = Sequence information from all samples is pooled together. Slow processing time, most sensitivity to rare taxa
    - pseudo = Sequence information is shared in a separate "prior" step. Intermediate processing time, intermediate sensitivity to rare taxa
- This step will also generate length_distribution_plot.png which shows the length distribution of the reads. You can remove reads that are too short or too long in the next step

### Download_silva_training_set:
- This will download SILVA database for the following step
- If you want to use different database change the link or location of your database in the config file

### Dada2_chimeras_taxonomy:
- At this step Dada2 will remove chimeras and assign taxonomy based on the SILVA database
- In the config file specify lower and upper cutoff for your read length. If you want to keep all of them set 0 and 1000 for lower and upper cutoff respectively

### Dada2_to_phyloseq:
- This step will generate phyloseq object for the downstream analysis

### Dada2_loss_plotting:
- This step will generate a plot in Dada2 folder which shows amount of reads filtered out at each step of the analysis
- Besides filtering, you should not lose too many reads at any step
    - If majority of reads were removed at merging step, truncLen at filtering step cut too much from the end of the reads
    - If majority of reads were removed as chimeric, revisit the Ns and primer removal steps

## Downstream analysis

### data_cleanup.ipynb
- Data filtering and quality control steps for microbiome data
- Removal of contaminants and outliers
- Validation of positive controls
- Rarefaction analysis and sequencing depth assessment
- Multiple data transformations (rarefaction, CLR, compositional)

### differential_abundance.ipynb
- Analysis of differential abundance of bacterial genera between different groups and timepoints
- Plotting of relative abundance changes over time
- Comparison of abundance patterns between treatment groups
- Statistical testing of abundance differences
- Visualization using various plot types (boxplots, violin plots)

### DESeq2.ipynb
- Differential abundance analysis using DESeq2 statistical framework
- Creation of volcano plots to visualize significant changes
- Multiple pairwise comparisons between treatment groups and timepoints
- Visualization of log2 fold changes

### beta_jaccard.ipynb
- Beta diversity analysis using Jaccard distance metric
- Ordination analysis (PCoA and NMDS)
- Statistical testing of group differences using PERMANOVA
- Visualization of sample clustering patterns

### alpha_beta.ipynb
- Combined alpha and beta diversity analyses
- Alpha diversity metrics (Shannon, Chao1)
- Beta diversity using Bray-Curtis dissimilarity
- Statistical testing of diversity differences
- Ordination analysis and visualization

### ancombc2.ipynb
- Differential abundance analysis using ANCOM-BC2 method
- Analysis of compositional data
- Multiple group comparisons
- Visualization using heatmaps
- Statistical testing with multiple comparison corrections

### clinical_correlations.ipynb
- Correlation analysis between microbiome composition and clinical parameters
- Spearman correlation calculations
- Multiple testing corrections
- Visualization using correlation heatmaps

## References

1. [DADA2 Ernakovich Lab](https://github.com/ErnakovichLab/dada2_ernakovichlab?tab=readme-ov-file)
2. [Zenodo Record](https://zenodo.org/records/4587955)
3. [DADA2 Tutorial](https://benjjneb.github.io/dada2/tutorial.html)
4. [DADA2 GitHub Issue #1307](https://github.com/benjjneb/dada2/issues/1307)