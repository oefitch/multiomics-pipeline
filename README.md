# multiomics-pipeline
Pipeline for processing RNA-seq and ATAC-seq data from fastq to bigwig

This pipeline is based on the ATAC-seq pipeline developed by John M. Gaspar and Aaron Kitzmiller https://github.com/harvardinformatics/ATAC-seq

## Overview

This pipeline is used to process raw reads in `.fastq` format to final `.bam` files that can be used for future analysis (e.g. read counting, etc.) as well as `.bigwig` files for visualization in IGV, and peak calling with MACS2 for `.narrowPeak` files. 

### Scripts 
So far I've written up the entire pipeline in one script, but I recommend you run this one step at a time for easier trouble shooting, I will update this soon. 

```
Fitch_full_ATAC_pipeline
```

## Step 1: Trim Reads

## Step 2: Map to Genome

## Step 3: Convert `.sam` to `.bam`

## Step 4: Remove Mitochondrial Reads

## Step 5: Remove PCR Duplicates

## Step 6: Remove non-unique reads

## Step 7: Index final `.bam` files

## Step 8: Convert `.bam` to `.bigwig`

## Step 9: Call Peaks with `MACS2`


