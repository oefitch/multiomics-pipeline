# multiomics-pipeline
Pipeline for processing RNA-seq and ATAC-seq data from fastq to bigwig

This file in **IN PROGRESS**

This pipeline is based on the ATAC-seq pipeline developed by John M. Gaspar and Aaron Kitzmiller https://github.com/harvardinformatics/ATAC-seq

## Overview

This pipeline is used to process raw reads in `.fastq` format to final `.bam` files that can be used for future analysis (e.g. read counting, etc.) as well as `.bigwig` files for visualization in IGV, and peak calling with MACS2 for `.narrowPeak` files. 

### Scripts 
So far I've written up the entire pipeline in one script, but I recommend you run this one step at a time for easier trouble shooting, I will update this soon. 

```
Fitch_full_ATAC_pipeline
```

## Step 1: Trim Reads

With next-generation sequencing (NGS), your samples contain fragments of genetic material (either DNA or RNA) and in order to sequence the fragments, you prepare a "lirary" by adding adapter sequences. These adapters have several puroses: 

1. Bind with primers - important for amplifying RNA/DNA fragments 
2. "Tag" your samples - unique adaparts are important for differenciating between sequence samples
3. Adapters contain regions that can bind to the sequence flow-cell - important for the sequencing process

Although adapters are important for the sequencing process, once you have your sequencing data back, you want to remove these adapters because they aren't biologically relevant. 
In some cases, you need to know your unique adapter sequences in order to trim them, but there are methods that can identify adapters without knowing the exact sequence. In my case, I knew the adapter sequences, so I trimmed those sequences specifically. 

For read trimming, I used to use `Trimmomatic` but I found that `bbduk` runs a lot faster for my purposes. 

To run `bbduk`, you will need to load/download the [latest version](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/)
If you know your adapter sequences, you can use the `ref` command to call the list of adapters. 

```
#For paired end reads
bbduk.sh in=READ_1 in2=READ_2 out=TRIMMED_1 out2=TRIMMED_2 outs=SINGLE_END ref=adapter_contamination_sequences.fasta ktrim=r minlength=25
```
After this importnat first step, you should have your trimmed reads! Congrats! Next let's map to the genome.

## Step 2: Map to Genome

## Step 3: Convert `.sam` to `.bam`

## Step 4: Remove Mitochondrial Reads

## Step 5: Remove PCR Duplicates

## Step 6: Remove non-unique reads

## Step 7: Index final `.bam` files

## Step 8: Convert `.bam` to `.bigwig`

## Step 9: Call Peaks with `MACS2`


