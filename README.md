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

With next-generation sequencing (NGS), your samples contain fragments of genetic material (either DNA or RNA) and in order to sequence the fragments, you prepare a "library" by adding adapter sequences. These adapters have several purposes: 

1. Bind with primers - important for amplifying RNA/DNA fragments 
2. "Tag" your samples - unique adaparts are important for differenciating between sequence samples
3. Adapters contain regions that can bind to the sequence flow-cell - important for the sequencing process

Although adapters are important for the sequencing process, once you have your sequencing data back, you want to remove these adapters because they aren't biologically relevant. 
In some cases, you need to know your unique adapter sequences in order to trim them, but there are methods that can identify adapters without knowing the exact sequence. In my case, I knew the adapter sequences, so I trimmed those sequences specifically. 

For read trimming, I used to use `Trimmomatic` but I found that `bbduk` runs a lot faster for my purposes. 

To run `bbduk`, you will need to load/download the [latest version](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/)
If you know your adapter sequences, you can use the `ref` parameter to call the list of adapters. The `ktrim` parameter allows you to set the trimming, from the left, right, or both directions, and `minlength` allows you to set the minimum length that sequences can be trimmed. I set the minimum length to 25 so anything 25bp or less would be removed. Consider the biological relevance of removing small reads and how stringent you would like to set this parameter for your system. 

```
#For paired end reads
bbduk.sh in=READ_1.fq in2=READ_2.fq out=TRIMMED_1.fq out2=TRIMMED_2.fq outs=SINGLE_END ref=adapter_contamination_sequences.fasta ktrim=r minlength=25
```
> **TIP:** Consider running [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) after this step to ensure that your adapter contamination is removed. 

After this important first step, you should have your trimmed reads! Congrats! Next let's map to the genome.

## Step 2: Map to Genome

Now that you have your trimmed reads, you can map them to your genome of interest. You may only have one option of a genome to work with, but if you're working with a common model organism such as human, mouse, zebrafish, etc. you'll have genome options. Consider if you want to map to the most recent genome, or to the genome with the most resources available, in some cases these may be the same genome, and if they are then you're in luck! Be sure to look into the sequencing methods used to build that genome, because not all genomes are created equally. Ultimately, you should consider for your purposes what the best genome version is, and in most cases the most recent genome will be best. If you'd like to know more about choosing the right genome, feel free to reach out to me!

To map to the genome, I used `bowtie2` which is a very widely used genome mapping software. 

Make sure you load/download the [latest version](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) of `bowtie2` before you get started. 

```
#For ATAC-seq
bowtie2  --very-sensitive  -x genome_index  -1 TRIMMED_1.fq -2 TRIMMED_2.fq -S {}.sam
```
For bowtie, this first parameter I have listed here `--very-sensitive` is important to consider when mapping RNA-seq vs. ATAC-seq data. 

There are 4 pre-set parameters in bowtie: `--very-sensitive`, `--sensitive`, `--fast`, `--very-fast`, but you can set all of the parameters individually if you'd like. But I've found that the pre-set parameters work well. When considering sensitivity vs. speed, more sensitive means it will search the genome for longer with more accuracy to find the best match, vs. increasing speed will search the genome for less time for the best match, meaning you may lose some accuracy. I think it is always better to be accurate than fast, but maybe with your system you know it will be accurate so you can be fast, but this is the trade-off to consider. 

With ATAC-seq, I've found that using the `--very-sensitive` parameter works to get around 50-70% mapped reads, which is acceptable for working with non-model organisms, but if you're working with common models, you may want a higher map percentage up to 80-90%. 

For RNA-seq, it's recommended to have 80-90% mapped reads, and you can increase map percentage using the preset `--local` parameter can help with this. Using the `--local` parameter will allow for non-exact sequences to allign. I find this is best for RNA samples because RNA sequences may be slightly different from the genome due to RNA processing, so using the `--local` parameter will find best mapped sequences that may ignore a few of the bases towards the end of the sequence known as "soft-clipping". For RNA-seq, I used the `--very-sensitve-local` preset. You could use this parameter for ATAC-seq too if you want to increase your map percentage. 

```
#For RNA-seq
bowtie2  --very-sensitive-local  -x genome_index  -1 TRIMMED_1.fq -2 TRIMMED_2.fq -S TRIMMED_MAPPED.sam
```
Alight, after this you should now have your mapped reads! Awesome!

## Step 3: Convert `.sam` to `.bam`

Up to now, you were working with `.fastq` sequence files, but now we're working with `.sam` files. And we're going to immediately convert those `.sam` files to `.bam` files. We're converting to `.bam` because, although `.sam` files are human readable, we actually want a file that is smaller and easier for the computer to work with, not for us. So, we will convert our `.sam` files to `.bam` files so that the computer can read them faster. 

To covert a `.sam` file to a `.bam` file, we are going to use an extremely useful software called `samtools`. I would say that `samtools` might be one of the most useful and versitle genomic tools, and learning and being familiar with `samtools` will be important for the rest of our pipeline. 

As always, make sure your running the [latest version](https://www.htslib.org/) of `samtools` and lets run a simple line of code 

```
#Make your sam file a bam file
samtools view -u TRIMMED_MAPPED.sam | samtools sort -n -o TRIMMED_MAPPED.bam
```
Now you have your mapped, trimmed, and sorted bam files! Congrats! You could stop here, but it is best to sort your files a bit further to use them for your downstream analysis. 

## Step 4: Remove Mitochondrial Reads

## Step 5: Remove PCR Duplicates

## Step 6: Remove non-unique reads

## Step 7: Index final `.bam` files

## Step 8: Convert `.bam` to `.bigwig`

## Step 9: Call Peaks with `MACS2`


