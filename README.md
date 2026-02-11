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

For ATAC-seq, it's important to remove mitochondrial reads during processing because, much like how rRNA reads can overwhelm your RNA-seq analysis, mitochondrial reads can make up a signifcant number of your ATAC reads, so it is best to remove them.  

For RNA-seq, you may not want to remove mitochondrial reads because you're interested in mitochondrial genes, you can also remove mitochondrail genes or rRNA later in your RNA analysis process. There is conflicting opinions on whether to remove mitochondrial reads from your `.bam` files in RNA-seq, but I think that treating your RNA and ATAC reads with the same processing is important for multiomics work, so I choose to remove mitochondrial reads at this step for RNA-seq, but this is up to your discretion. 

You can use the [python script](https://github.com/harvardinformatics/ATAC-seq/blob/master/atacseq/removeChrom.py) convieniently posted to Github by John M. Gaspar and Aaron Kitzmiller to remove mitochondrial reads. 

```
#For ATAC-seq
samtools view -h TRIMMED_MAPPED.bam  | python $PATH/removeChrom.py - - $NAME_OF_MITO_GENOME  |  samtools view -b | samtools sort - > TRIMMED_MAPPED_MITOUT.bam
```
> **NOTE:** Make sure you know the name of the mitochondrial genome in your species and update `$NAME_OF_MITO_GENOME`

You have now removed mitochondrial reads for your ATAC samples and you're ready to continue to the next processing steps. 

## Step 5: Remove PCR Duplicates

Next, we will remove PCR duplicates. These are identical reads that are a result of the sequencing process and should be removed so your samples are as biologically accurate as possible. 

To remove PCR duplicates, we can use `picard` tools, make sure you have the [lastest version](https://broadinstitute.github.io/picard/) downloaded. 

```
#Remove PCR duplicates
java -jar $PATH/picard.jar MarkDuplicates I=TRIMMED_MAPPED_MITOUT.bam O=TRIMMED_MAPPED_MITOUT_DUPOUT.bam M=duplicates.txt REMOVE_DUPLICATES=true
```
Great! Now we are ready to move on to the next processing step: removing non-unique reads. 

## Step 6: Remove non-uniquely mapped reads

This step is important after `bowtie2` processing, using another mapping software you may not have this issue. But `bowtie` may map a read to mutliple places in the genome and give it a low quality score. You can remove those with a low quality score using `samtools`. 

```
#Remove loq-qulaity mapped reads and sort again
# -n sorts reads alpha-numerically
samtools view -b  -q 10  TRIMMED_MAPPED_MITOUT_DUPOUT_UNIQ.bam | samtools sort -n -o  FINAL.bam
```

Alright, you've made it to the end! You now have your final `.bam` files. Congratulations! 

## Step 7: Index final `.bam` files

Before we move on to converting the `.bam` files to other file types for visualization, you'll want to index the `.bam` files. 

This is a simple `samtools` command 

```
#This will generate an index for your .bam files
samtools index FINAL.bam
```
Great, now you have your index! These next two steps can be performed in unison, they are not dependent on each other. 

## Step 8: Convert `.bam` to `.bigwig`

Now that we have our cleaned `.bam` files, we want to move on to visualizing our samples in IGV. The best way to visualize these samples is to make "bigwig" `.bw` files. That are mean for visualization so you don't have to haul around a giant `.bam` file to your genome viewers. 

We can use the `bamCoverage` tool to generate `.bw` files. Download the [latest version](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) of `bamCoverage` from `deepTools`

```
SIZE_in_BP= #size of your genome in base-pairs 
bamCoverage -b FINAL.bam -of bigwig --normalizeUsing BPM --ignoreForNormalization MT --effectiveGenomeSize $SIZE_in_BP -o FINAL.bw
```
Here, I'm normalizing using BPM, but there are different methods for normalizing ATAC-seq data for visualization. `bamCoverage` offers RPKM, CPM, BPM, RPGC, or none. It may be best to use the tool you use for RNA-seq analysis, or the same tool you use for downsteam ATAC analysis. My advice is to look into the best normalization practices for the type of biological system you're assessing. Is it a disease trial? Is it across development? 
> ** NOTE: ** I will post my guide for normalization for developmental data soon! Be on the lookout! 

## Step 9: Call Peaks with `macs2` 

In the previous step, we made `.bw` files to visualize our ATAC-seq read pile up in IGV. For ATAC-seq, we can determine if those piled-up reads are significant and call them "peaks". There are several software available for calling ATAC peaks, I use MACS2, but you could also use MACS3, HMMR-ATAC. The software listed here are the most popular and widely used.

```
SIZE_in_BP= #size of your genome in base-pairs 
macs2 callpeak -t FINAL.bam -f BAMPE -g $SIZE_in_BP  -n {}_23 -B -q 0.05 -s 75 --call-summits --outdir MACS2_FINAL
```
> **NOTE:** The `-f` parameter is set to `BAMPE` to recognize bam files with paired-end sequencing. 

##Final Thoughts 

Defintely check out the ATAC-seq pipeline developed by John M. Gaspar and Aaron Kitzmiller https://github.com/harvardinformatics/ATAC-seq for more tools. 

Next I will publish my analysis pipeline, so be on the look out for that!

Please reach out to me if you have any questions or comments. 
Thanks for reading and happy research! 
