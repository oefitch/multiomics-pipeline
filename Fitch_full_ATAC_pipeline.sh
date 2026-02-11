#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########
 
#SBATCH --time=150:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=40G                    # memory required per node - amount of memory (in bytes)
#SBATCH --job-name ATAC_QC      # you can give your job a name for easier identification (same as -J)
 
 # Author: Olivia Fitch 
 # Date: 3/11/24
 #Title: Full ATAC-seq Pipeline in Parallel
 #NOTE: You probably just want to do this one step at a time and not run the whole thing as one because it's likely to get stuck somewhere
#ADDITIONAL NOTE: Can I make this so that the file path is defined at the top so that's all I need to change? Maybe if cd to the path at the beginning I don't need to define the path? 
########## Command Lines to Run ##########

module purge #purge previous modules, not usually necessary

#bbduk
module load BBMap/39.01-GCC-12.3.0

#parallel
module load parallel/20230722-GCCcore-12.3.0

#bowtie
module load Bowtie2/2.5.1-GCC-12.3.0

#samtools
module load GCC/6.4.0-2.28  OpenMPI/2.1.1
module load SAMtools/1.9

#picard 
module load picard/2.18.1-Java-1.8.0_152

#MACS
module load GCC/8.3.0  OpenMPI/3.1.4
module load MACS2/2.2.5-Python-3.7.4

 #path to location of input files 

sample_names="CF2 CF5 CF6 CF8 CF9 CR10 CR17 NE5 NE6 NE8 NE9" #replace with your own sample names 

#### Trim Reads ####
#NOTE: with new package (faster than trimmomatic) in parallel (this is after re-naming read files from novogene)
parallel 'bbduk.sh in={}_1.fq.gz in2={}_2.fq.gz out={}_pe_bbtrimmed_1.fq out2={}_pe_bbtrimmmed_2.fq outs={}_se_bbtrimmed.fq ref=adapter_contamination_sequences.fasta ktrim=r minlength=25' ::: $sample_names #this list is what will go in {} one at a time

#### Map to Genome ####
#NOTE: need the genome file and index in your directory
parallel 'bowtie2  --very-sensitive  -x gar_index  -1 /mnt/home/fitcholi/gar_tails_23/ATAC_samples_24/QC/{}_pe_bbtrimmed_1.fq -2 /mnt/home/fitcholi/gar_tails_23/ATAC_samples_24/QC/{}_pe_bbtrimmmed_2.fq -S {}.sam' ::: $sample_names

#### Convert Sam to Bam ####
parallel 'samtools view -u {}.sam | samtools sort -n -o {}_sorted.bam' $sample_names

#after this you should have your cleaned and sorted bam files so you can now clean them up further 

#### Remove Mitochondrial Reads #### 
#NOTE: sometimes you have to sort again after this
parallel 'samtools view -h  {}_sorted.bam  |  python /mnt/home/fitcholi/gar_tails_23/ATAC_samples_24/QC/removeChrom.py - - MT  |  samtools view -b | samtools sort - > {}_sorted_mitout.bam' ::: $sample_names

#### Remove PCR Duplicates #### 
#NOTE: sometimes you have to sort again after this
parallel 'java -jar /opt/software/picard/2.18.1-Java-1.8.0_152/picard.jar MarkDuplicates I={}_sorted_mitout.bam O={}_sorted_mitout_dupout.bam M={}_dups.txt REMOVE_DUPLICATES=true' ::: $sample_names

#### Remove non-unique reads ####
#NOTE: sometimes you have to sort again after this, and sometimes samtools index doesnt like sort with "-n" so run w/o
parallel 'samtools view -b  -q 10  {}_sorted_mitout_dupout_sort.bam | samtools sort -n -o  {}_sorted_mitout_dupout_uniq.bam' ::: $sample_names

#### Index final bam files ####
parallel 'samtools index {}_sorted_mitout_dupout_uniq.bam' ::: $sample_names

#NOTE: These next two can be run at the same time

module load Conda/3
cd path #navigate to the environment where conda is installed 
conda activate atac #activate the conda environment where you installed bamcoverage 
#### Convert bam to bigwig ####
cd path #navigate to the directory where you want to work
parallel 'bamCoverage -b {}_sorted_mitout_dupout_uniq.bam -of bigwig --normalizeUsing BPM --ignoreForNormalization MT --effectiveGenomeSize 345577542 -o {}_23.bw' ::: $sample_names
conda deactivate


#### Call Peaks with MACS #### 
#NOTE: this will make a new folder, so make sure it's putting it where you want it
#NOTE: sometimes MACS doesn't want to run in parallel (idk why)

parallel 'macs2 callpeak -t {}_sorted_mitout_dupout_uniq.bam -f BAMPE -g 345577542 -n {}_23 -B -q 0.05 -s 75 --call-summits --outdir MACS2_{}_23' ::: $sample_names



scontrol show job $SLURM_JOB_ID     ### write job information to output file
