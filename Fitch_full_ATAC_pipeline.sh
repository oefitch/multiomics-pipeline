#bbduk
module load BBMap/39.01-GCC-12.3.0

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

sample_names="" #replace with your own sample names 

#### Trim Reads ####
#NOTE: with new package (faster than trimmomatic) in parallel (this is after re-naming read files from novogene)
bbduk.sh in={}_1.fq.gz in2={}_2.fq.gz out={}_pe_bbtrimmed_1.fq out2={}_pe_bbtrimmmed_2.fq outs={}_se_bbtrimmed.fq ref=adapter_contamination_sequences.fasta ktrim=r minlength=25

#### Map to Genome ####
#NOTE: need the genome file and index in your directory
bowtie2  --very-sensitive  -x gar_index  -1 {}_pe_bbtrimmed_1.fq -2 {}_pe_bbtrimmmed_2.fq -S {}.sam

#### Convert Sam to Bam ####

samtools view -u {}.sam | samtools sort -n -o {}_sorted.bam

#after this you should have your cleaned and sorted bam files so you can now clean them up further 

#### Remove Mitochondrial Reads #### 
#NOTE: sometimes you have to sort again after this
samtools view -h  {}_sorted.bam  |  python path/removeChrom.py - - MT  |  samtools view -b | samtools sort - > {}_sorted_mitout.bam

#### Remove PCR Duplicates #### 
#NOTE: sometimes you have to sort again after this
java -jar path/picard.jar MarkDuplicates I={}_sorted_mitout.bam O={}_sorted_mitout_dupout.bam M={}_dups.txt REMOVE_DUPLICATES=true

#### Remove non-unique reads ####
#NOTE: sometimes you have to sort again after this, and sometimes samtools index doesnt like sort with "-n" so run w/o
samtools view -b  -q 10  {}_sorted_mitout_dupout_sort.bam | samtools sort -n -o  {}_sorted_mitout_dupout_uniq.bam

#### Index final bam files ####
samtools index {}_sorted_mitout_dupout_uniq.bam

#NOTE: These next two can be run at the same time

module load Conda/3
cd path #navigate to the environment where conda is installed 
conda activate atac #activate the conda environment where you installed bamcoverage 
#### Convert bam to bigwig ####
cd path #navigate to the directory where you want to work
SIZE_in_BP= #size of your genome in base-pairs 
bamCoverage -b {}_sorted_mitout_dupout_uniq.bam -of bigwig --normalizeUsing BPM --ignoreForNormalization MT --effectiveGenomeSize $SIZE_in_BP -o {}_final.bw
conda deactivate


#### Call Peaks with MACS #### 
#NOTE: this will make a new folder, so make sure it's putting it where you want it
#NOTE: sometimes MACS doesn't want to run in parallel (idk why)
SIZE_in_BP= #size of your genome in base-pairs 
macs2 callpeak -t {}_sorted_mitout_dupout_uniq.bam -f BAMPE -g $SIZE_in_BP  -n {}_23 -B -q 0.05 -s 75 --call-summits --outdir MACS2_{}

