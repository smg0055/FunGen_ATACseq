#!/bin/bash

# STEP 1: Read mapping to reference genome using bwa-mem 
### Run on terminal window ###
module load bwa/0.7.12

# Make a directory and 
mkdir indexed_reference
cd indexed_reference

## Copy the reference genome .fna file 
cp -r /scratch/aubclsb0317/Atac_seq_raw/Reference/ncbi_dataset/data/GCF_019175285.1/GCF_019175285.1_SceUnd_v1.1_genomic.fna ./

# Indexing the reference on terminal window
bwa index GCF_019175285.1_SceUnd_v1.1_genomic.fna 
 
## In queue ##

#STEP 2 Mapping with BWA-MEM
## Load BWA module within the ASC
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bwa/0.7.12

# Input fastq files
for fq1 in /scratch/aubclsb0317/Atac_seq_raw/cat_raw/Adapter_trim/*_1_paired.fastq;
do
    echo "working with file $fq1"
    base=$(basename $fq1 _1_paired.fastq)
    echo "base name is $base"
    fq1=/scratch/aubclsb0317/Atac_seq_raw/cat_raw/Adapter_trim/${base}_1_paired.fastq
    fq2=/scratch/aubclsb0317/Atac_seq_raw/cat_raw/Adapter_trim/${base}_2_paired.fastq
    bwa mem /scratch/aubclsb0317/Atac_seq_raw/indexed_reference/GCF_019175285.1_SceUnd_v1.1_genomic.fna $fq1 $fq2 > ${base}.LH3.sam
done

## Run the queue with 25 cores and 25 gb memory space
### Increase run time if needed

#Step 3: Sorting with coardinate order using picards

# Input = output from Step 2
# Run the queue with 15 cores, default memory and time
#!/bin/bash

module load picard/1.79

for file in *.LH3.sam;
do
    tag=${file%.LH3.sam};
    java -Xmx2g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/SortSam.jar I="$file" O=$tag.sort.sam SORT_ORDER=coordinate
done


#Step 4: SAM to BAM converison using Samtools
# Input = output from Step 3
# Run the queue with 15 cores, default memory and time

#!/bin/bash

module load samtools/1.11 

for file in *.sort.sam;
do
    tag=${file%.sort.sam};
    samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -Sb "$file" > $tag.bam
done

#step5: PCR duplicates removal

#!/bin/bash
module load picard/2.24.0
for file in *.bam;
do
tag=${file%.bam};
picard MarkDuplicates -I "$file" -O $tag.sort.bam -M $tag.dupstat.txt --REMOVE_DUPLICATES true --RE
MOVE_SEQUENCING_DUPLICATES true 
done 

