#!/bin/bash

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bwa/0.7.12

bwa index GCF_019175285.1_SceUnd_v1.1_genomic.fna

#input fastq files
for fq1 in /scratch/aubclsb0334/ATACseq/index/*_1_paired.fastq;
do
    echo "working with file $fq1"
    base=$(basename $fq1 _1_paired.fastq)
    echo "base name is $base"
    fq1=/scratch/aubclsb0334/ATACseq/index/${base}_1_paired.fastq
    fq2=/scratch/aubclsb0334/ATACseq/index/${base}_2_paired.fastq
    bwa mem /scratch/aubclsb0334/ATACseq/index/GCF_019175285.1_SceUnd_v1.1_genomic.fna $fq1 $fq2 > ${base}.LH3.sam
done
