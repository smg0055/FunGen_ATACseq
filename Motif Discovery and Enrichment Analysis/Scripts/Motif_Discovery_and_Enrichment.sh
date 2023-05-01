#!/bin/bash

#loading the homer module
module load homer

#configuring a custom genome in Homer
    #This step would be especially important for other function of homer such as annotating peaks
./homer/bin/loadGenome.pl -name Fence_Lizard -fasta ./refgenome/GCF_019175285.1_SceUnd_v1.1_genomic.fna -gtf ./refgenome/GCF_019175285.1_SceUnd_v1.1_genomic.gff -promoters FenceLizard

#Performing motif discovery on the bed files
ls ./peakcalls/*.bed | while read id;
do
    output_file=$(basename ${id} .bed)
    echo ${output_file}

    mkdir -p ${output_file}

    ./homer/bin/findMotifsGenome.pl ${id} ./refgenome/GCF_019175285.1_SceUnd_v1.1_genomic.fna ${output_file} -p 12
done

