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

    ./homer/bin/findMotifsGenome.pl ${id} ./refgenome/GCF_019175285.1_SceUnd_v1.1_genomic.fna ${output_file} -p 8
done


#Motif Enrichment of control and heat stress treatments
    #you would want to configure a custom promoter set. It was attempted to do so with the reference genome using Homer's built-in
    #script for configuring a promoter set alongside the annotated reference genome as shown above. 
    #here, the promoter set configuration proved unsuccessful a number of times, so it required resorting to the closest default
    #promoter set to a lizard. 
./homer/bin/findMotifs.pl ./CatTsv/control_seq.autnorm.tsv chicken ./Motif_enrichment/control -p 12 -mset verbrates

./homer/bin/findMotifs.pl ./CatTsv/heatstress_seq.autonorm.tsv chicken ./Motif_enrichment/heatstress -p 12 -mset verbrates