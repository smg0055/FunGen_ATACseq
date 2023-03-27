#! /bin/bash

###############################################
## Purpose: Use FASTQC to evaluate the quality of the data: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
##      Input Data: NA
## 			Output: Downloaded read files, R1 and R2 files for each sample if paired-end data (FASTQ)
## FASTQC 	InPut: Downloade SRA files .fastq
##		Output: is a folder for each file. The last line of this script will make a tarball of the output directory to bring back to your computer
##
##  For running the script on the Alabama Super Computer.
##	For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
## 	After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
## 	then run the script by using "run_script [script name]"
## 	suggested paramenters are below to submit this script.
##		Core: 1
##		Time limit (HH:MM:SS): 04:00:00 
##		Memory: 4gb
##		Run on dmc
###############################################


## Load Modules
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load sra
module load fastqc/0.10.1

## Navigate to the directory containing the raw FASTQ files
cd /scratch/aubclsb0334/ATACSEQ/TS_01

## Run FastQC on all the files in the directory
fastqc *.fq.gz 
