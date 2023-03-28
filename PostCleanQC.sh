#! /bin/bash

##############################################
## Purpose: To evaluate the quality of the data after trimming using FastQC
## FASTQC 	Input: *.fastq
##		Output: is a folder for each file. The last line of this script will make a tarball of the output directory to bring back to your computer
##
##  For running the script on the Alabama Super Computer:
##	 For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
## 	After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
## 	then run the script by using "run_script [script name]"
## 	suggested paramenters are below to submit this script.
##		core: 1
##		time limit (HH:MM:SS): 04:00:00 
##		Memory: 4gb
##		run on dmc
###############################################


### Load Modules
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load fastqc/0.10.1

### Define variables and make directories
## Replace the numbers in the brackets with Your specific information
## make variable for your ASC ID so the directories are automatically made in YOUR directory (Example: MyID=aubtss)
MyID=aubclsb0317          

## Make variable that represents YOUR working directory(WD) in scratch, your Raw data directory (DD) and the pre or postcleaned status (CS).
## All the raw reads are present in the Atac_seq_raw directory
DD=/scratch/$MyID/Atac_seq_raw/cat_raw/Adapter_trim	    ## Example: DD=/scratch/$MyID/PracticeRNAseq/RawData
WD=/scratch/$MyID/Atac_seq_raw/cat_raw 		               ## Example: WD=/scratch/$MyID/PracticeRNAseq
CS=Postclean
 
## Make the directories in SCRATCH for holding the raw data 
## -p tells it to make any upper level directories that are not there. This will also make the WD.
mkdir -p $CS
## Move to the Data Directory
cd $DD

#### FASTQC to assess quality of the sequence data
## FastQC: run on each of the data files to check the quality of the data
## The output from this analysis is a folder of results and a zipped file of results
mkdir $DD/$CS
for file in ./*.fq.gz;
do 
fastqc $file --outdir=$DD/$CS
done
