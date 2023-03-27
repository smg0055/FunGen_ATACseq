#! /bin/bash

###############################################
##  FASTQC 	Input: Paired end tarball files for each sample (*fq.gz)
##		Output: A folder for each file. The last line of this script will make a tarball of the output directory to bring back to your computer
##  For running the script on the Alabama Super Computer.
##	For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
## 	After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
##  then run the script by using "run_script [script name]"
## 	suggested paramenters are below to submit this script.
##		core: 1
##		time limit (HH:MM:SS): 04:00:00 
##		Memory: 4gb
##		run on dmc
###############################################


########## Load Modules
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load fastqc/0.10.1

##########  Define variables and make directories
## Replace the ID with Your specific information
## Make variable for your ASC ID so the directories are automatically made in YOUR directory (Example: MyID=aubclsb0313)
MyID=aubclsb0317

## Make variable that represents your Raw data directory (DD) and the pre or postcleaned status (CS).
## All the raw reads are present in the Atac_seq_raw directory
DD=/scratch/$MyID/Atac_seq_raw/cat_raw
CS=preclean_cat
 
##  Make the directories in SCRATCH for holding the raw data 
## -p tells it to make any upper level directories that are not there. Notice how this will also make the WD.
mkdir -p $DD
## move to the Data Directory
cd $DD

## Concatenation of files
## TS_01
cat ../TS_01_L1_1.fq.gz ../TS_01_L2_1.fq.gz > TS_01_1.fq.gz
cat ../TS_01_L1_2.fq.gz ../TS_01_L2_2.fq.gz > TS_01_2.fq.gz
## TS_02  
cat ../TS_02_L1_1.fq.gz ../TS_02_L2_1.fq.gz > TS_02_1.fq.gz
cat ../TS_02_L1_2.fq.gz ../TS_02_L2_2.fq.gz > TS_02_2.fq.gz
## TS_03  
cat ../TS_03_L1_1.fq.gz ../TS_03_L2_1.fq.gz > TS_03_1.fq.gz
cat ../TS_03_L1_2.fq.gz ../TS_03_L2_2.fq.gz > TS_03_2.fq.gz
## TS_04  
cat ../TS_04_L1_1.fq.gz ../TS_04_L2_1.fq.gz > TS_04_1.fq.gz
cat ../TS_04_L1_2.fq.gz ../TS_04_L2_2.fq.gz > TS_04_2.fq.gz
## TS_05  
cat ../TS_05_L1_1.fq.gz ../TS_05_L2_1.fq.gz >  TS_05_1.fq.gz
cat ../TS_05_L1_2.fq.gz ../TS_05_L2_2.fq.gz > TS_05_2.fq.gz
## TS_09  
cat ../TS_09_L1_1.fq.gz ../TS_09_L2_1.fq.gz > TS_09_1.fq.gz  
cat ../TS_09_L1_2.fq.gz ../TS_09_L2_2.fq.gz > TS_09_2.fq.gz
## TS_10
cat ../TS_10_L1_1.fq.gz ../TS_10_L2_1.fq.gz   > TS_10_1.fq.gz
cat ../TS_10_L1_2.fq.gz   TS_10_L2_2.fq.gz > TS_10_2.fq.gz
## TS_11  
cat ../TS_11_L1_1.fq.gz ../TS_11_L2_1.fq.gz   > TS_11_1.fq.gz
cat ../TS_11_L1_2.fq.gz ../TS_11_L2_2.fq.gz > TS_11_2.fq.gz
## TS_12
cat ../TS_12_L1_1.fq.gz ../TS_12_L2_1.fq.gz > TS_12_1.fq.gz
cat ../TS_12_L1_2.fq.gz ../TS_12_L2_2.fq.gz > TS_12_2.fq.gz
## TS_13
cat ../TS_13_L1_1.fq.gz ../TS_13_L2_1.fq.gz > TS_13_1.fq.gz
cat ../TS_13_L1_2.fq.gz ../TS_13_L2_2.fq.gz > TS_13_2.fq.gz


#### FASTQC to assess quality of the sequence data
## The output from this analysis is a folder of results and a zipped file of results
mkdir $DD/$CS
cd $DD
for file in ./*.fq.gz;
do 
fastqc $file --outdir=$DD/$CS
done
