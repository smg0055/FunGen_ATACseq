#! /bin/bash

###############################################
## FASTQC   
## Input: Downloade SRA files .fastq
##      Output: is a folder for each file. The last line of this script will make a tarball of the output directory to bring back to your computer
##
##  For running the script on the Alabama Super Computer.
##  For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
##  After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
##  then run the script by using "run_script [script name]"
##  suggested paramenters are below to submit this script.
##      core: 1
##      Time limit (HH:MM:SS): 04:00:00 
##      Memory: 4gb
##      Run on dmc
###############################################


########## Load Modules
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load trimmomatic/0.39
module load fastqc/0.10.1

##########  Define variables and make directories
## Replace the numbers in the brackets with Your specific information
## make variable for your ASC ID so the directories are automatically made in YOUR directory (Example: MyID=aubclsb0313)
MyID=aubclsb0317

## Make variable that represents YOUR working directory(WD) in scratch, your Raw data directory (DD) and the pre or postcleaned status (CS)
DD=/scratch/aubclsb0317/Atac_seq_raw/cat_raw               ## Example: DD=/scratch/$MyID/PracticeRNAseq/RawData
WD=/scratch/$MyID/Atac_seq_raw/cat_raw/Adapter_trim        ## Example: WD=/scratch/$MyID/PracticeRNAseq

 
##  Make the directories in SCRATCH for holding the raw data 
## -p tells it to make any upper level directories that are not there. This will also make the WD.
mkdir -p $WD

## Move to the Data Directory
cd $WD

for fq1 in $DD/*_1.fq.gz
    do
    echo "working with file $fq1"

    base=$(basename $fq1 _1.fq.gz)
    echo "base name is $base"

    fq1=$DD/${base}_1.fq.gz
    fq2=$DD/${base}_2.fq.gz

        ### Run Trimmomatic in paired end (PE) mode with 6 threads using phred 33 quality score format.
        java -jar /mnt/beegfs/home/aubmxa/.conda/envs/BioInfo_Tools/share/trimmomatic-0.39-1/trimmomatic.jar  PE -threads 6 -phred33 \
        $fq1 $fq2 \
        ${base}_1_paired.fastq ${base}_1_unpaired.fastq  ${base}_2_paired.fastq ${base}_2_unpaired.fastq \
        ILLUMINACLIP:AdaptersToTrim_All.fa:2:35:10 HEADCROP:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36

done
