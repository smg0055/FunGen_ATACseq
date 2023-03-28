#! /bin/bash

###############################################
## Purpose: To trim sequencing adapters and low quality regions from the raw read data and using FASTQC to evaluate the quality of the data
## Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
##       Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
##                   Downloaded read files, R1 and R2 files for each sample if paired-end data (FASTQ)
##       Output: Trimmed R1 & R2 paired and unpaired reads (FASTQ)       
## FASTQC output is a folder for each file and a tarball of the output directory to bring back to your computer
##       Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
##                   Downloaded read files, R1 and R2 files for each sample if paired-end data (FASTQ)
##       Output: Trimmed R1 & R2 paired and unpaired reads (FASTQ)       
## After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
## then run the script by using "run_script [script name]"
## Suggested paramenters are below to submit this script to the ASC:
##      Core: 6
##      Time limit (HH:MM:SS): 02:00:00  (Increase if needed)
##      Memory: 12gb
##      Run on dmc
###############################################

# Load modules
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load trimmomatic/0.39
module load fastqc/0.10.1

## Make variables for your ASC ID so the directories are automatically made in YOUR directory (Example: aubclsb0313)
MyID=aubclsb0334

# Variables: raw data directory (DD), working directory(WD), cleaned status (CS), name of file containing the adpaters.
WD=/scratch/$MyID/ATACSEQ                                 ## Example: WD=/scratch/$MyID/PracticeRNAseq
DD=/scratch/$MyID/ATACSEQ/RawData                         ## Example: DD=/scratch/$MyID/PracticeRNAseq/RawData
CD=/scratch/$MyID/ATACSEQ/CleanData                       ## Example: CD=/scratch/$MyID/PracticeRNAseq/CleanData
CS=PostCleanQuality
adapters=AdaptersToTrim_All.fa  ## This is a fasta file that has a list of adapters commonly used in NGS sequencing. 
## You may need to edit this for other projects based on how your libraries were made in order to trim with the correct adapters.

## Make the directories to hold the Cleaned Data files, and the directory to hold the results for assessing quality of the cleaned data.
mkdir $CD
mkdir $WD/$CS

#### Trimmomatic
## Move to your Raw Data Directory
cd $DD

### Make list of file names to Trim
ls | grep ".fastq" |cut -d "_" -f 1 | sort | uniq > list

### Copy over the list of Sequencing Adapters that you want Trimmomatic to look for (including default adapters).
cp /home/$MyID/class_shared/AdaptersToTrim_All.fa . 

### Run a while loop to process through the names in the list and Trim them with the Trimmomatic Code
while read i
do

        ### Run Trimmomatic in paired end (PE) mode with 6 threads using phred 33 quality score format. 
        java -jar /mnt/beegfs/home/aubmxa/.conda/envs/BioInfo_Tools/share/trimmomatic-0.39-1/trimmomatic.jar  PE -threads 6 -phred33 \
        "$i"_1.fastq "$i"_2.fastq  \
        $CD/"$i"_1_paired.fastq $CD/"$i"_1_unpaired.fastq  $CD/"$i"_2_paired.fastq $CD/"$i"_2_unpaired.fastq \
        ILLUMINACLIP:AdaptersToTrim_All.fa:2:35:10 HEADCROP:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36
        
                ## Trim read for quality when quality drops below Q30 and remove sequences shorter than 36 bp
                ## PE for paired end phred-score-type  R1-Infile   R2-Infile  R1-Paired-outfile R1-unpaired-outfile R-Paired-outfile R2-unpaired-outfile  Trimming paramenter

#### FASTQC - assess quality of all cleaned sequence data files.
## The output from this analysis is both a folder of results and a zipped file of results

fastqc $CD/"$i"_1_paired.fastq --outdir=$WD/$CS
fastqc $CD/"$i"_2_paired.fastq --outdir=$WD/$CS

done<list			# This is the end of the loop

## Compress your results files from the Quality Assessment by FastQC and move to the directory with the cleaned data
cd $WD/$CS

#### Tarball the directory containing the FASTQC results so we can easily bring it back to your local computer.
## When finished, use scp or rsync to bring the .gz file to your computer and open the .html file to evaluate.
tar cvzf $CS.tar.gz $WD/$CS/*
