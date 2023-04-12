# STEP 5: DIfferential Peak calling 

## Input = output from step 4
## Run in in queue with 25 core, default memory and time

#!/bin/bash

#  Load the module for MACS
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load macs/3.0.0b1

for file in *.bam;
do
    tag=${file%.bam};
    macs3 callpeak -f BAMPE -t $tag.bam -g hs -n "$file" -B -q 0.01 --outdir Macs_result;
done
