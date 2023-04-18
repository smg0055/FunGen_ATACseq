## Kraken2 doesn’t estimate species abundances, for which we used Bracken:
## Bayesian Reestimation of Abundance with Kraken (https://ccb.jhu.edu/software/bracken/)
# We installed Bracken using conda similar to that of kraken.

conda create -n bracken

# Activate the environment

To activate this environment, use
#
#     $ conda activate bracken
#
# To deactivate an active environment, use
#
#     $ conda deactivate

# Install bracken

conda install -c bioconda bracken


# Then we build the database for braken for the first time from kraken. 
# We used the default kmer length of 35 and read length of 100 as couple of our base samples had lower read length

bracken-build -d kraken_db -t 10 -k 35 -l 100

#Then we ran Bracken on all the report files inside kraken_report using the following command.

#!/bin/bash
source activate bracken
for i in ./*_report.txt
do
  filename=$(basename "$i")
  fname="${filename%_report.txt}"
  bracken -d /scratch/aubrrb/kraken_db -i $i -r 100 -t 10 -l S -o ${fname}_bracken.txt -w ${fname}_bracken_report.txt
done

# Then we combined bracken result files into a single output using the python script:
# (https://github.com/jenniferlu717/Bracken/tree/master/analysis_scripts) and ran the following code 
### First copy the bracken_combine_outputs.py from the github page of braken and paste on hpc.
## Then run the command: 

bracken_combine_outputs.py –-files *_bracken.txt -o lizard_bracken.txt

# Finally we converted the bracken report files into .biom files for diversity analysis in phyloseq. 
# We installed a program kraken-biom  following the instruction and ran the following code to get a biom file.
# (https://github.com/smdabdoub/kraken-biom)
# Install the kraken-biom

pip install git+http://github.com/smdabdoub/kraken-biom.git

# Analyze OTU diversity

#!/bin/bash
kraken-biom *_bracken_report.txt -o bracken.biom --fmt json

#The bracken.biom file is then used as an input in phyloseq for further analysis in R.
