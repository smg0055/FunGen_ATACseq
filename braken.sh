## Kraken2 doesn’t estimate species abundances, so we used Bracken:
## (Bayesian Reestimation of Abundance with Kraken) (https://ccb.jhu.edu/software/bracken/).

# Install Bracken using conda - similar to kraken2
conda create -n bracken

# Activate the braken virtual environment
conda activate bracken

#Install bracken
conda install -c bioconda bracken

# Build the database for braken for the first time from kraken.
# We used the default kmer length of 35 and read length of 100 as couple of our base samples had lower read length
bracken-build -d kraken_db -t 10 -k 35 -l 100

# Run Bracken on all the report files inside of the kraken_report
#!/bin/bash
source activate bracken
for i in ./*_report.txt
do
  filename=$(basename "$i")
  fname="${filename%_report.txt}"
  bracken -d /scratch/aubrrb/kraken_db -i $i -r 100 -t 10 -l S -o ${fname}_bracken.txt -w ${fname}_bracken_report.txt
done

### Combine bracken result files into a single output using the python script
### (https://github.com/jenniferlu717/Bracken/tree/master/analysis_scripts) before running the rest of the script

# Copy the bracken_combine_outputs.py from the github page of braken and paste on hpc
bracken_combine_outputs.py –-files *_bracken.txt -o lizard_bracken.txt

# Convert the bracken report files into .biom files for diversity analysis in phyloseq
# Install kraken-biom (https://github.com/smdabdoub/kraken-biom) and ran the following code to get a biom file.

# Install the kraken-biom
pip install git+http://github.com/smdabdoub/kraken-biom.git

# OTU classification 
#!/bin/bash
kraken-biom *_bracken_report.txt -o bracken.biom --fmt json
