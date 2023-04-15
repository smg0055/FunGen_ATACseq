## Taxonomic classification of the kneaddata processed sequence files was done using Kraken2:
## (https://github.com/DerrickWood/Kraken22/wiki) using the standard database.                                                                                                             ##

# Install Kraken2 as a kraken2 environment using conda
conda create --yes -n kraken2

# Activate the kraken2 virtual environment
source activate kraken2

# Install kraken2 
conda install -c "bioconda/label/cf201901" kraken2

# Create a directory where all our database will be located.
mkdir kraken_db

# The latest version of the kraken database is downloaded from (https://benlangmead.github.io/aws-indexes/k2) and unzip it. 
# Download the standard database into the folder by running the command in queue, as it takes time and memory
# The database is more than 150GB so ownload the database in the working directory
#!/bin/bash
source activate kraken2
kraken2-build --standard –-threads 10 --db ./kraken_db

###############
### Kraken2 ###
###############

#!/bin/bash
for fq1 in ./*_paired_1.fastq
do
echo “working with file $fq1”
base=$(basename $fq1 _kneaddata_paired_1.fastq)
    echo "base name is $base"
    fq1=./${base}_paired_1.fastq
    fq2=./${base}_paired_2.fastq
  kraken2 --db /scratch/aubrrb/kraken_db \
  --confidence 0.1 \
  --threads 4 \
  --use-names \
  --output kraken_out/${base}_output.txt \
  --report kraken_reports/${base}_report.txt \
  --paired $fq1 $fq2 
done
