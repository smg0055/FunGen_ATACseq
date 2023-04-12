1.	Taxonomic classification of the kneaddata processed sequence files was done using Kraken2 (https://github.com/DerrickWood/Kraken22/wiki), using 
the standard database.
#First, we installed Kraken2 as a kraken2 environment using conda.

conda create --yes -n kraken2

#Then activate the environment with 
source activate kraken2

#Install kraken2 
conda install -c "bioconda/label/cf201901" kraken2


#For kraken database, we can downloaded the latest database from https://benlangmead.github.io/aws-indexes/k2 and unzip it. The database is more 
than 150GB so we can simply download the database in our working directory. We first created a directory where all our database will be located.

mkdir kraken_db

#Then we download the standard database into the folder by running the command in queue, as it takes time and memory

#!/bin/bash
source activate kraken2
kraken2-build --standard –-threads 10 --db ./kraken_db


##have an error 
ERROR:


Downloading nucleotide est accession to taxon map...rsync: link_stat "/taxonomy/accession2taxid/nucl_est.accession2taxid.gz" (in pub) failed: No such file 
or directory (2)
rsync error: some files/attrs were not transferred (see previous errors) (code 23) at main.c(1658) [Receiver=3.1.2]


###solve

1. Install the newest version with this command

wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.2.zip

2. unzip v2.1.2.zip 

 3. cd kraken2-2.1.2/

 mkdir ~/kraken2

4.  ./install_kraken2.sh ~/kraken2

5.export PATH="$PATH:/home/aubaclsb0317"


now download the database

#Then we download the standard database into the folder by running the command in queue, as it takes time and memory

#!/bin/bash

~/kraken2/scripts/kraken2-build --standard –-threads 10 --db ./kraken_db



#erorr2: Untarring taxonomy tree data... done.
rsync_from_ncbi.pl: unexpected FTP path (new server?) for https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/762/265/GC
F_000762265.1_ASM76226v1

change Within the 'rsync_from_ncbi.pl' script

if (! ($full_path =~ s#^ftp://${qm_server}${qm_server_path}/##)) {  die "$PROG: unexpected FTP path (new server?) for $ftp_path\n"; }

to:

if (! ($full_path =~ s#^https://${qm_server}${qm_server_path}/##)) { die "$PROG: unexpected FTP path (new server?) for $ftp_path\n"; }



#error3: All files processed, cleaning up extra sequence files... done, library complete.
Masking low-complexity regions of downloaded library...which: no dustmasker in (/scratch/aubclsb0317/kraken/kraken2-2
.1.2/scripts:/scratch/aubclsb0317/kraken/kraken2-2.1.2/scripts:/opt/asn/apps/lua_5.3.4/bin:/usr/lib64/qt-3.3/bin:/usr
/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/asn/bin:/apps/scripts:/opt/asn/apps/swami/bin:.:.)
Unable to find dustmasker in path, can't mask low-complexity sequences


#added --no masking
#to the command

#eror 4:
Creating sequence ID to taxonomy ID map (step 1)...
/scratch/aubclsb0317/kraken/kraken2-2.1.2/scripts/build_kraken2_db.sh: line 84: lookup_accession_numbers: command not
 found

#installed again to get all the scripts into the path.

###then conda installed into

#First, we installed Kraken2 as a kraken2 environment using conda.

conda create --yes -n kraken2

#Then activate the environment with 
source activate kraken2

#Install kraken2 
conda install -c "bioconda/label/cf201901" kraken2


#For kraken database, we can downloaded the latest database from https://benlangmead.github.io/aws-indexes/k2 and unzip it. The database is more 
than 150GB so we can simply download the database in our working directory. We first created a directory where all our database will be located.

mkdir kraken_db


#commadn
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


##got the database from rishi

