#### Downloading the reference genome for Sceloporus undulatus v1.1 (Fence Lizard) 
## https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_019175285.1/

curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_019175285.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_019175285.1.zip" -H "Accept: application/zip"

## Unzipping the tarball
unzip GCF_019175285.1.zip 
