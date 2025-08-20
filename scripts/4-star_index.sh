#!/bin/bash
#SBATCH --job-name=star_index # Job name
#SBATCH -n 1
#SBATCH -N 8
#SBATCH -c 14
#SBATCH --mem=120G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rogini.runghen@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
## assumes star

 start=`date +%s`
 echo $HOSTNAME

 outpath="References"
 mkdir -p ${outpath}

 cd ${outpath}

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/920/845/GCF_016920845.1_GAculeatus_UGA_version5/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna.gz
gunzip GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna.gz
FASTA="../GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna"


wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/920/845/GCF_016920845.1_GAculeatus_UGA_version5/GCF_016920845.1_GAculeatus_UGA_version5_genomic.gtf.gz
gunzip GCF_016920845.1_GAculeatus_UGA_version5_genomic.gtf.gz
GTF="../GCF_016920845.1_GAculeatus_UGA_version5_genomic.gtf"

 mkdir stickleback_v5

 cd stickleback_v5

 module load star

 call="STAR
     --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir . \
     --genomeFastaFiles ${FASTA} \
     --sjdbGTFfile ${GTF} \
     --sjdbOverhang 100
     --genomeSAindexNbases 14"

 echo $call
 eval $call

 end=`date +%s`
 runtime=$((end-start))
 echo $runtime

