#!/bin/bash
#SBATCH --job-name=gunzip_files
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=100G
#SBATCH --time 12:00:00
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.name@institute.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

for F in *.fastq.gz; do
    gunzip "$F"
done
