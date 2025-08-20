#!/bin/bash
#SBATCH --job-name=generate_sample_list
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@institution.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo "Generating samples list: $(date)"

# Change to trimmed files directory
cd 03_trim

# Create samples.txt by extracting sample names from .trim.gz files
ls -1 *.trim.gz | sed -e 's/\.trim\.gz$//' > ../04_read_count/samples.txt

echo "Generated samples list with $(wc -l < ../samples.txt) samples"
echo "Sample list saved as samples.txt"

echo "First 10 samples:"
head ../04_read_count/samples.txt
