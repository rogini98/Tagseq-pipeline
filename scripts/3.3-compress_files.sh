#!/bin/bash
#SBATCH --job-name=compress_trim
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo "Starting compression at $(date)"
## change to your directory here
cd /labs/Bolnick/ROL/test-samples/03_trim

echo "Files to compress: $(ls *.trim 2>/dev/null | wc -l)"
gzip *.trim
echo "Compressed files: $(ls *.trim.gz 2>/dev/null | wc -l)"

echo "Compression completed at $(date)"
