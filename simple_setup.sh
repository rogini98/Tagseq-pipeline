#!/bin/bash

# Simple TagSeq setup script
# Run this from the directory containing your FASTQ files

set -e

echo "=== Simple TagSeq Setup ==="
echo "Creating directory structure..."

# Create directories
mkdir -p 00_raw_fastqs
mkdir -p 01_merged_files  
mkdir -p 02_fastqc
mkdir -p 03_trim
mkdir -p 04_read_count
mkdir -p 04_read_count/References

echo "Directories created:"
ls -la | grep "^d" | grep -E "(00_|01_|02_|03_|04_)"

echo
echo "Finding FASTQ files..."

# Count files before moving
TOTAL_FILES=$(find . -name "*_R[12]_001.fastq.gz" -type f | wc -l)
echo "Found $TOTAL_FILES FASTQ files to organize"

if [ $TOTAL_FILES -eq 0 ]; then
    echo "ERROR: No FASTQ files found with pattern *_R[12]_001.fastq.gz"
    echo "Checking what files exist:"
    find . -name "*.fastq*" -type f | head -10
    exit 1
fi

echo
echo "Moving FASTQ files to 00_raw_fastqs/..."

# Move files
find . -name "*_R[12]_001.fastq.gz" -type f -exec mv {} 00_raw_fastqs/ \;

# Verify
MOVED_FILES=$(ls 00_raw_fastqs/*.fastq.gz 2>/dev/null | wc -l)
echo "Successfully moved $MOVED_FILES files to 00_raw_fastqs/"

echo
echo "Setup complete!"
echo "Next step: cd 00_raw_fastqs && sbatch ../concatenate_files.sh"
echo
echo "Directory structure:"
ls -la | grep "^d" | grep -E "(00_|01_|02_|03_|04_)"
echo
echo "Files in 00_raw_fastqs:"
ls 00_raw_fastqs/ | head -5
echo "..."
