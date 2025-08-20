#!/bin/bash
#SBATCH --job-name=merge_fastqs
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH --time 12:00:00
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rogini.runghen@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo "Starting lane merge on $(hostname) at $(date)"

# Create output directory 
mkdir -p ../01_merged_files

# Extract unique sample names (everything up to _L00)
# This gets the base sample name before the lane identifier
for file in *_L00*_R1_001.fastq.gz; do
    # Extract sample name by removing _L00X_R1_001.fastq.gz suffix
    sample=$(echo "$file" | sed 's/_L00[0-9]_R1_001\.fastq\.gz$//')
    echo "$sample"
done | sort | uniq > sample_list.txt

echo "Found samples to merge lanes for:"
cat sample_list.txt
echo ""

# Process each sample - merge lanes L001 and L002
while read sample; do
    echo "Processing sample: $sample"
    
    # Check what lane files exist for this sample
    lane_files=(${sample}_L00*_R1_001.fastq.gz)
    if [ -e "${lane_files[0]}" ]; then
        echo "  Found lane files: ${lane_files[@]}"
        echo "  Merging lanes L001 and L002..."
        cat ${sample}_L00*_R1_001.fastq.gz > "../01_merged_files/${sample}_merged_lanes_R1_001.fastq.gz"
        echo "    Created: ${sample}_merged_lanes_R1_001.fastq.gz"
        
        # Count how many lanes were merged
        num_files=$(ls ${sample}_L00*_R1_001.fastq.gz | wc -l)
        echo "    Merged $num_files lane files"
        
        # Verify the output file was created and has content
        if [ -s "../01_merged_files/${sample}_merged_lanes_R1_001.fastq.gz" ]; then
            echo "    Success: Output file created and has content"
        else
            echo "    Error: Output file is empty or not created"
        fi
    else
        echo "    Warning: No lane files found for $sample"
    fi
    echo ""
    
done < sample_list.txt

# Clean up
rm sample_list.txt

echo "Lane merge completed at $(date)"
echo ""
echo "Output files created:"
ls -la ../01_merged_files/

# Show file sizes to verify concatenation worked
echo ""
echo "File sizes (to verify lane concatenation worked):"
du -h ../01_merged_files/*
