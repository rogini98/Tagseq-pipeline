#!/bin/bash
#SBATCH --job-name=star_alignment_davisExtended
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 14
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rogini.runghen@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo "Starting STAR alignment and counting on $(hostname) at $(date)"

# Load required modules
module load STAR
module load samtools

# Set up directories - EXACTLY matching your workflow
WORK_DIR="/labs/Bolnick/ROL/test-samples/04_read_count"
TRIM_DIR="/labs/Bolnick/ROL/test-samples/03_trim"
REF_DIR="$WORK_DIR/References/stickleback_v5"
OUT_DIR="$WORK_DIR/02-STAR_alignment"  # Your exact directory name
SAMPLES_FILE="$WORK_DIR/samples.txt"

echo "Working directory: $WORK_DIR"
echo "Trimmed files: $TRIM_DIR"
echo "Reference index: $REF_DIR"
echo "Output directory: $OUT_DIR"
echo "Samples file: $SAMPLES_FILE"
echo ""

# Create output directory
mkdir -p $OUT_DIR

# Check if required files/directories exist
if [ ! -d "$REF_DIR" ]; then
    echo "❌ ERROR: STAR reference directory not found: $REF_DIR"
    echo "Please run star_index.sh first"
    exit 1
fi

if [ ! -f "$SAMPLES_FILE" ]; then
    echo "❌ ERROR: Samples file not found: $SAMPLES_FILE"
    echo "Please create samples.txt with sample names"
    exit 1
fi

# Count total samples
total_samples=$(wc -l < $SAMPLES_FILE)
echo "Found $total_samples samples to process"
echo ""

# Initialize counters
processed=0
successful=0
failed=0

# Process each sample
while read sample; do
    ((processed++))
    
    echo "[$processed/$total_samples] Processing sample: $sample"
    
    # Define input and output paths - matching your structure
    input_file="$TRIM_DIR/${sample}.trim.gz"
    sample_out_dir="$OUT_DIR/$sample"  # our subdirectory structure
    
    echo "  Input file: $input_file"
    echo "  Output directory: $sample_out_dir"
    
    # Check if input file exists
    if [ ! -f "$input_file" ]; then
        echo "  ❌ ERROR: Input file not found: $input_file"
        ((failed++))
        echo ""
        continue
    fi
    
    # Create sample output directory
    mkdir -p $sample_out_dir
    
    # Run STAR alignment with gene counting
    echo "  Running STAR alignment and counting..."
    
    if STAR --runThreadN $SLURM_CPUS_PER_TASK \
            --genomeDir $REF_DIR \
            --readFilesIn $input_file \
            --readFilesCommand zcat \
            --outFileNamePrefix $sample_out_dir/${sample}_ \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            --outSAMunmapped Within \
            --outSAMattributes Standard; then
        
        echo "  ✅ STAR completed successfully"
        
        # Check if output files were created - YOUR file structure
        bam_file="$sample_out_dir/${sample}_Aligned.sortedByCoord.out.bam"
        count_file="$sample_out_dir/${sample}_ReadsPerGene.out.tab"
        log_file="$sample_out_dir/${sample}_Log.final.out"
        
        if [ -f "$bam_file" ] && [ -f "$count_file" ]; then
            echo "  ✅ Output files created:"
            echo "    BAM: $(basename $bam_file)"
            echo "    Counts: $(basename $count_file)"
            echo "    Log: $(basename $log_file)"
            
            # Index the BAM file
            echo "  Indexing BAM file..."
            samtools index $bam_file
            
            # Show alignment rate (like you check)
            if [ -f "$log_file" ]; then
                unique_rate=$(grep "Uniquely mapped reads %" $log_file | awk '{print $6}')
                echo "  📊 Uniquely mapped reads: $unique_rate"
            fi
            
            ((successful++))
            
        else
            echo "  ❌ ERROR: Expected output files not found"
            ((failed++))
        fi
        
    else
        echo "  ❌ ERROR: STAR alignment failed"
        ((failed++))
    fi
    
    echo ""
    
done < $SAMPLES_FILE

# Summary - matching our approach
echo "=== STAR ALIGNMENT SUMMARY ==="
echo "Total samples processed: $processed"
echo "Successful alignments: $successful"
echo "Failed alignments: $failed"
echo ""

if [ $failed -eq 0 ]; then
    echo "✅ All samples aligned successfully!"
    echo ""
    echo "Your results structure:"
    echo "  📁 $OUT_DIR/"
    echo "    📁 sample1/"
    echo "      📄 sample1_Aligned.sortedByCoord.out.bam"
    echo "      📄 sample1_ReadsPerGene.out.tab (3 columns: unstranded, +strand, -strand)"
    echo "      📄 sample1_Log.final.out"
    echo "    📁 sample2/"
    echo "      ..."
    echo ""
    echo "Count file format (as you noted):"
    echo "  Column 1: Gene ID"
    echo "  Column 2: Unstranded counts"
    echo "  Column 3: Forward strand counts (USE THIS for TagSeq)"
    echo "  Column 4: Reverse strand counts"
    echo ""
    echo "Next steps (your workflow):"
    echo "  f) Check alignment rates (should be high 80s as you got before)"
    echo "  11) Collect count files: find . -name '*ReadsPerGene.out.tab' -exec cp {} all_counts/ \\;"
    echo "  12) Ready for DESeq2!"
    
else
    echo "⚠️ Some samples failed alignment"
    echo "Check the individual error messages above"
fi

echo ""
echo "Alignment completed at $(date)"

# Show alignment statistics preview
echo ""
echo "=== ALIGNMENT STATISTICS PREVIEW ==="
echo "Uniquely mapped reads % for first 5 samples:"
find $OUT_DIR -name "*Log.final.out" -exec grep "Uniquely mapped reads %" {} + | head -5
