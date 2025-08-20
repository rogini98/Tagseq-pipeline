#!/bin/bash

# Script to collect all gene count files from STAR output
# Usage: ./collect_counts.sh [output_directory_name]

echo "=== COLLECTING GENE COUNT FILES ==="
echo "Date: $(date)"
echo ""

# Set directories
WORK_DIR="/labs/Bolnick/ROL/test-samples/04_read_count"
STAR_DIR="$WORK_DIR/02-STAR_alignment"
COUNTS_DIR="$WORK_DIR/${1:-all_counts}"

echo "STAR alignment directory: $STAR_DIR"
echo "Output counts directory: $COUNTS_DIR"
echo ""

# Check if STAR output directory exists
if [ ! -d "$STAR_DIR" ]; then
    echo "❌ ERROR: STAR alignment directory not found: $STAR_DIR"
    echo "Please run davis_star.sh first"
    exit 1
fi

# Create output directory
mkdir -p "$COUNTS_DIR"

# Find and count gene count files
count_files=$(find "$STAR_DIR" -type f -name "*ReadsPerGene.out.tab" | wc -l)

echo "Found $count_files gene count files"

if [ $count_files -eq 0 ]; then
    echo "❌ ERROR: No gene count files found"
    echo "Check that STAR alignment completed successfully"
    exit 1
fi

echo ""
echo "Copying count files..."

# Copy all count files
copied=0
failed=0

while IFS= read -r -d '' file; do
    filename=$(basename "$file")
    echo "  Copying: $filename"
    
    if cp "$file" "$COUNTS_DIR/"; then
        ((copied++))
    else
        echo "    ❌ Failed to copy: $filename"
        ((failed++))
    fi
    
done < <(find "$STAR_DIR" -type f -name "*ReadsPerGene.out.tab" -print0)

echo ""
echo "=== COLLECTION SUMMARY ==="
echo "Total files found: $count_files"
echo "Successfully copied: $copied"
echo "Failed to copy: $failed"
echo ""

if [ $failed -eq 0 ]; then
    echo "✅ All count files collected successfully!"
    echo ""
    echo "Output directory: $COUNTS_DIR"
    echo "Files in directory:"
    ls -1 "$COUNTS_DIR" | head -10
    if [ $copied -gt 10 ]; then
        echo "... and $((copied - 10)) more files"
    fi
    echo ""
    echo "Count file format preview:"
    echo "Each file has 4 columns:"
    echo "  1. Gene ID"
    echo "  2. Unstranded counts"
    echo "  3. Forward strand counts" 
    echo "  4. Reverse strand counts"
    echo ""
    echo "For TagSeq, typically use column 3 (forward strand)"
    echo ""
    echo "Example from first file:"
    first_file=$(ls "$COUNTS_DIR"/*.tab | head -1)
    if [ -f "$first_file" ]; then
        echo "File: $(basename "$first_file")"
        head -10 "$first_file"
    fi
    echo ""
    echo "Ready for differential expression analysis!"
    
else
    echo "⚠️ Some files failed to copy - check permissions and disk space"
fi
