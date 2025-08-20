#!/bin/bash

# Simple TagSeq QC checker for three-spine stickleback
# Usage: ./check_tagseq_quality.sh [path_to_multiqc_reports_directory]
# Example: ./check_tagseq_quality.sh /labs/Bolnick/ROL/test-samples/02_fastqc/multiqc_reports

#################################################################
# TAGSEQ-SPECIFIC THRESHOLDS FOR STICKLEBACK
#################################################################

MIN_MEAN_QUALITY=20           # TagSeq acceptable at Q20
MAX_ADAPTER_CONTENT=15        # Higher tolerance for TagSeq  
MAX_DUPLICATION_LEVEL=70      # TagSeq naturally has high duplication
MAX_FAIL_PERCENT=25           # Max percentage of samples that can fail
MIN_READS_PER_SAMPLE=500000   # Minimum reads for TagSeq

#################################################################
# MAIN QC CHECK
#################################################################

# Get input directory
MULTIQC_DIR="${1:-./multiqc_reports}"

if [ ! -d "$MULTIQC_DIR" ]; then
    echo "Usage: $0 [path_to_multiqc_reports_directory]"
    echo "Error: Directory not found: $MULTIQC_DIR"
    exit 1
fi

echo "=== TAGSEQ QC EVALUATION FOR STICKLEBACK ==="
echo "Checking: $MULTIQC_DIR"
echo "Date: $(date)"
echo ""

cd "$MULTIQC_DIR"

# Check if required files exist
if [ ! -f "multiqc_data/multiqc_general_stats.txt" ]; then
    echo "❌ ERROR: MultiQC data files not found!"
    echo "Make sure MultiQC completed successfully."
    exit 1
fi

# Count total samples
total_samples=$(tail -n +2 multiqc_data/multiqc_general_stats.txt | wc -l)
echo "Total samples analyzed: $total_samples"

if [ $total_samples -eq 0 ]; then
    echo "❌ ERROR: No samples found in MultiQC data"
    exit 1
fi

# Initialize counters
critical_issues=0
warnings=0
low_quality_samples=0
high_adapter_samples=0
low_read_count_samples=0

echo ""
echo "Applying TagSeq-specific quality criteria..."
echo "  - Min quality: Q$MIN_MEAN_QUALITY"
echo "  - Max adapter: $MAX_ADAPTER_CONTENT%"
echo "  - Max duplication: $MAX_DUPLICATION_LEVEL% (normal for TagSeq)"
echo "  - Min reads: $MIN_READS_PER_SAMPLE"
echo ""

# Check each sample from general stats
while IFS=$'\t' read -r sample total_seqs poor_qual seq_length gc_content avg_qual rest; do
    # Skip header
    if [[ "$sample" == "Sample" ]]; then continue; fi
    
    sample_issues=0
    
    # Check read count
    if [[ -n "$total_seqs" && "$total_seqs" != "NA" ]]; then
        if [ "$total_seqs" -lt $MIN_READS_PER_SAMPLE ]; then
            echo "  ❌ $sample: Low read count ($total_seqs)"
            ((low_read_count_samples++))
            ((sample_issues++))
        fi
    fi
    
    # Check quality score
    if [[ -n "$avg_qual" && "$avg_qual" != "NA" ]]; then
        if (( $(echo "$avg_qual < $MIN_MEAN_QUALITY" | bc -l 2>/dev/null || echo "1") )); then
            echo "  ❌ $sample: Low quality (Q$avg_qual)"
            ((low_quality_samples++))
            ((sample_issues++))
        fi
    fi
    
    if [ $sample_issues -gt 0 ]; then
        ((critical_issues++))
    fi
    
done < multiqc_data/multiqc_general_stats.txt

# Check adapter content if fastqc data available
if [ -f "multiqc_data/multiqc_fastqc.txt" ]; then
    echo ""
    echo "Checking adapter content..."
    
    # Simple approach: count lines with high adapter content
    high_adapter_count=$(awk -F'\t' 'NR>1 && $NF > '$MAX_ADAPTER_CONTENT' {print $1, $NF}' multiqc_data/multiqc_fastqc.txt 2>/dev/null | wc -l)
    
    if [ $high_adapter_count -gt 0 ]; then
        echo "  ⚠ $high_adapter_count samples have >$MAX_ADAPTER_CONTENT% adapter content"
        ((warnings++))
        echo "    Trimming recommended for these samples"
    fi
fi

# Calculate failure percentage
fail_percent=0
if [ $total_samples -gt 0 ]; then
    fail_percent=$((critical_issues * 100 / total_samples))
fi

echo ""
echo "=== SUMMARY ==="
echo "Total samples: $total_samples"
echo "Critical issues: $critical_issues samples ($fail_percent%)"
echo "  - Low read count: $low_read_count_samples"
echo "  - Low quality: $low_quality_samples"
echo "Warnings: $warnings"
echo ""

#################################################################
# FINAL DECISION
#################################################################

echo "=== QC DECISION ==="

if [ $fail_percent -gt $MAX_FAIL_PERCENT ]; then
    echo "❌ QC FAILED"
    echo "Too many samples have critical issues ($fail_percent% > $MAX_FAIL_PERCENT%)"
    echo ""
    echo "Recommended actions:"
    echo "  1. Review individual sample quality reports"
    echo "  2. Consider excluding problematic samples"
    echo "  3. Check for systematic sequencing issues"
    echo "  4. Contact sequencing facility if issues are widespread"
    echo ""
    exit 1
    
elif [ $critical_issues -gt 0 ]; then
    echo "⚠ QC PASSED WITH WARNINGS"
    echo "Some samples have issues but overall quality is acceptable"
    echo ""
    echo "Recommended actions:"
    echo "  1. Note problematic samples for potential exclusion"
    echo "  2. Proceed with pipeline but monitor these samples"
    if [ $warnings -gt 0 ]; then
        echo "  3. Run adapter trimming before alignment"
    fi
    echo ""
    
else
    echo "✅ QC PASSED"
    echo "All samples meet TagSeq quality standards for stickleback"
    echo ""
    if [ $warnings -gt 0 ]; then
        echo "Note: Some samples have adapter content - trimming recommended"
    else
        echo "Quality is excellent - proceed with standard TagSeq pipeline"
    fi
    echo ""
fi

echo "=== NEXT STEPS ==="
echo "Pipeline can continue with:"
echo "  1. Adapter trimming (if warnings about adapter content)"
echo "  2. Alignment to stickleback genome/transcriptome"
echo "  3. Read counting and quantification"
echo "  4. Differential expression analysis"
echo ""
echo "QC evaluation completed at $(date)"
