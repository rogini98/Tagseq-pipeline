#!/bin/bash
#SBATCH --job-name=fastqc_multiqc_raw
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=100G
#SBATCH --time 12:00:00
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rogini.runghen@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo "Starting QC pipeline on $(hostname) at $(date)"

#################################################################
# STEP 1: Run FastQC on raw merged files
#################################################################

echo "STEP 1: Running FastQC on raw merged files..."

## Changing directories to where the fastq files are located
cd "/labs/Bolnick/ROL/test-samples/01_merged_files"

## Loading modules required for script commands
module load fastqc

# Output directory for FastQC results
FASTQC_OUTDIR="/labs/Bolnick/ROL/test-samples/02_fastqc"
mkdir -p $FASTQC_OUTDIR

## Running FASTQC
echo "  Running FastQC with 6 threads..."
fastqc -t 6 -o $FASTQC_OUTDIR *.fastq.gz

echo "  FastQC completed at $(date)"
echo "  FastQC results saved to: $FASTQC_OUTDIR"

#################################################################
# STEP 2: Run MultiQC on FastQC results
#################################################################

echo ""
echo "STEP 2: Running MultiQC on FastQC results..."

# Load MultiQC module
module load MultiQC/1.9

# Change to FastQC output directory
cd $FASTQC_OUTDIR

# Output directory for MultiQC reports
MULTIQC_OUTDIR="./multiqc_reports"
mkdir -p $MULTIQC_OUTDIR

echo "  Running MultiQC..."
multiqc --outdir $MULTIQC_OUTDIR ./

echo "  MultiQC completed at $(date)"
echo "  MultiQC report saved to: $MULTIQC_OUTDIR"

#################################################################
# STEP 3: Summary and next steps
#################################################################

echo ""
echo "=== QC PIPELINE COMPLETED ==="
echo "Completion time: $(date)"
echo ""
echo "Results locations:"
echo "  FastQC individual reports: $FASTQC_OUTDIR"
echo "  MultiQC summary report: $MULTIQC_OUTDIR/multiqc_report.html"
echo ""
echo "Next steps:"
echo "  1. Open the MultiQC report in a browser to examine quality"
echo "  2. Based on the results, determine trimming parameters"
echo "  3. Proceed with adapter trimming if needed"
echo ""
echo "To view the MultiQC report:"
echo "  firefox $MULTIQC_OUTDIR/multiqc_report.html"
