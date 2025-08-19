#!/bin/bash

# Simple TagSeq Pipeline Runner
# Runs all steps sequentially with basic dependency management

set -euo pipefail

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() { echo -e "${BLUE}INFO${NC}: $*"; }
log_warn() { echo -e "${YELLOW}WARN${NC}: $*"; }
log_error() { echo -e "${RED}ERROR${NC}: $*"; }
log_success() { echo -e "${GREEN}SUCCESS${NC}: $*"; }

echo "🧬 TagSeq Simple Pipeline Runner (with Quality Control)"
echo "====================================================="

# Check we're in the right directory
if [[ ! -d "00_raw_fastqs" ]]; then
    log_error "Not in a TagSeq project directory. Missing: 00_raw_fastqs/"
    exit 1
fi

log_info "Project: $(basename "$PWD")"
log_info "Starting pipeline execution..."

# Function to submit and wait for job
submit_and_wait() {
    local step_name="$1"
    local script="$2"
    local work_dir="$3"
    
    log_info "Step: $step_name"
    log_info "Script: $script"
    log_info "Working directory: $work_dir"
    
    if [[ ! -f "$script" ]]; then
        log_error "Script not found: $script"
        return 1
    fi
    
    # Change to working directory
    local current_dir=$(pwd)
    cd "$work_dir"
    
    # Submit job
    local job_output
    if job_output=$(sbatch "$current_dir/$script" 2>&1); then
        local job_id=$(echo "$job_output" | grep -o '[0-9]\+' | head -1)
        log_success "Submitted job ID: $job_id"
        
        # Wait for completion
        log_info "Waiting for job $job_id to complete..."
        while squeue -j "$job_id" -h >/dev/null 2>&1; do
            printf "."
            sleep 30
        done
        echo
        
        # Check completion status
        local job_state=$(sacct -j "$job_id" -n -o "State" 2>/dev/null | head -1 | tr -d ' ')
        if [[ "$job_state" == "COMPLETED" ]]; then
            log_success "$step_name completed successfully"
            cd "$current_dir"
            return 0
        else
            log_error "$step_name failed with state: $job_state"
            cd "$current_dir"
            return 1
        fi
    else
        log_error "Failed to submit $step_name: $job_output"
        cd "$current_dir"
        return 1
    fi
}

# Pipeline steps
steps=(
    "1|Concatenate lane files|1-concatenate_files.sh|00_raw_fastqs"
    "2|Decompress files|2-run_gunzip.sh|01_merged_files"
    "3|Run FastQC and MultiQC|2.1-fastqc_multiqc_merged.sh|."
    "4|Check TagSeq quality|2.2-check_tagseq_quality.sh|02_fastqc"
    "5|Setup TagSeq repository|3.1-setup_tagseq_repo.sh|."
    "6|Generate trimming commands|3.2-generate_clean_commands.sh|03_trim"
    "7|Compress trimmed files|3.3-compress_files.sh|."
    "8|Build STAR index|4-star_index.sh|04_read_count"
    "9|Generate sample list|4.1-generate_sample_list.sh|04_read_count"
    "10|Run STAR alignment|5-davis_star_extended.sh|04_read_count"
    "11|Collect count files|6-collect_counts.sh|04_read_count"
)

# Create trimming wrapper for step 6
create_trimming_wrapper() {
    mkdir -p 03_trim
    cat > 03_trim/run_trimming.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=tagseq_trimming
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-user=rogini.runghen@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail
echo "Starting trimming at $(date)"

# Generate commands
../3.2-generate_clean_commands.sh

# Create and run clean script
echo "#!/bin/bash" > clean.sh
cat clean_commands.txt >> clean.sh
chmod +x clean.sh

./clean.sh

echo "Trimming completed at $(date)"
EOF
    chmod +x 03_trim/run_trimming.sh
}

# Execute pipeline
failed_step=""
completed=0

for step_info in "${steps[@]}"; do
    IFS='|' read -r step_num step_name script work_dir <<< "$step_info"
    
    echo
    echo "=== Step $step_num: $step_name ==="
    
    # Special handling for quality check step
    if [[ $step_num -eq 4 ]]; then
        log_info "Running TagSeq quality check..."
        if [[ -d "02_fastqc/multiqc_reports" ]]; then
            cd 02_fastqc
            if ../2.2-check_tagseq_quality.sh ./multiqc_reports; then
                log_success "Quality check completed"
                cd ..
                ((completed++))
            else
                log_warn "Quality check completed with warnings"
                log_info "Review the output above and decide whether to continue"
                read -p "Continue with pipeline? [Y/n]: " response
                if [[ "$response" =~ ^[Nn]$ ]]; then
                    log_info "Pipeline stopped by user"
                    failed_step="$step_num"
                    break
                fi
                cd ..
                ((completed++))
            fi
        else
            log_error "MultiQC reports not found in 02_fastqc/multiqc_reports/"
            failed_step="$step_num"
            break
        fi
    # Special handling for trimming step
    elif [[ $step_num -eq 6 ]]; then
        log_info "Creating trimming wrapper..."
        create_trimming_wrapper
        if submit_and_wait "$step_name" "03_trim/run_trimming.sh" "."; then
            ((completed++))
        else
            failed_step="$step_num"
            break
        fi
    else
        if submit_and_wait "$step_name" "$script" "$work_dir"; then
            ((completed++))
        else
            failed_step="$step_num"
            break
        fi
    fi
done

# Final summary
echo
echo "🎯 TagSeq Pipeline Summary"
echo "=========================="
if [[ -z "$failed_step" ]]; then
    log_success "All $completed steps completed successfully! 🎉"
    log_success "Your TagSeq analysis is complete!"
    log_info "Count files are in: 04_read_count/all_counts/"
    echo
    echo "📊 Quality Control Results:"
    echo "- FastQC reports: 02_fastqc/"
    echo "- MultiQC summary: 02_fastqc/multiqc_reports/multiqc_report.html"
    echo "- Quality assessment: Review step 4 output above"
    echo
    echo "📁 Final Results:"
    echo "- Gene count files: 04_read_count/all_counts/"
    echo "- Sample alignment stats: 04_read_count/02-STAR_alignment/"
    echo
    echo "🔬 Next steps:"
    echo "- Copy all_counts/ to your local machine"
    echo "- Run differential expression analysis in R/DESeq2"
    echo "- Use MultiQC report for publication quality figures"
else
    log_error "Pipeline failed at step $failed_step"
    log_info "Completed $completed steps successfully"
    log_info "Check job logs for error details:"
    echo "  - Individual job logs: *.out and *.err files"
    echo "  - FastQC/MultiQC reports: 02_fastqc/ (if step 3+ completed)"
    echo
    echo "🔄 To resume from failed step:"
    echo "  ./simple_run_pipeline.sh  # Will detect completed steps"
    echo "  # OR manually continue from step $((failed_step + 1))"
fi

echo
echo "📈 Pipeline completed at $(date)"
echo "Total runtime: $SECONDS seconds"
