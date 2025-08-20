#!/bin/bash

#===============================================================================
# MASTER TAGSEQ PIPELINE RUNNER
# 
# Description: Runs the complete TagSeq pipeline from concatenation to count collection
# Author: Rogini Runghen
# Version: 2.0 (with Quality Control)
#
# Usage: 
#   ./run_tagseq_pipeline.sh [options]
#   
# Options:
#   --dry-run          Show what would be submitted without actually submitting
#   --start-at STEP    Start pipeline at specific step (1-10)
#   --stop-at STEP     Stop pipeline at specific step (1-10)
#   --interactive      Prompt before each step
#   --resume           Resume from failed jobs
#   --help             Show this help
#
# Steps:
#   1. Concatenate lane files (1-concatenate_files.sh)
#   2. Decompress files (2-run_gunzip.sh)  
#   3. Run FastQC and MultiQC (2.1-fastqc_multiqc_merged.sh)
#   4. Check TagSeq quality (2.2-check_tagseq_quality.sh)
#   5. Setup TagSeq repo (3.1-setup_tagseq_repo.sh)
#   6. Generate & run trimming (3.2-generate_clean_commands.sh + clean.sh)
#   7. Compress trimmed files (3.3-compress_files.sh)
#   8. Build STAR index (4-star_index.sh)
#   9. Generate samples & align (4.1-generate_sample_list.sh + 5-davis_star_extended.sh)
#   10. Collect counts (6-collect_counts.sh)
#===============================================================================

set -euo pipefail
IFS=$'\n\t'

# Script metadata
readonly SCRIPT_NAME="$(basename "$0")"
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly START_TIME=$(date +%s)
readonly TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Colors for output
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly PURPLE='\033[0;35m'
readonly CYAN='\033[0;36m'
readonly NC='\033[0m'

# Configuration
DRY_RUN=false
INTERACTIVE=false
START_AT=1
STOP_AT=10
RESUME_MODE=false
LOG_DIR="logs"
CHECKPOINT_FILE=".pipeline_checkpoints"

# Job tracking
declare -A JOB_IDS
declare -A JOB_STATUS
declare -A STEP_NAMES

#===============================================================================
# LOGGING AND UTILITY FUNCTIONS
#===============================================================================

log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo -e "$msg"
    [[ -f "$LOG_DIR/pipeline_${TIMESTAMP}.log" ]] && echo -e "$msg" >> "$LOG_DIR/pipeline_${TIMESTAMP}.log"
}

log_info() { log "${BLUE}INFO${NC}: $*"; }
log_warn() { log "${YELLOW}WARN${NC}: $*"; }
log_error() { log "${RED}ERROR${NC}: $*"; }
log_success() { log "${GREEN}SUCCESS${NC}: $*"; }
log_header() { log "${PURPLE}=== $* ===${NC}"; }
log_step() { log "${CYAN}STEP $1${NC}: $2"; }

# Progress tracking
show_pipeline_progress() {
    local current=$1
    local total=$2
    local step_name="$3"
    local percent=$((current * 100 / total))
    local filled=$((percent / 5))  # 20 chars max
    local empty=$((20 - filled))
    
    printf "\r${CYAN}Pipeline Progress: [%${filled}s%${empty}s] %d%% - %s${NC}" \
           "$(printf '█%.0s' $(seq 1 $filled))" \
           "$(printf '░%.0s' $(seq 1 $empty))" \
           "$percent" "$step_name"
}

#===============================================================================
# SLURM JOB MANAGEMENT
#===============================================================================

# Submit job with dependency
submit_job() {
    local step_num=$1
    local job_name="$2"
    local script_path="$3"
    local work_dir="$4"
    local dependency="${5:-}"
    
    local submit_cmd="sbatch"
    local full_job_name="tagseq_step${step_num}_${job_name}"
    
    # Add dependency if specified
    if [[ -n "$dependency" ]]; then
        submit_cmd="$submit_cmd --dependency=afterok:$dependency"
    fi
    
    # Add job name
    submit_cmd="$submit_cmd --job-name=$full_job_name"
    
    # Add script path
    submit_cmd="$submit_cmd $script_path"
    
    log_info "Submitting step $step_num: $job_name"
    log_info "Working directory: $work_dir"
    log_info "Script: $script_path"
    [[ -n "$dependency" ]] && log_info "Dependency: $dependency"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "DRY RUN: $submit_cmd"
        echo "fake_job_id_${step_num}"
        return 0
    fi
    
    # Change to working directory and submit
    local current_dir=$(pwd)
    cd "$work_dir"
    
    local job_output
    if job_output=$($submit_cmd 2>&1); then
        local job_id=$(echo "$job_output" | grep -o '[0-9]\+' | head -1)
        cd "$current_dir"
        
        if [[ -n "$job_id" ]]; then
            log_success "Submitted job ID: $job_id"
            echo "$job_id"
            return 0
        else
            log_error "Could not extract job ID from: $job_output"
            cd "$current_dir"
            return 1
        fi
    else
        log_error "Failed to submit job: $job_output"
        cd "$current_dir"
        return 1
    fi
}

# Wait for job completion
wait_for_job() {
    local job_id="$1"
    local step_name="$2"
    local timeout="${3:-3600}"  # Default 1 hour timeout
    
    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "DRY RUN: Would wait for job $job_id ($step_name)"
        return 0
    fi
    
    log_info "Waiting for job $job_id ($step_name) to complete..."
    
    local elapsed=0
    local check_interval=30
    
    while [[ $elapsed -lt $timeout ]]; do
        local job_state=$(squeue -j "$job_id" -h -o "%T" 2>/dev/null | head -1)
        
        if [[ -z "$job_state" ]]; then
            # Job not in queue, check if it completed successfully
            local job_info=$(sacct -j "$job_id" -n -o "State" 2>/dev/null | head -1 | tr -d ' ')
            
            case "$job_info" in
                "COMPLETED")
                    log_success "Job $job_id ($step_name) completed successfully"
                    return 0
                    ;;
                "FAILED"|"CANCELLED"|"TIMEOUT"|"NODE_FAIL")
                    log_error "Job $job_id ($step_name) failed with state: $job_info"
                    return 1
                    ;;
                *)
                    log_warn "Job $job_id ($step_name) finished with unknown state: $job_info"
                    return 1
                    ;;
            esac
        else
            case "$job_state" in
                "RUNNING"|"PENDING"|"CONFIGURING")
                    printf "\r${BLUE}Job $job_id ($step_name): $job_state for ${elapsed}s${NC}"
                    ;;
                "COMPLETED")
                    echo
                    log_success "Job $job_id ($step_name) completed"
                    return 0
                    ;;
                "FAILED"|"CANCELLED"|"TIMEOUT"|"NODE_FAIL")
                    echo
                    log_error "Job $job_id ($step_name) failed with state: $job_state"
                    return 1
                    ;;
            esac
        fi
        
        sleep $check_interval
        ((elapsed += check_interval))
    done
    
    echo
    log_error "Timeout waiting for job $job_id ($step_name)"
    return 1
}

#===============================================================================
# CHECKPOINT MANAGEMENT
#===============================================================================

save_checkpoint() {
    local step=$1
    local status="$2"
    local job_id="${3:-}"
    
    echo "$(date +%s):$step:$status:$job_id" >> "$CHECKPOINT_FILE"
}

load_checkpoints() {
    declare -g -A completed_steps
    
    if [[ -f "$CHECKPOINT_FILE" ]]; then
        while IFS=: read -r timestamp step status job_id; do
            if [[ "$status" == "completed" ]]; then
                completed_steps["$step"]="$job_id"
            fi
        done < "$CHECKPOINT_FILE"
        
        log_info "Loaded ${#completed_steps[@]} completed steps from checkpoint"
    fi
}

is_step_completed() {
    local step=$1
    [[ -n "${completed_steps[$step]:-}" ]]
}

#===============================================================================
# PIPELINE STEP DEFINITIONS
#===============================================================================

# Initialize step names
init_step_names() {
    STEP_NAMES[1]="Concatenate lane files"
    STEP_NAMES[2]="Decompress files"
    STEP_NAMES[3]="Run FastQC and MultiQC"
    STEP_NAMES[4]="Check TagSeq quality"
    STEP_NAMES[5]="Setup TagSeq repository"
    STEP_NAMES[6]="Generate and run trimming"
    STEP_NAMES[7]="Compress trimmed files"
    STEP_NAMES[8]="Build STAR index"
    STEP_NAMES[9]="Generate samples and align"
    STEP_NAMES[10]="Collect count files"
}

# Step 1: Concatenate lane files
run_step_1() {
    local step=1
    local step_name="${STEP_NAMES[$step]}"
    
    log_step $step "$step_name"
    
    # Check prerequisites
    if [[ ! -d "00_raw_fastqs" ]] || [[ $(ls 00_raw_fastqs/*.fastq* 2>/dev/null | wc -l || echo 0) -eq 0 ]]; then
        log_error "No FASTQ files found in 00_raw_fastqs/"
        return 1
    fi
    
    local script_path="1-concatenate_files.sh"
    if [[ ! -f "$script_path" ]]; then
        log_error "Script not found: $script_path"
        return 1
    fi
    
    local job_id
    if job_id=$(submit_job $step "concatenate" "$script_path" "00_raw_fastqs"); then
        JOB_IDS[$step]="$job_id"
        save_checkpoint $step "submitted" "$job_id"
        return 0
    else
        return 1
    fi
}

# Step 2: Decompress files
run_step_2() {
    local step=2
    local step_name="${STEP_NAMES[$step]}"
    local dependency="${JOB_IDS[1]:-}"
    
    log_step $step "$step_name"
    
    local script_path="2-run_gunzip.sh"
    if [[ ! -f "$script_path" ]]; then
        log_error "Script not found: $script_path"
        return 1
    fi
    
    local job_id
    if job_id=$(submit_job $step "gunzip" "$script_path" "01_merged_files" "$dependency"); then
        JOB_IDS[$step]="$job_id"
        save_checkpoint $step "submitted" "$job_id"
        return 0
    else
        return 1
    fi
}

# Step 3: Run FastQC and MultiQC
run_step_3() {
    local step=3
    local step_name="${STEP_NAMES[$step]}"
    local dependency="${JOB_IDS[2]:-}"
    
    log_step $step "$step_name"
    
    # Check prerequisites
    if [[ ! -d "01_merged_files" ]] || [[ $(ls 01_merged_files/*.fastq* 2>/dev/null | wc -l || echo 0) -eq 0 ]]; then
        log_error "No decompressed FASTQ files found in 01_merged_files/"
        return 1
    fi
    
    local script_path="2.1-fastqc_multiqc_merged.sh"
    if [[ ! -f "$script_path" ]]; then
        log_error "Script not found: $script_path"
        return 1
    fi
    
    local job_id
    if job_id=$(submit_job $step "fastqc_multiqc" "$script_path" "." "$dependency"); then
        JOB_IDS[$step]="$job_id"
        save_checkpoint $step "submitted" "$job_id"
        return 0
    else
        return 1
    fi
}

# Step 4: Check TagSeq quality
run_step_4() {
    local step=4
    local step_name="${STEP_NAMES[$step]}"
    local dependency="${JOB_IDS[3]:-}"
    
    log_step $step "$step_name"
    
    # Check prerequisites
    if [[ ! -d "02_fastqc/multiqc_reports" ]]; then
        log_error "MultiQC reports not found in 02_fastqc/multiqc_reports/"
        return 1
    fi
    
    local script_path="2.2-check_tagseq_quality.sh"
    if [[ ! -f "$script_path" ]]; then
        log_error "Script not found: $script_path"
        return 1
    fi
    
    # Create a wrapper script that runs the quality check
    local wrapper_script="02_fastqc/run_quality_check.sh"
    
    cat > "$wrapper_script" << 'EOF'
#!/bin/bash
#SBATCH --job-name=tagseq_qc_check
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.name@institution.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

echo "Starting TagSeq quality check at $(date)"

# Run quality check
if ../2.2-check_tagseq_quality.sh ./multiqc_reports; then
    echo "Quality check completed successfully at $(date)"
    echo "Review the output above for quality assessment"
else
    echo "Quality check indicated issues - review recommendations"
    echo "Pipeline can continue but consider the QC recommendations"
fi

echo "Quality check workflow completed at $(date)"
EOF
    
    chmod +x "$wrapper_script"
    
    local job_id
    if job_id=$(submit_job $step "qc_check" "$wrapper_script" "02_fastqc" "$dependency"); then
        JOB_IDS[$step]="$job_id"
        save_checkpoint $step "submitted" "$job_id"
        return 0
    else
        return 1
    fi
}

# Step 5: Setup TagSeq repository
run_step_5() {
    local step=5
    local step_name="${STEP_NAMES[$step]}"
    local dependency="${JOB_IDS[4]:-}"
    
    log_step $step "$step_name"
    
    local script_path="3.1-setup_tagseq_repo.sh"
    if [[ ! -f "$script_path" ]]; then
        log_error "Script not found: $script_path"
        return 1
    fi
    
    local job_id
    if job_id=$(submit_job $step "setup_tagseq" "$script_path" "." "$dependency"); then
        JOB_IDS[$step]="$job_id"
        save_checkpoint $step "submitted" "$job_id"
        return 0
    else
        return 1
    fi
}

# Step 6: Generate and run trimming
run_step_6() {
    local step=6
    local step_name="${STEP_NAMES[$step]}"
    local dependency="${JOB_IDS[5]:-}"
    
    log_step $step "$step_name"
    
    # This step involves multiple sub-steps
    local script_path="3.2-generate_clean_commands.sh"
    if [[ ! -f "$script_path" ]]; then
        log_error "Script not found: $script_path"
        return 1
    fi
    
    # Create a wrapper script that runs the generate command, then the trimming
    local wrapper_script="03_trim/run_trimming_wrapper.sh"
    mkdir -p 03_trim
    
    cat > "$wrapper_script" << 'EOF'
#!/bin/bash
#SBATCH --job-name=tagseq_trimming
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=username@institution.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

echo "Starting trimming workflow at $(date)"

# Generate clean commands
echo "Generating trimming commands..."
../3.2-generate_clean_commands.sh

# Check if clean_commands.txt was created
if [[ ! -f "clean_commands.txt" ]]; then
    echo "ERROR: clean_commands.txt not generated"
    exit 1
fi

# Create clean.sh from commands
echo "#!/bin/bash" > clean.sh
echo "set -euo pipefail" >> clean.sh
cat clean_commands.txt >> clean.sh
chmod +x clean.sh

# Run trimming
echo "Running trimming commands..."
if ./clean.sh; then
    echo "Trimming completed successfully at $(date)"
else
    echo "ERROR: Trimming failed"
    exit 1
fi

echo "Trimming workflow completed at $(date)"
EOF
    
    chmod +x "$wrapper_script"
    
    local job_id
    if job_id=$(submit_job $step "trimming" "$wrapper_script" "03_trim" "$dependency"); then
        JOB_IDS[$step]="$job_id"
        save_checkpoint $step "submitted" "$job_id"
        return 0
    else
        return 1
    fi
}

# Step 7: Compress trimmed files
run_step_7() {
    local step=7
    local step_name="${STEP_NAMES[$step]}"
    local dependency="${JOB_IDS[6]:-}"
    
    log_step $step "$step_name"
    
    local script_path="3.3-compress_files.sh"
    if [[ ! -f "$script_path" ]]; then
        log_error "Script not found: $script_path"
        return 1
    fi
    
    local job_id
    if job_id=$(submit_job $step "compress" "$script_path" "." "$dependency"); then
        JOB_IDS[$step]="$job_id"
        save_checkpoint $step "submitted" "$job_id"
        return 0
    else
        return 1
    fi
}

# Step 8: Build STAR index
run_step_8() {
    local step=8
    local step_name="${STEP_NAMES[$step]}"
    local dependency="${JOB_IDS[7]:-}"
    
    log_step $step "$step_name"
    
    local script_path="4-star_index.sh"
    if [[ ! -f "$script_path" ]]; then
        log_error "Script not found: $script_path"
        return 1
    fi
    
    local job_id
    if job_id=$(submit_job $step "star_index" "$script_path" "04_read_count" "$dependency"); then
        JOB_IDS[$step]="$job_id"
        save_checkpoint $step "submitted" "$job_id"
        return 0
    else
        return 1
    fi
}

# Step 9: Generate samples and align
run_step_9() {
    local step=9
    local step_name="${STEP_NAMES[$step]}"
    local dependency="${JOB_IDS[8]:-}"
    
    log_step $step "$step_name"
    
    # Create wrapper for sample generation + alignment
    local wrapper_script="04_read_count/run_alignment_wrapper.sh"
    
    cat > "$wrapper_script" << 'EOF'
#!/bin/bash
#SBATCH --job-name=tagseq_alignment
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 14
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=username@institution.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

echo "Starting alignment workflow at $(date)"

# Generate sample list
echo "Generating sample list..."
if ../4.1-generate_sample_list.sh; then
    echo "Sample list generated successfully"
else
    echo "ERROR: Failed to generate sample list"
    exit 1
fi

# Check if samples.txt was created
if [[ ! -f "samples.txt" ]]; then
    echo "ERROR: samples.txt not found"
    exit 1
fi

echo "Found $(wc -l < samples.txt) samples to process"

# Run STAR alignment
echo "Starting STAR alignment..."
if ../5-davis_star_extended.sh; then
    echo "STAR alignment completed successfully at $(date)"
else
    echo "ERROR: STAR alignment failed"
    exit 1
fi

echo "Alignment workflow completed at $(date)"
EOF
    
    chmod +x "$wrapper_script"
    
    local job_id
    if job_id=$(submit_job $step "alignment" "$wrapper_script" "04_read_count" "$dependency"); then
        JOB_IDS[$step]="$job_id"
        save_checkpoint $step "submitted" "$job_id"
        return 0
    else
        return 1
    fi
}

# Step 10: Collect counts
run_step_10() {
    local step=10
    local step_name="${STEP_NAMES[$step]}"
    local dependency="${JOB_IDS[9]:-}"
    
    log_step $step "$step_name"
    
    local script_path="6-collect_counts.sh"
    if [[ ! -f "$script_path" ]]; then
        log_error "Script not found: $script_path"
        return 1
    fi
    
    local job_id
    if job_id=$(submit_job $step "collect_counts" "$script_path" "04_read_count" "$dependency"); then
        JOB_IDS[$step]="$job_id"
        save_checkpoint $step "submitted" "$job_id"
        return 0
    else
        return 1
    fi
}

#===============================================================================
# MAIN PIPELINE EXECUTION
#===============================================================================

validate_environment() {
    log_info "Validating environment..."
    
    # Check required commands
    local required_commands=("sbatch" "squeue" "sacct")
    for cmd in "${required_commands[@]}"; do
        if ! command -v "$cmd" >/dev/null 2>&1; then
            log_error "Required command not found: $cmd"
            return 1
        fi
    done
    
    # Check directory structure
    if [[ ! -d "00_raw_fastqs" ]]; then
        log_error "Directory not found: 00_raw_fastqs"
        log_error "Are you in the correct project directory?"
        return 1
    fi
    
    # Check for pipeline scripts
    local scripts=("1-concatenate_files.sh" "2-run_gunzip.sh" "2.1-fastqc_multiqc_merged.sh"
                   "2.2-check_tagseq_quality.sh" "3.1-setup_tagseq_repo.sh" "3.2-generate_clean_commands.sh" 
                   "3.3-compress_files.sh" "4-star_index.sh" "4.1-generate_sample_list.sh" 
                   "5-davis_star_extended.sh" "6-collect_counts.sh")
    
    local missing_scripts=()
    for script in "${scripts[@]}"; do
        if [[ ! -f "$script" ]]; then
            missing_scripts+=("$script")
        fi
    done
    
    if [[ ${#missing_scripts[@]} -gt 0 ]]; then
        log_error "Missing pipeline scripts:"
        for script in "${missing_scripts[@]}"; do
            log_error "  - $script"
        done
        return 1
    fi
    
    log_success "Environment validation passed"
    return 0
}

run_pipeline() {
    log_header "Starting TagSeq Pipeline Execution"
    
    # Initialize
    init_step_names
    mkdir -p "$LOG_DIR"
    
    # Load checkpoints if resuming
    if [[ "$RESUME_MODE" == "true" ]]; then
        load_checkpoints
    fi
    
    # Define step functions
    local step_functions=(
        "" # Index 0 (unused)
        "run_step_1"
        "run_step_2" 
        "run_step_3"
        "run_step_4"
        "run_step_5"
        "run_step_6"
        "run_step_7"
        "run_step_8"
        "run_step_9"
        "run_step_10"
    )
    
    # Run pipeline steps
    local failed_steps=()
    local completed_steps=0
    local total_steps=$((STOP_AT - START_AT + 1))
    
    for ((step=START_AT; step<=STOP_AT; step++)); do
        local step_name="${STEP_NAMES[$step]}"
        
        # Check if step already completed
        if [[ "$RESUME_MODE" == "true" ]] && is_step_completed "$step"; then
            log_info "Step $step already completed: $step_name"
            ((completed_steps++))
            show_pipeline_progress $completed_steps $total_steps "$step_name (skipped)"
            continue
        fi
        
        # Interactive confirmation
        if [[ "$INTERACTIVE" == "true" ]]; then
            echo
            read -p "Run step $step ($step_name)? [Y/n]: " response
            if [[ "$response" =~ ^[Nn]$ ]]; then
                log_info "Skipping step $step: $step_name"
                continue
            fi
        fi
        
        # Submit step
        show_pipeline_progress $completed_steps $total_steps "Submitting: $step_name"
        
        if ${step_functions[$step]}; then
            log_success "Step $step submitted: $step_name"
            ((completed_steps++))
            show_pipeline_progress $completed_steps $total_steps "$step_name (submitted)"
        else
            log_error "Failed to submit step $step: $step_name"
            failed_steps+=($step)
            break
        fi
        
        echo  # New line for cleaner output
    done
    
    echo  # Final new line
    
    # Wait for completion if not dry run
    if [[ "$DRY_RUN" == "false" ]] && [[ ${#failed_steps[@]} -eq 0 ]]; then
        log_header "Waiting for Pipeline Completion"
        
        for ((step=START_AT; step<=STOP_AT; step++)); do
            local job_id="${JOB_IDS[$step]:-}"
            local step_name="${STEP_NAMES[$step]}"
            
            if [[ -n "$job_id" ]]; then
                if wait_for_job "$job_id" "$step_name" 7200; then  # 2 hour timeout per step
                    save_checkpoint $step "completed" "$job_id"
                    log_success "Step $step completed: $step_name"
                else
                    save_checkpoint $step "failed" "$job_id"
                    failed_steps+=($step)
                    log_error "Step $step failed: $step_name"
                    break
                fi
            fi
        done
    fi
    
    # Final report
    log_header "Pipeline Execution Summary"
    
    if [[ ${#failed_steps[@]} -eq 0 ]]; then
        log_success "All pipeline steps completed successfully! 🎉"
        log_info "Your TagSeq analysis is complete"
        log_info "Count files are ready in: 04_read_count/all_counts/"
        log_info "Quality reports available in: 02_fastqc/multiqc_reports/"
    else
        log_error "Pipeline failed at step(s): ${failed_steps[*]}"
        log_info "Check job logs for details"
        log_info "Use --resume to continue from failed step"
        return 1
    fi
}

#===============================================================================
# COMMAND LINE INTERFACE
#===============================================================================

show_help() {
    cat << EOF
TagSeq Master Pipeline Runner (with Quality Control)

DESCRIPTION:
    Runs the complete TagSeq pipeline from concatenation through count collection.
    Handles job dependencies, error recovery, and progress tracking.
    Includes integrated quality control steps with TagSeq-specific thresholds.

USAGE:
    $0 [OPTIONS]

OPTIONS:
    --dry-run              Show what would be submitted without actually submitting
    --start-at STEP        Start pipeline at specific step (1-10, default: 1)
    --stop-at STEP         Stop pipeline at specific step (1-10, default: 10)
    --interactive          Prompt before each step
    --resume               Resume from failed/incomplete jobs
    --help                 Show this help message

PIPELINE STEPS:
    1. Concatenate lane files         (1-concatenate_files.sh)
    2. Decompress files              (2-run_gunzip.sh)
    3. Run FastQC and MultiQC        (2.1-fastqc_multiqc_merged.sh)
    4. Check TagSeq quality          (2.2-check_tagseq_quality.sh)
    5. Setup TagSeq repository       (3.1-setup_tagseq_repo.sh)
    6. Generate and run trimming     (3.2-generate_clean_commands.sh + clean.sh)
    7. Compress trimmed files        (3.3-compress_files.sh)
    8. Build STAR index             (4-star_index.sh)
    9. Generate samples and align    (4.1-generate_sample_list.sh + 5-davis_star_extended.sh)
    10. Collect count files          (6-collect_counts.sh)

EXAMPLES:
    # Run complete pipeline
    $0

    # Dry run to see what would be submitted
    $0 --dry-run

    # Run interactively with confirmation prompts
    $0 --interactive

    # Run only steps 1-6 (through trimming)
    $0 --stop-at 6

    # Resume from step 7 after fixing an issue
    $0 --start-at 7 --resume

    # Run only the alignment steps
    $0 --start-at 8 --stop-at 9

    # Run only quality control steps
    $0 --start-at 3 --stop-at 4

QUALITY CONTROL FEATURES:
    - Automated FastQC and MultiQC analysis
    - TagSeq-specific quality thresholds for stickleback
    - Pass/fail decisions with clear recommendations
    - Integration with trimming workflow

PREREQUISITES:
    - Must be run from TagSeq project directory
    - All pipeline scripts must be present
    - Input data in 00_raw_fastqs/
    - SLURM cluster environment

MONITORING:
    - Check job status: squeue -u \$(whoami)
    - View logs in: logs/
    - Resume capability: $0 --resume
    - Quality reports: 02_fastqc/multiqc_reports/multiqc_report.html
EOF
}

parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --dry-run)
                DRY_RUN=true
                shift
                ;;
            --interactive)
                INTERACTIVE=true
                shift
                ;;
            --start-at)
                START_AT="$2"
                shift 2
                ;;
            --stop-at)
                STOP_AT="$2"
                shift 2
                ;;
            --resume)
                RESUME_MODE=true
                shift
                ;;
            --help|-h)
                show_help
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                show_help
                exit 1
                ;;
        esac
    done
    
    # Validate step ranges
    if [[ $START_AT -lt 1 ]] || [[ $START_AT -gt 10 ]]; then
        log_error "Invalid start step: $START_AT (must be 1-10)"
        exit 1
    fi
    
    if [[ $STOP_AT -lt 1 ]] || [[ $STOP_AT -gt 10 ]]; then
        log_error "Invalid stop step: $STOP_AT (must be 1-10)"
        exit 1
    fi
    
    if [[ $START_AT -gt $STOP_AT ]]; then
        log_error "Start step ($START_AT) cannot be greater than stop step ($STOP_AT)"
        exit 1
    fi
}

#===============================================================================
# MAIN EXECUTION
#===============================================================================

main() {
    # Parse command line arguments
    parse_arguments "$@"
    
    # Setup logging
    mkdir -p "$LOG_DIR"
    
    # Show configuration
    log_header "TagSeq Master Pipeline Runner (with Quality Control)"
    log_info "User: $(whoami)"
    log_info "Host: $(hostname)"
    log_info "Start time: $(date)"
    log_info "Project directory: $(pwd)"
    log_info "Steps to run: $START_AT to $STOP_AT"
    [[ "$DRY_RUN" == "true" ]] && log_info "Mode: DRY RUN"
    [[ "$INTERACTIVE" == "true" ]] && log_info "Mode: INTERACTIVE"
    [[ "$RESUME_MODE" == "true" ]] && log_info "Mode: RESUME"
    
    # Validate environment
    if ! validate_environment; then
        log_error "Environment validation failed"
        exit 1
    fi
    
    # Run pipeline
    if run_pipeline; then
        log_success "Pipeline execution completed successfully!"
        exit 0
    else
        log_error "Pipeline execution failed"
        exit 1
    fi
}

# Cleanup function
cleanup() {
    local exit_code=$?
    
    if [[ $exit_code -ne 0 ]]; then
        echo
        log_error "Pipeline runner interrupted or failed"
        log_info "Use --resume to continue from where it left off"
        log_info "Check logs in: $LOG_DIR/"
        log_info "Quality reports in: 02_fastqc/multiqc_reports/"
    fi
    
    exit $exit_code
}

trap cleanup EXIT INT TERM

# Execute main function
main "$@"
