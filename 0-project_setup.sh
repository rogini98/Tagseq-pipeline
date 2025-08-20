#!/bin/bash
#SBATCH --job-name=tagseq_setup_folders
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rogini.runghen@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

#===============================================================================
# TAGSEQ SETUP SCRIPT
#
# Description: Sets up project structure that works with your original scripts
# Author: Rogini Runghen
# Version: 1.0
#
# Usage:
#   ./setup_for_original_scripts.sh [project_path] [source_dir1] [source_dir2] ...
#   sbatch setup_for_original_scripts.sh /labs/Bolnick/ROL/my_new_project /data/source1 /data/source2
#
# Features:
# - Creates subfolders needed for subsequent steps of bioinformatics pipeline
# - Parallelized file copying for large datasets
# - Progress tracking and error handling
#===============================================================================

set -euo pipefail
IFS=$'\n\t'

# Script metadata
readonly SCRIPT_NAME="$(basename "$0")"
readonly START_TIME=$(date +%s)
readonly TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Colors for output
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly PURPLE='\033[0;35m'
readonly NC='\033[0m'

# Default configuration
PROJECT_PATH="${1:-/labs/Bolnick/ROL/test-samples}"
shift 2>/dev/null || true
SOURCE_DIRS=("$@")

# File copying parameters
MAX_PARALLEL="${MAX_PARALLEL:-4}"
VERIFY_COPIES="${VERIFY_COPIES:-true}"

#===============================================================================
# LOGGING FUNCTIONS
#===============================================================================

log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo -e "$msg"
}

log_info() { log "${BLUE}INFO${NC}: $*"; }
log_warn() { log "${YELLOW}WARN${NC}: $*"; }
log_error() { log "${RED}ERROR${NC}: $*"; }
log_success() { log "${GREEN}SUCCESS${NC}: $*"; }
log_header() { log "${PURPLE}=== $* ===${NC}"; }

# Progress tracking
show_progress() {
    local current=$1
    local total=$2
    local label="${3:-Processing}"
    local percent=$((current * 100 / total))
    local filled=$((percent / 2))
    local empty=$((50 - filled))

    printf "\r[%${filled}s%${empty}s] %d%% (%d/%d) %s" \
           "$(printf '#%.0s' $(seq 1 $filled))" \
           "$(printf ' %.0s' $(seq 1 $empty))" \
           "$percent" "$current" "$total" "$label"
}

# Cleanup function
cleanup() {
    local exit_code=$?
    log_info "Cleaning up..."

    # Kill background jobs
    jobs -p | xargs -r kill 2>/dev/null || true

    # Remove temp files
    rm -f /tmp/tagseq_setup_$$_*

    # Calculate runtime
    local end_time=$(date +%s)
    local runtime=$((end_time - START_TIME))
    local hours=$((runtime / 3600))
    local minutes=$(((runtime % 3600) / 60))
    local seconds=$((runtime % 60))

    log_info "Total runtime: ${hours}h ${minutes}m ${seconds}s"

    if [ $exit_code -eq 0 ]; then
        log_success "Setup completed successfully!"
    else
        log_error "Setup failed with exit code $exit_code"
    fi

    exit $exit_code
}

trap cleanup EXIT INT TERM

#===============================================================================
# DIRECTORY SETUP - EXACTLY AS YOUR SCRIPTS EXPECT
#===============================================================================

create_project_structure() {
    log_header "Creating Project Structure"
    log_info "Project location: $PROJECT_PATH"

    # Create main project directory
    mkdir -p "$PROJECT_PATH"
    cd "$PROJECT_PATH"

    # Create exactly the directories your original scripts expect
    local directories=(
        "00_raw_fastqs"
        "01_merged_files"
        "02_fastqc"
        "03_trim"
        "04_read_count"
        "04_read_count/References"
    )

    log_info "Creating directory structure..."
    for dir in "${directories[@]}"; do
        mkdir -p "$dir"
        log_info "Created: $dir"
    done

    log_success "Directory structure created"
    log_info "Current working directory: $(pwd)"
}

#===============================================================================
# FILE DISCOVERY AND COPYING
#===============================================================================

find_fastq_files() {
    local source_dir="$1"
    local temp_file="/tmp/tagseq_setup_$$_files_$(basename "$source_dir" | tr '/' '_')"

    log_info "Scanning for FASTQ files in: $source_dir"

    if [[ ! -d "$source_dir" ]]; then
        log_error "Source directory not found: $source_dir"
        return 1
    fi

    # Find all FASTQ files recursively (including Illumina naming convention)
    find "$source_dir" -type f \( \
        -name "*.fastq" -o \
        -name "*.fq" -o \
        -name "*.fastq.gz" -o \
        -name "*.fq.gz" -o \
        -name "*_R[12]_001.fastq.gz" -o \
        -name "*_R[12]_001.fq.gz" -o \
        -name "*_R[12]_[0-9][0-9][0-9].fastq.gz" -o \
        -name "*_R[12]_[0-9][0-9][0-9].fq.gz" \
    \) > "$temp_file" 2>/dev/null

    local file_count=$(wc -l < "$temp_file" 2>/dev/null || echo 0)

    if [[ $file_count -eq 0 ]]; then
        log_warn "No FASTQ files found in: $source_dir"
        rm -f "$temp_file"
        return 1
    fi

    log_success "Found $file_count FASTQ files in: $source_dir"
    echo "$temp_file"  # Return path to file list
}

# Parallel file copying function
copy_files_parallel() {
    local file_list="$1"
    local dest_dir="00_raw_fastqs"

    local total_files=$(wc -l < "$file_list")
    log_info "Copying $total_files files to $dest_dir..."

    if [[ $total_files -lt 50 ]]; then
        # For small datasets, use your original method (simple and fast)
        log_info "Using simple copy method for small dataset"
        cat "$file_list" | xargs -I {} cp {} "$dest_dir/"
    else
        # For large datasets, use parallel processing
        log_info "Using parallel copy method for large dataset ($MAX_PARALLEL processes)"

        # Split file list into chunks
        local chunk_size=$((total_files / MAX_PARALLEL + 1))
        split -l "$chunk_size" "$file_list" "/tmp/tagseq_setup_$$_chunk_"

        # Copy chunks in parallel
        local pids=()
        local chunk_files=($(ls /tmp/tagseq_setup_$$_chunk_* 2>/dev/null))

        for chunk_file in "${chunk_files[@]}"; do
            {
                while IFS= read -r source_file; do
                    [[ -n "$source_file" ]] && cp "$source_file" "$dest_dir/"
                done < "$chunk_file"
                rm -f "$chunk_file"
            } &
            pids+=($!)
        done

        # Wait for all copies to complete with progress tracking
        local completed=0
        local total_chunks=${#pids[@]}

        for pid in "${pids[@]}"; do
            wait "$pid"
            ((completed++))
            show_progress $completed $total_chunks "copying chunks"
        done
        echo  # New line after progress
    fi

    # Verify copy results
    local copied_count=$(ls "$dest_dir" | wc -l)
    log_success "Successfully copied $copied_count files"

    if [[ $copied_count -ne $total_files ]]; then
        log_warn "File count mismatch: expected $total_files, got $copied_count"
    fi
}

#===============================================================================
# MAIN PROCESSING
#===============================================================================

process_all_sources() {
    log_header "Processing Source Directories"

    if [[ ${#SOURCE_DIRS[@]} -eq 0 ]]; then
        log_error "No source directories specified"
        echo
        echo "Usage: $0 [project_path] [source_dir1] [source_dir2] ..."
        echo "Example: $0 /labs/Bolnick/ROL/my_project /data/JA23263_RoL_2019-2022_full_data"
        exit 1
    fi

    # Collect all FASTQ files from all sources
    local all_files_list="/tmp/tagseq_setup_$$_all_files"
    > "$all_files_list"  # Clear file

    local total_sources=0
    local successful_sources=0

    for source_dir in "${SOURCE_DIRS[@]}"; do
        ((total_sources++))
        log_info "Processing source $total_sources: $source_dir"

        if file_list=$(find_fastq_files "$source_dir"); then
            # Append to master list
            cat "$file_list" >> "$all_files_list"
            rm -f "$file_list"
            ((successful_sources++))
            log_success "Successfully processed: $source_dir"
        else
            log_warn "Skipped source directory: $source_dir"
        fi
    done

    if [[ $successful_sources -eq 0 ]]; then
        log_error "No FASTQ files found in any source directory"
        exit 1
    fi

    log_info "Found files from $successful_sources out of $total_sources source directories"

    # Remove duplicates and copy files
    sort -u "$all_files_list" > "${all_files_list}.unique"
    local unique_files=$(wc -l < "${all_files_list}.unique")
    log_info "Total unique FASTQ files to copy: $unique_files"

    copy_files_parallel "${all_files_list}.unique"

    # Cleanup temp files
    rm -f "$all_files_list" "${all_files_list}.unique"
}

#===============================================================================
# SETUP ORIGINAL SCRIPT ENVIRONMENT
#===============================================================================

setup_script_environment() {
    log_header "Setting Up Environment for Original Scripts"

    # Create any additional files/directories your scripts might need
    log_info "Preparing environment for original scripts..."

    # Your scripts work from the project root, so we're already in the right place
    log_info "Working directory: $(pwd)"
    log_info "This matches the structure your original scripts expect"

    # Create a simple samples list (this might be useful)
    if [[ -d "00_raw_fastqs" ]] && [[ $(ls 00_raw_fastqs/*.fastq* 2>/dev/null | wc -l || echo 0) -gt 0 ]]; then
        log_info "Generating initial samples list..."
        ls 00_raw_fastqs/*.fastq* 2>/dev/null | \
            sed 's|00_raw_fastqs/||; s/\.fastq.*$//' | \
            sort -u > samples_initial.txt || true

        local sample_count=$(wc -l < samples_initial.txt 2>/dev/null || echo 0)
        log_info "Identified $sample_count unique samples"
    fi

    log_success "Environment ready for original scripts"
}

#===============================================================================
# GENERATE USAGE INSTRUCTIONS
#===============================================================================

generate_instructions() {
    log_header "Generating Usage Instructions"

    local instruction_file="NEXT_STEPS.txt"

    cat > "$instruction_file" << EOF
TagSeq Pipeline - Ready to Use Your Original Scripts
===================================================

Generated: $(date)
Project: $(basename "$PROJECT_PATH")
Location: $PROJECT_PATH
User: $(whoami)

SETUP SUMMARY:
- Total FASTQ files organized: $(ls 00_raw_fastqs/*.fastq* 2>/dev/null | wc -l || echo 0)
- Directory structure: Ready for your original scripts
- Working directory: $(pwd)

YOUR ORIGINAL WORKFLOW (unchanged):
===================================

Step 1: ✅ COMPLETED - Files organized in 00_raw_fastqs/

Step 2: ✅ COMPLETED - Directory structure created

Step 3: Concatenate lane files
> cd 00_raw_fastqs
> sbatch ../concatenate_files.sh

Step 4: Quality control
> # Run your fastqc and multiqc scripts in 02_fastqc/

Step 5: Trim reads
> cd 01_merged_files
> # Run your gunzip, setup_tagseq_repo, and trimming scripts

Step 6: Alignment and counting
> cd 04_read_count
> # Run your STAR index, alignment, and counting scripts

IMPORTANT NOTES:
- Your original scripts will work without modification
- Run them from the same relative directories as before
- All paths and directory structures match your expectations
- SLURM job submission works exactly as before

MONITORING:
- Check job status: squeue -u $(whoami)
- Monitor progress in your script output files

NEXT STEP:
Run your concatenate_files.sh script to merge lane files:
> cd 00_raw_fastqs
> sbatch ../concatenate_files.sh

EOF

    log_success "Instructions saved to: $instruction_file"

    # Also display key information
    echo
    log_header "SETUP COMPLETE"
    log_success "Project location: $PROJECT_PATH"
    log_success "FASTQ files organized: $(ls 00_raw_fastqs/*.fastq* 2>/dev/null | wc -l || echo 0) files"
    log_info "Next step: cd 00_raw_fastqs && sbatch ../concatenate_files.sh"
    echo
    log_info "See $instruction_file for detailed next steps"
}

#===============================================================================
# MAIN EXECUTION
#===============================================================================

show_help() {
    cat << EOF
TagSeq Setup Script

DESCRIPTION:
    Sets up a TagSeq project with the exact structure of the various subfolders required for the subsequent steps

USAGE:
    $0 [PROJECT_PATH] [SOURCE_DIR1] [SOURCE_DIR2] ...

ARGUMENTS:
    PROJECT_PATH    Full path where to create the project (default: /labs/Bolnick/ROL/test-samples)
    SOURCE_DIR      One or more directories containing raw FASTQ files

EXAMPLES:
    # Create project at default location
    $0 /labs/Bolnick/ROL/test-samples /data/JA23263_RoL_2019-2022_full_data

    # Create project at custom location with multiple sources
    $0 /labs/Bolnick/ROL/my_new_project \\
        /data/JA23263_RoL_2019-2022_full_data \\
        /data/FASTQ_Generation_2024-01-26_15_04_37Z-715950346

    # Submit as SLURM job
    sbatch $0 /labs/Bolnick/ROL/RoL_2024 /data/source1 /data/source2

WHAT IT DOES:
    ✓ Creates exact directory structure your scripts expect
    ✓ Copies/organizes all FASTQ files into 00_raw_fastqs/
    ✓ Sets up working environment
    ✓ Generates usage instructions

AFTER SETUP:
    cd [PROJECT_PATH]
    # Continue with workflow
    cd 00_raw_fastqs
    sbatch ../concatenate_files.sh
EOF
}

# Add this near the top of main() function, before calling other functions
main() {
    # Handle help request
    if [[ $# -eq 0 ]] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
        show_help
        exit 0
    fi

    # Store original working directory BEFORE any directory changes
    readonly ORIGINAL_WD=$(pwd)

    # Convert source directories to absolute paths
    local abs_source_dirs=()
    for src_dir in "${SOURCE_DIRS[@]}"; do
        if [[ "$src_dir" == "." ]]; then
            abs_source_dirs+=("$ORIGINAL_WD")
        elif [[ "$src_dir" =~ ^/ ]]; then
            # Already absolute
            abs_source_dirs+=("$src_dir")
        else
            # Make relative path absolute
            abs_source_dirs+=("$ORIGINAL_WD/$src_dir")
        fi
    done

    # Update SOURCE_DIRS with absolute paths
    SOURCE_DIRS=("${abs_source_dirs[@]}")

    # Start setup
    log_header "TagSeq Project Setup - Original Script Compatible"
    log_info "Host: $(hostname)"
    log_info "User: $(whoami)"
    log_info "Start time: $(date)"
    log_info "Original working directory: $ORIGINAL_WD"
    log_info "Target location: $PROJECT_PATH"
    log_info "Source directories: ${SOURCE_DIRS[*]}"

    # Validate inputs
    if [[ ${#SOURCE_DIRS[@]} -eq 0 ]]; then
        log_error "No source directories specified"
        show_help
        exit 1
    fi

    for src_dir in "${SOURCE_DIRS[@]}"; do
        if [[ ! -d "$src_dir" ]]; then
            log_error "Source directory does not exist: $src_dir"
            exit 1
        fi
    done

    # Execute setup steps
    create_project_structure
    process_all_sources
    setup_script_environment
    generate_instructions

    log_success "Setup completed successfully!"
    log_info "Location: $PROJECT_PATH"
    log_info "Ready to start TagSeq Bioinformatics pipeline!"
}

# Execute main function
main "$@"
