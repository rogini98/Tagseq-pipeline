# TagSeq Analysis Pipeline

A comprehensive, production-ready pipeline for Tag-based RNA-seq (TagSeq) analysis from raw FASTQ files to gene expression count matrices. Optimized for three-spined stickleback (*Gasterosteus aculeatus*) but adaptable to other organisms.

## 🎯 Overview

This pipeline automates the complete TagSeq workflow with integrated quality control, error handling, and resume capabilities. Originally developed for stickleback transcriptomics research, it provides a robust framework for high-throughput TagSeq analysis on SLURM clusters.

### ✨ Key Features

- **🔄 Complete Automation**: Single-command execution from raw files to count matrices
- **🧪 Integrated Quality Control**: TagSeq-specific QC with automated pass/fail assessment
- **⚡ Parallel Processing**: Optimized for large datasets with parallel file handling
- **🔗 Smart Dependencies**: Automatic job dependency management in SLURM
- **📊 Resume Capability**: Continue from failed steps without restarting
- **📈 Progress Tracking**: Real-time monitoring with detailed logging
- **🎓 Student-Friendly**: Clear documentation and error messages
- **🔧 Configurable**: Easy adaptation to different organisms and cluster environments

## 📋 Pipeline Steps

| Step | Process | Script | Output |
|------|---------|--------|--------|
| **1** | File Organization | `setup_for_original_scripts.sh` | Organized project structure |
| **2** | Lane Concatenation | `1-concatenate_files.sh` | `01_merged_files/` |
| **3** | File Decompression | `2-run_gunzip.sh` | Uncompressed FASTQ files |
| **4** | Quality Control | `2.1-fastqc_multiqc_merged.sh` | FastQC/MultiQC reports |
| **5** | Quality Assessment | `2.2-check_tagseq_quality.sh` | QC pass/fail decision |
| **6** | TagSeq Setup | `3.1-setup_tagseq_repo.sh` | TagSeq clipper tools |
| **7** | Read Trimming | `3.2-generate_clean_commands.sh` | Adapter-trimmed reads |
| **8** | File Compression | `3.3-compress_files.sh` | Compressed trimmed files |
| **9** | Reference Index | `4-star_index.sh` | STAR genome index |
| **10** | Sample Preparation | `4.1-generate_sample_list.sh` | Sample manifest |
| **11** | Read Alignment | `5-davis_star_extended.sh` | BAM files + counts |
| **12** | Count Collection | `6-collect_counts.sh` | Gene expression matrix |

## 🚀 Quick Start

### Prerequisites

- **SLURM cluster environment**
- **Required software**: STAR, samtools, FastQC, MultiQC, cutadapt
- **Internet access** for reference genome download
- **Bash 4.0+** with standard UNIX utilities

### Installation

```bash
# Clone the repository
git clone https://github.com/rogini98/tagseq-pipeline.git
cd tagseq-pipeline

# Make scripts executable
chmod +x *.sh

# Check your environment
./utilities/check_dependencies.sh
```

### Basic Usage

```bash
# 1. Set up project and organize data
sbatch 0-project_setup.sh 

# 2. Navigate to project directory
cd /my_project_directory

# 3. Copy pipeline scripts
cp /path/to/tagseq-pipeline/*.sh .

# 4. Run complete pipeline
./simple_run_pipeline.sh
```

## 📖 Detailed Usage

### Project Setup

The pipeline begins with automated project setup and data organization:

```bash
# Basic setup with multiple data sources
./0-project_setup.sh PROJECT_PATH SOURCE_DIR1 SOURCE_DIR2

# Example for stickleback RoL project
sbatch 0-project_setup.sh /labs/Bolnick/ROL/RoL_2024_analysis \
    /data/JA23263_RoL_2019-2022_full_data \
    /data/FASTQ_Generation_2024-01-26_15_04_37Z-715950346
```

This creates a complete directory structure and organizes all FASTQ files from nested source directories.

### Pipeline Execution Options

#### Option 1: Simple Pipeline Runner (Recommended)
```bash
# Automatic execution with basic monitoring
./simple_run_pipeline.sh
```

#### Option 2: Master Pipeline Runner (Advanced)
```bash
# Complete pipeline with all features
./run_tagseq_pipeline.sh

# Interactive mode with step-by-step confirmation
./run_tagseq_pipeline.sh --interactive

# Preview jobs without submitting
./run_tagseq_pipeline.sh --dry-run

# Run specific steps only
./run_tagseq_pipeline.sh --start-at 4 --stop-at 6

# Resume from failed step
./run_tagseq_pipeline.sh --resume
```

### Quality Control Integration

The pipeline includes automated quality assessment:

```bash
# Run only QC steps
./run_tagseq_pipeline.sh --start-at 3 --stop-at 4

# View QC results
firefox 02_fastqc/multiqc_reports/multiqc_report.html
```

Quality thresholds optimized for TagSeq:
- **Minimum quality**: Q20 (suitable for TagSeq)
- **Adapter tolerance**: 15% (higher than standard RNA-seq)
- **Duplication tolerance**: 70% (TagSeq naturally high)
- **Minimum reads**: 500K per sample

## 📁 Repository Structure

```
tagseq-pipeline/
├── README.md                              # This file
├── LICENSE                                # MIT license
├── 0-project_setup.sh                     # Project setup and data organization
├── simple_run_pipeline.sh                 # Simple pipeline runner
├── run_tagseq_pipeline.sh                 # Advanced pipeline runner
├── scripts/                               # Individual pipeline steps
│   ├── 1-concatenate_files.sh            # Lane concatenation
│   ├── 2-run_gunzip.sh                   # File decompression
│   ├── 2.1-fastqc_multiqc_merged.sh      # Quality control analysis
│   ├── 2.2-check_tagseq_quality.sh       # Quality assessment
│   ├── 3.1-setup_tagseq_repo.sh          # TagSeq tools setup
│   ├── 3.2-generate_clean_commands.sh    # Trimming preparation
│   ├── 3.3-compress_files.sh             # File compression
│   ├── 4-star_index.sh                   # Reference indexing
│   ├── 4.1-generate_sample_list.sh       # Sample manifest
│   ├── 5-davis_star_extended.sh          # Read alignment
│   └── 6-collect_counts.sh               # Count collection
├── config/                                # Configuration templates
│   ├── pipeline.conf                     # Main configuration
│   └── slurm_template.sh                 # SLURM job template
├── utilities/                             # Helper scripts
│   ├── check_dependencies.sh             # Environment validation
│   ├── setup_project.sh                  # Quick project initialization
│   └── generate_report.sh                # Summary reporting
├── docs/                                  # Documentation
│   ├── installation.md                   # Detailed setup guide
│   ├── configuration.md                  # Customization options
│   ├── troubleshooting.md                # Common issues and solutions
│   └── examples/                         # Usage examples
└── examples/                              # Example datasets and configs
    ├── sample_data/                      # Small test dataset
    └── expected_outputs/                 # Reference results
```

## ⚙️ Configuration

### Basic Configuration

Edit `config/pipeline.conf` for your environment:

```bash
# Cluster settings
EMAIL="your.email@institution.edu"
PARTITION="general"
QOS="general"

# Resource allocation
STAR_THREADS=14
LARGE_MEMORY="120G"
DEFAULT_TIME="24:00:00"

# Module names (adjust for your cluster)
STAR_MODULE="STAR"
SAMTOOLS_MODULE="samtools"
FASTQC_MODULE="fastqc"
```

### Advanced Configuration

- **Reference genome**: Automatically downloads stickleback genome (customizable)
- **Quality thresholds**: TagSeq-optimized parameters (adjustable)
- **Resource limits**: Adaptive based on data size
- **Job dependencies**: Intelligent SLURM dependency chains

## 📊 Monitoring and Output

### Progress Monitoring

```bash
# Check job status
squeue -u $(whoami)

# Monitor pipeline logs
tail -f logs/pipeline_*.log

# View individual step outputs
ls *.out *.err
```

### Expected Outputs

**Quality Control:**
- `02_fastqc/multiqc_reports/multiqc_report.html` - Comprehensive QC report
- Individual FastQC reports for each sample

**Final Results:**
- `04_read_count/all_counts/` - Gene expression count files
- `04_read_count/02-STAR_alignment/` - Alignment BAM files and statistics

**Logs and Reports:**
- `logs/` - Detailed pipeline execution logs
- Job-specific `.out` and `.err` files

## 🕐 Runtime Expectations

| Dataset Size | Expected Runtime | Bottleneck Steps |
|--------------|------------------|------------------|
| **Small** (< 50 samples, < 100GB) | 4-8 hours | STAR alignment |
| **Medium** (50-200 samples, 100-500GB) | 8-24 hours | STAR alignment, trimming |
| **Large** (200+ samples, > 500GB) | 24-48 hours | File organization, alignment |

**Critical steps:**
- **Step 8** (STAR alignment): Usually 50-70% of total runtime
- **Step 4** (Quality control): 10-20% of total runtime
- **Step 7** (Trimming): 15-25% of total runtime

## 🐛 Troubleshooting

### Common Issues

**Pipeline fails at quality control:**
```bash
# Check FastQC results
firefox 02_fastqc/multiqc_reports/multiqc_report.html

# Review quality thresholds
./2.2-check_tagseq_quality.sh 02_fastqc/multiqc_reports/
```

**STAR alignment failures:**
```bash
# Check reference index
ls 04_read_count/References/stickleback_v5/

# Review memory usage
grep -i "memory\|killed" *.err

# Check sample file format
head 04_read_count/samples.txt
```

**Job dependency issues:**
```bash
# Check job status
squeue -u $(whoami)

# Review failed jobs
sacct -u $(whoami) --starttime=today --state=FAILED

# Resume from specific step
./run_tagseq_pipeline.sh --start-at 6 --resume
```

### Getting Help

1. **Check logs**: Always start with `.out` and `.err` files
2. **Review documentation**: See `docs/` for detailed guides
3. **Validate environment**: Run `utilities/check_dependencies.sh`
4. **Test with small dataset**: Use example data first
5. **Contact support**: Open an issue with log files

## 🔬 Adapting to Other Organisms

### Quick Adaptation

1. **Update reference URLs** in `4-star_index.sh`:
```bash
# Replace with your organism's genome and annotation
GENOME_URL="https://your-genome-url.fasta.gz"
GTF_URL="https://your-annotation-url.gtf.gz"
```

2. **Adjust quality thresholds** in `2.2-check_tagseq_quality.sh`:
```bash
# Modify based on your organism's typical quality metrics
MIN_MEAN_QUALITY=20
MIN_READS_PER_SAMPLE=500000
```

3. **Update configuration** in `config/pipeline.conf`:
```bash
# Organism-specific settings
SPECIES_NAME="Your organism"
GENOME_SA_INDEX_NBASES=14  # Adjust based on genome size
```

### Supported Organisms

- **Stickleback** (*Gasterosteus aculeatus*) - Default, fully optimized
- **Zebrafish** (*Danio rerio*) - Tested with minor modifications

## 📚 Citations and References

### Software Citations

If you use this pipeline, please cite the underlying tools:

- **STAR**: Dobin, A. et al. (2013) STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 29:15-21
- **FastQC**: Andrews, S. (2010) FastQC: a quality control tool for high throughput sequence data
- **MultiQC**: Ewels, P. et al. (2016) MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics* 32:3047-3048
- **TagSeq**: Meyer, E. et al. (2011) Profiling gene expression responses of coral larvae (Acropora millepora) to elevated temperature and settlement inducers using a novel RNA-Seq procedure. *Mol Ecol* 20:3599-3616


## 🤝 Contributing

This code was originally developed by Rogini Runghen and Grace Vaziri. We welcome contributions! Please see our contributing guidelines:

### How to Contribute

1. **Fork the repository**
2. **Create a feature branch**: `git checkout -b feature/amazing-feature`
3. **Make your changes**: Add features, fix bugs, improve documentation
4. **Test thoroughly**: Run on test data, check all scripts
5. **Commit changes**: `git commit -m 'Add amazing feature'`
6. **Push to branch**: `git push origin feature/amazing-feature`
7. **Open a Pull Request**

### Development Guidelines

- **Follow existing code style**: Consistent with current scripts
- **Add documentation**: Update README and docs/ as needed
- **Include tests**: Provide test cases for new features
- **Update examples**: Add usage examples if applicable

### Reporting Issues

Please include:
- **System information**: OS, SLURM version, cluster details
- **Error logs**: Relevant `.out` and `.err` files
- **Steps to reproduce**: Clear description of the issue
- **Expected behavior**: What should have happened

## 📄 License

This project is licensed under the MIT License.

## 🏛️ Acknowledgments

- **TagSeq method development**: Meyer et al. (2011)
- **STAR aligner development**: Dobin et al. (2013)  
- **Stickleback genome consortium**: Broad Institute
- **University of Connecticut HPC**: Computational resources
- **Bolnick Lab**: Method development and testing

## 📧 Contact

- **Primary Contact**: Rogini Runghen - rogini.runghen@gmail.com / Grace Vaziri
- **Lab Website**: [Lab URL]
- **Issues**: [GitHub Issues Page]

---

**⭐ If this pipeline helps your research, please star this repository and cite the relevant papers!**

*Last updated: 19 August 2025 *
