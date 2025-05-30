# Oligo Design Pipeline

A comprehensive Snakemake pipeline for designing PCR primers and hybrid capture probes with dual-reference specificity screening.

## 🚀 Quick Start

```bash
# Clone/download the pipeline
# Set up directory structure
./setup_pipeline.sh

# Configure your settings
vim config.yaml

# Run the complete pipeline
snakemake --cores 8

# Or run in test mode
snakemake test --cores 4
```

## 📋 Features

- **Dual Design Modes**: PCR primers and hybrid capture probes
- **Variant-Aware Design**: Uses population variants for enhanced specificity
- **Comprehensive Screening**: BLAST and Bowtie2 off-target analysis
- **Quality Control**: Built-in validation and regression testing
- **Rich Reporting**: HTML dashboards and detailed reports
- **Scalable**: Runs on local servers with configurable resources

## 🏗️ Architecture

### Design Modes

1. **PCR Mode**: Standard primer pairs for PCR amplification
   - Uses clean reference genome (hg38)
   - Optimized for primer pairs with product size constraints
   - Fast specificity screening

2. **Panel Mode**: Capture probes for hybrid capture panels
   - Uses variant-aware reference with IUPAC codes
   - Longer probe sequences (80-160bp)
   - Enhanced specificity against variant backgrounds

### Pipeline Phases

1. **Reference Preparation** (`workflow/rules/reference.smk`)
   - Download hg38 reference and gnomAD variants
   - Create variant-aware consensus sequences
   - Quality control and indexing

2. **Index Construction** (`workflow/rules/indexing.smk`)
   - Build BLAST and Bowtie2 indexes for both modes
   - Validation and performance monitoring

3. **Oligo Design** (`workflow/rules/design.smk`)
   - Primer3-based candidate generation
   - Dual-mode specificity screening
   - Thermodynamic and quality filtering

4. **Validation** (`workflow/rules/validate.smk`)
   - Test primer validation
   - Pipeline output verification
   - Performance benchmarking

5. **Reporting** (`workflow/rules/report.smk`)
   - Interactive HTML reports
   - Real-time dashboard
   - Export in multiple formats

## 📦 Installation

### Prerequisites

- **Snakemake** (≥7.0)
- **Python** (≥3.8) with pandas, biopython
- **Environment Modules** or direct tool installation

### Required Tools

- BLAST+ (≥2.12.0)
- Bowtie2 (≥2.4.5)
- BCFtools (≥1.17)
- SAMtools (≥1.17)
- Primer3 (≥2.6.1)

### Setup

1. **Download Pipeline**
   ```bash
   git clone <repository-url>
   cd oligo-pipeline
   ```

2. **Create Directory Structure**
   ```bash
   ./setup_pipeline.sh
   ```

3. **Configure Environment**
   ```bash
   # Edit config.yaml to match your system
   vim config.yaml
   
   # Update module names for your cluster
   # Adjust resource allocations
   # Set design parameters
   ```

## ⚙️ Configuration

### Basic Settings (`config.yaml`)

```yaml
# Genome and design modes
genome_build: "hg38"
design_modes: ["pcr", "panel"]

# Resource limits
max_threads: 8
max_memory: "32G"

# Tool modules (adjust for your system)
modules:
  blast: "blast/2.12.0"
  bowtie2: "bowtie2/2.4.5"
  # ... etc
```

### Design Parameters

```yaml
primer3:
  pcr:
    primer_opt_size: 20
    primer_opt_tm: 60.0
    product_size_range: "100-300"
  panel:
    primer_opt_size: 120
    primer_opt_tm: 65.0
```

### Quality Control

```yaml
screening:
  max_off_targets: 10
  min_specificity_score: 0.8

validation:
  max_acceptable_off_targets: 5
```

## 🚦 Usage

### Basic Execution

```bash
# Full pipeline
snakemake --cores 8

# Specific targets
snakemake results/design/pcr/filtered.tsv --cores 4
snakemake results/reports/pipeline_summary.html --cores 2

# Test mode (chr22 only)
snakemake test --cores 4
```

### Advanced Options

```bash
# Dry run (check what will be executed)
snakemake -n

# Force rebuild of specific components
snakemake --forcerun results/indexes/pcr/blast/done.txt

# With detailed logging
snakemake --cores 8 --verbose --reason

# Generate execution report
snakemake --cores 8 --report report.html
```

### Cluster Execution

```bash
# SLURM example
snakemake --cores 32 --cluster "sbatch -t 4:00:00 -c {threads} --mem={resources.mem_mb}"

# SGE example  
snakemake --cores 16 --cluster "qsub -pe smp {threads} -l h_vmem={resources.mem_mb}M"
```

## 📊 Outputs

### Primary Results

- `results/design/pcr/filtered.tsv` - Final PCR primer candidates
- `results/design/panel/filtered.tsv` - Final capture probe candidates
- `results/design/*/filtered.fa` - FASTA format sequences
- `results/design/*/filtered.bed` - Genomic coordinates

### Reports and Analysis

- `results/reports/pipeline_summary.html` - Main results dashboard
- `results/reports/dashboard.html` - Real-time monitoring
- `results/reports/*_detailed_report.html` - Mode-specific analysis
- `results/validation/validation_summary.html` - QC results

### Intermediate Files

- `resources/` - Reference genomes and variants
- `results/indexes/` - BLAST and Bowtie2 databases
- `results/stats/` - Performance and design statistics
- `logs/` - Execution logs for debugging

## 🔍 File Formats

### TSV Output Columns

| Column | Description |
|--------|-------------|
| primer_id | Unique identifier |
| type | forward/reverse/probe |
| sequence | Oligonucleotide sequence |
| length | Sequence length (bp) |
| tm | Melting temperature (°C) |
| gc_content | GC percentage |
| off_targets | Number of off-target hits |
| specificity_score | Specificity score (0-1) |

### FASTA Headers

```
>PCR_001_F length=20 tm=60.2 gc=55.0 off_targets=2
GCATACGTTGTATCCGGGCAT
```

## 🛠️ Customization

### Adding New Design Algorithms

1. Create new rule in `workflow/rules/design.smk`:
   ```python
   rule custom_design:
       input: "resources/{build}.fa"
       output: "results/design/{mode}/custom_candidates.tsv"
       shell: "custom_tool --input {input} --output {output}"
   ```

2. Update main workflow to include new rule

### Custom Quality Filters

Modify the filtering script in `filter_candidates` rule:
```python
# Add custom filtering logic
def custom_filter(row):
    return row['tm'] > 58 and row['gc_content'] < 65
```

### Additional Output Formats

Add new export rules in `workflow/rules/design.smk`:
```python
rule export_json:
    input: "results/design/{mode}/filtered.tsv"
    output: "results/design/{mode}/filtered.json"
    shell: "convert_to_json.py {input} {output}"
```

## 🧪 Testing

### Unit Tests

```bash
# Test individual components
snakemake results/qc/reference_qc_hg38.txt --cores 2
snakemake results/validation/test_results.txt --cores 2
```

### Integration Tests

```bash
# Full pipeline test with small dataset
snakemake test --cores 4

# Regression testing
snakemake results/validation/regression_test.txt --cores 2
```

### Performance Benchmarking

```bash
# Generate performance report
snakemake results/validation/performance_benchmark.txt --cores 4
```

## 🐛 Troubleshooting

### Common Issues

1. **Module Loading Errors**
   ```bash
   # Check available modules
   module avail blast
   
   # Update config.yaml with correct module names
   modules:
     blast: "blast/2.12.0"  # Adjust version
   ```

2. **Memory Errors**
   ```bash
   # Increase memory allocation in config.yaml
   memory:
     bowtie2_index: "32G"  # Increase from default
   ```

3. **Download Failures**
   ```bash
   # Check internet connectivity
   wget -O test.fa.gz https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
   
   # Alternative: manually download and place in resources/
   ```

4. **Empty Results**
   ```bash
   # Check log files
   tail -f logs/generate_candidates_pcr.log
   
   # Validate input files
   snakemake results/qc/reference_qc_hg38.txt
   ```

### Debug Mode

```bash
# Enable debug output
snakemake --cores 4 --verbose --printshellcmds

# Check individual rule execution
snakemake --cores 1 --debug-dag results/design/pcr/candidates.tsv
```

### Log Analysis

```bash
# Search for errors across all logs
grep -r "ERROR\|FAIL" logs/

# Check resource usage
grep -r "memory\|time" logs/

# Monitor progress
tail -f logs/*.log
```

## 📈 Performance Optimization

### Resource Tuning

1. **Thread Allocation**
   ```yaml
   threads:
     bowtie2_index: 16  # Increase for faster indexing
     screening: 8       # Balance with memory usage
   ```

2. **Memory Management**
   ```yaml
   memory:
     blast_index: "16G"    # Adjust based on genome size
     consensus: "12G"      # For variant processing
   ```

### Workflow Optimization

1. **Parallel Execution**
   ```bash
   # Run multiple modes simultaneously
   snakemake --cores 16 --jobs 4
   ```

2. **Intermediate Cleanup**
   ```bash
   # Remove temporary files automatically
   snakemake --cores 8 --delete-temp-output
   ```

3. **Checkpoint Optimization**
   ```bash
   # Use workflow profiles for different scenarios
   snakemake --profile config/slurm/
   ```

## 📚 Advanced Usage

### Custom Reference Genomes

1. **Alternative Builds**
   ```yaml
   genome_build: "hg19"  # Change in config.yaml
   urls:
     reference_fasta: "https://path/to/hg19.fa.gz"
   ```

2. **Non-Human Species**
   ```yaml
   genome_build: "mm10"
   urls:
     reference_fasta: "https://path/to/mouse/genome.fa.gz"
     # Update variant sources accordingly
   ```

### Batch Processing

```python
# Create config files for multiple projects
import yaml

projects = ["BRCA_panel", "cardiac_genes", "cancer_hotspots"]
for project in projects:
    config = load_base_config()
    config["project_name"] = project
    # Customize parameters per project
    save_config(f"config_{project}.yaml", config)
```

### Integration with LIMS

```python
# Example REST API integration
rule upload_results:
    input: "results/design/{mode}/filtered.tsv"
    shell: """
    curl -X POST -F "file=@{input}" \
         https://lims.example.com/api/upload_oligos
    """
```

## 🤝 Contributing

### Development Setup

```bash
# Fork repository and create development branch
git checkout -b feature/new-algorithm

# Install development dependencies
pip install snakemake pandas biopython pytest

# Run tests
pytest tests/
```

### Code Style

- Follow PEP 8 for Python code
- Use descriptive variable names
- Add comments for complex logic
- Update documentation for new features

### Submitting Changes

1. Create feature branch
2. Add tests for new functionality
3. Update documentation
4. Submit pull request with clear description

## 📄 License

This pipeline is released under the MIT License. See LICENSE file for details.

## 🆘 Support

- **Issues**: Report bugs and feature requests via GitHub Issues
- **Documentation**: Check the wiki for detailed examples
- **Community**: Join the discussion forum for questions

## 📖 Citation

If you use this pipeline in your research, please cite:

```
Oligo Design Pipeline v1.0
A Snakemake workflow for PCR primer and capture probe design
https://github.com/your-repo/oligo-pipeline
```

## 🔄 Version History

- **v1.0** - Initial release with PCR and panel modes
- **v1.1** - Added variant-aware design and enhanced reporting
- **v1.2** - Performance optimizations and cluster support

---

*Generated with ❤️ by the Oligo Design Pipeline*