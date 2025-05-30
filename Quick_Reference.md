# Oligo Design Pipeline - Quick Reference

## 🚀 Getting Started

```bash
# 1. Setup pipeline
./setup_pipeline.sh

# 2. Configure settings
vim config.yaml

# 3. Run pipeline
python3 run_pipeline.py full --cores 8
```

## 📋 Common Commands

### Pipeline Execution
```bash
# Test mode (chr22 only)
python3 run_pipeline.py test --cores 4

# PCR primers only
python3 run_pipeline.py pcr --cores 8

# Capture probes only
python3 run_pipeline.py panel --cores 8

# Full pipeline
python3 run_pipeline.py full --cores 16

# Dry run (check without executing)
python3 run_pipeline.py full --dryrun
```

### Direct Snakemake Commands
```bash
# Basic execution
snakemake --cores 8

# With profile
snakemake --profile config/profile_local

# Cluster execution
snakemake --profile config/profile_slurm

# Force rebuild
snakemake --cores 8 --forceall

# Generate DAG
snakemake --dag | dot -Tpng > dag.png
```

### Status and Monitoring
```bash
# Check status
python3 run_pipeline.py status

# View results
python3 run_pipeline.py view --report-type summary

# List targets
python3 run_pipeline.py targets

# Monitor logs
tail -f logs/*.log
```

## 📊 Key Outputs

| File | Description |
|------|-------------|
| `results/design/pcr/filtered.tsv` | Final PCR primers |
| `results/design/panel/filtered.tsv` | Final capture probes |
| `results/reports/pipeline_summary.html` | Main results report |
| `results/reports/dashboard.html` | Real-time dashboard |
| `results/validation/test_results.txt` | Validation results |

## ⚙️ Configuration Essentials

### Basic Settings (`config.yaml`)
```yaml
genome_build: "hg38"
design_modes: ["pcr", "panel"]
max_threads: 8
max_memory: "32G"

# Adjust module names for your system
modules:
  blast: "blast/2.12.0"
  bowtie2: "bowtie2/2.4.5"
```

### Design Parameters
```yaml
primer3:
  pcr:
    primer_opt_tm: 60.0
    product_size_range: "100-300"
  panel:
    primer_opt_tm: 65.0
    primer_opt_size: 120
```

### Quality Control
```yaml
screening:
  max_off_targets: 10
  min_specificity_score: 0.8
```

## 🔧 Troubleshooting

### Common Issues
1. **Module loading errors**: Update module names in `config.yaml`
2. **Memory errors**: Increase memory allocation
3. **Download failures**: Check internet connectivity
4. **Empty results**: Check log files for errors

### Debug Commands
```bash
# Check syntax
snakemake --lint

# Verbose output
snakemake --cores 4 --verbose --printshellcmds

# Debug specific rule
snakemake --cores 1 --debug results/design/pcr/candidates.tsv
```

### Log Locations
- `logs/` - Rule execution logs
- `logs/slurm/` - SLURM job logs (if using cluster)
- `logs/pipeline_execution.log` - Main pipeline log

## 📈 Performance Tuning

### Resource Allocation
```yaml
# High-memory rules
threads:
  bowtie2_index: 16
  consensus: 8

memory:
  bowtie2_index: "32G"
  blast_index: "16G"
```

### Parallel Execution
```bash
# Multiple jobs
snakemake --cores 16 --jobs 4

# Cluster with job limits
snakemake --profile config/profile_slurm --jobs 100
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
# Full test
snakemake test --cores 4

# Validation
snakemake results/validation/integration_test.txt --cores 2
```

## 📁 Directory Structure

```
oligo-pipeline/
├── Snakefile                 # Main workflow
├── config.yaml              # Configuration
├── setup_pipeline.sh        # Setup script
├── run_pipeline.py          # Launcher script
├── workflow/
│   └── rules/               # Rule definitions
├── config/
│   ├── profile_local/       # Local execution profile
│   └── profile_slurm/       # SLURM cluster profile
├── resources/               # Downloaded references
├── results/                 # All outputs
│   ├── design/             # Design results
│   ├── reports/            # HTML reports
│   └── validation/         # QC results
└── logs/                   # Execution logs
```

## 🎯 Pipeline Targets

### Main Targets
- `all` - Complete pipeline
- `test` - Test mode (chr22)
- `results/design/pcr/filtered.tsv` - PCR primers
- `results/design/panel/filtered.tsv` - Capture probes

### Intermediate Targets
- `resources/hg38.fa` - Reference genome
- `results/indexes/pcr/blast/done.txt` - PCR BLAST index
- `results/validation/test_results.txt` - Validation

### Reports
- `results/reports/pipeline_summary.html` - Main report
- `results/reports/dashboard.html` - Dashboard
- `results/validation/validation_summary.html` - QC report

## 🔄 Workflow Phases

1. **Reference Preparation**
   - Download hg38 and variants
   - Create variant-aware sequences
   - Quality control

2. **Index Construction**
   - Build BLAST databases
   - Build Bowtie2 indexes
   - Validate indexes

3. **Oligo Design**
   - Generate candidates
   - Screen for specificity
   - Filter and rank

4. **Validation**
   - Test known primers
   - Validate outputs
   - Performance metrics

5. **Reporting**
   - Generate HTML reports
   - Create dashboards
   - Export results

## 🛠️ Customization

### Adding New Rules
```python
# In workflow/rules/design.smk
rule custom_filter:
    input: "results/design/{mode}/candidates.tsv"
    output: "results/design/{mode}/custom_filtered.tsv"
    shell: "custom_script.py {input} {output}"
```

### Custom Parameters
```yaml
# In config.yaml
custom_params:
  min_tm: 58.0
  max_gc: 65.0
```

### Additional Outputs
```python
# Export in new format
rule export_json:
    input: "results/design/{mode}/filtered.tsv"
    output: "results/design/{mode}/filtered.json"
    shell: "convert_to_json.py {input} {output}"
```

## 💡 Tips

1. **Start small**: Use test mode first
2. **Monitor resources**: Check memory and disk usage
3. **Save configurations**: Keep project-specific configs
4. **Regular backups**: Archive important results
5. **Version control**: Track configuration changes

## 🆘 Getting Help

- Check log files in `logs/`
- Use `--verbose` flag for detailed output
- Run `snakemake --list` to see available targets
- Validate syntax with `snakemake --lint`
- Generate report with `snakemake --report`

---

*For detailed documentation, see README.md*