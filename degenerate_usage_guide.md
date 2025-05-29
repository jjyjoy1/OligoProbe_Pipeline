# Degenerate Primer Design Workflow Usage Guide

## Overview

This workflow extends your Snakemake pipeline with advanced degenerate primer design capabilities using:

- **HYDEN/DegePrime**: Design degenerate primers for sequence variants
- **ThermoAlign**: Validate thermodynamic properties of degenerate primers  
- **PriMux**: Optimize multiplex primer combinations

## Workflow Steps

```
Input Sequences → MSA → HYDEN/DegePrime → ThermoAlign → PriMux → Final Design
```

## Setup Instructions

### 1. Add Rules to Your Pipeline

Save the degenerate design rules as `workflow/rules/degenerate.smk` and update your main Snakefile to include:

```python
include: "workflow/rules/degenerate.smk"
```

### 2. Update Configuration

Add the degenerate primer configuration to your `config.yaml` (see configuration template).

### 3. Prepare Input Files

#### Option A: Multiple Sequences per Target

Create FASTA files with sequence variants for each target:

```bash
# Example: BRCA1_variants.fa
>BRCA1_variant_1
ATGGATTTCCGTCTGAGGGCTTCGTGGCCTGGCCCCGCCGGGCCGACCTGGAAGCTGTCCTACAGCCACACC
>BRCA1_variant_2  
ATGGATTTCCGTCTGAGGGCTTCGTGGCCTGGCCCCGCCGGGCCGACCTGGAAGCTGTCCTACAGCCACACT
>BRCA1_variant_3
ATGGATTTCCGTCTGAGGGCTTCGTGGCCTGGCCCCGCCGGGCCGACCTGGAAGCTGTCCTACAGCCACGCC
```

#### Option B: Extract from VCF/BED Files

Use the provided extraction scripts:

```bash
# Extract variant sequences from ClinVar VCF
python vcf_to_fasta.py --vcf clinvar_BRCA1.vcf --reference hg38.fa --output BRCA1_variants.fa --flank 150

# Extract sequences from BED intervals
python bed_to_fasta.py --bed cancer_genes.bed --reference hg38.fa --output cancer_targets.fa --extend 100
```

## Running the Workflow

### Run Complete Pipeline (Standard + Degenerate)

```bash
snakemake --cores 8 all
```

### Run Only Degenerate Primer Design

```bash
snakemake --cores 8 degenerate_only
```

### Run Specific Steps

```bash
# Run only HYDEN design
snakemake --cores 4 results/degenerate/{target}/hyden_candidates.fa

# Run ThermoAlign validation
snakemake --cores 2 results/degenerate/{target}/thermoalign_filtered.fa

# Run PriMux multiplex design
snakemake --cores 4 results/degenerate/multiplex/primux_design.txt
```

## Input Requirements

### File Structure

```
resources/
├── targets/
│   ├── BRCA1_variants.fa    # Multiple sequences for BRCA1
│   ├── BRCA2_variants.fa    # Multiple sequences for BRCA2
│   └── TP53_variants.fa     # Multiple sequences for TP53
└── reference/
    └── hg38.fa              # Reference genome
```

### Sequence Requirements

- **Minimum 3 sequences** per target for meaningful degenerate design
- **Similar length sequences** (±20% length variation recommended)
- **Pre-aligned sequences** (optional - workflow will align if needed)
- **Quality sequences** without excessive Ns or ambiguous bases

## Output Files

### Individual Target Results

```
results/degenerate/{target}/
├── alignment.aln                    # Multiple sequence alignment
├── hyden_candidates.fa             # HYDEN primer candidates
├── hyden_report.txt                # HYDEN design report
├── thermoalign_results.txt         # Thermodynamic validation
├── thermoalign_filtered.fa         # Validated primers
└── thermoalign_report.txt          # Validation report
```

### Multiplex Design Results

```
results/degenerate/multiplex/
├── primux_design.txt               # Multiplex primer sets
├── primux_sets.fa                  # Primer sequences organized by set
├── primux_interactions.txt         # Primer-dimer analysis
└── primux_report.txt               # Multiplex design report
```

### Export Files

```
results/exports/
├── degenerate_primers_all.fa       # All primers in FASTA format
├── multiplex_sets.csv              # Multiplex sets in CSV format
└── primer_ordering_sheet.csv       # Ready for primer synthesis ordering
```

## Understanding the Results

### HYDEN Output

- **Degeneracy**: Number of possible sequences (e.g., degeneracy=4 means 4 variants)
- **Conservation**: Fraction of sequences covered by the primer
- **Tm estimate**: Melting temperature prediction

### ThermoAlign Validation

- **Efficiency**: Fraction of targets successfully amplified
- **Tm**: Refined melting temperature calculation
- **Passed**: Whether primer meets quality criteria

### PriMux Multiplex Design

- **Set ID**: Unique identifier for each multiplex group
- **Tm range**: Melting temperature range within set
- **Interaction score**: Primer-dimer formation potential

## Troubleshooting

### Common Issues

1. **No primers generated**
   - Check sequence alignment quality
   - Reduce degeneracy limits in config
   - Increase conservation threshold tolerance

2. **Poor primer quality**
   - Adjust Tm ranges in configuration
   - Modify GC content constraints
   - Check input sequence quality

3. **No multiplex sets found**
   - Reduce Tm difference limits
   - Increase interaction threshold (less negative)
   - Design more primers per target

### Quality Control Checks

```bash
# Check alignment quality
head -n 20 results/degenerate/{target}/alignment.aln

# Review primer degeneracy
grep "degeneracy" results/degenerate/{target}/hyden_report.txt

# Check multiplex compatibility
cat results/degenerate/multiplex/primux_design.txt
```

## Advanced Usage

### Custom Tool Parameters

Modify `config.yaml` to fine-tune tool behavior:

```yaml
degenerate:
  hyden:
    max_degeneracy: 32      # Reduce for less complex primers
    coverage: 90            # Increase for higher specificity
  
  thermoalign:
    min_efficiency: 0.9     # Increase for stricter validation
    tm_tolerance: 2.0       # Reduce for tighter Tm matching
```

### Integration with Existing Workflows

The degenerate workflow complements standard primer design:

1. **Standard primers** for single targets with Primer3
2. **Degenerate primers** for targets with known variants
3. **Multiplex optimization** for high-throughput applications

### Experimental Validation Workflow

1. **Order primers** using the generated ordering sheet
2. **Test individual primers** with gradient PCR
3. **Optimize multiplex conditions** starting with predicted Tm ranges
4. **Validate specificity** using the provided interaction analysis

## Example Commands

### Complete Workflow for Cancer Gene Panel

```bash
# 1. Extract sequences for cancer genes
python vcf_to_fasta.py --vcf clinvar_cancer_variants.vcf --reference hg38.fa --output cancer_variants.fa --flank 200

# 2. Split by gene and run degenerate design
snakemake --cores 8 degenerate_only

# 3. Review results
cat results/degenerate/multiplex/primux_design.txt
cat results/exports/primer_ordering_sheet.csv
```

### Quick Test with Sample Data

```bash
# Test with provided examples
snakemake --cores 2 results/degenerate/test_target/hyden_candidates.fa --dry-run
```

## Best Practices

1. **Start small**: Test with 2-3 targets before scaling up
2. **Validate computationally**: Review all QC metrics before ordering primers
3. **Plan experiments**: Use multiplex predictions to design validation experiments
4. **Document parameters**: Keep track of configuration settings for reproducibility
5. **Iterative design**: Adjust parameters based on experimental results

For questions or issues, check the log files in `logs/` directory for detailed error messages and processing information.