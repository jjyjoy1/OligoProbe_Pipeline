# =============================================================================
# Oligo Design Pipeline - Main Snakefile
# =============================================================================

import os
from pathlib import Path

# Load configuration
configfile: "config.yaml"

# Set up paths
WORKDIR = config.get("workdir", os.getcwd())
RESOURCES = os.path.join(WORKDIR, "resources")
RESULTS = os.path.join(WORKDIR, "results")
LOGS = os.path.join(WORKDIR, "logs")

# Ensure directories exist
os.makedirs(RESOURCES, exist_ok=True)
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(LOGS, exist_ok=True)

# Global variables
GENOME_BUILD = config["genome_build"]
DESIGN_MODES = config["design_modes"]

# =============================================================================
# Target Rules
# =============================================================================

rule all:
    input:
        # Reference files
        expand("resources/{build}.fa", build=GENOME_BUILD),
        expand("resources/{build}_panel.fa", build=GENOME_BUILD),
        
        # Index files
        expand("results/indexes/{mode}/blast/done.txt", mode=DESIGN_MODES),
        expand("results/indexes/{mode}/bowtie2/done.txt", mode=DESIGN_MODES),
        
        # Design outputs
        expand("results/design/{mode}/candidates.tsv", mode=DESIGN_MODES),
        expand("results/design/{mode}/filtered.tsv", mode=DESIGN_MODES),
        
        # Validation
        "results/validation/test_results.txt",
        
        # Database integration
        "results/database/database_summary.html",
        
        # Final report
        "results/reports/pipeline_summary.html"

# Database-only target
rule database:
    input:
        "results/database/database_summary.html",
        "results/database/exports/all_oligos.fasta"

# Test rule for development
rule test:
    input:
        "resources/test_chr22.fa",
        "results/indexes/pcr/blast/test_done.txt",
        "results/design/pcr/test_candidates.tsv"

# =============================================================================
# Include rule modules
# =============================================================================

include: "workflow/rules/reference.smk"
include: "workflow/rules/indexing.smk" 
include: "workflow/rules/design.smk"
include: "workflow/rules/validate.smk"
include: "workflow/rules/report.smk"
include: "workflow/rules/database.smk"

# =============================================================================
# Utility rules
# =============================================================================

rule clean:
    shell:
        """
        rm -rf results/design/*
        rm -rf logs/*
        echo "Cleaned design results and logs"
        """

rule clean_all:
    shell:
        """
        rm -rf results/*
        rm -rf resources/*.fa*
        rm -rf resources/*.vcf*
        rm -rf logs/*
        echo "Cleaned all generated files"
        """

rule info:
    run:
        print(f"Genome build: {GENOME_BUILD}")
        print(f"Design modes: {', '.join(DESIGN_MODES)}")
        print(f"Working directory: {WORKDIR}")
        print(f"Resources directory: {RESOURCES}")
        print(f"Results directory: {RESULTS}")

# =============================================================================
# Configuration validation
# =============================================================================

def validate_config():
    """Validate configuration parameters"""
    required_keys = ["genome_build", "design_modes", "urls", "primer3"]
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required config key: {key}")
    
    if not isinstance(config["design_modes"], list):
        raise ValueError("design_modes must be a list")
    
    valid_modes = ["pcr", "panel"]
    for mode in config["design_modes"]:
        if mode not in valid_modes:
            raise ValueError(f"Invalid design mode: {mode}. Must be one of {valid_modes}")

# Validate config on load
validate_config()

# =============================================================================
# Helper functions
# =============================================================================

def get_fasta_for_mode(mode):
    """Get the appropriate FASTA file for design mode"""
    if mode == "pcr":
        return f"resources/{GENOME_BUILD}.fa"
    elif mode == "panel":
        return f"resources/{GENOME_BUILD}_panel.fa"
    else:
        raise ValueError(f"Unknown design mode: {mode}")

def get_threads(rule_name):
    """Get thread count for specific rules"""
    thread_config = config.get("threads", {})
    return thread_config.get(rule_name, 1)

def get_memory(rule_name):
    """Get memory allocation for specific rules"""
    memory_config = config.get("memory", {})
    return memory_config.get(rule_name, "4G")

# =============================================================================
# Workflow configuration
# =============================================================================

# Global workflow settings
shell.prefix("set -euo pipefail; ")

# Logging configuration
logging_level = config.get("logging_level", "INFO")

# Resource constraints
MAX_THREADS = config.get("max_threads", 8)
MAX_MEMORY = config.get("max_memory", "32G")

print(f"Oligo Design Pipeline initialized")
print(f"Target genome: {GENOME_BUILD}")
print(f"Design modes: {', '.join(DESIGN_MODES)}")
print(f"Max threads: {MAX_THREADS}")
print(f"Max memory: {MAX_MEMORY}")
