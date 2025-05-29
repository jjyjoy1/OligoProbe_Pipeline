# =============================================================================
# Pathogen Detection Oligo Design Pipeline - Main Snakefile
# =============================================================================

import os
from pathlib import Path

# Load configuration - use pathogen-specific config
configfile: "config_pathogen.yaml"

# Set up paths
WORKDIR = config.get("workdir", os.getcwd())
RESOURCES = os.path.join(WORKDIR, "resources")
RESULTS = os.path.join(WORKDIR, "results")
LOGS = os.path.join(WORKDIR, "logs")

# Ensure directories exist
os.makedirs(RESOURCES, exist_ok=True)
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(LOGS, exist_ok=True)

# Pathogen-specific variables
PIPELINE_TYPE = config.get("pipeline_type", "pathogen_detection")
TARGET_PATHOGENS = config["target_pathogens"]
DESIGN_MODES = config["design_modes"]

# Get enabled pathogen types
ENABLED_PATHOGEN_TYPES = [ptype for ptype, settings in TARGET_PATHOGENS.items() 
                         if settings.get("enabled", False)]

# =============================================================================
# Target Rules
# =============================================================================

rule all:
    input:
        # Pathogen reference databases
        expand("resources/pathogens/{ptype}/reference_complete.txt", ptype=ENABLED_PATHOGEN_TYPES),
        
        # Host exclusion databases
        "resources/hosts/exclusion_databases_ready.txt",
        
        # Taxonomy integration
        "resources/taxonomy/ncbi_taxonomy_loaded.txt",
        
        # Pathogen-specific indexes
        expand("results/indexes/{mode}/pathogen_blast/done.txt", mode=DESIGN_MODES),
        expand("results/indexes/{mode}/host_exclusion/done.txt", mode=DESIGN_MODES),
        
        # Design outputs with pathogen specificity
        expand("results/design/{mode}/pathogen_candidates.tsv", mode=DESIGN_MODES),
        expand("results/design/{mode}/host_screened.tsv", mode=DESIGN_MODES),
        expand("results/design/{mode}/cross_reactivity_filtered.tsv", mode=DESIGN_MODES),
        
        # Pathogen-specific validation
        "results/validation/pathogen_specificity_test.txt",
        "results/validation/host_exclusion_test.txt",
        "results/validation/cross_reactivity_analysis.txt",
        
        # Enhanced database with pathogen schema
        "results/database/pathogen_database_ready.txt",
        "results/database/taxonomy_integrated.txt",
        "results/database/pathogen_summary.html",
        
        # Clinical validation tracking
        "results/clinical/validation_summary.txt",
        
        # Final pathogen detection report
        "results/reports/pathogen_detection_summary.html"

# Pathogen panel targets
rule pathogen_panels:
    input:
        # Multi-pathogen panels
        "results/panels/respiratory/panel_design.tsv",
        "results/panels/foodborne/panel_design.tsv", 
        "results/panels/antimicrobial_resistance/amr_panel.tsv",
        "results/panels/virulence_factors/virulence_panel.tsv"

# Clinical validation target
rule clinical_validation:
    input:
        "results/clinical/performance_validation.txt",
        "results/clinical/cross_reactivity_validation.txt",
        "results/clinical/regulatory_documentation.txt"

# Test rule for pathogen detection development
rule test_pathogen:
    input:
        "resources/test_pathogens/test_panel.fa",
        "results/indexes/species_specific/pathogen_blast/test_done.txt",
        "results/design/species_specific/test_pathogen_candidates.tsv"

# =============================================================================
# Include rule modules
# =============================================================================

include: "workflow/rules/pathogen_reference.smk"      # Pathogen database acquisition
include: "workflow/rules/host_exclusion.smk"         # Host genome exclusion  
include: "workflow/rules/taxonomy_integration.smk"   # NCBI taxonomy integration
include: "workflow/rules/pathogen_indexing.smk"      # Pathogen-specific indexing
include: "workflow/rules/pathogen_design.smk"        # Pathogen oligo design
include: "workflow/rules/cross_reactivity.smk"       # Cross-reactivity analysis
include: "workflow/rules/pathogen_validation.smk"    # Pathogen-specific validation
include: "workflow/rules/pathogen_database.smk"      # Enhanced database integration
include: "workflow/rules/clinical_validation.smk"    # Clinical validation tracking
include: "workflow/rules/pathogen_panels.smk"        # Multi-pathogen panel design
include: "workflow/rules/pathogen_reports.smk"       # Pathogen-specific reporting

# =============================================================================
# Utility rules for pathogen detection
# =============================================================================

rule clean_pathogen:
    shell:
        """
        rm -rf results/design/*/pathogen_*
        rm -rf results/design/*/host_*
        rm -rf results/design/*/cross_reactivity_*
        rm -rf logs/pathogen_*
        echo "Cleaned pathogen design results and logs"
        """

rule clean_all_pathogen:
    shell:
        """
        rm -rf results/*
        rm -rf resources/pathogens/*
        rm -rf resources/hosts/*
        rm -rf resources/taxonomy/*
        rm -rf logs/*
        echo "Cleaned all pathogen detection files"
        """

rule pathogen_info:
    run:
        print(f"Pipeline Type: {PIPELINE_TYPE}")
        print(f"Enabled Pathogen Types: {', '.join(ENABLED_PATHOGEN_TYPES)}")
        print(f"Design Modes: {', '.join(DESIGN_MODES)}")
        print(f"Working Directory: {WORKDIR}")
        
        # Print target pathogens
        for ptype in ENABLED_PATHOGEN_TYPES:
            pathogens = TARGET_PATHOGENS[ptype].get("priority_pathogens", [])
            print(f"{ptype.title()} Targets ({len(pathogens)}): {', '.join(pathogens[:3])}{'...' if len(pathogens) > 3 else ''}")

# =============================================================================
# Configuration validation for pathogen detection
# =============================================================================

def validate_pathogen_config():
    """Validate pathogen-specific configuration parameters"""
    required_keys = ["pipeline_type", "target_pathogens", "design_modes", "databases"]
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required config key: {key}")
    
    if config["pipeline_type"] != "pathogen_detection":
        print("Warning: Pipeline type is not set to 'pathogen_detection'")
    
    # Validate pathogen types
    valid_pathogen_types = ["bacteria", "viruses", "fungi", "parasites"]
    for ptype in TARGET_PATHOGENS.keys():
        if ptype not in valid_pathogen_types:
            raise ValueError(f"Invalid pathogen type: {ptype}. Must be one of {valid_pathogen_types}")
    
    # Check that at least one pathogen type is enabled
    if not any(TARGET_PATHOGENS[ptype].get("enabled", False) for ptype in TARGET_PATHOGENS):
        raise ValueError("At least one pathogen type must be enabled")
    
    # Validate design modes
    valid_modes = ["species_specific", "genus_specific", "multi_pathogen"]
    for mode in DESIGN_MODES:
        if mode not in valid_modes:
            raise ValueError(f"Invalid design mode: {mode}. Must be one of {valid_modes}")

# Validate config on load
validate_pathogen_config()

# =============================================================================
# Pathogen-specific helper functions
# =============================================================================

def get_pathogen_fasta_for_mode(mode):
    """Get the appropriate pathogen FASTA file for design mode"""
    if mode == "species_specific":
        return "resources/pathogens/species_specific_combined.fa"
    elif mode == "genus_specific":
        return "resources/pathogens/genus_representative.fa"
    elif mode == "multi_pathogen":
        return "resources/pathogens/multi_pathogen_panel.fa"
    else:
        raise ValueError(f"Unknown pathogen design mode: {mode}")

def get_host_exclusion_db(host_type="human"):
    """Get host exclusion database path"""
    return f"resources/hosts/{host_type}_exclusion.fa"

def get_pathogen_threads(rule_name):
    """Get thread count for pathogen-specific rules"""
    pathogen_thread_config = {
        "pathogen_download": 4,
        "taxonomy_processing": 4,
        "pathogen_indexing": 8,
        "host_exclusion": 6,
        "cross_reactivity_screening": 12,
        "phylogenetic_analysis": 4
    }
    
    thread_config = config.get("threads", {})
    return thread_config.get(rule_name, pathogen_thread_config.get(rule_name, 2))

def get_pathogen_memory(rule_name):
    """Get memory allocation for pathogen-specific rules"""
    pathogen_memory_config = {
        "pathogen_download": "4G",
        "taxonomy_processing": "16G", 
        "pathogen_indexing": "32G",
        "host_exclusion": "16G",
        "cross_reactivity_screening": "24G",
        "phylogenetic_analysis": "8G"
    }
    
    memory_config = config.get("memory", {})
    return memory_config.get(rule_name, pathogen_memory_config.get(rule_name, "4G"))

def get_clinical_requirements():
    """Get clinical validation requirements"""
    return config.get("clinical", {
        "performance_targets": {
            "analytical_sensitivity": "10 CFU/ml",
            "analytical_specificity": ">99%",
            "clinical_sensitivity": ">95%",
            "clinical_specificity": ">98%"
        }
    })

def get_pathogen_list_for_type(pathogen_type):
    """Get list of target pathogens for a specific type"""
    if pathogen_type in TARGET_PATHOGENS and TARGET_PATHOGENS[pathogen_type].get("enabled", False):
        return TARGET_PATHOGENS[pathogen_type].get("priority_pathogens", [])
    return []

# =============================================================================
# Pathogen detection workflow configuration
# =============================================================================

# Global workflow settings for pathogen detection
shell.prefix("set -euo pipefail; ")

# Pathogen-specific logging
logging_level = config.get("logging_level", "INFO")

# Resource constraints for pathogen databases
MAX_THREADS = config.get("max_threads", 16)  # Higher for pathogen screening
MAX_MEMORY = config.get("max_memory", "64G")  # More memory for large pathogen databases

# Clinical application settings
CLINICAL_APPLICATION = config.get("clinical", {}).get("intended_use", "research_use_only")
REGULATORY_LEVEL = "basic" if CLINICAL_APPLICATION == "research_use_only" else "clinical"

print(f"Pathogen Detection Pipeline initialized")
print(f"Target pathogen types: {', '.join(ENABLED_PATHOGEN_TYPES)}")
print(f"Design modes: {', '.join(DESIGN_MODES)}")
print(f"Clinical application: {CLINICAL_APPLICATION}")
print(f"Max threads: {MAX_THREADS}")
print(f"Max memory: {MAX_MEMORY}")

# =============================================================================
# Pathogen detection quality control
# =============================================================================

def get_pathogen_qc_thresholds():
    """Get pathogen-specific quality control thresholds"""
    return config.get("qc", {
        "specificity_thresholds": {
            "species_level": 0.95,
            "genus_level": 0.90,
            "family_level": 0.85
        },
        "cross_reactivity": {
            "max_off_targets": 3,
            "max_host_similarity": 0.75,
            "max_related_pathogen_similarity": 0.85
        }
    })

# =============================================================================
# Workflow execution summary
# =============================================================================

def print_pathogen_workflow_summary():
    """Print summary of pathogen detection workflow"""
    print("\n" + "="*80)
    print("PATHOGEN DETECTION WORKFLOW SUMMARY")
    print("="*80)
    
    print(f"Pipeline Type: {PIPELINE_TYPE}")
    print(f"Application: {config.get('application', 'diagnostic')}")
    
    print(f"\nTarget Pathogens:")
    for ptype in ENABLED_PATHOGEN_TYPES:
        pathogens = get_pathogen_list_for_type(ptype)
        print(f"  {ptype.title()}: {len(pathogens)} species")
    
    print(f"\nDesign Modes: {', '.join(DESIGN_MODES)}")
    
    qc = get_pathogen_qc_thresholds()
    print(f"\nQuality Thresholds:")
    print(f"  Species specificity: {qc['specificity_thresholds']['species_level']}")
    print(f"  Max host similarity: {qc['cross_reactivity']['max_host_similarity']}")
    
    clinical = get_clinical_requirements()
    print(f"\nClinical Requirements:")
    print(f"  Analytical sensitivity: {clinical['performance_targets']['analytical_sensitivity']}")
    print(f"  Analytical specificity: {clinical['performance_targets']['analytical_specificity']}")
    
    print("="*80)

# Print workflow summary when pipeline loads
print_pathogen_workflow_summary()

