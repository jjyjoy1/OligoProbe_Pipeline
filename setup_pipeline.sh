#!/bin/bash

# =============================================================================
# Oligo Design Pipeline Setup Script
# =============================================================================

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Pipeline information
PIPELINE_NAME="Oligo Design Pipeline"
PIPELINE_VERSION="1.0"
REQUIRED_SNAKEMAKE_VERSION="7.0"

# Print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Display header
display_header() {
    echo "============================================================================="
    echo "                      $PIPELINE_NAME v$PIPELINE_VERSION"
    echo "============================================================================="
    echo ""
    echo "This script will set up the directory structure and check dependencies"
    echo "for the oligo design pipeline."
    echo ""
}

# Check dependencies
check_dependencies() {
    print_status "Checking dependencies..."
    
    local missing_deps=0
    
    # Check Python
    if command -v python3 &> /dev/null; then
        PYTHON_VERSION=$(python3 --version | cut -d' ' -f2)
        print_success "Python 3 found: $PYTHON_VERSION"
    else
        print_error "Python 3 not found"
        missing_deps=$((missing_deps + 1))
    fi
    
    # Check Snakemake
    if command -v snakemake &> /dev/null; then
        SNAKEMAKE_VERSION=$(snakemake --version)
        print_success "Snakemake found: $SNAKEMAKE_VERSION"
        
        # Check version compatibility
        if python3 -c "import sys; sys.exit(0 if '$SNAKEMAKE_VERSION' >= '$REQUIRED_SNAKEMAKE_VERSION' else 1)"; then
            print_success "Snakemake version is compatible"
        else
            print_warning "Snakemake version $SNAKEMAKE_VERSION may be incompatible (required: >=$REQUIRED_SNAKEMAKE_VERSION)"
        fi
    else
        print_error "Snakemake not found"
        print_status "Install with: pip install snakemake"
        missing_deps=$((missing_deps + 1))
    fi
    
    # Check wget
    if command -v wget &> /dev/null; then
        print_success "wget found"
    else
        print_error "wget not found (required for downloading reference data)"
        missing_deps=$((missing_deps + 1))
    fi
    
    # Check bc (for calculations)
    if command -v bc &> /dev/null; then
        print_success "bc (calculator) found"
    else
        print_error "bc not found (required for calculations)"
        missing_deps=$((missing_deps + 1))
    fi
    
    # Check environment modules (optional)
    if command -v module &> /dev/null; then
        print_success "Environment modules found"
    else
        print_warning "Environment modules not found (you'll need to adjust tool paths)"
    fi
    
    return $missing_deps
}

# Create directory structure
create_directories() {
    print_status "Creating directory structure..."
    
    local dirs=(
        "workflow/rules"
        "workflow/scripts"
        "workflow/envs"
        "resources"
        "results/design/pcr"
        "results/design/panel"
        "results/indexes/pcr/blast"
        "results/indexes/pcr/bowtie2"
        "results/indexes/panel/blast"
        "results/indexes/panel/bowtie2"
        "results/stats"
        "results/qc"
        "results/validation"
        "results/reports"
        "results/archives"
        "results/maintenance"
        "logs"
        "test/baseline_results"
        "config"
        "scripts"
    )
    
    for dir in "${dirs[@]}"; do
        if [ ! -d "$dir" ]; then
            mkdir -p "$dir"
            print_success "Created directory: $dir"
        else
            print_status "Directory already exists: $dir"
        fi
    done
}

# Check disk space
check_disk_space() {
    print_status "Checking available disk space..."
    
    local available_gb=$(df . --output=avail --block-size=1G | tail -n1 | tr -d ' ')
    local required_gb=50
    
    if [ "$available_gb" -ge "$required_gb" ]; then
        print_success "Sufficient disk space available: ${available_gb}GB (required: ${required_gb}GB)"
    else
        print_warning "Limited disk space: ${available_gb}GB available (recommended: ${required_gb}GB)"
        print_warning "Consider cleaning up space or using a different location"
    fi
}

# Check/create test data
setup_test_data() {
    print_status "Setting up test data..."
    
    # Create minimal test primer list
    if [ ! -f "test/test_primers.tsv" ]; then
        cat > test/test_primers.tsv << EOF
primer_id	type	sequence	expected_hits
GAPDH_F	forward	GTCAACGGATTTGGTCGTATTG	1
GAPDH_R	reverse	CATGGGTGGAATCATATTGGAA	1
ACTB_F	forward	CATGTACGTTGCTATCCAGGC	1
ACTB_R	reverse	CTCCTTAATGTCACGCACGAT	1
EOF
        print_success "Created test primer file: test/test_primers.tsv"
    fi
    
    # Create example input BED file
    if [ ! -f "test/example_targets.bed" ]; then
        cat > test/example_targets.bed << EOF
chr22	16000000	16001000	GENE1_exon1
chr22	16050000	16051000	GENE1_exon2
chr22	17000000	17001000	GENE2_exon1
chr22	17100000	17101000	GENE2_exon2
EOF
        print_success "Created example targets file: test/example_targets.bed"
    fi
}

# Check tool availability
check_tools() {
    print_status "Checking tool availability..."
    
    local tools=(
        "blast:blastn"
        "bowtie2:bowtie2"
        "bcftools:bcftools"
        "samtools:samtools"
        "primer3:primer3_core"
    )
    
    local missing_tools=0
    
    for tool in "${tools[@]}"; do
        local tool_name=$(echo "$tool" | cut -d':' -f1)
        local tool_cmd=$(echo "$tool" | cut -d':' -f2)
        
        if command -v "$tool_cmd" &> /dev/null; then
            print_success "$tool_name found: $tool_cmd"
        else
            # Try with module load
            if command -v module &> /dev/null; then
                if module load "$tool_name" 2>/dev/null && command -v "$tool_cmd" &> /dev/null; then
                    print_success "$tool_name found via module: $tool_cmd"
                    module unload "$tool_name" 2>/dev/null || true
                else
                    print_warning "$tool_name not found (checked command and module)"
                    missing_tools=$((missing_tools + 1))
                fi
            else
                print_warning "$tool_name not found"
                missing_tools=$((missing_tools + 1))
            fi
        fi
    done
    
    if [ $missing_tools -eq 0 ]; then
        print_success "All required tools are available"
    else
        print_warning "$missing_tools tool(s) not found. You may need to:"
        print_warning "  1. Install missing tools"
        print_warning "  2. Update module names in config.yaml"
        print_warning "  3. Adjust PATH environment variable"
    fi
    
    return $missing_tools
}

# Generate example configuration
generate_example_config() {
    print_status "Checking configuration..."
    
    if [ ! -f "config.yaml" ]; then
        print_warning "config.yaml not found, this should be created by the main pipeline"
        return 1
    else
        print_success "config.yaml found"
    fi
    
    # Create a local config override example
    if [ ! -f "config_local.yaml" ]; then
        cat > config_local.yaml << EOF
# Local configuration overrides
# Copy settings from config.yaml and modify as needed

# Example: Adjust for your cluster
#modules:
#  blast: "BLAST+/2.12.0"
#  bowtie2: "Bowtie2/2.4.5"

# Example: Adjust resource limits
#max_threads: 16
#max_memory: "64G"

# Example: Adjust for testing
#chromosomes:
#  - "22"  # Test with chr22 only
EOF
        print_success "Created example local config: config_local.yaml"
    fi
}

# Perform basic validation
validate_setup() {
    print_status "Validating setup..."
    
    # Check for required files
    local required_files=(
        "Snakefile"
        "config.yaml"
        "workflow/rules/reference.smk"
        "workflow/rules/indexing.smk"
        "workflow/rules/design.smk"
        "workflow/rules/validate.smk"
        "workflow/rules/report.smk"
    )
    
    local missing_files=0
    
    for file in "${required_files[@]}"; do
        if [ -f "$file" ]; then
            print_success "Found: $file"
        else
            print_error "Missing: $file"
            missing_files=$((missing_files + 1))
        fi
    done
    
    if [ $missing_files -eq 0 ]; then
        print_success "All required pipeline files found"
    else
        print_error "$missing_files required file(s) missing"
        return 1
    fi
    
    # Test Snakemake syntax
    print_status "Testing Snakemake syntax..."
    if snakemake --dryrun --quiet &> /dev/null; then
        print_success "Snakemake syntax validation passed"
    else
        print_error "Snakemake syntax validation failed"
        print_status "Run 'snakemake --dryrun' for detailed error messages"
        return 1
    fi
    
    return 0
}

# Generate helpful commands
generate_commands() {
    print_status "Generating helpful commands..."
    
    cat > pipeline_commands.sh << 'EOF'
#!/bin/bash
# Helpful commands for the Oligo Design Pipeline

# Basic pipeline execution
alias run_pipeline="snakemake --cores 8"
alias run_test="snakemake test --cores 4"
alias run_pcr="snakemake results/design/pcr/filtered.tsv --cores 4"
alias run_panel="snakemake results/design/panel/filtered.tsv --cores 4"

# Debugging and monitoring
alias check_pipeline="snakemake --dryrun"
alias check_syntax="snakemake --lint"
alias show_dag="snakemake --dag"
alias show_rules="snakemake --list"

# Cleanup commands
alias clean_results="snakemake clean"
alias clean_all="snakemake clean_all"

# Monitoring
alias watch_logs="tail -f logs/*.log"
alias pipeline_status="ls -la results/design/*/filtered.tsv"

# Reports
alias view_report="firefox results/reports/pipeline_summary.html"
alias view_dashboard="firefox results/reports/dashboard.html"

echo "Pipeline commands loaded. Available aliases:"
echo "  run_pipeline, run_test, run_pcr, run_panel"
echo "  check_pipeline, check_syntax, show_dag, show_rules"
echo "  clean_results, clean_all"
echo "  watch_logs, pipeline_status"
echo "  view_report, view_dashboard"
EOF
    
    chmod +x pipeline_commands.sh
    print_success "Created pipeline_commands.sh with helpful aliases"
}

# Main setup function
main() {
    display_header
    
    local overall_status=0
    
    # Run checks
    if ! check_dependencies; then
        print_error "Dependency check failed"
        overall_status=1
    fi
    
    # Create directories
    create_directories
    
    # Check disk space
    check_disk_space
    
    # Set up test data
    setup_test_data
    
    # Check tools
    if ! check_tools; then
        print_warning "Some tools are missing, but setup can continue"
    fi
    
    # Generate configurations
    generate_example_config
    
    # Validate setup
    if ! validate_setup; then
        print_error "Setup validation failed"
        overall_status=1
    fi
    
    # Generate helper commands
    generate_commands
    
    # Final status
    echo ""
    echo "============================================================================="
    if [ $overall_status -eq 0 ]; then
        print_success "Pipeline setup completed successfully!"
        echo ""
        echo "Next steps:"
        echo "  1. Review and customize config.yaml"
        echo "  2. Run: snakemake test --cores 4  (to test with chr22)"
        echo "  3. Run: snakemake --cores 8  (for full pipeline)"
        echo "  4. View results: firefox results/reports/pipeline_summary.html"
        echo ""
        echo "For help: source pipeline_commands.sh"
    else
        print_error "Pipeline setup completed with warnings/errors"
        echo ""
        echo "Please address the issues above before running the pipeline."
        echo "Check the documentation for troubleshooting guidance."
    fi
    echo "============================================================================="
    
    return $overall_status
}

# Handle command line arguments
case "${1:-setup}" in
    "setup"|"")
        main
        ;;
    "check")
        print_status "Running dependency check only..."
        check_dependencies
        check_tools
        ;;
    "clean")
        print_status "Cleaning up setup files..."
        rm -rf results/ resources/ logs/
        print_success "Cleanup completed"
        ;;
    "validate")
        print_status "Running validation only..."
        validate_setup
        ;;
    "help"|"-h"|"--help")
        echo "Oligo Design Pipeline Setup Script"
        echo ""
        echo "Usage: $0 [command]"
        echo ""
        echo "Commands:"
        echo "  setup     - Full setup (default)"
        echo "  check     - Check dependencies only"
        echo "  clean     - Clean up generated files"
        echo "  validate  - Validate setup only"
        echo "  help      - Show this help message"
        echo ""
        ;;
    *)
        print_error "Unknown command: $1"
        echo "Use '$0 help' for usage information"
        exit 1
        ;;
esac
