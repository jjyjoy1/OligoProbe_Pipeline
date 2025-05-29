# =============================================================================
# Pathogen Reference Database Acquisition Rules
# =============================================================================

# Download and prepare pathogen reference databases
rule download_pathogen_genomes:
    output:
        bacteria="resources/pathogens/bacteria/genomes_downloaded.txt",
        viruses="resources/pathogens/viruses/genomes_downloaded.txt" if "viruses" in ENABLED_PATHOGEN_TYPES else [],
        fungi="resources/pathogens/fungi/genomes_downloaded.txt" if "fungi" in ENABLED_PATHOGEN_TYPES else []
    params:
        target_pathogens=TARGET_PATHOGENS,
        databases=config["databases"]["pathogens"]
    log:
        "logs/download_pathogen_genomes.log"
    threads: get_pathogen_threads("pathogen_download")
    resources:
        mem_mb=lambda w: get_pathogen_memory("pathogen_download").replace("G", "000").replace("M", "")
    script:
        "scripts/download_pathogen_genomes.py"

# Download NCBI Taxonomy database
rule download_taxonomy:
    output:
        taxonomy_dir=directory("resources/taxonomy"),
        nodes="resources/taxonomy/nodes.dmp",
        names="resources/taxonomy/names.dmp",
        complete="resources/taxonomy/ncbi_taxonomy_loaded.txt"
    params:
        taxonomy_url=config["databases"]["taxonomy"]["ncbi_taxonomy"]
    log:
        "logs/download_taxonomy.log"
    threads: get_pathogen_threads("taxonomy_processing")
    shell:
        """
        echo "Downloading NCBI Taxonomy database..." > {log}
        
        mkdir -p {output.taxonomy_dir}
        cd {output.taxonomy_dir}
        
        # Download taxonomy dump
        wget -O taxdump.tar.gz {params.taxonomy_url} 2>> ../{log}
        
        # Extract taxonomy files
        tar -xzf taxdump.tar.gz 2>> ../{log}
        
        # Verify required files exist
        if [ -f nodes.dmp ] && [ -f names.dmp ]; then
            echo "Taxonomy download completed successfully" > {output.complete}
            echo "Files: nodes.dmp ($(wc -l < nodes.dmp) lines), names.dmp ($(wc -l < names.dmp) lines)" >> {output.complete}
        else
            echo "Error: Required taxonomy files not found" >> ../{log}
            exit 1
        fi
        
        echo "Taxonomy download completed: $(date)" >> ../{log}
        """

# Process bacterial 16S/23S databases
rule download_bacterial_markers:
    output:
        silva_16s="resources/pathogens/bacteria/silva_16s.fasta",
        silva_processed="resources/pathogens/bacteria/silva_16s_processed.fasta"
    params:
        silva_url=config["databases"]["pathogens"]["bacteria"]["silva_16s"]
    log:
        "logs/download_bacterial_markers.log"
    threads: 4
    shell:
        """
        echo "Downloading bacterial marker gene databases..." > {log}
        
        mkdir -p resources/pathogens/bacteria
        
        # Download SILVA 16S database
        wget -O {output.silva_16s}.gz {params.silva_url} 2>> {log}
        gunzip {output.silva_16s}.gz 2>> {log}
        
        # Process SILVA sequences (clean headers, filter by length)
        python3 << 'PYTHON_EOF' 2>> {log}
import re
from Bio import SeqIO

def clean_silva_header(header):
    """Clean SILVA header to extract organism name and taxonomy"""
    # SILVA format: >XXXXX.1.1373 Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;
    parts = header.split(' ', 1)
    accession = parts[0].replace('>', '')
    
    if len(parts) > 1:
        taxonomy = parts[1].strip()
        # Extract species name (last two taxonomic levels)
        tax_parts = taxonomy.rstrip(';').split(';')
        if len(tax_parts) >= 2:
            species = f"{tax_parts[-2]} {tax_parts[-1]}"
        else:
            species = tax_parts[-1] if tax_parts else "Unknown"
    else:
        species = "Unknown"
    
    return f"{accession}|{species}"

# Process SILVA 16S sequences
processed_count = 0
with open("{output.silva_processed}", "w") as out_file:
    for record in SeqIO.parse("{output.silva_16s}", "fasta"):
        # Filter by length (16S should be ~1500bp)
        if 1200 <= len(record.seq) <= 1800:
            # Clean header
            clean_header = clean_silva_header(record.description)
            out_file.write(f">{clean_header}\\n{record.seq}\\n")
            processed_count += 1

print(f"Processed {processed_count} 16S sequences")
PYTHON_EOF
        
        echo "Bacterial marker download completed: $(date)" >> {log}
        """

# Create pathogen-specific reference compilations
rule compile_pathogen_references:
    input:
        bacteria_genomes="resources/pathogens/bacteria/genomes_downloaded.txt",
        silva_16s="resources/pathogens/bacteria/silva_16s_processed.fasta",
        taxonomy_complete="resources/taxonomy/ncbi_taxonomy_loaded.txt"
    output:
        species_specific="resources/pathogens/species_specific_combined.fa",
        genus_representative="resources/pathogens/genus_representative.fa",
        multi_pathogen="resources/pathogens/multi_pathogen_panel.fa",
        metadata="resources/pathogens/pathogen_metadata.tsv"
    params:
        target_pathogens=TARGET_PATHOGENS,
        panels=config.get("panels", {})
    log:
        "logs/compile_pathogen_references.log"
    threads: get_pathogen_threads("pathogen_indexing")
    script:
        "scripts/compile_pathogen_references.py"

# Download antimicrobial resistance gene databases
rule download_amr_databases:
    output:
        card_db="resources/amr/card_database.fasta",
        resfinder_db="resources/amr/resfinder_database.fasta",
        amr_metadata="resources/amr/amr_metadata.tsv"
    params:
        amr_databases=config["databases"]["amr"]
    log:
        "logs/download_amr_databases.log"
    threads: 4
    shell:
        """
        echo "Downloading AMR gene databases..." > {log}
        
        mkdir -p resources/amr
        
        # Download CARD database
        echo "Downloading CARD database..." >> {log}
        wget -O resources/amr/card_data.tar.bz2 {params.amr_databases[card_database]} 2>> {log}
        cd resources/amr
        tar -xjf card_data.tar.bz2 2>> ../{log}
        
        # Extract sequences from CARD
        if [ -f card.json ]; then
            python3 << 'PYTHON_EOF' 2>> ../{log}
import json
import re

# Parse CARD JSON and extract sequences
with open('card.json', 'r') as f:
    card_data = json.load(f)

with open('card_database.fasta', 'w') as fasta_out, open('amr_metadata.tsv', 'w') as meta_out:
    # Write metadata header
    meta_out.write("gene_id\\tgene_name\\tresistance_mechanism\\tantibiotic_class\\torganism\\n")
    
    for aro_id, entry in card_data.items():
        if aro_id.isdigit() and 'model_sequences' in entry:
            gene_name = entry.get('model_name', 'Unknown')
            resistance_mechanism = entry.get('ARO_description', 'Unknown')
            
            # Extract sequences
            for seq_id, seq_data in entry['model_sequences'].items():
                if 'dna_sequence' in seq_data:
                    sequence = seq_data['dna_sequence']['sequence']
                    organism = seq_data.get('NCBI_taxonomy', {}).get('NCBI_taxonomy_name', 'Unknown')
                    
                    # Write FASTA
                    fasta_out.write(f">CARD:{aro_id}|{gene_name}|{organism}\\n{sequence}\\n")
                    
                    # Write metadata
                    meta_out.write(f"{aro_id}\\t{gene_name}\\t{resistance_mechanism}\\tbeta-lactam\\t{organism}\\n")

print("CARD database processing completed")
PYTHON_EOF
        else
            echo "Warning: CARD JSON file not found, creating empty database" >> ../{log}
            touch card_database.fasta
            echo -e "gene_id\\tgene_name\\tresistance_mechanism\\tantibiotic_class\\torganism" > amr_metadata.tsv
        fi
        
        cd ../..
        
        # Create placeholder for ResFinder (requires separate download)
        touch {output.resfinder_db}
        
        echo "AMR database download completed: $(date)" >> {log}
        """

# Download virulence factor databases
rule download_virulence_databases:
    output:
        vfdb="resources/virulence/vfdb_sequences.fasta",
        virulence_metadata="resources/virulence/virulence_metadata.tsv"
    params:
        virulence_databases=config["databases"]["virulence"]
    log:
        "logs/download_virulence_databases.log"
    threads: 4
    shell:
        """
        echo "Downloading virulence factor databases..." > {log}
        
        mkdir -p resources/virulence
        
        # Download VFDB
        echo "Downloading VFDB..." >> {log}
        wget -O {output.vfdb}.gz {params.virulence_databases[vfdb]} 2>> {log}
        gunzip {output.vfdb}.gz 2>> {log}
        
        # Process VFDB headers and create metadata
        python3 << 'PYTHON_EOF' 2>> {log}
import re
from Bio import SeqIO

# Parse VFDB and create metadata
with open("{output.virulence_metadata}", "w") as meta_out:
    meta_out.write("vf_id\\tfactor_name\\tfactor_type\\torganism\\tfunction\\n")
    
    for record in SeqIO.parse("{output.vfdb}", "fasta"):
        # VFDB header format: >VFG000001(gb|AAC73113) (sdiA) von Willebrand factor type A domain protein [Adherence (VF0044)] [Escherichia coli CFT073]
        header = record.description
        
        # Extract VF ID
        vf_match = re.search(r'(VFG\\d+)', header)
        vf_id = vf_match.group(1) if vf_match else "Unknown"
        
        # Extract factor name
        name_match = re.search(r'\\(([^)]+)\\)', header)
        factor_name = name_match.group(1) if name_match else "Unknown"
        
        # Extract organism
        org_match = re.search(r'\\[([^\\]]+)\\]$', header)
        organism = org_match.group(1) if org_match else "Unknown"
        
        # Extract function/type
        func_match = re.search(r'\\[([^\\]]+)\\]\\s*\\[', header)
        function = func_match.group(1) if func_match else "Unknown"
        
        # Determine factor type based on function
        if "toxin" in function.lower():
            factor_type = "toxin"
        elif "adherence" in function.lower():
            factor_type = "adhesin"
        elif "invasion" in function.lower():
            factor_type = "invasin"
        else:
            factor_type = "other"
        
        meta_out.write(f"{vf_id}\\t{factor_name}\\t{factor_type}\\t{organism}\\t{function}\\n")

print("Virulence factor database processing completed")
PYTHON_EOF
        
        echo "Virulence database download completed: $(date)" >> {log}
        """

# Create comprehensive pathogen database status
rule pathogen_reference_status:
    input:
        expand("resources/pathogens/{ptype}/genomes_downloaded.txt", ptype=ENABLED_PATHOGEN_TYPES),
        "resources/pathogens/species_specific_combined.fa",
        "resources/taxonomy/ncbi_taxonomy_loaded.txt",
        "resources/amr/amr_metadata.tsv",
        "resources/virulence/virulence_metadata.tsv"
    output:
        expand("resources/pathogens/{ptype}/reference_complete.txt", ptype=ENABLED_PATHOGEN_TYPES),
        "resources/pathogen_databases_ready.txt"
    log:
        "logs/pathogen_reference_status.log"
    shell:
        """
        echo "Checking pathogen database status..." > {log}
        
        # Create completion markers for each pathogen type
        for ptype in {ENABLED_PATHOGEN_TYPES}; do
            echo "Pathogen databases ready for $ptype" > resources/pathogens/$ptype/reference_complete.txt
            echo "Generated: $(date)" >> resources/pathogens/$ptype/reference_complete.txt
            
            # Count sequences
            if [ -f "resources/pathogens/species_specific_combined.fa" ]; then
                SEQ_COUNT=$(grep -c "^>" resources/pathogens/species_specific_combined.fa)
                echo "Total sequences: $SEQ_COUNT" >> resources/pathogens/$ptype/reference_complete.txt
            fi
        done
        
        # Create overall status file
        cat > {output[1]} << EOF
Pathogen Reference Databases Status
===================================
Generated: $(date)

Enabled Pathogen Types: {ENABLED_PATHOGEN_TYPES}

Database Files:
- Species-specific: $([ -f resources/pathogens/species_specific_combined.fa ] && echo "Ready" || echo "Missing")
- Genus representative: $([ -f resources/pathogens/genus_representative.fa ] && echo "Ready" || echo "Missing")
- Multi-pathogen panel: $([ -f resources/pathogens/multi_pathogen_panel.fa ] && echo "Ready" || echo "Missing")
- NCBI Taxonomy: $([ -f resources/taxonomy/nodes.dmp ] && echo "Ready" || echo "Missing")
- AMR Database: $([ -f resources/amr/amr_metadata.tsv ] && echo "Ready" || echo "Missing")
- Virulence Database: $([ -f resources/virulence/virulence_metadata.tsv ] && echo "Ready" || echo "Missing")

Sequence Counts:
EOF
        
        # Add sequence counts
        for fa_file in resources/pathogens/*.fa; do
            if [ -f "$fa_file" ]; then
                COUNT=$(grep -c "^>" "$fa_file")
                echo "- $(basename $fa_file): $COUNT sequences" >> {output[1]}
            fi
        done
        
        echo "Pathogen reference status check completed: $(date)" >> {log}
        """

# Quality control for pathogen databases
rule qc_pathogen_databases:
    input:
        "resources/pathogen_databases_ready.txt"
    output:
        "results/qc/pathogen_database_qc.txt"
    log:
        "logs/qc_pathogen_databases.log"
    shell:
        """
        echo "Running quality control on pathogen databases..." > {log}
        
        mkdir -p results/qc
        
        cat > {output} << EOF
Pathogen Database Quality Control Report
========================================
Generated: $(date)

Database Integrity Checks:
EOF
        
        # Check file integrity
        echo "File Integrity:" >> {output}
        
        for db_file in resources/pathogens/*.fa resources/amr/*.fasta resources/virulence/*.fasta; do
            if [ -f "$db_file" ]; then
                # Check if file is valid FASTA
                HEADER_COUNT=$(grep -c "^>" "$db_file" 2>/dev/null || echo "0")
                FILE_SIZE=$(stat -c%s "$db_file" 2>/dev/null || echo "0")
                
                if [ "$HEADER_COUNT" -gt 0 ] && [ "$FILE_SIZE" -gt 1000 ]; then
                    echo "✓ $(basename $db_file): $HEADER_COUNT sequences, $FILE_SIZE bytes" >> {output}
                else
                    echo "✗ $(basename $db_file): Invalid or empty file" >> {output}
                fi
            fi
        done
        
        # Check for duplicates
        echo "" >> {output}
        echo "Duplicate Analysis:" >> {output}
        
        if [ -f "resources/pathogens/species_specific_combined.fa" ]; then
            TOTAL_SEQS=$(grep -c "^>" resources/pathogens/species_specific_combined.fa)
            UNIQUE_SEQS=$(grep "^>" resources/pathogens/species_specific_combined.fa | sort | uniq | wc -l)
            DUPLICATES=$((TOTAL_SEQS - UNIQUE_SEQS))
            
            echo "Total sequences: $TOTAL_SEQS" >> {output}
            echo "Unique headers: $UNIQUE_SEQS" >> {output}
            echo "Potential duplicates: $DUPLICATES" >> {output}
        fi
        
        echo "" >> {output}
        echo "Quality Control Summary:" >> {output}
        if grep -q "✗" {output}; then
            echo "Status: ISSUES FOUND - Review failures above" >> {output}
        else
            echo "Status: ALL CHECKS PASSED" >> {output}
        fi
        
        echo "Pathogen database QC completed: $(date)" >> {log}
        """
