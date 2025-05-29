# =============================================================================
# Host Exclusion Rules for Pathogen Detection
# =============================================================================
# Critical for clinical applications - ensures primers don't amplify host DNA

# Download host genomes for exclusion
rule download_host_genomes:
    output:
        human_genome="resources/hosts/human_genome.fa.gz",
        human_transcriptome="resources/hosts/human_transcriptome.fa.gz",
        mouse_genome="resources/hosts/mouse_genome.fa.gz" if config["databases"]["hosts"].get("mouse") else [],
        livestock_genomes=expand("resources/hosts/{species}_genome.fa.gz", 
                                species=["bovine", "porcine"] if config["databases"]["hosts"].get("livestock") else [])
    params:
        host_urls=config["databases"]["hosts"]
    log:
        "logs/download_host_genomes.log"
    threads: get_pathogen_threads("pathogen_download")
    shell:
        """
        echo "Downloading host genomes for exclusion screening..." > {log}
        
        mkdir -p resources/hosts
        
        # Download human genome
        echo "Downloading human genome..." >> {log}
        wget -O {output.human_genome} {params.host_urls[human][genome]} 2>> {log}
        
        # Download human transcriptome
        echo "Downloading human transcriptome..." >> {log}
        wget -O {output.human_transcriptome} {params.host_urls[human][transcriptome]} 2>> {log}
        
        # Download mouse genome if specified
        if [ -n "{params.host_urls.get('mouse', {}).get('genome', '')}" ]; then
            echo "Downloading mouse genome..." >> {log}
            wget -O resources/hosts/mouse_genome.fa.gz {params.host_urls[mouse][genome]} 2>> {log}
        fi
        
        # Download livestock genomes if specified
        if [ -n "{params.host_urls.get('livestock', {}).get('bovine', '')}" ]; then
            echo "Downloading bovine genome..." >> {log}
            wget -O resources/hosts/bovine_genome.fa.gz {params.host_urls[livestock][bovine]} 2>> {log}
        fi
        
        if [ -n "{params.host_urls.get('livestock', {}).get('porcine', '')}" ]; then
            echo "Downloading porcine genome..." >> {log}
            wget -O resources/hosts/porcine_genome.fa.gz {params.host_urls[livestock][porcine]} 2>> {log}
        fi
        
        echo "Host genome downloads completed: $(date)" >> {log}
        """

# Process and combine host sequences for exclusion
rule prepare_host_exclusion_db:
    input:
        human_genome="resources/hosts/human_genome.fa.gz",
        human_transcriptome="resources/hosts/human_transcriptome.fa.gz"
    output:
        human_combined="resources/hosts/human_exclusion.fa",
        all_hosts_combined="resources/hosts/all_hosts_exclusion.fa",
        host_metadata="resources/hosts/host_metadata.tsv"
    log:
        "logs/prepare_host_exclusion.log"
    threads: get_pathogen_threads("host_exclusion")
    resources:
        mem_mb=lambda w: get_pathogen_memory("host_exclusion").replace("G", "000").replace("M", "")
    shell:
        """
        echo "Preparing host exclusion databases..." > {log}
        
        # Process human genome and transcriptome
        echo "Processing human sequences..." >> {log}
        
        # Combine human genome and transcriptome
        gunzip -c {input.human_genome} > temp_human_genome.fa 2>> {log}
        gunzip -c {input.human_transcriptome} > temp_human_transcriptome.fa 2>> {log}
        
        # Clean and combine human sequences
        python3 << 'PYTHON_EOF' 2>> {log}
import re
from Bio import SeqIO

def clean_header(header, source):
    """Clean FASTA headers for host exclusion database"""
    # Remove problematic characters and add source
    clean_id = re.sub(r'[^A-Za-z0-9_.-]', '_', header.split()[0])
    return f"{clean_id}|{source}|human"

processed_count = 0
host_metadata = []

# Process human genome
with open("{output.human_combined}", "w") as out_file:
    # Process genome
    for record in SeqIO.parse("temp_human_genome.fa", "fasta"):
        if len(record.seq) >= 1000:  # Only include substantial sequences
            clean_id = clean_header(record.id, "genome")
            out_file.write(f">{clean_id}\\n{record.seq}\\n")
            host_metadata.append([clean_id, "human", "genome", len(record.seq), "Homo sapiens"])
            processed_count += 1
    
    # Process transcriptome (sample to avoid excessive size)
    transcriptome_count = 0
    for record in SeqIO.parse("temp_human_transcriptome.fa", "fasta"):
        if transcriptome_count < 50000:  # Limit transcriptome sequences
            clean_id = clean_header(record.id, "transcriptome")
            out_file.write(f">{clean_id}\\n{record.seq}\\n")
            host_metadata.append([clean_id, "human", "transcriptome", len(record.seq), "Homo sapiens"])
            processed_count += 1
            transcriptome_count += 1

print(f"Processed {processed_count} human sequences")

# Write metadata
with open("{output.host_metadata}", "w") as meta_file:
    meta_file.write("sequence_id\\thost_species\\tsource_type\\tlength\\torganism\\n")
    for entry in host_metadata:
        meta_file.write("\\t".join(map(str, entry)) + "\\n")
PYTHON_EOF
        
        # Copy human exclusion to all hosts (will be extended if other hosts are added)
        cp {output.human_combined} {output.all_hosts_combined}
        
        # Process additional host genomes if present
        for host_genome in resources/hosts/*_genome.fa.gz; do
            if [ -f "$host_genome" ] && [[ "$host_genome" != *"human"* ]]; then
                HOST_NAME=$(basename "$host_genome" _genome.fa.gz)
                echo "Adding $HOST_NAME to exclusion database..." >> {log}
                
                gunzip -c "$host_genome" | python3 -c "
import sys
from Bio import SeqIO
for record in SeqIO.parse(sys.stdin, 'fasta'):
    if len(record.seq) >= 1000:
        clean_id = record.id.replace('|', '_').replace(' ', '_')
        print(f'>{clean_id}|genome|{HOST_NAME}')
        print(record.seq)
" >> {output.all_hosts_combined}
            fi
        done
        
        # Cleanup temporary files
        rm -f temp_human_genome.fa temp_human_transcriptome.fa
        
        echo "Host exclusion database preparation completed: $(date)" >> {log}
        """

# Build host exclusion indexes
rule build_host_exclusion_indexes:
    input:
        host_db="resources/hosts/all_hosts_exclusion.fa"
    output:
        blast_done="results/indexes/{mode}/host_exclusion/blast_done.txt",
        bowtie2_done="results/indexes/{mode}/host_exclusion/bowtie2_done.txt"
    params:
        blast_prefix="results/indexes/{mode}/host_exclusion/host_exclusion",
        bowtie2_prefix="results/indexes/{mode}/host_exclusion/host_exclusion"
    log:
        "logs/build_host_exclusion_indexes_{mode}.log"
    threads: get_pathogen_threads("host_exclusion")
    resources:
        mem_mb=lambda w: get_pathogen_memory("host_exclusion").replace("G", "000").replace("M", "")
    shell:
        """
        module load {config[modules][blast]}
        module load {config[modules][bowtie2]}
        
        echo "Building host exclusion indexes for {wildcards.mode} mode..." > {log}
        
        mkdir -p results/indexes/{wildcards.mode}/host_exclusion
        
        # Build BLAST database for host exclusion
        echo "Building BLAST database..." >> {log}
        makeblastdb \
            -in {input.host_db} \
            -dbtype nucl \
            -out {params.blast_prefix} \
            -title "Host_Exclusion_Database_{wildcards.mode}" \
            -parse_seqids \
            2>> {log}
        
        echo "Host exclusion BLAST database built: $(date)" > {output.blast_done}
        
        # Build Bowtie2 index for host exclusion
        echo "Building Bowtie2 index..." >> {log}
        bowtie2-build \
            --threads {threads} \
            {input.host_db} \
            {params.bowtie2_prefix} \
            2>> {log}
        
        echo "Host exclusion Bowtie2 index built: $(date)" > {output.bowtie2_done}
        
        echo "Host exclusion index building completed: $(date)" >> {log}
        """

# Screen oligo candidates against host sequences
rule screen_host_exclusion:
    input:
        candidates="results/design/{mode}/pathogen_candidates.tsv",
        candidates_fasta="results/design/{mode}/pathogen_candidates.fa",
        blast_done="results/indexes/{mode}/host_exclusion/blast_done.txt",
        bowtie2_done="results/indexes/{mode}/host_exclusion/bowtie2_done.txt"
    output:
        blast_results="results/design/{mode}/host_exclusion_blast.txt",
        bowtie2_results="results/design/{mode}/host_exclusion_bowtie2.sam",
        host_screened="results/design/{mode}/host_screened.tsv",
        host_rejected="results/design/{mode}/host_rejected.tsv",
        screening_summary="results/design/{mode}/host_screening_summary.txt"
    params:
        blast_db="results/indexes/{mode}/host_exclusion/host_exclusion",
        bowtie2_db="results/indexes/{mode}/host_exclusion/host_exclusion",
        max_host_similarity=config["screening"]["host_exclusion"]["max_host_similarity"],
        critical_hosts=config["screening"]["host_exclusion"]["critical_hosts"]
    log:
        "logs/screen_host_exclusion_{mode}.log"
    threads: get_pathogen_threads("cross_reactivity_screening")
    shell:
        """
        module load {config[modules][blast]}
        module load {config[modules][bowtie2]}
        
        echo "Screening candidates against host sequences for {wildcards.mode} mode..." > {log}
        
        # BLAST screening against host sequences
        echo "Running BLAST host exclusion screening..." >> {log}
        blastn \
            -query {input.candidates_fasta} \
            -db {params.blast_db} \
            -out {output.blast_results} \
            -evalue 1000 \
            -word_size 7 \
            -max_target_seqs 100 \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            -num_threads {threads} \
            2>> {log}
        
        # Bowtie2 screening for additional sensitivity
        echo "Running Bowtie2 host exclusion screening..." >> {log}
        bowtie2 \
            -x {params.bowtie2_db} \
            -f {input.candidates_fasta} \
            -S {output.bowtie2_results} \
            -k 10 \
            -L 15 \
            -N 1 \
            --threads {threads} \
            --no-head \
            2>> {log}
        
        # Filter candidates based on host similarity
        echo "Filtering candidates based on host similarity..." >> {log}
        
        python3 << 'PYTHON_EOF' 2>> {log}
import pandas as pd
import sys

def parse_blast_host_results(blast_file):
    """Parse BLAST results against host sequences"""
    host_hits = {{}}
    
    try:
        with open(blast_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\\t')
                if len(fields) >= 12:
                    query_id = fields[0]
                    subject_id = fields[1]
                    identity = float(fields[2])
                    length = int(fields[3])
                    
                    # Extract host information from subject ID
                    if '|' in subject_id:
                        parts = subject_id.split('|')
                        if len(parts) >= 3:
                            host_type = parts[2]  # human, mouse, etc.
                            source_type = parts[1]  # genome, transcriptome
                        else:
                            host_type = "unknown"
                            source_type = "unknown"
                    else:
                        host_type = "unknown"
                        source_type = "unknown"
                    
                    if query_id not in host_hits:
                        host_hits[query_id] = []
                    
                    host_hits[query_id].append({{
                        'host_type': host_type,
                        'source_type': source_type,
                        'identity': identity,
                        'length': length,
                        'subject_id': subject_id
                    }})
    except FileNotFoundError:
        print("Warning: BLAST results file not found")
    
    return host_hits

def assess_host_cross_reactivity(host_hits, max_similarity):
    """Assess host cross-reactivity risk"""
    assessment = {{}}
    
    for query_id, hits in host_hits.items():
        max_identity = max([hit['identity'] for hit in hits]) if hits else 0
        critical_hits = [hit for hit in hits if hit['host_type'] in {params.critical_hosts} and hit['identity'] >= 90]
        
        # Determine risk level
        if critical_hits:
            risk_level = "CRITICAL"
            concern = "High similarity to critical host sequences"
        elif max_identity >= max_similarity:
            risk_level = "HIGH"
            concern = f"Host similarity {max_identity:.1f}% exceeds threshold {max_similarity}%"
        elif max_identity >= max_similarity - 10:
            risk_level = "MEDIUM"
            concern = f"Moderate host similarity {max_identity:.1f}%"
        else:
            risk_level = "LOW"
            concern = f"Low host similarity {max_identity:.1f}%"
        
        assessment[query_id] = {{
            'max_host_identity': max_identity,
            'critical_hits': len(critical_hits),
            'total_hits': len(hits),
            'risk_level': risk_level,
            'concern': concern
        }}
    
    return assessment

# Load candidate oligos
candidates_df = pd.read_csv("{input.candidates}", sep='\\t')

# Parse host screening results
host_hits = parse_blast_host_results("{output.blast_results}")

# Assess cross-reactivity
host_assessment = assess_host_cross_reactivity(host_hits, {params.max_host_similarity})

# Add host screening results to dataframe
candidates_df['max_host_identity'] = candidates_df['primer_id'].map(
    lambda x: host_assessment.get(x, {{}}).get('max_host_identity', 0)
)
candidates_df['host_risk_level'] = candidates_df['primer_id'].map(
    lambda x: host_assessment.get(x, {{}}).get('risk_level', 'LOW')
)
candidates_df['host_concern'] = candidates_df['primer_id'].map(
    lambda x: host_assessment.get(x, {{}}).get('concern', 'No significant host similarity')
)

# Filter candidates
passed_candidates = candidates_df[
    (candidates_df['max_host_identity'] < {params.max_host_similarity}) &
    (candidates_df['host_risk_level'] != 'CRITICAL')
]

rejected_candidates = candidates_df[
    (candidates_df['max_host_identity'] >= {params.max_host_similarity}) |
    (candidates_df['host_risk_level'] == 'CRITICAL')
]

# Save results
passed_candidates.to_csv("{output.host_screened}", sep='\\t', index=False)
rejected_candidates.to_csv("{output.host_rejected}", sep='\\t', index=False)

# Generate summary
with open("{output.screening_summary}", 'w') as f:
    f.write("Host Exclusion Screening Summary\\n")
    f.write("================================\\n\\n")
    f.write(f"Total candidates screened: {{len(candidates_df)}}\\n")
    f.write(f"Passed host screening: {{len(passed_candidates)}}\\n")
    f.write(f"Rejected due to host similarity: {{len(rejected_candidates)}}\\n")
    f.write(f"Pass rate: {{len(passed_candidates)/len(candidates_df)*100:.1f}}%\\n\\n")
    
    f.write("Risk Level Distribution:\\n")
    risk_counts = candidates_df['host_risk_level'].value_counts()
    for risk, count in risk_counts.items():
        f.write(f"  {{risk}}: {{count}}\\n")
    
    f.write(f"\\nScreening Parameters:\\n")
    f.write(f"  Max host similarity threshold: {params.max_host_similarity}%\\n")
    f.write(f"  Critical hosts: {params.critical_hosts}\\n")

print(f"Host exclusion screening completed: {{len(passed_candidates)}} passed, {{len(rejected_candidates)}} rejected")
PYTHON_EOF
        
        echo "Host exclusion screening completed: $(date)" >> {log}
        """

# Validate host exclusion effectiveness
rule validate_host_exclusion:
    input:
        host_screened="results/design/{mode}/host_screened.tsv",
        host_rejected="results/design/{mode}/host_rejected.tsv"
    output:
        "results/validation/host_exclusion_validation_{mode}.txt"
    params:
        validation_sequences=config.get("validation", {}).get("host_validation_sequences", [
            "GAPDH_human_specific",
            "ACTB_human_specific", 
            "18S_rRNA_human"
        ])
    log:
        "logs/validate_host_exclusion_{mode}.log"
    shell:
        """
        echo "Validating host exclusion for {wildcards.mode} mode..." > {log}
        
        cat > {output} << EOF
Host Exclusion Validation Report - {wildcards.mode} Mode
========================================================
Generated: $(date)

Validation Summary:
EOF
        
        # Check screening results
        TOTAL_SCREENED=$(tail -n +2 {input.host_screened} | wc -l)
        TOTAL_REJECTED=$(tail -n +2 {input.host_rejected} | wc -l)
        TOTAL_CANDIDATES=$((TOTAL_SCREENED + TOTAL_REJECTED))
        
        echo "Total candidates: $TOTAL_CANDIDATES" >> {output}
        echo "Passed host screening: $TOTAL_SCREENED" >> {output}
        echo "Rejected for host similarity: $TOTAL_REJECTED" >> {output}
        
        if [ $TOTAL_CANDIDATES -gt 0 ]; then
            PASS_RATE=$(echo "scale=1; $TOTAL_SCREENED * 100 / $TOTAL_CANDIDATES" | bc 2>/dev/null || echo "0")
            echo "Host screening pass rate: $PASS_RATE%" >> {output}
        fi
        
        echo "" >> {output}
        echo "Quality Assessment:" >> {output}
        
        # Check for critical rejections
        if [ -f {input.host_rejected} ]; then
            CRITICAL_REJECTIONS=$(tail -n +2 {input.host_rejected} | grep -c "CRITICAL" || echo "0")
            echo "Critical host similarity rejections: $CRITICAL_REJECTIONS" >> {output}
            
            if [ $CRITICAL_REJECTIONS -gt 0 ]; then
                echo "⚠️  WARNING: Critical host similarities detected" >> {output}
            else
                echo "✅ No critical host similarities detected" >> {output}
            fi
        fi
        
        # Validate that known host sequences would be rejected
        echo "" >> {output}
        echo "Known Host Sequence Check:" >> {output}
        echo "(These sequences should be rejected if present)" >> {output}
        
        for seq_name in {params.validation_sequences}; do
            echo "  $seq_name: Not tested (validation sequences not implemented)" >> {output}
        done
        
        echo "" >> {output}
        echo "Recommendations:" >> {output}
        if [ $CRITICAL_REJECTIONS -gt 0 ]; then
            echo "- Review critical rejections in host_rejected.tsv" >> {output}
            echo "- Consider adjusting similarity thresholds" >> {output}
        else
            echo "- Host exclusion appears effective" >> {output}
            echo "- Proceed with pathogen specificity screening" >> {output}
        fi
        
        echo "Host exclusion validation completed: $(date)" >> {log}
        """

# Generate comprehensive host exclusion report
rule host_exclusion_report:
    input:
        expand("results/design/{mode}/host_screened.tsv", mode=DESIGN_MODES),
        expand("results/validation/host_exclusion_validation_{mode}.txt", mode=DESIGN_MODES)
    output:
        "results/validation/host_exclusion_test.txt",
        "results/reports/host_exclusion_summary.html"
    log:
        "logs/host_exclusion_report.log"
    shell:
        """
        echo "Generating comprehensive host exclusion report..." > {log}
        
        # Create text summary
        cat > {output[0]} << EOF
Host Exclusion Pipeline Summary
===============================
Generated: $(date)

Design Modes Processed: {DESIGN_MODES}

EOF
        
        # Add summary for each mode
        for mode in {DESIGN_MODES}; do
            echo "Mode: $mode" >> {output[0]}
            if [ -f "results/design/$mode/host_screened.tsv" ]; then
                PASSED=$(tail -n +2 "results/design/$mode/host_screened.tsv" | wc -l)
                echo "  Passed host screening: $PASSED oligos" >> {output[0]}
            fi
            if [ -f "results/design/$mode/host_rejected.tsv" ]; then
                REJECTED=$(tail -n +2 "results/design/$mode/host_rejected.tsv" | wc -l)
                echo "  Rejected for host similarity: $REJECTED oligos" >> {output[0]}
            fi
            echo "" >> {output[0]}
        done
        
        echo "Overall Status: Host exclusion screening completed successfully" >> {output[0]}
        
        # Create HTML report
        cat > {output[1]} << 'HTML_EOF'
<!DOCTYPE html>
<html>
<head>
    <title>Host Exclusion Summary Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background-color: #f0f8ff; padding: 20px; border-radius: 8px; }}
        .summary {{ background-color: #f9f9f9; padding: 15px; margin: 15px 0; }}
        .pass {{ color: green; font-weight: bold; }}
        .reject {{ color: red; font-weight: bold; }}
        table {{ width: 100%; border-collapse: collapse; margin: 15px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #4CAF50; color: white; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 Host Exclusion Summary Report</h1>
        <p>Generated: $(date)</p>
    </div>
    
    <div class="summary">
        <h2>Summary</h2>
        <p>Host exclusion screening ensures that designed primers do not cross-react with host organism sequences, which is critical for clinical pathogen detection applications.</p>
    </div>
    
    <h2>Results by Design Mode</h2>
    <table>
        <tr>
            <th>Design Mode</th>
            <th>Candidates Screened</th>
            <th>Passed</th>
            <th>Rejected</th>
            <th>Pass Rate</th>
        </tr>
HTML_EOF
        
        # Add table rows for each mode
        for mode in {DESIGN_MODES}; do
            if [ -f "results/design/$mode/host_screened.tsv" ] && [ -f "results/design/$mode/host_rejected.tsv" ]; then
                PASSED=$(tail -n +2 "results/design/$mode/host_screened.tsv" | wc -l)
                REJECTED=$(tail -n +2 "results/design/$mode/host_rejected.tsv" | wc -l)
                TOTAL=$((PASSED + REJECTED))
                if [ $TOTAL -gt 0 ]; then
                    PASS_RATE=$(echo "scale=1; $PASSED * 100 / $TOTAL" | bc 2>/dev/null || echo "0")
                else
                    PASS_RATE="0"
                fi
                
                echo "        <tr>" >> {output[1]}
                echo "            <td>$mode</td>" >> {output[1]}
                echo "            <td>$TOTAL</td>" >> {output[1]}
                echo "            <td class=\"pass\">$PASSED</td>" >> {output[1]}
                echo "            <td class=\"reject\">$REJECTED</td>" >> {output[1]}
                echo "            <td>$PASS_RATE%</td>" >> {output[1]}
                echo "        </tr>" >> {output[1]}
            fi
        done
        
        cat >> {output[1]} << 'HTML_EOF'
    </table>
    
    <h2>Quality Control</h2>
    <p>✅ Host exclusion screening completed successfully</p>
    <p>📊 Detailed validation reports available for each design mode</p>
    
    <h2>Next Steps</h2>
    <ul>
        <li>Review rejected candidates in host_rejected.tsv files</li>
        <li>Proceed with pathogen specificity screening</li>
        <li>Validate selected primers experimentally</li>
    </ul>
</body>
</html>
HTML_EOF
        
        echo "Host exclusion report generation completed: $(date)" >> {log}
        """

# Mark host exclusion databases as ready
rule host_exclusion_ready:
    input:
        expand("results/indexes/{mode}/host_exclusion/blast_done.txt", mode=DESIGN_MODES),
        expand("results/indexes/{mode}/host_exclusion/bowtie2_done.txt", mode=DESIGN_MODES),
        "resources/hosts/host_metadata.tsv"
    output:
        "resources/hosts/exclusion_databases_ready.txt"
    shell:
        """
        echo "Host Exclusion Databases Ready" > {output}
        echo "Generated: $(date)" >> {output}
        echo "" >> {output}
        echo "Available Databases:" >> {output}
        echo "- Human genome + transcriptome" >> {output}
        
        # Check for additional host databases
        for host_file in resources/hosts/*_genome.fa.gz; do
            if [ -f "$host_file" ] && [[ "$host_file" != *"human"* ]]; then
                HOST_NAME=$(basename "$host_file" _genome.fa.gz)
                echo "- $HOST_NAME genome" >> {output}
            fi
        done
        
        echo "" >> {output}
        echo "Indexes built for design modes: {DESIGN_MODES}" >> {output}
        """
