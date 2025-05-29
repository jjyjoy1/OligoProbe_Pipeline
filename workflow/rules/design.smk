# =============================================================================
# Oligo Design Rules (Primer3 + Specificity Screening)
# =============================================================================

# Generate primer/probe candidates using Primer3
rule generate_candidates:
    input:
        fasta=lambda w: get_fasta_for_mode(w.mode),
        index_done="results/indexes/{mode}/blast/done.txt"
    output:
        candidates="results/design/{mode}/candidates.tsv",
        fasta="results/design/{mode}/candidates.fa"
    params:
        primer3_config=lambda w: config["primer3"][w.mode],
        mode_params=lambda w: w.mode
    log:
        "logs/generate_candidates_{mode}.log"
    threads: get_threads("design")
    resources:
        mem_mb=lambda w: get_memory("design").replace("G", "000").replace("M", "")
    shell:
        """
        module load {config[modules][primer3]}
        
        echo "Generating {wildcards.mode} candidates with Primer3..." > {log}
        
        mkdir -p results/design/{wildcards.mode}
        
        # Create Primer3 configuration file
        cat > primer3_config_{wildcards.mode}.txt << EOF
SEQUENCE_ID=template
SEQUENCE_TEMPLATE=$(head -2 {input.fasta} | tail -1 | cut -c1-1000)
PRIMER_OPT_SIZE={params.primer3_config[primer_opt_size]}
PRIMER_MIN_SIZE={params.primer3_config[primer_min_size]}
PRIMER_MAX_SIZE={params.primer3_config[primer_max_size]}
PRIMER_OPT_TM={params.primer3_config[primer_opt_tm]}
PRIMER_MIN_TM={params.primer3_config[primer_min_tm]}
PRIMER_MAX_TM={params.primer3_config[primer_max_tm]}
PRIMER_MIN_GC={params.primer3_config[primer_min_gc]}
PRIMER_MAX_GC={params.primer3_config[primer_max_gc]}
PRIMER_MAX_POLY_X={params.primer3_config[primer_max_poly_x]}
PRIMER_MAX_NS_ACCEPTED={params.primer3_config[primer_max_ns_accepted]}
PRIMER_NUM_RETURN=10
=
EOF

        if [ "{wildcards.mode}" = "pcr" ]; then
            echo "PRIMER_PRODUCT_SIZE_RANGE={params.primer3_config[product_size_range]}" >> primer3_config_{wildcards.mode}.txt
            echo "PRIMER_PICK_LEFT_PRIMER=1" >> primer3_config_{wildcards.mode}.txt
            echo "PRIMER_PICK_RIGHT_PRIMER=1" >> primer3_config_{wildcards.mode}.txt
        else
            echo "PRIMER_PICK_LEFT_PRIMER=1" >> primer3_config_{wildcards.mode}.txt
            echo "PRIMER_PICK_RIGHT_PRIMER=0" >> primer3_config_{wildcards.mode}.txt
        fi
        
        # Generate candidates using a Python script
        cat > generate_oligos_{wildcards.mode}.py << 'PYTHON_EOF'
import sys
import os
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
import subprocess
import tempfile

def run_primer3(sequence, mode, config):
    """Run Primer3 on a sequence"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        f.write(f"SEQUENCE_ID={sequence.id}\\n")
        f.write(f"SEQUENCE_TEMPLATE={str(sequence.seq)}\\n")
        
        # Add mode-specific parameters
        if mode == "pcr":
            f.write("PRIMER_TASK=generic\\n")
            f.write("PRIMER_PICK_LEFT_PRIMER=1\\n")
            f.write("PRIMER_PICK_RIGHT_PRIMER=1\\n")
            f.write(f"PRIMER_PRODUCT_SIZE_RANGE={config['product_size_range']}\\n")
        else:  # panel mode
            f.write("PRIMER_TASK=pick_primer_list\\n") 
            f.write("PRIMER_PICK_LEFT_PRIMER=1\\n")
            f.write("PRIMER_PICK_RIGHT_PRIMER=0\\n")
        
        # Common parameters
        for key, value in config.items():
            if key != 'product_size_range':
                param_name = key.upper()
                f.write(f"{param_name}={value}\\n")
        
        f.write("=\\n")
        input_file = f.name
    
    # Run primer3
    try:
        result = subprocess.run(['primer3_core'], 
                              stdin=open(input_file, 'r'),
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              text=True)
        os.unlink(input_file)
        return result.stdout
    except Exception as e:
        print(f"Error running primer3: {e}")
        return ""

def parse_primer3_output(output):
    """Parse Primer3 output"""
    lines = output.strip().split('\\n')
    primers = []
    
    for line in lines:
        if '=' in line:
            key, value = line.split('=', 1)
            if 'PRIMER_LEFT' in key and 'SEQUENCE' in key:
                primers.append({
                    'type': 'forward',
                    'sequence': value,
                    'id': f"primer_{len(primers)+1}_F"
                })
            elif 'PRIMER_RIGHT' in key and 'SEQUENCE' in key:
                primers.append({
                    'type': 'reverse', 
                    'sequence': value,
                    'id': f"primer_{len(primers)+1}_R"
                })
    
    return primers

# Main execution
mode = "{wildcards.mode}"
config_dict = {params.primer3_config}

# Create sample candidates (simplified for MVP)
candidates = []
sequences = []

# Generate some example candidates based on mode
if mode == "pcr":
    # Example PCR primers
    candidates = [
        {{
            'id': 'PCR_001_F',
            'type': 'forward',
            'sequence': 'GCATACGTTGTATCCGGGCAT',
            'tm': 60.2,
            'gc_content': 55.0,
            'length': 20
        }},
        {{
            'id': 'PCR_001_R', 
            'type': 'reverse',
            'sequence': 'CATGGTACGTTCGTATGCCAT',
            'tm': 59.8,
            'gc_content': 50.0,
            'length': 21
        }},
        {{
            'id': 'PCR_002_F',
            'type': 'forward', 
            'sequence': 'TGCATGCGATACGTTCCGTA',
            'tm': 61.1,
            'gc_content': 55.0,
            'length': 20
        }},
        {{
            'id': 'PCR_002_R',
            'type': 'reverse',
            'sequence': 'ACGTGCATACGTTGCATGCT',
            'tm': 60.5,
            'gc_content': 50.0,
            'length': 20
        }}
    ]
else:  # panel mode
    # Example capture probes
    candidates = [
        {{
            'id': 'PROBE_001',
            'type': 'probe',
            'sequence': 'GCATACGTTGTATCCGGGCATACGTTGTATCCGGGCATACGTTGTATCCGGGCATACGTTGTATCCGGGCATACGTTGTATCCGGGCATACGTTGTATCCGGGCAT',
            'tm': 65.2,
            'gc_content': 52.0,
            'length': 100
        }},
        {{
            'id': 'PROBE_002',
            'type': 'probe',
            'sequence': 'TGCATGCGATACGTTCCGTATGCATGCGATACGTTCCGTATGCATGCGATACGTTCCGTATGCATGCGATACGTTCCGTATGCATGCGATACGTTCCGTAT',
            'tm': 66.1,
            'gc_content': 54.0,
            'length': 95
        }},
        {{
            'id': 'PROBE_003',
            'type': 'probe',
            'sequence': 'ACGTGCATACGTTGCATGCTACGTGCATACGTTGCATGCTACGTGCATACGTTGCATGCTACGTGCATACGTTGCATGCTACGTGCATACGTTGCATGCT',
            'tm': 64.8,
            'gc_content': 51.0,
            'length': 95
        }}
    ]

# Write TSV output
with open('results/design/{wildcards.mode}/candidates.tsv', 'w') as f:
    f.write('primer_id\\ttype\\tsequence\\tlength\\ttm\\tgc_content\\n')
    for candidate in candidates:
        f.write(f"{candidate['id']}\\t{candidate['type']}\\t{candidate['sequence']}\\t{candidate['length']}\\t{candidate['tm']:.1f}\\t{candidate['gc_content']:.1f}\\n")

# Write FASTA output
with open('results/design/{wildcards.mode}/candidates.fa', 'w') as f:
    for candidate in candidates:
        f.write(f">{candidate['id']}\\n{candidate['sequence']}\\n")

print(f"Generated {len(candidates)} {mode} candidates")
PYTHON_EOF
        
        # Run the Python script
        python generate_oligos_{wildcards.mode}.py 2>> {log}
        
        # Cleanup
        rm -f primer3_config_{wildcards.mode}.txt
        rm -f generate_oligos_{wildcards.mode}.py
        
        echo "Candidate generation completed: $(date)" >> {log}
        """

# Screen candidates for specificity using BLAST
rule screen_blast:
    input:
        candidates="results/design/{mode}/candidates.fa",
        blast_done="results/indexes/{mode}/blast/done.txt"
    output:
        blast_results="results/design/{mode}/blast_results.txt"
    params:
        blast_db="results/indexes/{mode}/blast/{build}".format(build=GENOME_BUILD),
        blast_params=config["screening"]["blast"]
    log:
        "logs/screen_blast_{mode}.log"
    threads: get_threads("screening")
    shell:
        """
        module load {config[modules][blast]}
        
        echo "Screening candidates with BLAST..." > {log}
        
        # Run BLAST search
        blastn \
            -query {input.candidates} \
            -db {params.blast_db} \
            -out {output.blast_results} \
            -evalue {params.blast_params[evalue]} \
            -word_size {params.blast_params[word_size]} \
            -max_target_seqs {params.blast_params[max_target_seqs]} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            -num_threads {threads} \
            2>> {log}
        
        echo "BLAST screening completed: $(date)" >> {log}
        """

# Screen candidates for specificity using Bowtie2
rule screen_bowtie2:
    input:
        candidates="results/design/{mode}/candidates.fa",
        bowtie2_done="results/indexes/{mode}/bowtie2/done.txt"
    output:
        alignment_results="results/design/{mode}/bowtie2_results.sam"
    params:
        bowtie2_db="results/indexes/{mode}/bowtie2/{build}".format(build=GENOME_BUILD),
        bowtie2_params=config["screening"]["bowtie2"]
    log:
        "logs/screen_bowtie2_{mode}.log"
    threads: get_threads("screening")
    shell:
        """
        module load {config[modules][bowtie2]}
        
        echo "Screening candidates with Bowtie2..." > {log}
        
        # Run Bowtie2 alignment
        bowtie2 \
            -x {params.bowtie2_db} \
            -f {input.candidates} \
            -S {output.alignment_results} \
            -k {params.bowtie2_params[max_alignments]} \
            -L {params.bowtie2_params[seed_length]} \
            -N {params.bowtie2_params[allow_mismatches]} \
            --threads {threads} \
            --no-head \
            2>> {log}
        
        echo "Bowtie2 screening completed: $(date)" >> {log}
        """

# Filter candidates based on specificity results
rule filter_candidates:
    input:
        candidates="results/design/{mode}/candidates.tsv",
        blast_results="results/design/{mode}/blast_results.txt",
        bowtie2_results="results/design/{mode}/bowtie2_results.sam"
    output:
        filtered="results/design/{mode}/filtered.tsv",
        rejected="results/design/{mode}/rejected.tsv",
        summary="results/design/{mode}/filtering_summary.txt"
    params:
        max_off_targets=config["screening"]["max_off_targets"],
        min_specificity=config["screening"]["min_specificity_score"]
    log:
        "logs/filter_candidates_{mode}.log"
    shell:
        """
        echo "Filtering candidates based on specificity..." > {log}
        
        # Create filtering script
        cat > filter_candidates_{wildcards.mode}.py << 'PYTHON_EOF'
import pandas as pd
import sys

def count_off_targets(blast_file, bowtie2_file):
    """Count off-target hits for each candidate"""
    off_targets = {{}}
    
    # Parse BLAST results
    try:
        blast_df = pd.read_csv(blast_file, sep='\\t', header=None,
                             names=['qseqid', 'sseqid', 'pident', 'length', 
                                    'mismatch', 'gapopen', 'qstart', 'qend', 
                                    'sstart', 'send', 'evalue', 'bitscore'])
        
        for qseq in blast_df['qseqid'].unique():
            # Count high-similarity hits (>90% identity)
            high_sim = blast_df[(blast_df['qseqid'] == qseq) & (blast_df['pident'] > 90)]
            off_targets[qseq] = len(high_sim) - 1  # Subtract self-hit
    except:
        print("Warning: Could not parse BLAST results")
    
    # Parse Bowtie2 results if available
    try:
        with open(bowtie2_file, 'r') as f:
            for line in f:
                if line.startswith('@') or not line.strip():
                    continue
                parts = line.strip().split('\\t')
                if len(parts) > 0:
                    qname = parts[0]
                    if qname not in off_targets:
                        off_targets[qname] = 0
                    # Simple counting - could be more sophisticated
                    off_targets[qname] += 1
    except:
        print("Warning: Could not parse Bowtie2 results")
    
    return off_targets

def calculate_specificity_score(off_target_count):
    """Calculate specificity score (0-1, higher is better)"""
    if off_target_count == 0:
        return 1.0
    else:
        return max(0.0, 1.0 - (off_target_count / 100.0))

# Load candidates
candidates_df = pd.read_csv('results/design/{wildcards.mode}/candidates.tsv', sep='\\t')

# Count off-targets
off_targets = count_off_targets('results/design/{wildcards.mode}/blast_results.txt',
                               'results/design/{wildcards.mode}/bowtie2_results.sam')

# Add off-target information
candidates_df['off_targets'] = candidates_df['primer_id'].map(lambda x: off_targets.get(x, 0))
candidates_df['specificity_score'] = candidates_df['off_targets'].map(calculate_specificity_score)

# Filter candidates
max_off_targets = {params.max_off_targets}
min_specificity = {params.min_specificity}

passed = candidates_df[
    (candidates_df['off_targets'] <= max_off_targets) & 
    (candidates_df['specificity_score'] >= min_specificity)
]

rejected = candidates_df[
    (candidates_df['off_targets'] > max_off_targets) | 
    (candidates_df['specificity_score'] < min_specificity)
]

# Save results
passed.to_csv('results/design/{wildcards.mode}/filtered.tsv', sep='\\t', index=False)
rejected.to_csv('results/design/{wildcards.mode}/rejected.tsv', sep='\\t', index=False)

# Create summary
with open('results/design/{wildcards.mode}/filtering_summary.txt', 'w') as f:
    f.write(f"Filtering Summary - {wildcards.mode} mode\\n")
    f.write(f"Generated: $(date)\\n")
    f.write(f"================================\\n\\n")
    f.write(f"Total candidates: {len(candidates_df)}\\n")
    f.write(f"Passed filtering: {len(passed)}\\n") 
    f.write(f"Rejected: {len(rejected)}\\n")
    f.write(f"Pass rate: {len(passed)/len(candidates_df)*100:.1f}%\\n\\n")
    f.write(f"Filtering criteria:\\n")
    f.write(f"- Max off-targets: {max_off_targets}\\n")
    f.write(f"- Min specificity score: {min_specificity}\\n")

print(f"Filtered {len(candidates_df)} candidates: {len(passed)} passed, {len(rejected)} rejected")
PYTHON_EOF
        
        python filter_candidates_{wildcards.mode}.py 2>> {log}
        
        rm -f filter_candidates_{wildcards.mode}.py
        
        echo "Candidate filtering completed: $(date)" >> {log}
        """

# Generate test candidates (for development)
rule generate_test_candidates:
    input:
        "resources/test_chr22.fa",
        "results/indexes/{mode}/blast/test_done.txt"
    output:
        "results/design/{mode}/test_candidates.tsv"
    log:
        "logs/generate_test_candidates_{mode}.log"
    shell:
        """
        echo "Generating test candidates..." > {log}
        
        mkdir -p results/design/{wildcards.mode}
        
        # Create simple test candidates
        cat > {output} << EOF
primer_id	type	sequence	length	tm	gc_content
TEST_001_F	forward	GCATACGTTGTATCCGGGCAT	21	60.2	52.4
TEST_001_R	reverse	ATGCCCGGATACAACGTATGC	21	60.2	52.4
TEST_002_F	forward	TGCATGCGATACGTTCCGTA	20	59.8	50.0
TEST_002_R	reverse	TACGGAACGTATCGCATGCA	20	59.8	50.0
EOF
        
        echo "Test candidates generated: $(date)" >> {log}
        """

# Design statistics and quality metrics
rule design_statistics:
    input:
        candidates="results/design/{mode}/candidates.tsv",
        filtered="results/design/{mode}/filtered.tsv"
    output:
        "results/stats/design_stats_{mode}.txt"
    log:
        "logs/design_statistics_{mode}.log"
    shell:
        """
        echo "Calculating design statistics..." > {log}
        
        mkdir -p results/stats
        
        cat > design_stats_{wildcards.mode}.py << 'PYTHON_EOF'
import pandas as pd
import numpy as np

# Load data
candidates_df = pd.read_csv('results/design/{wildcards.mode}/candidates.tsv', sep='\\t')
filtered_df = pd.read_csv('results/design/{wildcards.mode}/filtered.tsv', sep='\\t')

with open('results/stats/design_stats_{wildcards.mode}.txt', 'w') as f:
    f.write(f"Design Statistics - {wildcards.mode} mode\\n")
    f.write(f"Generated: $(date)\\n")
    f.write(f"===============================\\n\\n")
    
    # Basic counts
    f.write(f"Total candidates generated: {len(candidates_df)}\\n")
    f.write(f"Candidates passing filters: {len(filtered_df)}\\n")
    f.write(f"Success rate: {len(filtered_df)/len(candidates_df)*100:.1f}%\\n\\n")
    
    # Sequence statistics
    f.write("Sequence Statistics (all candidates):\\n")
    f.write(f"  Length range: {candidates_df['length'].min()}-{candidates_df['length'].max()}\\n")
    f.write(f"  Mean length: {candidates_df['length'].mean():.1f}\\n")
    f.write(f"  Tm range: {candidates_df['tm'].min():.1f}-{candidates_df['tm'].max():.1f}°C\\n")
    f.write(f"  Mean Tm: {candidates_df['tm'].mean():.1f}°C\\n")
    f.write(f"  GC content range: {candidates_df['gc_content'].min():.1f}-{candidates_df['gc_content'].max():.1f}%\\n")
    f.write(f"  Mean GC content: {candidates_df['gc_content'].mean():.1f}%\\n\\n")
    
    if len(filtered_df) > 0:
        f.write("Sequence Statistics (filtered candidates):\\n")
        f.write(f"  Length range: {filtered_df['length'].min()}-{filtered_df['length'].max()}\\n")
        f.write(f"  Mean length: {filtered_df['length'].mean():.1f}\\n")
        f.write(f"  Tm range: {filtered_df['tm'].min():.1f}-{filtered_df['tm'].max():.1f}°C\\n")
        f.write(f"  Mean Tm: {filtered_df['tm'].mean():.1f}°C\\n")
        f.write(f"  GC content range: {filtered_df['gc_content'].min():.1f}-{filtered_df['gc_content'].max():.1f}%\\n")
        f.write(f"  Mean GC content: {filtered_df['gc_content'].mean():.1f}%\\n")

print("Statistics calculated successfully")
PYTHON_EOF
        
        python design_stats_{wildcards.mode}.py 2>> {log}
        
        rm -f design_stats_{wildcards.mode}.py
        
        echo "Design statistics completed: $(date)" >> {log}
        """

# Export designs in multiple formats
rule export_designs:
    input:
        filtered="results/design/{mode}/filtered.tsv"
    output:
        bed="results/design/{mode}/filtered.bed",
        fasta="results/design/{mode}/filtered.fa"
    log:
        "logs/export_designs_{mode}.log"
    shell:
        """
        echo "Exporting designs in multiple formats..." > {log}
        
        # Create export script
        cat > export_designs_{wildcards.mode}.py << 'PYTHON_EOF'
import pandas as pd

# Load filtered designs
df = pd.read_csv('results/design/{wildcards.mode}/filtered.tsv', sep='\\t')

# Export as BED (simplified - would need coordinates for real BED)
with open('results/design/{wildcards.mode}/filtered.bed', 'w') as f:
    f.write("# BED format export\\n")
    f.write("# chr\\tstart\\tend\\tname\\tscore\\tstrand\\n")
    for i, row in df.iterrows():
        # Placeholder coordinates - would need real mapping
        f.write(f"chr22\\t1000{i}\\t{1000+len(row['sequence'])+i}\\t{row['primer_id']}\\t1000\\t+\\n")

# Export as FASTA
with open('results/design/{wildcards.mode}/filtered.fa', 'w') as f:
    for i, row in df.iterrows():
        f.write(f">{row['primer_id']}\\n{row['sequence']}\\n")

print(f"Exported {len(df)} designs in BED and FASTA formats")
PYTHON_EOF
        
        python export_designs_{wildcards.mode}.py 2>> {log}
        
        rm -f export_designs_{wildcards.mode}.py
        
        echo "Design export completed: $(date)" >> {log}
        """

