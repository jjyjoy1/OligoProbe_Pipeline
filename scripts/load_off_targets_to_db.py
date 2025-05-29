#!/usr/bin/env python3
"""
Load off-target analysis results into PostgreSQL database
"""

import pandas as pd
import psycopg2
import sys
import os
import re
from datetime import datetime

def connect_to_database(db_config):
    """Connect to PostgreSQL database"""
    return psycopg2.connect(
        host=db_config.get('host', 'localhost'),
        port=db_config.get('port', 5432),
        database=db_config.get('database', 'oligodb'),
        user=db_config.get('user', 'postgres'),
        password=db_config.get('password', '')
    )

def parse_blast_results(blast_file):
    """Parse BLAST results file"""
    if not os.path.exists(blast_file):
        print(f"BLAST file not found: {blast_file}")
        return []
    
    blast_hits = []
    
    with open(blast_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # BLAST tabular format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
            fields = line.split('\t')
            if len(fields) >= 12:
                hit = {
                    'query_id': fields[0],
                    'subject_id': fields[1],
                    'identity_percent': float(fields[2]),
                    'alignment_length': int(fields[3]),
                    'mismatch_count': int(fields[4]),
                    'gap_count': int(fields[5]),
                    'query_start': int(fields[6]),
                    'query_end': int(fields[7]),
                    'subject_start': int(fields[8]),
                    'subject_end': int(fields[9]),
                    'e_value': float(fields[10]),
                    'bit_score': float(fields[11])
                }
                blast_hits.append(hit)
    
    return blast_hits

def parse_chromosome_location(subject_id):
    """Parse chromosome and position from subject ID"""
    # Try to extract chromosome and position from subject_id
    # Format examples: "chr1:1000-2000", "chr22_123456_123556"
    
    chromosome = None
    start_pos = None
    end_pos = None
    
    # Pattern 1: chr1:1000-2000
    match = re.match(r'chr(\w+):(\d+)-(\d+)', subject_id)
    if match:
        chromosome = f"chr{match.group(1)}"
        start_pos = int(match.group(2))
        end_pos = int(match.group(3))
        return chromosome, start_pos, end_pos
    
    # Pattern 2: chr22_123456_123556
    match = re.match(r'chr(\w+)_(\d+)_(\d+)', subject_id)
    if match:
        chromosome = f"chr{match.group(1)}"
        start_pos = int(match.group(2))
        end_pos = int(match.group(3))
        return chromosome, start_pos, end_pos
    
    # Pattern 3: Just chromosome
    match = re.match(r'chr(\w+)', subject_id)
    if match:
        chromosome = f"chr{match.group(1)}"
        return chromosome, start_pos, end_pos
    
    return chromosome, start_pos, end_pos

def get_oligo_design_id(cursor, query_id):
    """Get oligo_design_id from primer name"""
    cursor.execute("""
        SELECT od.oligo_design_id 
        FROM oligo_design od
        JOIN bioentry be ON od.bioentry_id = be.bioentry_id
        WHERE be.name = %s
    """, (query_id,))
    
    result = cursor.fetchone()
    return result[0] if result else None

def insert_off_target_hit(cursor, oligo_design_id, hit):
    """Insert off-target hit into database"""
    
    # Parse location information
    chromosome, start_pos, end_pos = parse_chromosome_location(hit['subject_id'])
    
    # Determine strand (simplified)
    strand = '+' if hit['subject_start'] < hit['subject_end'] else '-'
    
    # Calculate actual positions on chromosome
    if start_pos and end_pos:
        actual_start = start_pos + min(hit['subject_start'], hit['subject_end']) - 1
        actual_end = start_pos + max(hit['subject_start'], hit['subject_end']) - 1
    else:
        actual_start = hit['subject_start']
        actual_end = hit['subject_end']
    
    cursor.execute("""
        INSERT INTO off_target_hit (
            oligo_design_id, chromosome, start_pos, end_pos, strand,
            identity_percent, mismatch_count, gap_count, alignment_length,
            e_value, bit_score, hit_sequence
        ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
    """, (
        oligo_design_id, chromosome, actual_start, actual_end, strand,
        hit['identity_percent'], hit['mismatch_count'], hit['gap_count'],
        hit['alignment_length'], hit['e_value'], hit['bit_score'],
        hit['subject_id']  # Store subject_id as hit_sequence for now
    ))

def load_off_targets(blast_files, db_config):
    """Load off-target hits from BLAST results"""
    
    conn = connect_to_database(db_config)
    cursor = conn.cursor()
    
    total_hits = 0
    loaded_hits = 0
    
    try:
        for blast_file in blast_files:
            if not os.path.exists(blast_file):
                print(f"Skipping missing file: {blast_file}")
                continue
                
            print(f"Processing {blast_file}")
            blast_hits = parse_blast_results(blast_file)
            total_hits += len(blast_hits)
            
            for hit in blast_hits:
                # Skip self-hits (100% identity over full length)
                if hit['identity_percent'] == 100.0 and hit['mismatch_count'] == 0:
                    continue
                
                oligo_design_id = get_oligo_design_id(cursor, hit['query_id'])
                if oligo_design_id:
                    try:
                        insert_off_target_hit(cursor, oligo_design_id, hit)
                        loaded_hits += 1
                        
                        if loaded_hits % 1000 == 0:
                            print(f"Loaded {loaded_hits} off-target hits...")
                            conn.commit()
                            
                    except Exception as e:
                        print(f"Error inserting hit for {hit['query_id']}: {e}")
                        continue
                else:
                    print(f"Oligo not found for query: {hit['query_id']}")
        
        # Update off-target counts in oligo_design table
        print("Updating off-target counts...")
        cursor.execute("""
            UPDATE oligo_design 
            SET off_target_count = (
                SELECT COUNT(*) 
                FROM off_target_hit 
                WHERE off_target_hit.oligo_design_id = oligo_design.oligo_design_id
            )
        """)
        
        conn.commit()
        print(f"Successfully loaded {loaded_hits} off-target hits from {total_hits} total hits")
        
        return loaded_hits, total_hits
        
    except Exception as e:
        conn.rollback()
        raise e
    finally:
        cursor.close()
        conn.close()

def main():
    """Main function for Snakemake script"""
    
    # Get parameters from snakemake object
    blast_files = [snakemake.input.pcr_blast, snakemake.input.panel_blast]
    db_config = snakemake.params.db_config
    output_file = snakemake.output[0]
    log_file = snakemake.log[0]
    
    try:
        # Redirect output to log file
        with open(log_file, 'w') as log:
            sys.stdout = log
            sys.stderr = log
            
            print("Loading off-target analysis results")
            print(f"BLAST files: {blast_files}")
            print(f"Database: {db_config.get('database', 'oligodb')}")
            print("-" * 50)
            
            loaded_hits, total_hits = load_off_targets(blast_files, db_config)
            
            # Write completion marker
            with open(output_file, 'w') as f:
                f.write(f"Loaded {loaded_hits} off-target hits from {total_hits} total\n")
                f.write(f"Completed at {datetime.now()}\n")
            
            print("Off-target loading completed successfully")
        
    except Exception as e:
        with open(log_file, 'w') as log:
            print(f"Error loading off-targets: {e}", file=log)
        
        # Write error marker
        with open(output_file, 'w') as f:
            f.write(f"ERROR: Failed to load off-targets - {e}\n")
        
        raise

if __name__ == "__main__":
    main()
