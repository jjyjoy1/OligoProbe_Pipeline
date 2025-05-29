#!/usr/bin/env python3
"""
Export oligo data from PostgreSQL database in various formats
"""

import psycopg2
import json
import csv
import sys
import os
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

def export_fasta(cursor, output_file):
    """Export all oligos in FASTA format"""
    
    cursor.execute("""
        SELECT 
            be.name,
            bs.seq,
            ot.name as oligo_type,
            od.tm_calculated,
            od.gc_content,
            od.design_mode,
            od.specificity_score,
            bd.name as project_name
        FROM oligo_design od
        JOIN bioentry be ON od.bioentry_id = be.bioentry_id
        JOIN biosequence bs ON be.bioentry_id = bs.bioentry_id
        JOIN oligo_type ot ON od.oligo_type_id = ot.oligo_type_id
        JOIN biodatabase bd ON be.biodatabase_id = bd.biodatabase_id
        ORDER BY be.name
    """)
    
    with open(output_file, 'w') as f:
        for row in cursor.fetchall():
            name, sequence, oligo_type, tm, gc, design_mode, specificity, project = row
            
            # Create FASTA header with metadata
            header = f">{name}"
            
            # Add metadata to header
            metadata = []
            if tm:
                metadata.append(f"tm={tm:.1f}")
            if gc:
                metadata.append(f"gc={gc:.1f}")
            if specificity:
                metadata.append(f"specificity={specificity:.3f}")
            
            metadata.extend([
                f"type={oligo_type}",
                f"mode={design_mode}",
                f"project={project}"
            ])
            
            if metadata:
                header += " " + "|".join(metadata)
            
            f.write(header + "\n")
            f.write(sequence + "\n")
    
    print(f"FASTA export completed: {output_file}")

def export_csv(cursor, output_file):
    """Export all oligos in CSV format"""
    
    cursor.execute("""
        SELECT 
            be.name as oligo_name,
            be.accession,
            bs.seq as sequence,
            bs.length,
            ot.name as oligo_type,
            od.design_mode,
            od.tm_calculated,
            od.gc_content,
            od.specificity_score,
            od.off_target_count,
            od.validation_status,
            od.target_region,
            bd.name as project_name,
            dm.name as design_method,
            od.created_date
        FROM oligo_design od
        JOIN bioentry be ON od.bioentry_id = be.bioentry_id
        JOIN biosequence bs ON be.bioentry_id = bs.bioentry_id
        JOIN oligo_type ot ON od.oligo_type_id = ot.oligo_type_id
        JOIN design_method dm ON od.design_method_id = dm.design_method_id
        JOIN biodatabase bd ON be.biodatabase_id = bd.biodatabase_id
        ORDER BY be.name
    """)
    
    # Get column names
    columns = [desc[0] for desc in cursor.description]
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        
        for row in cursor.fetchall():
            # Convert datetime to string
            row_data = []
            for item in row:
                if isinstance(item, datetime):
                    row_data.append(item.isoformat())
                else:
                    row_data.append(item)
            writer.writerow(row_data)
    
    print(f"CSV export completed: {output_file}")

def export_json(cursor, output_file):
    """Export all oligos in JSON format"""
    
    cursor.execute("""
        SELECT 
            od.oligo_design_id,
            be.name as oligo_name,
            be.accession,
            be.description,
            bs.seq as sequence,
            bs.length,
            ot.name as oligo_type,
            od.design_mode,
            od.target_region,
            od.tm_calculated,
            od.gc_content,
            od.specificity_score,
            od.off_target_count,
            od.validation_status,
            bd.name as project_name,
            dm.name as design_method,
            dm.version as method_version,
            od.created_date,
            od.modified_date
        FROM oligo_design od
        JOIN bioentry be ON od.bioentry_id = be.bioentry_id
        JOIN biosequence bs ON be.bioentry_id = bs.bioentry_id
        JOIN oligo_type ot ON od.oligo_type_id = ot.oligo_type_id
        JOIN design_method dm ON od.design_method_id = dm.design_method_id
        JOIN biodatabase bd ON be.biodatabase_id = bd.biodatabase_id
        ORDER BY be.name
    """)
    
    # Convert to list of dictionaries
    columns = [desc[0] for desc in cursor.description]
    oligos = []
    
    for row in cursor.fetchall():
        oligo_dict = {}
        for i, value in enumerate(row):
            if isinstance(value, datetime):
                oligo_dict[columns[i]] = value.isoformat()
            else:
                oligo_dict[columns[i]] = value
        oligos.append(oligo_dict)
    
    # Get summary statistics
    cursor.execute("""
        SELECT 
            COUNT(*) as total_oligos,
            COUNT(DISTINCT bd.biodatabase_id) as total_projects,
            AVG(od.tm_calculated) as avg_tm,
            AVG(od.gc_content) as avg_gc,
            AVG(od.specificity_score) as avg_specificity
        FROM oligo_design od
        JOIN bioentry be ON od.bioentry_id = be.bioentry_id
        JOIN biodatabase bd ON be.biodatabase_id = bd.biodatabase_id
        WHERE od.tm_calculated IS NOT NULL
    """)
    
    stats_row = cursor.fetchone()
    statistics = {
        'total_oligos': stats_row[0],
        'total_projects': stats_row[1],
        'avg_tm': float(stats_row[2]) if stats_row[2] else None,
        'avg_gc': float(stats_row[3]) if stats_row[3] else None,
        'avg_specificity': float(stats_row[4]) if stats_row[4] else None
    }
    
    # Create final JSON structure
    export_data = {
        'metadata': {
            'export_date': datetime.now().isoformat(),
            'total_records': len(oligos),
            'source': 'OligoDesignPipeline Database',
            'format_version': '1.0'
        },
        'statistics': statistics,
        'oligos': oligos
    }
    
    with open(output_file, 'w') as f:
        json.dump(export_data, f, indent=2, default=str)
    
    print(f"JSON export completed: {output_file}")

def export_primer_pairs_csv(cursor, output_file):
    """Export primer pairs in CSV format"""
    
    cursor.execute("""
        SELECT 
            pp.primer_pair_id,
            pp.product_size,
            pp.pair_specificity,
            pp.validation_status as pair_status,
            
            -- Forward primer
            f_be.name as forward_name,
            f_bs.seq as forward_sequence,
            f_od.tm_calculated as forward_tm,
            f_od.gc_content as forward_gc,
            f_od.specificity_score as forward_specificity,
            
            -- Reverse primer
            r_be.name as reverse_name,
            r_bs.seq as reverse_sequence,
            r_od.tm_calculated as reverse_tm,
            r_od.gc_content as reverse_gc,
            r_od.specificity_score as reverse_specificity,
            
            -- Additional info
            f_od.target_region,
            ABS(f_od.tm_calculated - r_od.tm_calculated) as tm_difference,
            f_bd.name as project_name,
            pp.created_date
            
        FROM primer_pair pp
        
        -- Forward primer join
        JOIN oligo_design f_od ON pp.forward_oligo_id = f_od.oligo_design_id
        JOIN bioentry f_be ON f_od.bioentry_id = f_be.bioentry_id
        JOIN biosequence f_bs ON f_be.bioentry_id = f_bs.bioentry_id
        JOIN biodatabase f_bd ON f_be.biodatabase_id = f_bd.biodatabase_id
        
        -- Reverse primer join
        JOIN oligo_design r_od ON pp.reverse_oligo_id = r_od.oligo_design_id
        JOIN bioentry r_be ON r_od.bioentry_id = r_be.bioentry_id
        JOIN biosequence r_bs ON r_be.bioentry_id = r_bs.bioentry_id
        
        ORDER BY pp.primer_pair_id
    """)
    
    columns = [desc[0] for desc in cursor.description]
    
    primer_pairs_file = output_file.replace('.csv', '_primer_pairs.csv')
    
    with open(primer_pairs_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        
        for row in cursor.fetchall():
            # Convert datetime to string
            row_data = []
            for item in row:
                if isinstance(item, datetime):
                    row_data.append(item.isoformat())
                else:
                    row_data.append(item)
            writer.writerow(row_data)
    
    print(f"Primer pairs CSV export completed: {primer_pairs_file}")
    return primer_pairs_file

def export_off_targets_csv(cursor, output_file):
    """Export off-target analysis results"""
    
    cursor.execute("""
        SELECT 
            be.name as oligo_name,
            oth.chromosome,
            oth.start_pos,
            oth.end_pos,
            oth.strand,
            oth.identity_percent,
            oth.mismatch_count,
            oth.gap_count,
            oth.alignment_length,
            oth.e_value,
            oth.bit_score,
            oth.gene_symbol,
            oth.feature_type,
            od.design_mode,
            ot.name as oligo_type
        FROM off_target_hit oth
        JOIN oligo_design od ON oth.oligo_design_id = od.oligo_design_id
        JOIN bioentry be ON od.bioentry_id = be.bioentry_id
        JOIN oligo_type ot ON od.oligo_type_id = ot.oligo_type_id
        WHERE oth.identity_percent >= 80  -- Only significant hits
        ORDER BY be.name, oth.identity_percent DESC
    """)
    
    columns = [desc[0] for desc in cursor.description]
    off_targets_file = output_file.replace('.csv', '_off_targets.csv')
    
    with open(off_targets_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(cursor.fetchall())
    
    print(f"Off-targets CSV export completed: {off_targets_file}")
    return off_targets_file

def main():
    """Main function for Snakemake script"""
    
    # Get parameters from snakemake object
    db_config = snakemake.params.db_config
    fasta_output = snakemake.output.fasta
    csv_output = snakemake.output.csv
    json_output = snakemake.output.json
    log_file = snakemake.log[0]
    
    try:
        # Redirect output to log file
        with open(log_file, 'w') as log:
            original_stdout = sys.stdout
            sys.stdout = log
            sys.stderr = log
            
            print("Starting database export")
            print(f"Database: {db_config.get('database', 'oligodb')}")
            print("-" * 50)
            
            # Create output directory
            os.makedirs(os.path.dirname(fasta_output), exist_ok=True)
            
            # Connect to database
            conn = connect_to_database(db_config)
            cursor = conn.cursor()
            
            # Export in different formats
            print("Exporting FASTA format...")
            export_fasta(cursor, fasta_output)
            
            print("Exporting CSV format...")
            export_csv(cursor, csv_output)
            
            print("Exporting JSON format...")
            export_json(cursor, json_output)
            
            # Export additional files
            print("Exporting primer pairs...")
            export_primer_pairs_csv(cursor, csv_output)
            
            print("Exporting off-target analysis...")
            export_off_targets_csv(cursor, csv_output)
            
            cursor.close()
            conn.close()
            
            print("Database export completed successfully")
            
            # Restore stdout for any final output
            sys.stdout = original_stdout
        
    except Exception as e:
        with open(log_file, 'w') as log:
            print(f"Error during database export: {e}", file=log)
        raise

if __name__ == "__main__":
    main()
