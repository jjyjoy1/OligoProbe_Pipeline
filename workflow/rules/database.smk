# =============================================================================
# Database Integration Rules - PostgreSQL + BioSQL
# =============================================================================

# Load database credentials and settings
rule setup_database:
    output:
        "results/database/setup_complete.txt"
    params:
        db_config=config.get("database", {}),
        setup_sql="database/setup_biosql_db.sql"
    log:
        "logs/setup_database.log"
    run:
        import psycopg2
        import os
        
        # Database connection parameters
        db_params = {
            'host': params.db_config.get('host', 'localhost'),
            'port': params.db_config.get('port', 5432),
            'database': params.db_config.get('database', 'oligodb'),
            'user': params.db_config.get('user', 'postgres'),
            'password': params.db_config.get('password', '')
        }
        
        print("Setting up BioSQL database...", file=open(log[0], 'w'))
        
        try:
            # Connect to PostgreSQL
            conn = psycopg2.connect(**db_params)
            conn.autocommit = True
            cur = conn.cursor()
            
            # Read and execute setup SQL
            if os.path.exists(params.setup_sql):
                with open(params.setup_sql, 'r') as f:
                    setup_commands = f.read()
                
                # Execute setup commands
                cur.execute(setup_commands)
                print("Database setup completed successfully", file=open(log[0], 'a'))
            else:
                print(f"Setup SQL file not found: {params.setup_sql}", file=open(log[0], 'a'))
                raise FileNotFoundError(f"Setup SQL file not found: {params.setup_sql}")
            
            # Create completion marker
            with open(output[0], 'w') as f:
                f.write(f"Database setup completed at {datetime.datetime.now()}\n")
                f.write(f"Host: {db_params['host']}\n")
                f.write(f"Database: {db_params['database']}\n")
            
            cur.close()
            conn.close()
            
        except Exception as e:
            print(f"Database setup failed: {e}", file=open(log[0], 'a'))
            raise

# Create project entry for this pipeline run
rule create_project:
    input:
        "results/database/setup_complete.txt"
    output:
        "results/database/project_created.txt"
    params:
        db_config=config.get("database", {}),
        project_name=config.get("project_name", f"OligoDesign_{GENOME_BUILD}"),
        description=config.get("project_description", "Automated oligo design pipeline results")
    log:
        "logs/create_project.log"
    run:
        import psycopg2
        from datetime import datetime
        
        db_params = {
            'host': params.db_config.get('host', 'localhost'),
            'port': params.db_config.get('port', 5432),
            'database': params.db_config.get('database', 'oligodb'),
            'user': params.db_config.get('user', 'postgres'),
            'password': params.db_config.get('password', '')
        }
        
        try:
            conn = psycopg2.connect(**db_params)
            cur = conn.cursor()
            
            # Insert or update project
            cur.execute("""
                INSERT INTO biodatabase (name, authority, description)
                VALUES (%s, %s, %s)
                ON CONFLICT (name) DO UPDATE SET
                    description = EXCLUDED.description,
                    modified_date = CURRENT_TIMESTAMP
                RETURNING biodatabase_id
            """, (params.project_name, "OligoDesignPipeline", params.description))
            
            project_id = cur.fetchone()[0]
            conn.commit()
            
            # Record pipeline run
            cur.execute("""
                INSERT INTO pipeline_run (
                    run_name, pipeline_version, genome_build, 
                    start_time, status, design_parameters
                ) VALUES (%s, %s, %s, %s, %s, %s)
                RETURNING run_id
            """, (
                params.project_name,
                "1.0",
                GENOME_BUILD,
                datetime.now(),
                "running",
                psycopg2.extras.Json(config)
            ))
            
            run_id = cur.fetchone()[0]
            conn.commit()
            
            with open(output[0], 'w') as f:
                f.write(f"Project created: {params.project_name}\n")
                f.write(f"Project ID: {project_id}\n")
                f.write(f"Run ID: {run_id}\n")
            
            cur.close()
            conn.close()
            
        except Exception as e:
            print(f"Project creation failed: {e}", file=open(log[0], 'w'))
            raise

# Load PCR primer results into database
rule load_pcr_results:
    input:
        project="results/database/project_created.txt",
        pcr_results="results/design/pcr/filtered.tsv",
        blast_results="results/design/pcr/blast_results.txt"
    output:
        "results/database/pcr_loaded.txt"
    params:
        db_config=config.get("database", {}),
        project_name=config.get("project_name", f"OligoDesign_{GENOME_BUILD}")
    log:
        "logs/load_pcr_results.log"
    script:
        "scripts/load_oligos_to_db.py"

# Load panel probe results into database  
rule load_panel_results:
    input:
        project="results/database/project_created.txt",
        panel_results="results/design/panel/filtered.tsv",
        blast_results="results/design/panel/blast_results.txt"
    output:
        "results/database/panel_loaded.txt"
    params:
        db_config=config.get("database", {}),
        project_name=config.get("project_name", f"OligoDesign_{GENOME_BUILD}")
    log:
        "logs/load_panel_results.log"
    script:
        "scripts/load_oligos_to_db.py"

# Create primer pairs from individual primers
rule create_primer_pairs:
    input:
        "results/database/pcr_loaded.txt"
    output:
        "results/database/primer_pairs_created.txt"
    params:
        db_config=config.get("database", {})
    log:
        "logs/create_primer_pairs.log"
    run:
        import psycopg2
        import re
        
        db_params = {
            'host': params.db_config.get('host', 'localhost'),
            'port': params.db_config.get('port', 5432),
            'database': params.db_config.get('database', 'oligodb'),
            'user': params.db_config.get('user', 'postgres'),
            'password': params.db_config.get('password', '')
        }
        
        try:
            conn = psycopg2.connect(**db_params)
            cur = conn.cursor()
            
            # Find matching forward and reverse primers
            cur.execute("""
                WITH primer_info AS (
                    SELECT 
                        od.oligo_design_id,
                        be.name,
                        ot.name as primer_type,
                        bs.seq,
                        od.target_region
                    FROM oligo_design od
                    JOIN bioentry be ON od.bioentry_id = be.bioentry_id
                    JOIN biosequence bs ON be.bioentry_id = bs.bioentry_id
                    JOIN oligo_type ot ON od.oligo_type_id = ot.oligo_type_id
                    WHERE ot.name IN ('forward_primer', 'reverse_primer')
                    AND od.design_mode = 'pcr'
                )
                SELECT 
                    f.oligo_design_id as forward_id,
                    r.oligo_design_id as reverse_id,
                    f.name as forward_name,
                    r.name as reverse_name,
                    LENGTH(f.seq) + LENGTH(r.seq) as estimated_product_size
                FROM primer_info f
                JOIN primer_info r ON f.target_region = r.target_region
                WHERE f.primer_type = 'forward_primer'
                AND r.primer_type = 'reverse_primer'
                AND f.name ~ '^(.+)_F$'
                AND r.name ~ '^(.+)_R$'
                AND SUBSTRING(f.name FROM '^(.+)_F$') = SUBSTRING(r.name FROM '^(.+)_R$')
            """)
            
            primer_pairs = cur.fetchall()
            pair_count = 0
            
            for forward_id, reverse_id, forward_name, reverse_name, product_size in primer_pairs:
                # Insert primer pair
                cur.execute("""
                    INSERT INTO primer_pair (
                        forward_oligo_id, reverse_oligo_id, product_size, validation_status
                    ) VALUES (%s, %s, %s, %s)
                    ON CONFLICT DO NOTHING
                """, (forward_id, reverse_id, product_size, 'pending'))
                
                if cur.rowcount > 0:
                    pair_count += 1
            
            conn.commit()
            
            with open(output[0], 'w') as f:
                f.write(f"Created {pair_count} primer pairs from {len(primer_pairs)} candidates\n")
            
            cur.close()
            conn.close()
            
        except Exception as e:
            print(f"Primer pair creation failed: {e}", file=open(log[0], 'w'))
            raise

# Load off-target analysis results
rule load_off_targets:
    input:
        pcr_loaded="results/database/pcr_loaded.txt",
        panel_loaded="results/database/panel_loaded.txt",
        pcr_blast="results/design/pcr/blast_results.txt",
        panel_blast="results/design/panel/blast_results.txt"
    output:
        "results/database/off_targets_loaded.txt"
    params:
        db_config=config.get("database", {})
    log:
        "logs/load_off_targets.log"
    script:
        "scripts/load_off_targets_to_db.py"

# Generate database summary report
rule database_summary:
    input:
        pcr="results/database/pcr_loaded.txt",
        panel="results/database/panel_loaded.txt",
        pairs="results/database/primer_pairs_created.txt",
        off_targets="results/database/off_targets_loaded.txt"
    output:
        summary="results/database/database_summary.html",
        stats="results/database/database_stats.txt"
    params:
        db_config=config.get("database", {})
    log:
        "logs/database_summary.log"
    script:
        "scripts/generate_db_summary.py"

# Database backup and export
rule backup_database:
    input:
        "results/database/database_summary.html"
    output:
        "results/database/backup_{timestamp}.sql"
    params:
        db_config=config.get("database", {}),
        timestamp=datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log:
        "logs/backup_database.log"
    shell:
        """
        pg_dump \
            --host={params.db_config[host]} \
            --port={params.db_config[port]} \
            --username={params.db_config[user]} \
            --dbname={params.db_config[database]} \
            --no-password \
            --verbose \
            --clean \
            --create \
            --file={output} \
        2> {log}
        
        echo "Database backup completed: $(date)" >> {log}
        """

# Database maintenance and optimization
rule optimize_database:
    output:
        "results/database/optimization_complete.txt"
    params:
        db_config=config.get("database", {})
    log:
        "logs/optimize_database.log"
    run:
        import psycopg2
        
        db_params = {
            'host': params.db_config.get('host', 'localhost'),
            'port': params.db_config.get('port', 5432),
            'database': params.db_config.get('database', 'oligodb'),
            'user': params.db_config.get('user', 'postgres'),
            'password': params.db_config.get('password', '')
        }
        
        try:
            conn = psycopg2.connect(**db_params)
            conn.autocommit = True
            cur = conn.cursor()
            
            print("Starting database optimization...", file=open(log[0], 'w'))
            
            # Update table statistics
            cur.execute("ANALYZE;")
            print("Table statistics updated", file=open(log[0], 'a'))
            
            # Vacuum and reindex
            cur.execute("VACUUM ANALYZE;")
            print("Vacuum completed", file=open(log[0], 'a'))
            
            # Reindex trigram indexes
            cur.execute("REINDEX INDEX CONCURRENTLY idx_biosequence_seq_gin;")
            cur.execute("REINDEX INDEX CONCURRENTLY idx_biosequence_seq_gist;")
            print("Trigram indexes rebuilt", file=open(log[0], 'a'))
            
            with open(output[0], 'w') as f:
                f.write(f"Database optimization completed at {datetime.datetime.now()}\n")
            
            cur.close()
            conn.close()
            
        except Exception as e:
            print(f"Database optimization failed: {e}", file=open(log[0], 'a'))
            raise

# Export oligos in various formats from database
rule export_from_database:
    input:
        "results/database/database_summary.html"
    output:
        fasta="results/database/exports/all_oligos.fasta",
        csv="results/database/exports/all_oligos.csv",
        json="results/database/exports/all_oligos.json"
    params:
        db_config=config.get("database", {})
    log:
        "logs/export_from_database.log"
    script:
        "scripts/export_from_db.py"

# Database quality control and validation
rule validate_database:
    input:
        "results/database/off_targets_loaded.txt"
    output:
        "results/database/validation_report.txt"
    params:
        db_config=config.get("database", {})
    log:
        "logs/validate_database.log"
    run:
        import psycopg2
        
        db_params = {
            'host': params.db_config.get('host', 'localhost'),
            'port': params.db_config.get('port', 5432),
            'database': params.db_config.get('database', 'oligodb'),
            'user': params.db_config.get('user', 'postgres'),
            'password': params.db_config.get('password', '')
        }
        
        validation_queries = [
            ("Sequence integrity", "SELECT COUNT(*) FROM biosequence WHERE seq IS NULL OR LENGTH(seq) = 0"),
            ("GC content validity", "SELECT COUNT(*) FROM oligo_design WHERE gc_content < 0 OR gc_content > 100"),
            ("Tm validity", "SELECT COUNT(*) FROM oligo_design WHERE tm_calculated < 0 OR tm_calculated > 100"),
            ("Orphaned sequences", "SELECT COUNT(*) FROM biosequence bs LEFT JOIN bioentry be ON bs.bioentry_id = be.bioentry_id WHERE be.bioentry_id IS NULL"),
            ("Duplicate sequences", "SELECT seq, COUNT(*) as count FROM biosequence GROUP BY seq HAVING COUNT(*) > 1"),
        ]
        
        try:
            conn = psycopg2.connect(**db_params)
            cur = conn.cursor()
            
            with open(output[0], 'w') as report:
                report.write("Database Validation Report\n")
                report.write("=" * 50 + "\n\n")
                
                for check_name, query in validation_queries:
                    cur.execute(query)
                    result = cur.fetchall()
                    
                    report.write(f"{check_name}:\n")
                    if check_name == "Duplicate sequences":
                        if result:
                            report.write(f"  Found {len(result)} duplicate sequences\n")
                            for seq, count in result[:5]:  # Show first 5
                                report.write(f"    {seq[:30]}... (count: {count})\n")
                        else:
                            report.write("  No duplicates found\n")
                    else:
                        count = result[0][0] if result else 0
                        status = "PASS" if count == 0 else "FAIL"
                        report.write(f"  {status}: {count} issues found\n")
                    report.write("\n")
            
            cur.close()
            conn.close()
            
        except Exception as e:
            with open(output[0], 'w') as report:
                report.write(f"Database validation failed: {e}\n")
            raise
