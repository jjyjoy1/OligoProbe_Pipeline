# =============================================================================
# Index Construction Rules (BLAST and Bowtie2)
# =============================================================================

# Build BLAST database for specificity screening
rule build_blast_index:
    input:
        fasta=lambda w: get_fasta_for_mode(w.mode)
    output:
        done="results/indexes/{mode}/blast/done.txt",
        db="results/indexes/{mode}/blast/{build}.ndb".format(build=GENOME_BUILD)
    params:
        prefix="results/indexes/{mode}/blast/{build}".format(build=GENOME_BUILD),
        title=lambda w: f"{GENOME_BUILD}_{w.mode}_database"
    log:
        "logs/blast_index_{mode}.log"
    threads: get_threads("blast_index")
    resources:
        mem_mb=lambda w: get_memory("blast_index").replace("G", "000").replace("M", "")
    shell:
        """
        module load {config[modules][blast]}
        
        echo "Building BLAST database for {wildcards.mode} mode..." > {log}
        
        # Create output directory
        mkdir -p results/indexes/{wildcards.mode}/blast
        
        # Build BLAST database
        makeblastdb \
            -in {input.fasta} \
            -dbtype nucl \
            -out {params.prefix} \
            -title "{params.title}" \
            -parse_seqids \
            -logfile {log}
        
        # Create completion marker
        echo "BLAST database built: $(date)" > {output.done}
        echo "Database files:" >> {output.done}
        ls -la results/indexes/{wildcards.mode}/blast/ >> {output.done}
        
        echo "BLAST index completed: $(date)" >> {log}
        """

# Build Bowtie2 index for alignment-based screening
rule build_bowtie2_index:
    input:
        fasta=lambda w: get_fasta_for_mode(w.mode)
    output:
        done="results/indexes/{mode}/bowtie2/done.txt",
        index="results/indexes/{mode}/bowtie2/{build}.1.bt2".format(build=GENOME_BUILD)
    params:
        prefix="results/indexes/{mode}/bowtie2/{build}".format(build=GENOME_BUILD)
    log:
        "logs/bowtie2_index_{mode}.log"
    threads: get_threads("bowtie2_index")
    resources:
        mem_mb=lambda w: get_memory("bowtie2_index").replace("G", "000").replace("M", "")
    shell:
        """
        module load {config[modules][bowtie2]}
        
        echo "Building Bowtie2 index for {wildcards.mode} mode..." > {log}
        
        # Create output directory
        mkdir -p results/indexes/{wildcards.mode}/bowtie2
        
        # Build Bowtie2 index
        bowtie2-build \
            --threads {threads} \
            {input.fasta} \
            {params.prefix} \
            2>> {log}
        
        # Create completion marker
        echo "Bowtie2 index built: $(date)" > {output.done}
        echo "Index files:" >> {output.done}
        ls -la results/indexes/{wildcards.mode}/bowtie2/ >> {output.done}
        
        echo "Bowtie2 index completed: $(date)" >> {log}
        """

# Build test indexes (smaller, for development)
rule build_test_blast_index:
    input:
        "resources/test_chr22.fa"
    output:
        done="results/indexes/{mode}/blast/test_done.txt"
    params:
        prefix="results/indexes/{mode}/blast/test_chr22"
    log:
        "logs/test_blast_index_{mode}.log"
    threads: 2
    shell:
        """
        module load {config[modules][blast]}
        
        echo "Building test BLAST database..." > {log}
        
        mkdir -p results/indexes/{wildcards.mode}/blast
        
        makeblastdb \
            -in {input} \
            -dbtype nucl \
            -out {params.prefix} \
            -title "Test_chr22_{wildcards.mode}" \
            -parse_seqids \
            -logfile {log}
        
        echo "Test BLAST database built: $(date)" > {output.done}
        echo "Test BLAST index completed: $(date)" >> {log}
        """

rule build_test_bowtie2_index:
    input:
        "resources/test_chr22.fa"
    output:
        done="results/indexes/{mode}/bowtie2/test_done.txt"
    params:
        prefix="results/indexes/{mode}/bowtie2/test_chr22"
    log:
        "logs/test_bowtie2_index_{mode}.log"
    threads: 4
    shell:
        """
        module load {config[modules][bowtie2]}
        
        echo "Building test Bowtie2 index..." > {log}
        
        mkdir -p results/indexes/{wildcards.mode}/bowtie2
        
        bowtie2-build \
            --threads {threads} \
            {input} \
            {params.prefix} \
            2>> {log}
        
        echo "Test Bowtie2 index built: $(date)" > {output.done}
        echo "Test Bowtie2 index completed: $(date)" >> {log}
        """

# Validate index integrity
rule validate_indexes:
    input:
        blast="results/indexes/{mode}/blast/done.txt",
        bowtie2="results/indexes/{mode}/bowtie2/done.txt"
    output:
        "results/qc/index_validation_{mode}.txt"
    params:
        blast_db="results/indexes/{mode}/blast/{build}".format(build=GENOME_BUILD),
        bowtie2_db="results/indexes/{mode}/bowtie2/{build}".format(build=GENOME_BUILD)
    log:
        "logs/validate_indexes_{mode}.log"
    shell:
        """
        module load {config[modules][blast]}
        module load {config[modules][bowtie2]}
        
        echo "Validating indexes for {wildcards.mode} mode..." > {log}
        
        mkdir -p results/qc
        
        echo "Index Validation Report - {wildcards.mode} mode" > {output}
        echo "Generated: $(date)" >> {output}
        echo "==========================================" >> {output}
        echo "" >> {output}
        
        # Test BLAST database
        echo "Testing BLAST database..." >> {output}
        blastdbcmd -db {params.blast_db} -info >> {output} 2>> {log}
        
        if [ $? -eq 0 ]; then
            echo "✓ BLAST database is valid" >> {output}
        else
            echo "✗ BLAST database validation failed" >> {output}
            exit 1
        fi
        
        echo "" >> {output}
        
        # Test Bowtie2 index
        echo "Testing Bowtie2 index..." >> {output}
        bowtie2-inspect -s {params.bowtie2_db} | head -5 >> {output} 2>> {log}
        
        if [ $? -eq 0 ]; then
            echo "✓ Bowtie2 index is valid" >> {output}
        else
            echo "✗ Bowtie2 index validation failed" >> {output}
            exit 1
        fi
        
        echo "" >> {output}
        echo "All indexes validated successfully" >> {output}
        echo "Index validation completed: $(date)" >> {log}
        """

# Index statistics and information
rule index_stats:
    input:
        blast="results/indexes/{mode}/blast/done.txt",
        bowtie2="results/indexes/{mode}/bowtie2/done.txt"
    output:
        "results/stats/index_stats_{mode}.txt"
    params:
        blast_db="results/indexes/{mode}/blast/{build}".format(build=GENOME_BUILD),
        bowtie2_db="results/indexes/{mode}/bowtie2/{build}".format(build=GENOME_BUILD)
    log:
        "logs/index_stats_{mode}.log"
    shell:
        """
        module load {config[modules][blast]}
        module load {config[modules][bowtie2]}
        
        echo "Collecting index statistics..." > {log}
        
        mkdir -p results/stats
        
        echo "Index Statistics - {wildcards.mode} mode" > {output}
        echo "Generated: $(date)" >> {output}
        echo "===============================" >> {output}
        echo "" >> {output}
        
        # BLAST database info
        echo "BLAST Database Information:" >> {output}
        blastdbcmd -db {params.blast_db} -info >> {output} 2>> {log}
        echo "" >> {output}
        
        # File sizes
        echo "File Sizes:" >> {output}
        echo "BLAST files:" >> {output}
        du -h results/indexes/{wildcards.mode}/blast/* >> {output}
        echo "" >> {output}
        echo "Bowtie2 files:" >> {output}
        du -h results/indexes/{wildcards.mode}/bowtie2/* >> {output}
        echo "" >> {output}
        
        # Total space usage
        TOTAL_SIZE=$(du -sh results/indexes/{wildcards.mode}/ | cut -f1)
        echo "Total index size for {wildcards.mode} mode: $TOTAL_SIZE" >> {output}
        
        echo "Index statistics completed: $(date)" >> {log}
        """

# Clean up old indexes (utility rule)
rule clean_indexes:
    shell:
        """
        echo "Cleaning old indexes..."
        rm -rf results/indexes/*/blast/*
        rm -rf results/indexes/*/bowtie2/*
        echo "Index cleanup completed"
        """

# Archive indexes for backup
rule archive_indexes:
    input:
        expand("results/indexes/{mode}/blast/done.txt", mode=DESIGN_MODES),
        expand("results/indexes/{mode}/bowtie2/done.txt", mode=DESIGN_MODES)
    output:
        "results/archives/indexes_{build}.tar.gz".format(build=GENOME_BUILD)
    log:
        "logs/archive_indexes.log"
    shell:
        """
        echo "Creating index archive..." > {log}
        
        mkdir -p results/archives
        
        tar -czf {output} results/indexes/ 2>> {log}
        
        echo "Archive created: {output}" >> {log}
        echo "Archive size: $(du -h {output} | cut -f1)" >> {log}
        echo "Index archiving completed: $(date)" >> {log}
        """

# Index maintenance check
rule index_maintenance:
    input:
        expand("results/indexes/{mode}/blast/done.txt", mode=DESIGN_MODES),
        expand("results/indexes/{mode}/bowtie2/done.txt", mode=DESIGN_MODES)
    output:
        "results/maintenance/index_check.txt"
    log:
        "logs/index_maintenance.log"
    shell:
        """
        echo "Performing index maintenance check..." > {log}
        
        mkdir -p results/maintenance
        
        echo "Index Maintenance Report" > {output}
        echo "Generated: $(date)" >> {output}
        echo "========================" >> {output}
        echo "" >> {output}
        
        # Check index ages
        echo "Index Creation Times:" >> {output}
        for mode in {DESIGN_MODES}; do
            echo "Mode: $mode" >> {output}
            if [ -f "results/indexes/$mode/blast/done.txt" ]; then
                echo "  BLAST: $(cat results/indexes/$mode/blast/done.txt | head -1)" >> {output}
            fi
            if [ -f "results/indexes/$mode/bowtie2/done.txt" ]; then
                echo "  Bowtie2: $(cat results/indexes/$mode/bowtie2/done.txt | head -1)" >> {output}
            fi
            echo "" >> {output}
        done
        
        # Check disk usage
        echo "Disk Usage:" >> {output}
        du -sh results/indexes/*/ >> {output}
        echo "" >> {output}
        
        # Recommendations
        DAYS_OLD=$(find results/indexes -name "done.txt" -mtime +30 | wc -l)
        if [ $DAYS_OLD -gt 0 ]; then
            echo "RECOMMENDATION: Some indexes are >30 days old. Consider rebuilding." >> {output}
        else
            echo "STATUS: All indexes are current." >> {output}
        fi
        
        echo "Maintenance check completed: $(date)" >> {log}
        """


