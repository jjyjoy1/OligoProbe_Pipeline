# =============================================================================
# Validation and Quality Control Rules
# =============================================================================

# Validate test primers against the database
rule validate_test_primers:
    input:
        blast_done=expand("results/indexes/{mode}/blast/done.txt", mode=DESIGN_MODES),
        bowtie2_done=expand("results/indexes/{mode}/bowtie2/done.txt", mode=DESIGN_MODES)
    output:
        "results/validation/test_results.txt"
    params:
        test_primers=config["validation"]["test_primers"],
        blast_db_pcr="results/indexes/pcr/blast/{build}".format(build=GENOME_BUILD),
        blast_db_panel="results/indexes/panel/blast/{build}".format(build=GENOME_BUILD),
        max_off_targets=config["validation"]["max_acceptable_off_targets"]
    log:
        "logs/validate_test_primers.log"
    threads: 4
    shell:
        """
        module load {config[modules][blast]}
        
        echo "Validating test primers..." > {log}
        
        mkdir -p results/validation
        
        # Create test primer FASTA
        cat > test_primers.fa << EOF
>{params.test_primers[0][name]}
{params.test_primers[0][sequence]}
>{params.test_primers[1][name]}
{params.test_primers[1][sequence]}
>{params.test_primers[2][name]}
{params.test_primers[2][sequence]}
>{params.test_primers[3][name]}
{params.test_primers[3][sequence]}
EOF
        
        echo "Test Primer Validation Report" > {output}
        echo "Generated: $(date)" >> {output}
        echo "==============================" >> {output}
        echo "" >> {output}
        
        # Test against PCR database
        echo "Testing against PCR database..." >> {output}
        blastn -query test_primers.fa -db {params.blast_db_pcr} \
               -outfmt "6 qseqid sseqid pident length mismatch" \
               -max_target_seqs 50 -evalue 1000 \
               > pcr_blast_results.txt 2>> {log}
        
        PCR_HITS=$(cat pcr_blast_results.txt | wc -l)
        echo "PCR database hits: $PCR_HITS" >> {output}
        
        # Test against Panel database  
        echo "Testing against Panel database..." >> {output}
        blastn -query test_primers.fa -db {params.blast_db_panel} \
               -outfmt "6 qseqid sseqid pident length mismatch" \
               -max_target_seqs 50 -evalue 1000 \
               > panel_blast_results.txt 2>> {log}
        
        PANEL_HITS=$(cat panel_blast_results.txt | wc -l)
        echo "Panel database hits: $PANEL_HITS" >> {output}
        echo "" >> {output}
        
        # Analyze results
        echo "Detailed Results:" >> {output}
        echo "=================" >> {output}
        
        for primer in GAPDH_F GAPDH_R ACTB_F ACTB_R; do
            echo "Primer: $primer" >> {output}
            
            PCR_COUNT=$(grep "^$primer" pcr_blast_results.txt | wc -l)
            PANEL_COUNT=$(grep "^$primer" panel_blast_results.txt | wc -l)
            
            echo "  PCR hits: $PCR_COUNT" >> {output}
            echo "  Panel hits: $PANEL_COUNT" >> {output}
            
            # Check if within acceptable range
            if [ $PCR_COUNT -le {params.max_off_targets} ] && [ $PANEL_COUNT -le {params.max_off_targets} ]; then
                echo "  Status: PASS" >> {output}
            else
                echo "  Status: FAIL (too many off-targets)" >> {output}
            fi
            echo "" >> {output}
        done
        
        # Overall validation status
        FAILED_PRIMERS=$(grep "Status: FAIL" {output} | wc -l)
        
        echo "Overall Validation:" >> {output}
        echo "==================" >> {output}
        
        if [ $FAILED_PRIMERS -eq 0 ]; then
            echo "✓ All test primers PASSED validation" >> {output}
            echo "✓ Databases are functioning correctly" >> {output}
        else
            echo "✗ $FAILED_PRIMERS test primers FAILED validation" >> {output}
            echo "✗ Check database integrity or primer sequences" >> {output}
        fi
        
        # Cleanup
        rm -f test_primers.fa pcr_blast_results.txt panel_blast_results.txt
        
        echo "Test primer validation completed: $(date)" >> {log}
        """

# Validate pipeline outputs
rule validate_pipeline_outputs:
    input:
        pcr_candidates="results/design/pcr/candidates.tsv",
        pcr_filtered="results/design/pcr/filtered.tsv",
        panel_candidates="results/design/panel/candidates.tsv", 
        panel_filtered="results/design/panel/filtered.tsv"
    output:
        "results/validation/pipeline_validation.txt"
    params:
        min_candidates=config["qc"]["min_records"]["candidates"],
        min_filtered=config["qc"]["min_records"]["filtered"]
    log:
        "logs/validate_pipeline_outputs.log"
    shell:
        """
        echo "Validating pipeline outputs..." > {log}
        
        echo "Pipeline Output Validation Report" > {output}
        echo "Generated: $(date)" >> {output}
        echo "==================================" >> {output}
        echo "" >> {output}
        
        # Check file existence and sizes
        echo "File Existence and Size Checks:" >> {output}
        echo "===============================" >> {output}
        
        for file in {input}; do
            if [ -f "$file" ]; then
                SIZE=$(stat -c%s "$file")
                LINES=$(wc -l < "$file")
                echo "✓ $file: $SIZE bytes, $LINES lines" >> {output}
            else
                echo "✗ $file: MISSING" >> {output}
            fi
        done
        
        echo "" >> {output}
        
        # Check record counts
        echo "Record Count Validation:" >> {output}
        echo "========================" >> {output}
        
        # PCR candidates
        PCR_CANDIDATES=$(tail -n +2 {input.pcr_candidates} | wc -l)
        PCR_FILTERED=$(tail -n +2 {input.pcr_filtered} | wc -l)
        
        echo "PCR Mode:" >> {output}
        echo "  Candidates: $PCR_CANDIDATES (min: {params.min_candidates})" >> {output}
        echo "  Filtered: $PCR_FILTERED (min: {params.min_filtered})" >> {output}
        
        if [ $PCR_CANDIDATES -ge {params.min_candidates} ]; then
            echo "  ✓ PCR candidates count OK" >> {output}
        else
            echo "  ✗ PCR candidates count too low" >> {output}
        fi
        
        if [ $PCR_FILTERED -ge {params.min_filtered} ]; then
            echo "  ✓ PCR filtered count OK" >> {output}
        else
            echo "  ✗ PCR filtered count too low" >> {output}
        fi
        
        echo "" >> {output}
        
        # Panel candidates
        PANEL_CANDIDATES=$(tail -n +2 {input.panel_candidates} | wc -l)
        PANEL_FILTERED=$(tail -n +2 {input.panel_filtered} | wc -l)
        
        echo "Panel Mode:" >> {output}
        echo "  Candidates: $PANEL_CANDIDATES (min: {params.min_candidates})" >> {output}
        echo "  Filtered: $PANEL_FILTERED (min: {params.min_filtered})" >> {output}
        
        if [ $PANEL_CANDIDATES -ge {params.min_candidates} ]; then
            echo "  ✓ Panel candidates count OK" >> {output}
        else
            echo "  ✗ Panel candidates count too low" >> {output}
        fi
        
        if [ $PANEL_FILTERED -ge {params.min_filtered} ]; then
            echo "  ✓ Panel filtered count OK" >> {output}
        else
            echo "  ✗ Panel filtered count too low" >> {output}
        fi
        
        echo "" >> {output}
        
        # Content validation
        echo "Content Validation:" >> {output}
        echo "==================" >> {output}
        
        # Check for required columns
        REQUIRED_COLS="primer_id type sequence length tm gc_content"
        
        for file in {input.pcr_candidates} {input.pcr_filtered} {input.panel_candidates} {input.panel_filtered}; do
            echo "Checking columns in $file:" >> {output}
            HEADER=$(head -1 "$file")
            
            for col in $REQUIRED_COLS; do
                if echo "$HEADER" | grep -q "$col"; then
                    echo "  ✓ $col" >> {output}
                else
                    echo "  ✗ $col (missing)" >> {output}
                fi
            done
            echo "" >> {output}
        done
        
        # Overall validation status
        ERRORS=$(grep "✗" {output} | wc -l)
        
        echo "Overall Status:" >> {output}
        echo "===============" >> {output}
        
        if [ $ERRORS -eq 0 ]; then
            echo "✓ ALL VALIDATIONS PASSED" >> {output}
            echo "✓ Pipeline outputs are valid and complete" >> {output}
        else
            echo "✗ $ERRORS VALIDATION ERRORS FOUND" >> {output}
            echo "✗ Check individual validations above" >> {output}
        fi
        
        echo "Pipeline output validation completed: $(date)" >> {log}
        """

# Performance benchmarking
rule benchmark_pipeline:
    input:
        "results/validation/test_results.txt",
        "results/validation/pipeline_validation.txt"
    output:
        "results/validation/performance_benchmark.txt"
    log:
        "logs/benchmark_pipeline.log"
    shell:
        """
        echo "Running performance benchmark..." > {log}
        
        echo "Pipeline Performance Benchmark" > {output}
        echo "Generated: $(date)" >> {output}
        echo "==============================" >> {output}
        echo "" >> {output}
        
        # Collect timing information from logs
        echo "Timing Information:" >> {output}
        echo "==================" >> {output}
        
        for logfile in logs/*.log; do
            if [ -f "$logfile" ]; then
                BASENAME=$(basename "$logfile" .log)
                START_TIME=$(head -1 "$logfile" | grep -o '[0-9][0-9]:[0-9][0-9]:[0-9][0-9]' | head -1)
                END_TIME=$(tail -5 "$logfile" | grep -o '[0-9][0-9]:[0-9][0-9]:[0-9][0-9]' | tail -1)
                
                if [ -n "$START_TIME" ] && [ -n "$END_TIME" ]; then
                    echo "$BASENAME: $START_TIME - $END_TIME" >> {output}
                fi
            fi
        done
        
        echo "" >> {output}
        
        # Resource usage
        echo "Resource Usage:" >> {output}
        echo "===============" >> {output}
        
        # Disk usage
        echo "Disk Usage:" >> {output}
        du -sh resources/ results/ 2>/dev/null >> {output} || echo "Could not calculate disk usage" >> {output}
        echo "" >> {output}
        
        # File counts
        echo "File Counts:" >> {output}
        echo "Resources: $(find resources/ -type f 2>/dev/null | wc -l) files" >> {output}
        echo "Results: $(find results/ -type f 2>/dev/null | wc -l) files" >> {output}
        echo "Logs: $(find logs/ -type f 2>/dev/null | wc -l) files" >> {output}
        echo "" >> {output}
        
        # Performance recommendations
        echo "Performance Recommendations:" >> {output}
        echo "============================" >> {output}
        
        TOTAL_SIZE=$(du -sm results/ 2>/dev/null | cut -f1 || echo "0")
        
        if [ "$TOTAL_SIZE" -gt 1000 ]; then
            echo "• Results directory >1GB - consider cleanup" >> {output}
        else
            echo "• Results directory size reasonable" >> {output}
        fi
        
        LOG_COUNT=$(find logs/ -type f 2>/dev/null | wc -l)
        if [ "$LOG_COUNT" -gt 50 ]; then
            echo "• Many log files ($LOG_COUNT) - consider log rotation" >> {output}
        else
            echo "• Log file count reasonable" >> {output}
        fi
        
        echo "• For faster runs, increase thread counts in config.yaml" >> {output}
        echo "• Consider using test mode for development" >> {output}
        
        echo "Performance benchmark completed: $(date)" >> {log}
        """

# Integration test with known sequences
rule integration_test:
    input:
        indexes_done=expand("results/indexes/{mode}/blast/done.txt", mode=DESIGN_MODES)
    output:
        "results/validation/integration_test.txt"
    log:
        "logs/integration_test.log"
    shell:
        """
        echo "Running integration tests..." > {log}
        
        echo "Integration Test Report" > {output}
        echo "Generated: $(date)" >> {output}
        echo "=======================" >> {output}
        echo "" >> {output}
        
        # Test 1: Database connectivity
        echo "Test 1: Database Connectivity" >> {output}
        echo "==============================" >> {output}
        
        module load {config[modules][blast]}
        
        for mode in pcr panel; do
            DB="results/indexes/$mode/blast/{GENOME_BUILD}"
            echo "Testing $mode database..." >> {output}
            
            if blastdbcmd -db "$DB" -info >/dev/null 2>&1; then
                echo "✓ $mode database accessible" >> {output}
            else
                echo "✗ $mode database not accessible" >> {output}
            fi
        done
        
        echo "" >> {output}
        
        # Test 2: Sample sequence search
        echo "Test 2: Sample Sequence Search" >> {output}
        echo "===============================" >> {output}
        
        # Create a simple test sequence
        echo ">test_seq" > test_seq.fa
        echo "GCATACGTTGTATCCGGGCAT" >> test_seq.fa
        
        for mode in pcr panel; do
            DB="results/indexes/$mode/blast/{GENOME_BUILD}"
            echo "Searching $mode database..." >> {output}
            
            HITS=$(blastn -query test_seq.fa -db "$DB" -outfmt 6 | wc -l)
            echo "  Hits found: $HITS" >> {output}
            
            if [ "$HITS" -gt 0 ]; then
                echo "  ✓ Search successful" >> {output}
            else
                echo "  ✗ No hits found (potential issue)" >> {output}
            fi
        done
        
        echo "" >> {output}
        
        # Test 3: File integrity
        echo "Test 3: File Integrity" >> {output}
        echo "=======================" >> {output}
        
        # Check critical files exist and are non-empty
        CRITICAL_FILES="resources/{GENOME_BUILD}.fa resources/{GENOME_BUILD}_panel.fa"
        
        for file in $CRITICAL_FILES; do
            if [ -f "$file" ] && [ -s "$file" ]; then
                echo "✓ $file: OK" >> {output}
            else
                echo "✗ $file: Missing or empty" >> {output}
            fi
        done
        
        echo "" >> {output}
        
        # Overall test result
        FAILURES=$(grep "✗" {output} | wc -l)
        
        echo "Integration Test Summary:" >> {output}
        echo "========================" >> {output}
        
        if [ "$FAILURES" -eq 0 ]; then
            echo "✓ ALL TESTS PASSED" >> {output}
            echo "✓ Pipeline is ready for production use" >> {output}
        else
            echo "✗ $FAILURES TEST(S) FAILED" >> {output}
            echo "✗ Address failures before production use" >> {output}
        fi
        
        # Cleanup
        rm -f test_seq.fa
        
        echo "Integration test completed: $(date)" >> {log}
        """

# Regression test with version comparison
rule regression_test:
    input:
        current_results=expand("results/design/{mode}/filtered.tsv", mode=DESIGN_MODES)
    output:
        "results/validation/regression_test.txt"
    log:
        "logs/regression_test.log"
    shell:
        """
        echo "Running regression tests..." > {log}
        
        echo "Regression Test Report" > {output}
        echo "Generated: $(date)" >> {output}
        echo "======================" >> {output}
        echo "" >> {output}
        
        # Check if baseline results exist
        if [ -d "test/baseline_results" ]; then
            echo "Comparing against baseline results..." >> {output}
            echo "====================================" >> {output}
            
            for mode in pcr panel; do
                CURRENT="results/design/$mode/filtered.tsv"
                BASELINE="test/baseline_results/design_$mode\_filtered.tsv"
                
                if [ -f "$BASELINE" ]; then
                    CURRENT_COUNT=$(tail -n +2 "$CURRENT" | wc -l)
                    BASELINE_COUNT=$(tail -n +2 "$BASELINE" | wc -l)
                    
                    echo "$mode mode comparison:" >> {output}
                    echo "  Current: $CURRENT_COUNT candidates" >> {output}
                    echo "  Baseline: $BASELINE_COUNT candidates" >> {output}
                    
                    DIFF=$(echo "$CURRENT_COUNT - $BASELINE_COUNT" | bc 2>/dev/null || echo "0")
                    echo "  Difference: $DIFF candidates" >> {output}
                    
                    # Allow 10% variation
                    THRESHOLD=$(echo "$BASELINE_COUNT * 0.1" | bc 2>/dev/null || echo "1")
                    
                    if [ "${DIFF#-}" -le "${THRESHOLD%.*}" ]; then
                        echo "  ✓ Within acceptable range" >> {output}
                    else
                        echo "  ⚠ Outside acceptable range (>10% change)" >> {output}
                    fi
                    echo "" >> {output}
                else
                    echo "$mode mode: No baseline available" >> {output}
                fi
            done
        else
            echo "No baseline results found." >> {output}
            echo "To enable regression testing:" >> {output}
            echo "1. Create test/baseline_results/ directory" >> {output}
            echo "2. Copy current results as baseline" >> {output}
            echo "3. Re-run regression test" >> {output}
            echo "" >> {output}
            
            echo "Creating baseline from current results..." >> {output}
            mkdir -p test/baseline_results
            
            for mode in pcr panel; do
                cp "results/design/$mode/filtered.tsv" "test/baseline_results/design_$mode\_filtered.tsv"
                echo "✓ Baseline created for $mode mode" >> {output}
            done
        fi
        
        echo "Regression test completed: $(date)" >> {log}
        """

# Comprehensive validation report
rule validation_summary:
    input:
        test_results="results/validation/test_results.txt",
        pipeline_validation="results/validation/pipeline_validation.txt", 
        integration_test="results/validation/integration_test.txt",
        benchmark="results/validation/performance_benchmark.txt"
    output:
        "results/validation/validation_summary.html"
    log:
        "logs/validation_summary.log"
    shell:
        """
        echo "Creating validation summary..." > {log}
        
        cat > {output} << 'HTML_EOF'
<!DOCTYPE html>
<html>
<head>
    <title>Oligo Pipeline Validation Summary</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
        .section {{ margin: 20px 0; }}
        .pass {{ color: green; font-weight: bold; }}
        .fail {{ color: red; font-weight: bold; }}
        .warn {{ color: orange; font-weight: bold; }}
        pre {{ background-color: #f8f8f8; padding: 10px; border-left: 3px solid #ccc; }}
        .summary {{ background-color: #e8f4f8; padding: 15px; border-radius: 5px; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Oligo Design Pipeline - Validation Summary</h1>
        <p>Generated: $(date)</p>
        <p>Pipeline Version: {GENOME_BUILD}</p>
    </div>
    
    <div class="summary">
        <h2>Overall Status</h2>
HTML_EOF
        
        # Count total passes and failures
        TOTAL_PASS=$(grep -h "✓" {input} | wc -l)
        TOTAL_FAIL=$(grep -h "✗" {input} | wc -l)
        
        if [ "$TOTAL_FAIL" -eq 0 ]; then
            echo "        <p class=\"pass\">✓ ALL VALIDATIONS PASSED ($TOTAL_PASS checks)</p>" >> {output}
        else
            echo "        <p class=\"fail\">✗ $TOTAL_FAIL validation(s) failed out of $(($TOTAL_PASS + $TOTAL_FAIL)) total checks</p>" >> {output}
        fi
        
        cat >> {output} << 'HTML_EOF'
    </div>
    
    <div class="section">
        <h2>Test Results</h2>
        <pre>
HTML_EOF
        
        cat {input.test_results} >> {output}
        
        cat >> {output} << 'HTML_EOF'
        </pre>
    </div>
    
    <div class="section">
        <h2>Pipeline Validation</h2>
        <pre>
HTML_EOF
        
        cat {input.pipeline_validation} >> {output}
        
        cat >> {output} << 'HTML_EOF'
        </pre>
    </div>
    
    <div class="section">
        <h2>Integration Test</h2>
        <pre>
HTML_EOF
        
        cat {input.integration_test} >> {output}
        
        cat >> {output} << 'HTML_EOF'
        </pre>
    </div>
    
    <div class="section">
        <h2>Performance Benchmark</h2>
        <pre>
HTML_EOF
        
        cat {input.benchmark} >> {output}
        
        cat >> {output} << 'HTML_EOF'
        </pre>
    </div>
    
    <div class="section">
        <h2>Recommendations</h2>
        <ul>
            <li>Review any failed validations above</li>
            <li>Monitor performance metrics for optimization opportunities</li>
            <li>Update baseline results after significant pipeline changes</li>
            <li>Run full validation before production deployment</li>
        </ul>
    </div>
</body>
</html>
HTML_EOF
        
        echo "Validation summary created: $(date)" >> {log}
        """
