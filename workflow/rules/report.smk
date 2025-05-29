# =============================================================================
# Reporting and Visualization Rules
# =============================================================================

# Generate comprehensive pipeline summary report
rule pipeline_summary:
    input:
        pcr_filtered="results/design/pcr/filtered.tsv",
        panel_filtered="results/design/panel/filtered.tsv",
        pcr_stats="results/stats/design_stats_pcr.txt",
        panel_stats="results/stats/design_stats_panel.txt",
        validation="results/validation/test_results.txt"
    output:
        html="results/reports/pipeline_summary.html",
        txt="results/reports/pipeline_summary.txt"
    log:
        "logs/pipeline_summary.log"
    shell:
        """
        echo "Generating pipeline summary report..." > {log}
        
        mkdir -p results/reports
        
        # Create text summary first
        cat > {output.txt} << EOF
OLIGO DESIGN PIPELINE SUMMARY REPORT
====================================
Generated: $(date)
Pipeline Version: {GENOME_BUILD}
Design Modes: {DESIGN_MODES}

RESULTS OVERVIEW
================
EOF
        
        # Add PCR results
        echo "" >> {output.txt}
        echo "PCR Mode Results:" >> {output.txt}
        PCR_COUNT=$(tail -n +2 {input.pcr_filtered} | wc -l)
        echo "  Final candidates: $PCR_COUNT" >> {output.txt}
        
        if [ $PCR_COUNT -gt 0 ]; then
            echo "  Sample candidates:" >> {output.txt}
            head -6 {input.pcr_filtered} | tail -n +2 | while read line; do
                echo "    $line" >> {output.txt}
            done
        fi
        
        # Add Panel results
        echo "" >> {output.txt}
        echo "Panel Mode Results:" >> {output.txt}
        PANEL_COUNT=$(tail -n +2 {input.panel_filtered} | wc -l)
        echo "  Final candidates: $PANEL_COUNT" >> {output.txt}
        
        if [ $PANEL_COUNT -gt 0 ]; then
            echo "  Sample candidates:" >> {output.txt}
            head -6 {input.panel_filtered} | tail -n +2 | while read line; do
                echo "    $line" >> {output.txt}
            done
        fi
        
        # Add statistics
        echo "" >> {output.txt}
        echo "DETAILED STATISTICS" >> {output.txt}
        echo "===================" >> {output.txt}
        echo "" >> {output.txt}
        echo "PCR Statistics:" >> {output.txt}
        cat {input.pcr_stats} >> {output.txt}
        echo "" >> {output.txt}
        echo "Panel Statistics:" >> {output.txt}
        cat {input.panel_stats} >> {output.txt}
        
        # Add validation results
        echo "" >> {output.txt}
        echo "VALIDATION RESULTS" >> {output.txt}
        echo "==================" >> {output.txt}
        cat {input.validation} >> {output.txt}
        
        # Now create HTML version
        cat > {output.html} << 'HTML_START'
<!DOCTYPE html>
<html>
<head>
    <title>Oligo Design Pipeline Summary</title>
    <style>
        body {{ 
            font-family: 'Segoe UI', Arial, sans-serif; 
            line-height: 1.6; 
            margin: 0; 
            padding: 20px; 
            background-color: #f5f5f5; 
        }}
        .container {{ 
            max-width: 1200px; 
            margin: 0 auto; 
            background-color: white; 
            padding: 30px; 
            border-radius: 10px; 
            box-shadow: 0 0 20px rgba(0,0,0,0.1); 
        }}
        .header {{ 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
            color: white; 
            padding: 30px; 
            border-radius: 10px; 
            margin-bottom: 30px; 
            text-align: center; 
        }}
        .header h1 {{ margin: 0; font-size: 2.5em; }}
        .header p {{ margin: 10px 0 0 0; opacity: 0.9; }}
        .section {{ 
            margin: 30px 0; 
            padding: 20px; 
            border-left: 4px solid #667eea; 
            background-color: #f8f9ff; 
        }}
        .section h2 {{ 
            color: #333; 
            margin-top: 0; 
            border-bottom: 2px solid #eee; 
            padding-bottom: 10px; 
        }}
        .results-grid {{ 
            display: grid; 
            grid-template-columns: 1fr 1fr; 
            gap: 20px; 
            margin: 20px 0; 
        }}
        .result-card {{ 
            background: white; 
            padding: 20px; 
            border-radius: 8px; 
            border: 1px solid #ddd; 
            box-shadow: 0 2px 4px rgba(0,0,0,0.1); 
        }}
        .result-card h3 {{ 
            color: #667eea; 
            margin-top: 0; 
        }}
        .metric {{ 
            font-size: 2em; 
            font-weight: bold; 
            color: #333; 
            text-align: center; 
            margin: 10px 0; 
        }}
        .candidates-table {{ 
            width: 100%; 
            border-collapse: collapse; 
            margin: 15px 0; 
        }}
        .candidates-table th, .candidates-table td {{ 
            border: 1px solid #ddd; 
            padding: 8px; 
            text-align: left; 
        }}
        .candidates-table th {{ 
            background-color: #667eea; 
            color: white; 
        }}
        .candidates-table tr:nth-child(even) {{ 
            background-color: #f2f2f2; 
        }}
        .status {{ 
            padding: 5px 10px; 
            border-radius: 20px; 
            font-weight: bold; 
            text-align: center; 
            display: inline-block; 
        }}
        .status.success {{ background-color: #d4edda; color: #155724; }}
        .status.warning {{ background-color: #fff3cd; color: #856404; }}
        .status.error {{ background-color: #f8d7da; color: #721c24; }}
        .download-links {{ 
            background-color: #e8f4f8; 
            padding: 20px; 
            border-radius: 8px; 
            margin: 20px 0; 
        }}
        .download-links h3 {{ margin-top: 0; }}
        .download-links a {{ 
            display: inline-block; 
            margin: 5px 10px 5px 0; 
            padding: 8px 15px; 
            background-color: #667eea; 
            color: white; 
            text-decoration: none; 
            border-radius: 5px; 
        }}
        .download-links a:hover {{ background-color: #5568d3; }}
        pre {{ 
            background-color: #f4f4f4; 
            padding: 15px; 
            border-radius: 5px; 
            overflow-x: auto; 
            font-size: 0.9em; 
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 Oligo Design Pipeline</h1>
            <p>Comprehensive Summary Report</p>
            <p>Generated: $(date) | Build: {GENOME_BUILD}</p>
        </div>
HTML_START
        
        # Add results overview section
        cat >> {output.html} << 'HTML_OVERVIEW'
        <div class="section">
            <h2>📊 Results Overview</h2>
            <div class="results-grid">
                <div class="result-card">
                    <h3>PCR Primers</h3>
                    <div class="metric">PCR_COUNT_PLACEHOLDER</div>
                    <p>High-quality primer pairs for PCR amplification</p>
                </div>
                <div class="result-card">
                    <h3>Capture Probes</h3>
                    <div class="metric">PANEL_COUNT_PLACEHOLDER</div>
                    <p>Specific probes for hybrid capture panels</p>
                </div>
            </div>
        </div>
HTML_OVERVIEW
        
        # Replace placeholders with actual counts
        sed -i "s/PCR_COUNT_PLACEHOLDER/$PCR_COUNT/g" {output.html}
        sed -i "s/PANEL_COUNT_PLACEHOLDER/$PANEL_COUNT/g" {output.html}
        
        # Add sample results tables
        cat >> {output.html} << 'HTML_SAMPLES'
        <div class="section">
            <h2>🔬 Sample Results</h2>
            
            <h3>PCR Primers (Top 5)</h3>
            <table class="candidates-table">
                <tr>
                    <th>Primer ID</th>
                    <th>Type</th>
                    <th>Sequence</th>
                    <th>Length</th>
                    <th>Tm (°C)</th>
                    <th>GC%</th>
                </tr>
HTML_SAMPLES
        
        # Add PCR primer rows
        head -6 {input.pcr_filtered} | tail -n +2 | while IFS=$'\t' read -r id type seq len tm gc rest; do
            cat >> {output.html} << HTML_ROW
                <tr>
                    <td>$id</td>
                    <td>$type</td>
                    <td style="font-family: monospace; font-size: 0.8em;">${{seq:0:30}}...</td>
                    <td>$len</td>
                    <td>$tm</td>
                    <td>$gc</td>
                </tr>
HTML_ROW
        done
        
        cat >> {output.html} << 'HTML_PANEL_TABLE'
            </table>
            
            <h3>Capture Probes (Top 5)</h3>
            <table class="candidates-table">
                <tr>
                    <th>Probe ID</th>
                    <th>Type</th>
                    <th>Sequence</th>
                    <th>Length</th>
                    <th>Tm (°C)</th>
                    <th>GC%</th>
                </tr>
HTML_PANEL_TABLE
        
        # Add panel probe rows
        head -6 {input.panel_filtered} | tail -n +2 | while IFS=$'\t' read -r id type seq len tm gc rest; do
            cat >> {output.html} << HTML_ROW
                <tr>
                    <td>$id</td>
                    <td>$type</td>
                    <td style="font-family: monospace; font-size: 0.8em;">${{seq:0:30}}...</td>
                    <td>$len</td>
                    <td>$tm</td>
                    <td>$gc</td>
                </tr>
HTML_ROW
        done
        
        cat >> {output.html} << 'HTML_DOWNLOADS'
            </table>
        </div>
        
        <div class="section">
            <h2>📥 Download Results</h2>
            <div class="download-links">
                <h3>Available Files:</h3>
                <a href="../design/pcr/filtered.tsv">PCR Primers (TSV)</a>
                <a href="../design/panel/filtered.tsv">Capture Probes (TSV)</a>
                <a href="../design/pcr/filtered.fa">PCR Primers (FASTA)</a>
                <a href="../design/panel/filtered.fa">Capture Probes (FASTA)</a>
                <a href="../design/pcr/filtered.bed">PCR Primers (BED)</a>
                <a href="../design/panel/filtered.bed">Capture Probes (BED)</a>
                <a href="pipeline_summary.txt">Text Summary</a>
            </div>
        </div>
        
        <div class="section">
            <h2>📈 Statistics</h2>
            <h3>PCR Mode</h3>
            <pre>
HTML_DOWNLOADS
        
        cat {input.pcr_stats} >> {output.html}
        
        cat >> {output.html} << 'HTML_STATS2'
            </pre>
            
            <h3>Panel Mode</h3>
            <pre>
HTML_STATS2
        
        cat {input.panel_stats} >> {output.html}
        
        cat >> {output.html} << 'HTML_VALIDATION'
            </pre>
        </div>
        
        <div class="section">
            <h2>✅ Validation Results</h2>
            <pre>
HTML_VALIDATION
        
        cat {input.validation} >> {output.html}
        
        cat >> {output.html} << 'HTML_FOOTER'
            </pre>
        </div>
        
        <div class="section">
            <h2>ℹ️ Pipeline Information</h2>
            <p><strong>Workflow Engine:</strong> Snakemake</p>
            <p><strong>Reference Build:</strong> {GENOME_BUILD}</p>
            <p><strong>Design Modes:</strong> PCR, Hybrid Capture Panel</p>
            <p><strong>Generated:</strong> $(date)</p>
        </div>
        
    </div>
</body>
</html>
HTML_FOOTER
        
        echo "Pipeline summary report completed: $(date)" >> {log}
        """

# Generate design-specific reports for each mode
rule mode_specific_report:
    input:
        candidates="results/design/{mode}/candidates.tsv",
        filtered="results/design/{mode}/filtered.tsv",
        rejected="results/design/{mode}/rejected.tsv",
        summary="results/design/{mode}/filtering_summary.txt"
    output:
        "results/reports/{mode}_detailed_report.html"
    log:
        "logs/mode_report_{mode}.log"
    shell:
        """
        echo "Generating {wildcards.mode} mode detailed report..." > {log}
        
        cat > {output} << 'HTML_START'
<!DOCTYPE html>
<html>
<head>
    <title>Oligo Design Report - MODE_PLACEHOLDER Mode</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; }}
        .header {{ background-color: #f0f8ff; padding: 20px; border-radius: 8px; margin-bottom: 20px; }}
        .section {{ margin: 20px 0; }}
        .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; }}
        .stat-card {{ background: #f9f9f9; padding: 15px; border-radius: 5px; text-align: center; }}
        .stat-number {{ font-size: 2em; font-weight: bold; color: #2c5aa0; }}
        table {{ width: 100%; border-collapse: collapse; margin: 15px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 10px; text-align: left; }}
        th {{ background-color: #2c5aa0; color: white; }}
        tr:nth-child(even) {{ background-color: #f2f2f2; }}
        .sequence {{ font-family: monospace; background-color: #f0f0f0; padding: 2px 4px; }}
        .pass {{ color: green; font-weight: bold; }}
        .fail {{ color: red; font-weight: bold; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 Oligo Design Report</h1>
        <h2>MODE_PLACEHOLDER Mode Analysis</h2>
        <p>Generated: $(date)</p>
    </div>
HTML_START
        
        # Replace mode placeholder
        sed -i "s/MODE_PLACEHOLDER/{wildcards.mode}/g" {output}
        
        # Add statistics
        TOTAL_CANDIDATES=$(tail -n +2 {input.candidates} | wc -l)
        FILTERED_CANDIDATES=$(tail -n +2 {input.filtered} | wc -l)
        REJECTED_CANDIDATES=$(tail -n +2 {input.rejected} | wc -l)
        SUCCESS_RATE=$(echo "scale=1; $FILTERED_CANDIDATES * 100 / $TOTAL_CANDIDATES" | bc 2>/dev/null || echo "0")
        
        cat >> {output} << HTML_STATS
    <div class="section">
        <h2>📊 Summary Statistics</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-number">$TOTAL_CANDIDATES</div>
                <div>Total Candidates</div>
            </div>
            <div class="stat-card">
                <div class="stat-number">$FILTERED_CANDIDATES</div>
                <div>Passed Filtering</div>
            </div>
            <div class="stat-card">
                <div class="stat-number">$REJECTED_CANDIDATES</div>
                <div>Rejected</div>
            </div>
            <div class="stat-card">
                <div class="stat-number">$SUCCESS_RATE%</div>
                <div>Success Rate</div>
            </div>
        </div>
    </div>
HTML_STATS
        
        # Add filtered candidates table
        cat >> {output} << 'HTML_TABLE'
    <div class="section">
        <h2>✅ Accepted Candidates</h2>
        <table>
            <thead>
                <tr>
                    <th>ID</th>
                    <th>Type</th>
                    <th>Sequence</th>
                    <th>Length</th>
                    <th>Tm (°C)</th>
                    <th>GC%</th>
                    <th>Off-targets</th>
                    <th>Specificity</th>
                </tr>
            </thead>
            <tbody>
HTML_TABLE
        
        # Add table rows (limit to first 20 for readability)
        head -21 {input.filtered} | tail -n +2 | while IFS=$'\t' read -r id type seq len tm gc off_targets specificity rest; do
            # Truncate long sequences
            display_seq="${{seq:0:40}}"
            if [ ${{#seq}} -gt 40 ]; then
                display_seq="$display_seq..."
            fi
            
            cat >> {output} << HTML_ROW
                <tr>
                    <td>$id</td>
                    <td>$type</td>
                    <td class="sequence">$display_seq</td>
                    <td>$len</td>
                    <td>$tm</td>
                    <td>$gc</td>
                    <td>${{off_targets:-0}}</td>
                    <td>${{specificity:-1.0}}</td>
                </tr>
HTML_ROW
        done
        
        cat >> {output} << 'HTML_FILTERING'
            </tbody>
        </table>
    </div>
    
    <div class="section">
        <h2>🔍 Filtering Summary</h2>
        <pre>
HTML_FILTERING
        
        cat {input.summary} >> {output}
        
        cat >> {output} << 'HTML_END'
        </pre>
    </div>
    
    <div class="section">
        <h2>📥 Download Options</h2>
        <ul>
            <li><a href="../design/{wildcards.mode}/filtered.tsv">Accepted candidates (TSV)</a></li>
            <li><a href="../design/{wildcards.mode}/filtered.fa">Accepted candidates (FASTA)</a></li>
            <li><a href="../design/{wildcards.mode}/filtered.bed">Accepted candidates (BED)</a></li>
            <li><a href="../design/{wildcards.mode}/rejected.tsv">Rejected candidates (TSV)</a></li>
            <li><a href="../design/{wildcards.mode}/candidates.tsv">All candidates (TSV)</a></li>
        </ul>
    </div>
    
</body>
</html>
HTML_END
        
        echo "{wildcards.mode} mode detailed report completed: $(date)" >> {log}
        """

# Create a simple dashboard with key metrics
rule create_dashboard:
    input:
        expand("results/design/{mode}/filtered.tsv", mode=DESIGN_MODES),
        expand("results/stats/design_stats_{mode}.txt", mode=DESIGN_MODES)
    output:
        "results/reports/dashboard.html"
    log:
        "logs/create_dashboard.log"
    shell:
        """
        echo "Creating pipeline dashboard..." > {log}
        
        cat > {output} << 'DASHBOARD_START'
<!DOCTYPE html>
<html>
<head>
    <title>Oligo Pipeline Dashboard</title>
    <meta http-equiv="refresh" content="30">
    <style>
        body {{ 
            font-family: 'Segoe UI', Arial, sans-serif; 
            margin: 0; 
            padding: 20px; 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }}
        .dashboard {{ 
            max-width: 1200px; 
            margin: 0 auto; 
        }}
        .header {{ 
            background: rgba(255,255,255,0.95); 
            padding: 30px; 
            border-radius: 15px; 
            text-align: center; 
            margin-bottom: 30px;
            box-shadow: 0 8px 32px rgba(0,0,0,0.1);
        }}
        .header h1 {{ 
            margin: 0; 
            color: #333; 
            font-size: 2.5em;
        }}
        .metrics-grid {{ 
            display: grid; 
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); 
            gap: 20px; 
            margin-bottom: 30px;
        }}
        .metric-card {{ 
            background: rgba(255,255,255,0.95); 
            padding: 30px; 
            border-radius: 15px; 
            text-align: center;
            box-shadow: 0 8px 32px rgba(0,0,0,0.1);
            transition: transform 0.3s ease;
        }}
        .metric-card:hover {{ 
            transform: translateY(-5px);
        }}
        .metric-number {{ 
            font-size: 3em; 
            font-weight: bold; 
            color: #667eea; 
            margin: 10px 0;
        }}
        .metric-label {{ 
            color: #666; 
            font-size: 1.1em;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        .status-indicator {{ 
            display: inline-block; 
            width: 12px; 
            height: 12px; 
            border-radius: 50%; 
            margin-right: 8px;
        }}
        .status-running {{ background-color: #ffc107; }}
        .status-complete {{ background-color: #28a745; }}
        .status-error {{ background-color: #dc3545; }}
        .quick-links {{ 
            background: rgba(255,255,255,0.95); 
            padding: 25px; 
            border-radius: 15px;
            box-shadow: 0 8px 32px rgba(0,0,0,0.1);
        }}
        .quick-links h2 {{ 
            margin-top: 0; 
            color: #333;
        }}
        .link-grid {{ 
            display: grid; 
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); 
            gap: 15px;
        }}
        .link-card {{ 
            background: #f8f9fa; 
            padding: 15px; 
            border-radius: 8px; 
            text-align: center;
            transition: background-color 0.3s ease;
        }}
        .link-card:hover {{ 
            background: #e9ecef;
        }}
        .link-card a {{ 
            text-decoration: none; 
            color: #667eea; 
            font-weight: bold;
        }}
        .timestamp {{ 
            color: #888; 
            font-size: 0.9em; 
            margin-top: 10px;
        }}
    </style>
</head>
<body>
    <div class="dashboard">
        <div class="header">
            <h1>🧬 Oligo Design Pipeline</h1>
            <p>Real-time Dashboard</p>
            <div class="timestamp">Last updated: $(date)</div>
        </div>
        
        <div class="metrics-grid">
DASHBOARD_START
        
        # Add metrics for each mode
        for mode in pcr panel; do
            if [ -f "results/design/$mode/filtered.tsv" ]; then
                COUNT=$(tail -n +2 "results/design/$mode/filtered.tsv" | wc -l)
                STATUS="complete"
                MODE_DISPLAY=$(echo "$mode" | tr '[:lower:]' '[:upper:]')
            else
                COUNT="0"
                STATUS="running"
                MODE_DISPLAY=$(echo "$mode" | tr '[:lower:]' '[:upper:]')
            fi
            
            cat >> {output} << METRIC_CARD
            <div class="metric-card">
                <div class="status-indicator status-$STATUS"></div>
                <div class="metric-number">$COUNT</div>
                <div class="metric-label">$MODE_DISPLAY Candidates</div>
            </div>
METRIC_CARD
        done
        
        # Add pipeline status metrics
        TOTAL_FILES=$(find results/ -name "*.tsv" -o -name "*.fa" -o -name "*.bed" 2>/dev/null | wc -l)
        TOTAL_SIZE=$(du -sh results/ 2>/dev/null | cut -f1 || echo "0")
        
        cat >> {output} << 'ADDITIONAL_METRICS'
            <div class="metric-card">
                <div class="status-indicator status-complete"></div>
                <div class="metric-number">TOTAL_FILES_PLACEHOLDER</div>
                <div class="metric-label">Output Files</div>
            </div>
            <div class="metric-card">
                <div class="status-indicator status-complete"></div>
                <div class="metric-number">TOTAL_SIZE_PLACEHOLDER</div>
                <div class="metric-label">Total Size</div>
            </div>
        </div>
        
        <div class="quick-links">
            <h2>📋 Quick Links</h2>
            <div class="link-grid">
                <div class="link-card">
                    <a href="pipeline_summary.html">📊 Full Report</a>
                </div>
                <div class="link-card">
                    <a href="pcr_detailed_report.html">🔬 PCR Details</a>
                </div>
                <div class="link-card">
                    <a href="panel_detailed_report.html">🧪 Panel Details</a>
                </div>
                <div class="link-card">
                    <a href="../design/pcr/filtered.tsv">📄 PCR Results</a>
                </div>
                <div class="link-card">
                    <a href="../design/panel/filtered.tsv">📄 Panel Results</a>
                </div>
                <div class="link-card">
                    <a href="../validation/validation_summary.html">✅ Validation</a>
                </div>
            </div>
        </div>
    </div>
</body>
</html>
ADDITIONAL_METRICS
        
        # Replace placeholders
        sed -i "s/TOTAL_FILES_PLACEHOLDER/$TOTAL_FILES/g" {output}
        sed -i "s/TOTAL_SIZE_PLACEHOLDER/$TOTAL_SIZE/g" {output}
        
        echo "Dashboard created: $(date)" >> {log}
        """

# Generate configuration report
rule config_report:
    output:
        "results/reports/configuration.html"
    log:
        "logs/config_report.log"
    shell:
        """
        echo "Generating configuration report..." > {log}
        
        cat > {output} << 'CONFIG_START'
<!DOCTYPE html>
<html>
<head>
    <title>Pipeline Configuration Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; }}
        .header {{ background-color: #f8f9fa; padding: 20px; border-radius: 8px; }}
        .section {{ margin: 20px 0; }}
        .config-block {{ background-color: #f4f4f4; padding: 15px; border-radius: 5px; font-family: monospace; }}
        table {{ width: 100%; border-collapse: collapse; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #343a40; color: white; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🔧 Pipeline Configuration Report</h1>
        <p>Generated: $(date)</p>
    </div>
    
    <div class="section">
        <h2>Environment Settings</h2>
        <table>
            <tr><th>Parameter</th><th>Value</th></tr>
            <tr><td>Genome Build</td><td>{GENOME_BUILD}</td></tr>
            <tr><td>Design Modes</td><td>{DESIGN_MODES}</td></tr>
            <tr><td>Max Threads</td><td>{config[max_threads]}</td></tr>
            <tr><td>Max Memory</td><td>{config[max_memory]}</td></tr>
        </table>
    </div>
    
    <div class="section">
        <h2>Primer3 Configuration</h2>
        <h3>PCR Mode</h3>
        <div class="config-block">
CONFIG_START
        
        # Add Primer3 PCR config
        for key in primer_opt_size primer_min_size primer_max_size primer_opt_tm primer_min_tm primer_max_tm; do
            value=$(echo "{config[primer3][pcr][$key]}" | sed 's/[{}]//g')
            echo "$key: $value<br>" >> {output}
        done
        
        cat >> {output} << 'CONFIG_PANEL'
        </div>
        
        <h3>Panel Mode</h3>
        <div class="config-block">
CONFIG_PANEL
        
        # Add Primer3 Panel config
        for key in primer_opt_size primer_min_size primer_max_size primer_opt_tm primer_min_tm primer_max_tm; do
            value=$(echo "{config[primer3][panel][$key]}" | sed 's/[{}]//g')
            echo "$key: $value<br>" >> {output}
        done
        
        cat >> {output} << 'CONFIG_END'
        </div>
    </div>
    
    <div class="section">
        <h2>Tool Versions</h2>
        <table>
            <tr><th>Tool</th><th>Module</th></tr>
            <tr><td>BLAST</td><td>{config[modules][blast]}</td></tr>
            <tr><td>Bowtie2</td><td>{config[modules][bowtie2]}</td></tr>
            <tr><td>BCFtools</td><td>{config[modules][bcftools]}</td></tr>
            <tr><td>SAMtools</td><td>{config[modules][samtools]}</td></tr>
            <tr><td>Primer3</td><td>{config[modules][primer3]}</td></tr>
        </table>
    </div>
    
    <div class="section">
        <h2>Quality Control Thresholds</h2>
        <table>
            <tr><th>Parameter</th><th>Value</th></tr>
            <tr><td>Max Off-targets</td><td>{config[screening][max_off_targets]}</td></tr>
            <tr><td>Min Specificity Score</td><td>{config[screening][min_specificity_score]}</td></tr>
            <tr><td>Min Candidates</td><td>{config[qc][min_records][candidates]}</td></tr>
            <tr><td>Min Filtered</td><td>{config[qc][min_records][filtered]}</td></tr>
        </table>
    </div>
    
</body>
</html>
CONFIG_END
        
        echo "Configuration report completed: $(date)" >> {log}
        """

# Archive all reports
rule archive_reports:
    input:
        "results/reports/pipeline_summary.html",
        "results/reports/dashboard.html",
        expand("results/reports/{mode}_detailed_report.html", mode=DESIGN_MODES),
        "results/reports/configuration.html"
    output:
        "results/archives/reports_{build}.tar.gz".format(build=GENOME_BUILD)
    log:
        "logs/archive_reports.log"
    shell:
        """
        echo "Archiving reports..." > {log}
        
        mkdir -p results/archives
        
        tar -czf {output} results/reports/ 2>> {log}
        
        echo "Report archive created: {output}" >> {log}
        echo "Archive size: $(du -h {output} | cut -f1)" >> {log}
        echo "Report archiving completed: $(date)" >> {log}
        """
