#!/usr/bin/env python3
"""
Generate database summary report for oligo design results
"""

import psycopg2
import sys
from datetime import datetime
import json

def connect_to_database(db_config):
    """Connect to PostgreSQL database"""
    return psycopg2.connect(
        host=db_config.get('host', 'localhost'),
        port=db_config.get('port', 5432),
        database=db_config.get('database', 'oligodb'),
        user=db_config.get('user', 'postgres'),
        password=db_config.get('password', '')
    )

def get_database_statistics(cursor):
    """Get comprehensive database statistics"""
    stats = {}
    
    # Basic counts
    queries = {
        'total_projects': "SELECT COUNT(*) FROM biodatabase",
        'total_oligos': "SELECT COUNT(*) FROM oligo_design",
        'total_sequences': "SELECT COUNT(*) FROM biosequence",
        'pcr_primers': """
            SELECT COUNT(*) FROM oligo_design od 
            JOIN oligo_type ot ON od.oligo_type_id = ot.oligo_type_id 
            WHERE ot.name IN ('forward_primer', 'reverse_primer')
        """,
        'capture_probes': """
            SELECT COUNT(*) FROM oligo_design od 
            JOIN oligo_type ot ON od.oligo_type_id = ot.oligo_type_id 
            WHERE ot.name = 'capture_probe'
        """,
        'primer_pairs': "SELECT COUNT(*) FROM primer_pair",
        'off_target_hits': "SELECT COUNT(*) FROM off_target_hit",
        'validated_oligos': "SELECT COUNT(*) FROM oligo_design WHERE validation_status = 'validated'",
    }
    
    for key, query in queries.items():
        cursor.execute(query)
        stats[key] = cursor.fetchone()[0]
    
    # Design mode breakdown
    cursor.execute("""
        SELECT design_mode, COUNT(*) 
        FROM oligo_design 
        GROUP BY design_mode
    """)
    stats['by_design_mode'] = dict(cursor.fetchall())
    
    # Validation status breakdown
    cursor.execute("""
        SELECT validation_status, COUNT(*) 
        FROM oligo_design 
        GROUP BY validation_status
    """)
    stats['by_validation_status'] = dict(cursor.fetchall())
    
    # Quality metrics
    cursor.execute("""
        SELECT 
            AVG(tm_calculated) as avg_tm,
            MIN(tm_calculated) as min_tm,
            MAX(tm_calculated) as max_tm,
            AVG(gc_content) as avg_gc,
            MIN(gc_content) as min_gc,
            MAX(gc_content) as max_gc,
            AVG(specificity_score) as avg_specificity,
            AVG(off_target_count) as avg_off_targets
        FROM oligo_design 
        WHERE tm_calculated IS NOT NULL
    """)
    quality_metrics = cursor.fetchone()
    stats['quality_metrics'] = {
        'avg_tm': float(quality_metrics[0]) if quality_metrics[0] else 0,
        'min_tm': float(quality_metrics[1]) if quality_metrics[1] else 0,
        'max_tm': float(quality_metrics[2]) if quality_metrics[2] else 0,
        'avg_gc': float(quality_metrics[3]) if quality_metrics[3] else 0,
        'min_gc': float(quality_metrics[4]) if quality_metrics[4] else 0,
        'max_gc': float(quality_metrics[5]) if quality_metrics[5] else 0,
        'avg_specificity': float(quality_metrics[6]) if quality_metrics[6] else 0,
        'avg_off_targets': float(quality_metrics[7]) if quality_metrics[7] else 0,
    }
    
    # Recent activity
    cursor.execute("""
        SELECT COUNT(*) FROM oligo_design 
        WHERE created_date >= CURRENT_DATE - INTERVAL '7 days'
    """)
    stats['recent_oligos'] = cursor.fetchone()[0]
    
    return stats

def generate_html_report(stats, output_file):
    """Generate HTML summary report"""
    
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Oligo Database Summary</title>
    <style>
        body {{ 
            font-family: 'Segoe UI', Arial, sans-serif; 
            margin: 20px; 
            background-color: #f5f5f5; 
        }}
        .container {{ 
            max-width: 1200px; 
            margin: 0 auto; 
            background: white; 
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
        .stats-grid {{ 
            display: grid; 
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); 
            gap: 20px; 
            margin: 20px 0; 
        }}
        .stat-card {{ 
            background: #f8f9ff; 
            padding: 20px; 
            border-radius: 8px; 
            border-left: 4px solid #667eea; 
        }}
        .stat-number {{ 
            font-size: 2.5em; 
            font-weight: bold; 
            color: #667eea; 
            margin: 10px 0; 
        }}
        .stat-label {{ 
            color: #666; 
            font-size: 1.1em; 
        }}
        .quality-metrics {{ 
            background: #e8f5e8; 
            padding: 20px; 
            border-radius: 8px; 
            margin: 20px 0; 
        }}
        .breakdown {{ 
            display: grid; 
            grid-template-columns: 1fr 1fr; 
            gap: 20px; 
            margin: 20px 0; 
        }}
        .breakdown-card {{ 
            background: white; 
            padding: 20px; 
            border-radius: 8px; 
            border: 1px solid #ddd; 
        }}
        table {{ 
            width: 100%; 
            border-collapse: collapse; 
            margin: 10px 0; 
        }}
        th, td {{ 
            border: 1px solid #ddd; 
            padding: 8px; 
            text-align: left; 
        }}
        th {{ 
            background-color: #667eea; 
            color: white; 
        }}
        .metric-value {{ 
            font-weight: bold; 
            color: #333; 
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 Oligo Database Summary</h1>
            <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-number">{stats['total_oligos']:,}</div>
                <div class="stat-label">Total Oligos</div>
            </div>
            <div class="stat-card">
                <div class="stat-number">{stats['pcr_primers']:,}</div>
                <div class="stat-label">PCR Primers</div>
            </div>
            <div class="stat-card">
                <div class="stat-number">{stats['capture_probes']:,}</div>
                <div class="stat-label">Capture Probes</div>
            </div>
            <div class="stat-card">
                <div class="stat-number">{stats['primer_pairs']:,}</div>
                <div class="stat-label">Primer Pairs</div>
            </div>
            <div class="stat-card">
                <div class="stat-number">{stats['off_target_hits']:,}</div>
                <div class="stat-label">Off-target Hits</div>
            </div>
            <div class="stat-card">
                <div class="stat-number">{stats['validated_oligos']:,}</div>
                <div class="stat-label">Validated Oligos</div>
            </div>
        </div>
        
        <div class="quality-metrics">
            <h2>📊 Quality Metrics</h2>
            <table>
                <tr>
                    <th>Metric</th>
                    <th>Average</th>
                    <th>Min</th>
                    <th>Max</th>
                </tr>
                <tr>
                    <td>Melting Temperature (°C)</td>
                    <td class="metric-value">{stats['quality_metrics']['avg_tm']:.1f}</td>
                    <td>{stats['quality_metrics']['min_tm']:.1f}</td>
                    <td>{stats['quality_metrics']['max_tm']:.1f}</td>
                </tr>
                <tr>
                    <td>GC Content (%)</td>
                    <td class="metric-value">{stats['quality_metrics']['avg_gc']:.1f}</td>
                    <td>{stats['quality_metrics']['min_gc']:.1f}</td>
                    <td>{stats['quality_metrics']['max_gc']:.1f}</td>
                </tr>
                <tr>
                    <td>Specificity Score</td>
                    <td class="metric-value">{stats['quality_metrics']['avg_specificity']:.3f}</td>
                    <td>-</td>
                    <td>-</td>
                </tr>
                <tr>
                    <td>Off-target Count</td>
                    <td class="metric-value">{stats['quality_metrics']['avg_off_targets']:.1f}</td>
                    <td>-</td>
                    <td>-</td>
                </tr>
            </table>
        </div>
        
        <div class="breakdown">
            <div class="breakdown-card">
                <h3>By Design Mode</h3>
                <table>
                    <tr><th>Mode</th><th>Count</th></tr>"""
    
    for mode, count in stats['by_design_mode'].items():
        html_content += f"<tr><td>{mode}</td><td>{count:,}</td></tr>"
    
    html_content += """
                </table>
            </div>
            <div class="breakdown-card">
                <h3>By Validation Status</h3>
                <table>
                    <tr><th>Status</th><th>Count</th></tr>"""
    
    for status, count in stats['by_validation_status'].items():
        html_content += f"<tr><td>{status}</td><td>{count:,}</td></tr>"
    
    html_content += f"""
                </table>
            </div>
        </div>
        
        <div style="margin-top: 30px; padding: 20px; background: #f0f8ff; border-radius: 8px;">
            <h3>📈 Recent Activity</h3>
            <p><strong>{stats['recent_oligos']:,}</strong> oligos added in the last 7 days</p>
        </div>
        
        <div style="margin-top: 20px; font-size: 0.9em; color: #666; text-align: center;">
            <p>Generated by Oligo Design Pipeline Database Module</p>
        </div>
    </div>
</body>
</html>
    """
    
    with open(output_file, 'w') as f:
        f.write(html_content)

def generate_text_stats(stats, output_file):
    """Generate text-based statistics file"""
    
    with open(output_file, 'w') as f:
        f.write("OLIGO DATABASE STATISTICS\n")
        f.write("=" * 50 + "\n")
        f.write(f"Generated: {datetime.now()}\n\n")
        
        f.write("OVERVIEW\n")
        f.write("-" * 20 + "\n")
        f.write(f"Total Projects: {stats['total_projects']:,}\n")
        f.write(f"Total Oligos: {stats['total_oligos']:,}\n")
        f.write(f"Total Sequences: {stats['total_sequences']:,}\n")
        f.write(f"PCR Primers: {stats['pcr_primers']:,}\n")
        f.write(f"Capture Probes: {stats['capture_probes']:,}\n")
        f.write(f"Primer Pairs: {stats['primer_pairs']:,}\n")
        f.write(f"Off-target Hits: {stats['off_target_hits']:,}\n")
        f.write(f"Validated Oligos: {stats['validated_oligos']:,}\n\n")
        
        f.write("DESIGN MODE BREAKDOWN\n")
        f.write("-" * 20 + "\n")
        for mode, count in stats['by_design_mode'].items():
            f.write(f"{mode}: {count:,}\n")
        f.write("\n")
        
        f.write("VALIDATION STATUS\n")
        f.write("-" * 20 + "\n")
        for status, count in stats['by_validation_status'].items():
            f.write(f"{status}: {count:,}\n")
        f.write("\n")
        
        f.write("QUALITY METRICS\n")
        f.write("-" * 20 + "\n")
        qm = stats['quality_metrics']
        f.write(f"Average Tm: {qm['avg_tm']:.1f}°C (range: {qm['min_tm']:.1f}-{qm['max_tm']:.1f})\n")
        f.write(f"Average GC: {qm['avg_gc']:.1f}% (range: {qm['min_gc']:.1f}-{qm['max_gc']:.1f})\n")
        f.write(f"Average Specificity: {qm['avg_specificity']:.3f}\n")
        f.write(f"Average Off-targets: {qm['avg_off_targets']:.1f}\n\n")
        
        f.write("RECENT ACTIVITY\n")
        f.write("-" * 20 + "\n")
        f.write(f"Oligos added (last 7 days): {stats['recent_oligos']:,}\n")

def main():
    """Main function for Snakemake script"""
    
    # Get parameters from snakemake object
    db_config = snakemake.params.db_config
    html_output = snakemake.output.summary
    stats_output = snakemake.output.stats
    log_file = snakemake.log[0]
    
    try:
        # Redirect output to log file
        with open(log_file, 'w') as log:
            sys.stdout = log
            sys.stderr = log
            
            print("Generating database summary report")
            print(f"Database: {db_config.get('database', 'oligodb')}")
            print("-" * 50)
            
            # Connect to database and get statistics
            conn = connect_to_database(db_config)
            cursor = conn.cursor()
            
            print("Collecting database statistics...")
            stats = get_database_statistics(cursor)
            
            cursor.close()
            conn.close()
            
            print("Generating HTML report...")
            generate_html_report(stats, html_output)
            
            print("Generating text statistics...")
            generate_text_stats(stats, stats_output)
            
            print("Database summary generation completed successfully")
        
    except Exception as e:
        with open(log_file, 'w') as log:
            print(f"Error generating database summary: {e}", file=log)
        raise

if __name__ == "__main__":
    main()

