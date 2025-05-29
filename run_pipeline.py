#!/usr/bin/env python3
"""
Oligo Design Pipeline Launcher

A user-friendly launcher for the Snakemake oligo design pipeline.
Provides simple command-line interface with common options.
"""

import argparse
import os
import sys
import subprocess
import yaml
from pathlib import Path
import time
import webbrowser

def print_banner():
    """Print the pipeline banner"""
    banner = """
╔══════════════════════════════════════════════════════════════════════════════╗
║                           🧬 Oligo Design Pipeline                           ║
║                              Version 1.0                                    ║
╚══════════════════════════════════════════════════════════════════════════════╝
    """
    print(banner)

def check_dependencies():
    """Check if required dependencies are available"""
    deps = {
        'snakemake': 'Snakemake workflow engine',
        'python3': 'Python 3 interpreter', 
        'wget': 'File download utility'
    }
    
    missing = []
    for dep, desc in deps.items():
        if subprocess.run(['which', dep], capture_output=True).returncode != 0:
            missing.append(f"{dep} ({desc})")
    
    if missing:
        print("❌ Missing dependencies:")
        for dep in missing:
            print(f"   - {dep}")
        print("\nPlease install missing dependencies before continuing.")
        return False
    
    print("✅ All dependencies found")
    return True

def load_config():
    """Load and validate configuration"""
    if not os.path.exists('config.yaml'):
        print("❌ config.yaml not found in current directory")
        print("Please run this script from the pipeline root directory")
        return None
    
    try:
        with open('config.yaml', 'r') as f:
            config = yaml.safe_load(f)
        print("✅ Configuration loaded successfully")
        return config
    except Exception as e:
        print(f"❌ Error loading config.yaml: {e}")
        return None

def estimate_runtime(mode, cores):
    """Estimate runtime based on mode and resources"""
    estimates = {
        'test': {'base': 15, 'factor': 0.5},
        'pcr': {'base': 120, 'factor': 0.7},
        'panel': {'base': 180, 'factor': 0.8},
        'full': {'base': 240, 'factor': 0.9}
    }
    
    if mode in estimates:
        base_time = estimates[mode]['base']
        factor = estimates[mode]['factor']
        estimated = int(base_time * factor * (8 / cores))  # Scale by cores
        return estimated
    return 60  # Default

def estimate_disk_usage(mode):
    """Estimate disk usage for different modes"""
    usage = {
        'test': '2-5 GB',
        'pcr': '10-20 GB', 
        'panel': '15-30 GB',
        'full': '20-40 GB'
    }
    return usage.get(mode, '10-30 GB')

def run_snakemake(target, cores, extra_args=None):
    """Run snakemake with specified parameters"""
    cmd = ['snakemake', target, '--cores', str(cores)]
    
    if extra_args:
        cmd.extend(extra_args)
    
    print(f"🚀 Executing: {' '.join(cmd)}")
    print("=" * 80)
    
    try:
        result = subprocess.run(cmd, check=False)
        return result.returncode == 0
    except KeyboardInterrupt:
        print("\n⏹️  Pipeline interrupted by user")
        return False
    except Exception as e:
        print(f"❌ Error running pipeline: {e}")
        return False

def open_results(report_type='summary'):
    """Open results in web browser"""
    report_files = {
        'summary': 'results/reports/pipeline_summary.html',
        'dashboard': 'results/reports/dashboard.html',
        'pcr': 'results/reports/pcr_detailed_report.html',
        'panel': 'results/reports/panel_detailed_report.html',
        'validation': 'results/validation/validation_summary.html'
    }
    
    report_file = report_files.get(report_type, report_files['summary'])
    
    if os.path.exists(report_file):
        print(f"📊 Opening {report_type} report...")
        webbrowser.open(f'file://{os.path.abspath(report_file)}')
    else:
        print(f"❌ Report not found: {report_file}")
        print("Run the pipeline first to generate reports")

def show_status():
    """Show current pipeline status"""
    print("📊 Pipeline Status")
    print("=" * 50)
    
    # Check key output files
    outputs = {
        'PCR Results': 'results/design/pcr/filtered.tsv',
        'Panel Results': 'results/design/panel/filtered.tsv', 
        'Validation': 'results/validation/test_results.txt',
        'Main Report': 'results/reports/pipeline_summary.html'
    }
    
    for name, path in outputs.items():
        if os.path.exists(path):
            size = os.path.getsize(path)
            mtime = time.ctime(os.path.getmtime(path))
            print(f"✅ {name}: {size:,} bytes ({mtime})")
        else:
            print(f"⭕ {name}: Not found")
    
    # Check disk usage
    if os.path.exists('results'):
        try:
            result = subprocess.run(['du', '-sh', 'results'], capture_output=True, text=True)
            if result.returncode == 0:
                size = result.stdout.strip().split()[0]
                print(f"\n💾 Total results size: {size}")
        except:
            pass

def list_targets():
    """List available Snakemake targets"""
    print("🎯 Available Pipeline Targets")
    print("=" * 50)
    
    try:
        result = subprocess.run(['snakemake', '--list'], capture_output=True, text=True)
        if result.returncode == 0:
            print(result.stdout)
        else:
            print("Could not list targets")
    except Exception as e:
        print(f"Error listing targets: {e}")

def main():
    """Main launcher function"""
    parser = argparse.ArgumentParser(
        description='Oligo Design Pipeline Launcher',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s test                    # Run test mode (chr22 only)
  %(prog)s full --cores 16        # Run full pipeline with 16 cores
  %(prog)s pcr --cores 8          # Run PCR mode only
  %(prog)s status                 # Show pipeline status
  %(prog)s view                   # Open results in browser
        """
    )
    
    parser.add_argument('mode', choices=['test', 'pcr', 'panel', 'full', 'status', 'view', 'targets', 'clean'],
                       help='Pipeline mode to run')
    parser.add_argument('--cores', type=int, default=4,
                       help='Number of CPU cores to use (default: 4)')
    parser.add_argument('--dryrun', action='store_true',
                       help='Show what would be executed without running')
    parser.add_argument('--force', action='store_true',
                       help='Force re-execution of all rules')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose output')
    parser.add_argument('--report-type', choices=['summary', 'dashboard', 'pcr', 'panel', 'validation'],
                       default='summary', help='Type of report to open')
    
    args = parser.parse_args()
    
    print_banner()
    
    # Handle special modes
    if args.mode == 'status':
        show_status()
        return
    elif args.mode == 'view':
        open_results(args.report_type)
        return
    elif args.mode == 'targets':
        list_targets()
        return
    elif args.mode == 'clean':
        if input("⚠️  This will delete all results. Continue? (y/N): ").lower() == 'y':
            subprocess.run(['snakemake', 'clean_all'])
            print("🧹 Cleanup completed")
        return
    
    # Check dependencies and config
    if not check_dependencies():
        sys.exit(1)
    
    config = load_config()
    if not config:
        sys.exit(1)
    
    # Build snakemake command
    extra_args = []
    
    if args.dryrun:
        extra_args.append('--dryrun')
    if args.force:
        extra_args.append('--forceall')
    if args.verbose:
        extra_args.extend(['--verbose', '--printshellcmds'])
    
    # Determine target
    targets = {
        'test': 'test',
        'pcr': 'results/design/pcr/filtered.tsv',
        'panel': 'results/design/panel/filtered.tsv', 
        'full': 'all'
    }
    
    target = targets[args.mode]
    
    # Show execution info
    print(f"🎯 Target: {args.mode}")
    print(f"⚙️  Cores: {args.cores}")
    print(f"⏱️  Estimated runtime: {estimate_runtime(args.mode, args.cores)} minutes")
    print(f"💾 Estimated disk usage: {estimate_disk_usage(args.mode)}")
    
    if not args.dryrun:
        if input("\n▶️  Start pipeline? (Y/n): ").lower() not in ['', 'y', 'yes']:
            print("🛑 Pipeline cancelled")
            return
    
    print("\n" + "=" * 80)
    start_time = time.time()
    
    # Run pipeline
    success = run_snakemake(target, args.cores, extra_args)
    
    end_time = time.time()
    runtime = int((end_time - start_time) / 60)
    
    print("=" * 80)
    
    if success:
        print(f"✅ Pipeline completed successfully in {runtime} minutes!")
        
        if not args.dryrun and args.mode != 'test':
            if input("📊 Open results report? (Y/n): ").lower() in ['', 'y', 'yes']:
                open_results()
    else:
        print(f"❌ Pipeline failed after {runtime} minutes")
        print("Check log files in logs/ directory for details")
        sys.exit(1)

if __name__ == '__main__':
    main()
