#!/usr/bin/env python3
"""
Extract FASTA sequences from reference genome based on BED file intervals.
Supports standard BED format with optional strand information.

Usage:
    python bed_to_fasta.py --bed regions.bed --reference genome.fa --output sequences.fa
    python bed_to_fasta.py --bed targets.bed --reference hg38.fa --output probes.fa --extend 50
"""

import argparse
import sys
import os
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("Error: BioPython is required. Install with: pip install biopython")
    sys.exit(1)

def parse_bed_line(line, line_num):
    """
    Parse a BED file line and extract interval information
    
    BED format: chrom start end [name] [score] [strand] [other fields...]
    Note: BED uses 0-based start, 1-based end coordinates
    """
    if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
        return None
    
    fields = line.strip().split('\t')
    if len(fields) < 3:
        print(f"Warning: Line {line_num} has fewer than 3 required BED fields, skipping")
        return None
    
    try:
        interval = {
            'chrom': fields[0],
            'start': int(fields[1]),  # 0-based
            'end': int(fields[2]),    # 1-based (exclusive)
            'name': fields[3] if len(fields) > 3 and fields[3] != '.' else f"region_{line_num}",
            'score': fields[4] if len(fields) > 4 and fields[4] != '.' else '0',
            'strand': fields[5] if len(fields) > 5 and fields[5] in ['+', '-'] else '+',
            'line_num': line_num
        }
        
        # Validate coordinates
        if interval['start'] >= interval['end']:
            print(f"Warning: Invalid coordinates at line {line_num}: start >= end")
            return None
            
        if interval['start'] < 0:
            print(f"Warning: Negative start coordinate at line {line_num}")
            return None
        
        return interval
        
    except ValueError as e:
        print(f"Warning: Error parsing coordinates at line {line_num}: {e}")
        return None

def normalize_chromosome_name(chrom, reference_chroms):
    """
    Normalize chromosome names to match reference genome
    Handles chr1 vs 1, chrM vs chrMT, etc.
    """
    # Try exact match first
    if chrom in reference_chroms:
        return chrom
    
    # Try adding 'chr' prefix
    if f"chr{chrom}" in reference_chroms:
        return f"chr{chrom}"
    
    # Try removing 'chr' prefix
    if chrom.startswith('chr') and chrom[3:] in reference_chroms:
        return chrom[3:]
    
    # Handle mitochondrial chromosome variations
    mt_variants = {
        'M': ['chrM', 'chrMT', 'MT'],
        'chrM': ['M', 'chrMT', 'MT'], 
        'chrMT': ['M', 'chrM', 'MT'],
        'MT': ['M', 'chrM', 'chrMT']
    }
    
    if chrom in mt_variants:
        for variant in mt_variants[chrom]:
            if variant in reference_chroms:
                return variant
    
    return None

def extract_bed_sequence(reference_dict, interval, extend=0):
    """
    Extract sequence for a BED interval with optional extension
    
    Args:
        reference_dict: Dictionary of chromosome sequences
        interval: BED interval dictionary
        extend: Number of bases to extend on both sides
    
    Returns:
        tuple: (sequence, actual_start, actual_end, success)
    """
    # Normalize chromosome name
    reference_chroms = list(reference_dict.keys())
    normalized_chrom = normalize_chromosome_name(interval['chrom'], reference_chroms)
    
    if not normalized_chrom:
        print(f"Warning: Chromosome '{interval['chrom']}' not found in reference genome")
        return None, None, None, False
    
    # Get chromosome sequence
    chrom_seq = reference_dict[normalized_chrom]
    chrom_length = len(chrom_seq)
    
    # Calculate extraction boundaries with extension
    extract_start = max(0, interval['start'] - extend)
    extract_end = min(chrom_length, interval['end'] + extend)
    
    # Extract sequence
    extracted_seq = chrom_seq[extract_start:extract_end]
    
    # Handle strand
    if interval['strand'] == '-':
        extracted_seq = extracted_seq.reverse_complement()
    
    return str(extracted_seq), extract_start, extract_end, True

def load_reference_genome(reference_file):
    """Load reference genome into memory"""
    print(f"Loading reference genome: {reference_file}")
    
    reference_dict = {}
    
    try:
        for record in SeqIO.parse(reference_file, "fasta"):
            reference_dict[record.id] = record.seq
            print(f"  Loaded chromosome: {record.id} ({len(record.seq):,} bp)")
    except Exception as e:
        print(f"Error loading reference genome: {e}")
        sys.exit(1)
    
    print(f"Loaded {len(reference_dict)} chromosomes/contigs")
    return reference_dict

def process_bed_file(bed_file, reference_dict, extend, output_file, max_intervals=None):
    """Process BED file and extract sequences for each interval"""
    
    sequences = []
    processed_count = 0
    success_count = 0
    
    print(f"Processing BED file: {bed_file}")
    if extend > 0:
        print(f"Extending intervals by {extend} bp on each side")
    
    try:
        with open(bed_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.strip() == '':
                    continue
                
                interval = parse_bed_line(line, line_num)
                if not interval:
                    continue
                
                processed_count += 1
                
                # Check max intervals limit
                if max_intervals and processed_count > max_intervals:
                    print(f"Reached maximum interval limit ({max_intervals}), stopping...")
                    break
                
                # Extract sequence
                seq, actual_start, actual_end, success = extract_bed_sequence(
                    reference_dict, 
                    interval, 
                    extend
                )
                
                if success:
                    # Create sequence ID
                    original_coords = f"{interval['chrom']}:{interval['start']}-{interval['end']}"
                    if extend > 0:
                        extended_coords = f"{interval['chrom']}:{actual_start}-{actual_end}"
                        seq_id = f"{interval['name']}_{original_coords}_ext{extend}bp"
                        description = (f"name={interval['name']} "
                                     f"original={original_coords} "
                                     f"extended={extended_coords} "
                                     f"strand={interval['strand']} "
                                     f"length={len(seq)}bp")
                    else:
                        seq_id = f"{interval['name']}_{original_coords}"
                        description = (f"name={interval['name']} "
                                     f"coords={original_coords} "
                                     f"strand={interval['strand']} "
                                     f"length={len(seq)}bp")
                    
                    # Create SeqRecord
                    seq_record = SeqRecord(
                        Seq(seq),
                        id=seq_id,
                        description=description
                    )
                    
                    sequences.append(seq_record)
                    success_count += 1
                    
                    if success_count % 100 == 0:
                        print(f"  Processed {success_count} intervals successfully...")
                
                else:
                    print(f"  Failed to extract sequence for interval at line {line_num}")
    
    except Exception as e:
        print(f"Error processing BED file: {e}")
        sys.exit(1)
    
    # Write output FASTA
    print(f"Writing {len(sequences)} sequences to: {output_file}")
    
    try:
        with open(output_file, 'w') as f:
            SeqIO.write(sequences, f, "fasta")
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)
    
    print(f"Summary:")
    print(f"  Total intervals processed: {processed_count}")
    print(f"  Successful extractions: {success_count}")
    print(f"  Failed extractions: {processed_count - success_count}")
    
    return success_count

def validate_bed_file(bed_file):
    """Basic validation of BED file format"""
    print(f"Validating BED file format...")
    
    line_count = 0
    valid_lines = 0
    
    try:
        with open(bed_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if line == '' or line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                    continue
                
                line_count += 1
                fields = line.split('\t')
                
                if len(fields) >= 3:
                    try:
                        start = int(fields[1])
                        end = int(fields[2])
                        if start >= 0 and end > start:
                            valid_lines += 1
                    except ValueError:
                        print(f"  Warning: Invalid coordinates at line {line_num}")
                else:
                    print(f"  Warning: Insufficient fields at line {line_num}")
    
    except Exception as e:
        print(f"Error validating BED file: {e}")
        return False
    
    print(f"  Total data lines: {line_count}")
    print(f"  Valid intervals: {valid_lines}")
    print(f"  Invalid intervals: {line_count - valid_lines}")
    
    return valid_lines > 0

def main():
    parser = argparse.ArgumentParser(
        description="Extract FASTA sequences from reference genome based on BED file intervals",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
BED File Format:
  Column 1: Chromosome name
  Column 2: Start position (0-based)
  Column 3: End position (1-based, exclusive)
  Column 4: Name (optional)
  Column 5: Score (optional)
  Column 6: Strand (optional, + or -)

Examples:
  # Basic usage
  python bed_to_fasta.py --bed regions.bed --reference hg38.fa --output sequences.fa
  
  # Extend intervals by 50bp on each side
  python bed_to_fasta.py --bed targets.bed --reference genome.fa --output extended.fa --extend 50
  
  # Process only first 500 intervals
  python bed_to_fasta.py --bed large.bed --reference genome.fa --output subset.fa --max-intervals 500
  
  # Validate BED file format
  python bed_to_fasta.py --bed regions.bed --reference genome.fa --output out.fa --validate
        """
    )
    
    parser.add_argument('--bed', '-b', required=True,
                       help='Input BED file with genomic intervals')
    parser.add_argument('--reference', '-r', required=True,
                       help='Reference genome FASTA file')
    parser.add_argument('--output', '-o', required=True,
                       help='Output FASTA file')
    parser.add_argument('--extend', '-e', type=int, default=0,
                       help='Extend intervals by N bp on both sides (default: 0)')
    parser.add_argument('--max-intervals', type=int,
                       help='Maximum number of intervals to process (for testing)')
    parser.add_argument('--validate', action='store_true',
                       help='Validate BED file format before processing')
    
    args = parser.parse_args()
    
    # Validate input files
    if not os.path.exists(args.bed):
        print(f"Error: BED file not found: {args.bed}")
        sys.exit(1)
    
    if not os.path.exists(args.reference):
        print(f"Error: Reference genome file not found: {args.reference}")
        sys.exit(1)
    
    # Validate BED file if requested
    if args.validate:
        if not validate_bed_file(args.bed):
            print("BED file validation failed. Please check the format.")
            sys.exit(1)
    
    # Create output directory if needed
    output_dir = Path(args.output).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load reference genome
    reference_dict = load_reference_genome(args.reference)
    
    # Process BED file
    success_count = process_bed_file(
        args.bed, 
        reference_dict, 
        args.extend, 
        args.output,
        args.max_intervals
    )
    
    print(f"\nCompleted! Extracted {success_count} sequences from BED intervals.")
    print(f"Output written to: {args.output}")

if __name__ == "__main__":
    main()

