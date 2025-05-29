#!/usr/bin/env python3
"""
Extract FASTA sequences from reference genome based on ClinVar VCF file.
Extracts flanking sequences (default 150bp on each side) around variant positions.

Usage:
    python vcf_to_fasta.py --vcf variants.vcf --reference genome.fa --output targets.fa
    python vcf_to_fasta.py --vcf variants.vcf --reference genome.fa --output targets.fa --flank 200
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

def parse_vcf_line(line):
    """Parse a VCF line and extract relevant information"""
    if line.startswith('#'):
        return None
    
    fields = line.strip().split('\t')
    if len(fields) < 8:
        return None
    
    return {
        'chrom': fields[0],
        'pos': int(fields[1]),  # 1-based position
        'id': fields[2] if fields[2] != '.' else None,
        'ref': fields[3],
        'alt': fields[4],
        'qual': fields[5],
        'filter': fields[6],
        'info': fields[7]
    }

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

def extract_flanking_sequence(reference_dict, chrom, pos, ref_allele, flank_size=150):
    """
    Extract flanking sequence around a variant position
    
    Args:
        reference_dict: Dictionary of chromosome sequences
        chrom: Chromosome name
        pos: 1-based position
        ref_allele: Reference allele
        flank_size: Number of bases to extract on each side
    
    Returns:
        tuple: (sequence, actual_start, actual_end, success)
    """
    # Normalize chromosome name
    reference_chroms = list(reference_dict.keys())
    normalized_chrom = normalize_chromosome_name(chrom, reference_chroms)
    
    if not normalized_chrom:
        print(f"Warning: Chromosome '{chrom}' not found in reference genome")
        return None, None, None, False
    
    # Get chromosome sequence
    chrom_seq = reference_dict[normalized_chrom]
    chrom_length = len(chrom_seq)
    
    # Convert to 0-based coordinates
    start_pos = pos - 1
    
    # Verify reference allele matches (optional check)
    if len(ref_allele) == 1:  # SNV
        ref_in_genome = str(chrom_seq[start_pos]).upper()
        if ref_in_genome != ref_allele.upper():
            print(f"Warning: Reference allele mismatch at {chrom}:{pos}. "
                  f"VCF: {ref_allele}, Genome: {ref_in_genome}")
    
    # Calculate extraction boundaries
    extract_start = max(0, start_pos - flank_size)
    extract_end = min(chrom_length, start_pos + len(ref_allele) + flank_size)
    
    # Extract sequence
    extracted_seq = chrom_seq[extract_start:extract_end]
    
    return str(extracted_seq), extract_start + 1, extract_end, True  # Convert back to 1-based

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

def process_vcf_file(vcf_file, reference_dict, flank_size, output_file):
    """Process VCF file and extract flanking sequences"""
    
    sequences = []
    processed_count = 0
    success_count = 0
    
    print(f"Processing VCF file: {vcf_file}")
    print(f"Flanking region size: {flank_size} bp on each side")
    
    try:
        with open(vcf_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                variant = parse_vcf_line(line)
                if not variant:
                    continue
                
                processed_count += 1
                
                # Extract flanking sequence
                seq, start, end, success = extract_flanking_sequence(
                    reference_dict, 
                    variant['chrom'], 
                    variant['pos'], 
                    variant['ref'], 
                    flank_size
                )
                
                if success:
                    # Create sequence ID
                    var_id = variant['id'] if variant['id'] else f"var_{processed_count}"
                    seq_id = f"{var_id}_{variant['chrom']}_{variant['pos']}_{variant['ref']}_{variant['alt']}"
                    
                    # Create description
                    description = (f"chr={variant['chrom']} pos={variant['pos']} "
                                 f"ref={variant['ref']} alt={variant['alt']} "
                                 f"extracted_region={variant['chrom']}:{start}-{end} "
                                 f"flank_size={flank_size}bp")
                    
                    # Create SeqRecord
                    seq_record = SeqRecord(
                        Seq(seq),
                        id=seq_id,
                        description=description
                    )
                    
                    sequences.append(seq_record)
                    success_count += 1
                    
                    if success_count % 100 == 0:
                        print(f"  Processed {success_count} variants successfully...")
                
                else:
                    print(f"  Failed to extract sequence for variant at {variant['chrom']}:{variant['pos']}")
    
    except Exception as e:
        print(f"Error processing VCF file: {e}")
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
    print(f"  Total variants processed: {processed_count}")
    print(f"  Successful extractions: {success_count}")
    print(f"  Failed extractions: {processed_count - success_count}")
    
    return success_count

def main():
    parser = argparse.ArgumentParser(
        description="Extract FASTA sequences from reference genome based on ClinVar VCF file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with 150bp flanking regions
  python vcf_to_fasta.py --vcf clinvar.vcf --reference hg38.fa --output targets.fa
  
  # Custom flanking region size
  python vcf_to_fasta.py --vcf variants.vcf --reference genome.fa --output targets.fa --flank 200
  
  # Process only first 1000 variants
  python vcf_to_fasta.py --vcf large_file.vcf --reference genome.fa --output subset.fa --max-variants 1000
        """
    )
    
    parser.add_argument('--vcf', '-v', required=True,
                       help='Input VCF file (ClinVar or any standard VCF)')
    parser.add_argument('--reference', '-r', required=True,
                       help='Reference genome FASTA file')
    parser.add_argument('--output', '-o', required=True,
                       help='Output FASTA file')
    parser.add_argument('--flank', '-f', type=int, default=150,
                       help='Flanking region size in bp (default: 150)')
    parser.add_argument('--max-variants', type=int,
                       help='Maximum number of variants to process (for testing)')
    
    args = parser.parse_args()
    
    # Validate input files
    if not os.path.exists(args.vcf):
        print(f"Error: VCF file not found: {args.vcf}")
        sys.exit(1)
    
    if not os.path.exists(args.reference):
        print(f"Error: Reference genome file not found: {args.reference}")
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = Path(args.output).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load reference genome
    reference_dict = load_reference_genome(args.reference)
    
    # Process VCF file
    success_count = process_vcf_file(
        args.vcf, 
        reference_dict, 
        args.flank, 
        args.output
    )
    
    print(f"\nCompleted! Extracted {success_count} target sequences.")
    print(f"Output written to: {args.output}")

if __name__ == "__main__":
    main()


