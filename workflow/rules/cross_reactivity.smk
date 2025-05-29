#!/usr/bin/env python3
"""
Pathogen-Specific Oligo Design Engine

Specialized oligo design engine for pathogen detection applications
with enhanced specificity screening and taxonomic considerations.
"""

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
import tempfile
import subprocess
import os
import sys
import json
from typing import List, Dict, Tuple, Optional
import logging
from dataclasses import dataclass
from collections import defaultdict

@dataclass
class PathogenTarget:
    """Pathogen target information"""
    organism_name: str
    tax_id: Optional[int]
    target_gene: str
    sequence: str
    specificity_level: str  # 'species', 'genus', 'family'
    pathogen_type: str  # 'bacteria', 'virus', 'fungi', 'parasite'
    clinical_significance: str

@dataclass
class DesignParameters:
    """Design parameters for pathogen-specific oligos"""
    primer_opt_size: int = 22
    primer_min_size: int = 18
    primer_max_size: int = 28
    primer_opt_tm: float = 60.0
    primer_min_tm: float = 58.0
    primer_max_tm: float = 62.0
    primer_min_gc: float = 40.0
    primer_max_gc: float = 65.0
    product_size_range: str = "80-200"
    max_poly_x: int = 3
    max_ns_accepted: int = 0
    specificity_threshold: float = 95.0
    host_exclusion_threshold: float = 75.0

class PathogenOligoDesigner:
    """Main pathogen oligo design engine"""
    
    def __init__(self, config: Dict, mode: str = "species_specific"):
        self.config = config
        self.mode = mode
        self.logger = logging.getLogger(__name__)
        
        # Load design parameters based on mode
        self.design_params = self._load_design_parameters()
        
        # Initialize target database
        self.pathogen_targets = []
        self.exclusion_sequences = []
        
    def _load_design_parameters(self) -> DesignParameters:
        """Load design parameters from config"""
        params_config = self.config.get("primer3", {}).get(self.mode, {})
        
        return DesignParameters(
            primer_opt_size=params_config.get("primer_opt_size", 22),
            primer_min_size=params_config.get("primer_min_size", 18),
            primer_max_size=params_config.get("primer_max_size", 28),
            primer_opt_tm=params_config.get("primer_opt_tm", 60.0),
            primer_min_tm=params_config.get("primer_min_tm", 58.0),
            primer_max_tm=params_config.get("primer_max_tm", 62.0),
            primer_min_gc=params_config.get("primer_min_gc", 40.0),
            primer_max_gc=params_config.get("primer_max_gc", 65.0),
            product_size_range=params_config.get("product_size_range", "80-200"),
            max_poly_x=params_config.get("primer_max_poly_x", 3),
            max_ns_accepted=params_config.get("primer_max_ns_accepted", 0)
        )
    
    def load_pathogen_targets(self, target_file: str) -> None:
        """Load pathogen target sequences"""
        self.logger.info(f"Loading pathogen targets from {target_file}")
        
        target_count = 0
        for record in SeqIO.parse(target_file, "fasta"):
            # Parse header for pathogen information
            header_parts = record.description.split("|")
            
            if len(header_parts) >= 3:
                organism_name = header_parts[1] if len(header_parts) > 1 else "Unknown"
                target_gene = header_parts[2] if len(header_parts) > 2 else "Unknown"
                pathogen_type = header_parts[3] if len(header_parts) > 3 else "bacteria"
            else:
                organism_name = record.id
                target_gene = "Unknown"
                pathogen_type = "bacteria"
            
            target = PathogenTarget(
                organism_name=organism_name,
                tax_id=None,  # Could be parsed from header if available
                target_gene=target_gene,
                sequence=str(record.seq),
                specificity_level=self.mode.replace("_", " "),
                pathogen_type=pathogen_type,
                clinical_significance="diagnostic"
            )
            
            self.pathogen_targets.append(target)
            target_count += 1
        
        self.logger.info(f"Loaded {target_count} pathogen targets")
    
    def load_exclusion_sequences(self, exclusion_files: List[str]) -> None:
        """Load host/contaminant exclusion sequences"""
        self.logger.info(f"Loading exclusion sequences from {len(exclusion_files)} files")
        
        exclusion_count = 0
        for file_path in exclusion_files:
            if os.path.exists(file_path):
                for record in SeqIO.parse(file_path, "fasta"):
                    self.exclusion_sequences.append(str(record.seq))
                    exclusion_count += 1
        
        self.logger.info(f"Loaded {exclusion_count} exclusion sequences")
    
    def design_primers_for_target(self, target: PathogenTarget) -> List[Dict]:
        """Design primers for a specific pathogen target"""
        self.logger.debug(f"Designing primers for {target.organism_name}")
        
        candidates = []
        
        # Use sliding window approach for primer candidate generation
        sequence = target.sequence.upper()
        min_size = self.design_params.primer_min_size
        max_size = self.design_params.primer_max_size
        
        # Generate forward primer candidates
        for i in range(len(sequence) - min_size + 1):
            for length in range(min_size, min(max_size + 1, len(sequence) - i + 1)):
                primer_seq = sequence[i:i + length]
                
                if self._is_valid_primer_sequence(primer_seq):
                    primer_data = self._analyze_primer_sequence(primer_seq, target, "forward", i)
                    if primer_data:
                        candidates.append(primer_data)
        
        # For PCR mode, also generate reverse primers
        if self.mode == "species_specific" and len(candidates) > 0:
            # Generate reverse primers from the 3' end
            for i in range(len(sequence) - max_size, len(sequence) - min_size):
                for length in range(min_size, min(max_size + 1, len(sequence) - i + 1)):
                    if i + length <= len(sequence):
                        primer_seq = str(Seq(sequence[i:i + length]).reverse_complement())
                        
                        if self._is_valid_primer_sequence(primer_seq):
                            primer_data = self._analyze_primer_sequence(primer_seq, target, "reverse", i)
                            if primer_data:
                                candidates.append(primer_data)
        
        # Rank candidates by quality score
        candidates.sort(key=lambda x: x['quality_score'], reverse=True)
        
        # Return top candidates
        max_candidates = 20 if self.mode == "multi_pathogen" else 50
        return candidates[:max_candidates]
    
    def _is_valid_primer_sequence(self, sequence: str) -> bool:
        """Check if sequence meets basic primer criteria"""
        # Check for poly-nucleotide runs
        for nucleotide in 'ATGC':
            if nucleotide * (self.design_params.max_poly_x + 1) in sequence:
                return False
        
        # Check for ambiguous nucleotides
        if self.design_params.max_ns_accepted == 0 and 'N' in sequence:
            return False
        
        # Check GC content
        gc_content = GC(sequence)
        if not (self.design_params.primer_min_gc <= gc_content <= self.design_params.primer_max_gc):
            return False
        
        # Check for excessive secondary structure potential
        if self._has_excessive_secondary_structure(sequence):
            return False
        
        return True
    
    def _analyze_primer_sequence(self, sequence: str, target: PathogenTarget, 
                                primer_type: str, position: int) -> Optional[Dict]:
        """Analyze primer sequence and calculate quality metrics"""
        
        # Calculate melting temperature
        try:
            tm = mt.Tm_NN(sequence, strict=False)
        except:
            tm = mt.Tm_GC(sequence)  # Fallback method
        
        # Check Tm range
        if not (self.design_params.primer_min_tm <= tm <= self.design_params.primer_max_tm):
            return None
        
        # Calculate GC content
        gc_content = GC(sequence)
        
        # Calculate quality score based on multiple factors
        quality_score = self._calculate_quality_score(sequence, tm, gc_content)
        
        # Generate unique primer ID
        primer_id = f"{target.organism_name.replace(' ', '_')}_{target.target_gene}_{primer_type}_{position}"
        primer_id = primer_id.replace('|', '_').replace('/', '_')[:50]  # Limit length
        
        return {
            'primer_id': primer_id,
            'sequence': sequence,
            'type': primer_type,
            'length': len(sequence),
            'tm_calculated': round(tm, 2),
            'gc_content': round(gc_content, 2),
            'quality_score': round(quality_score, 3),
            'target_organism': target.organism_name,
            'target_gene': target.target_gene,
            'target_position': position,
            'pathogen_type': target.pathogen_type,
            'specificity_level': target.specificity_level,
            'clinical_significance': target.clinical_significance,
            'design_mode': self.mode
        }
    
    def _calculate_quality_score(self, sequence: str, tm: float, gc_content: float) -> float:
        """Calculate overall quality score for primer"""
        score = 1.0
        
        # Tm score (prefer values close to optimal)
        tm_optimal = self.design_params.primer_opt_tm
        tm_deviation = abs(tm - tm_optimal)
        tm_score = max(0, 1 - (tm_deviation / 5.0))  # Penalty for deviation > 5°C
        score *= tm_score
        
        # GC content score (prefer 45-55%)
        gc_optimal = 50.0
        gc_deviation = abs(gc_content - gc_optimal)
        gc_score = max(0, 1 - (gc_deviation / 25.0))  # Penalty for deviation > 25%
        score *= gc_score
        
        # Length score (prefer optimal size)
        length_optimal = self.design_params.primer_opt_size
        length_deviation = abs(len(sequence) - length_optimal)
        length_score = max(0, 1 - (length_deviation / 10.0))
        score *= length_score
        
        # Sequence complexity score
        complexity_score = self._calculate_sequence_complexity(sequence)
        score *= complexity_score
        
        return score
    
    def _calculate_sequence_complexity(self, sequence: str) -> float:
        """Calculate sequence complexity (avoid repeats and low complexity)"""
        # Penalize di-nucleotide repeats
        complexity = 1.0
        
        for i in range(len(sequence) - 1):
            dinucleotide = sequence[i:i+2]
            count = sequence.count(dinucleotide)
            if count > len(sequence) / 4:  # More than 25% of sequence
                complexity *= 0.8
        
        # Penalize mono-nucleotide runs
        for nucleotide in 'ATGC':
            max_run = 0
            current_run = 0
            for base in sequence:
                if base == nucleotide:
                    current_run += 1
                    max_run = max(max_run, current_run)
                else:
                    current_run = 0
            
            if max_run > 3:
                complexity *= (0.9 ** (max_run - 3))
        
        return complexity
    
    def _has_excessive_secondary_structure(self, sequence: str) -> bool:
        """Check for potential secondary structure issues"""
        # Simple hairpin detection
        min_stem = 4
        min_loop = 3
        
        for i in range(len(sequence) - min_stem - min_loop - min_stem):
            for j in range(i + min_stem + min_loop, len(sequence) - min_stem + 1):
                stem5 = sequence[i:i + min_stem]
                stem3 = sequence[j:j + min_stem]
                
                # Check for complementarity
                matches = 0
                for k in range(min_stem):
                    if self._is_complement(stem5[k], stem3[min_stem - 1 - k]):
                        matches += 1
                
                if matches >= min_stem - 1:  # Allow one mismatch
                    return True
        
        return False
    
    def _is_complement(self, base1: str, base2: str) -> bool:
        """Check if two bases are complementary"""
        complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return complements.get(base1) == base2
    
    def screen_for_specificity(self, candidates: List[Dict], 
                             blast_db: str, max_off_targets: int = 10) -> List[Dict]:
        """Screen candidates for specificity using BLAST"""
        self.logger.info(f"Screening {len(candidates)} candidates for specificity")
        
        # Create temporary FASTA file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_fasta:
            for candidate in candidates:
                tmp_fasta.write(f">{candidate['primer_id']}\n{candidate['sequence']}\n")
            tmp_fasta_path = tmp_fasta.name
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp_output:
            tmp_output_path = tmp_output.name
        
        try:
            # Run BLAST
            blast_cmd = [
                'blastn',
                '-query', tmp_fasta_path,
                '-db', blast_db,
                '-out', tmp_output_path,
                '-evalue', '1000',
                '-word_size', '7',
                '-max_target_seqs', '100',
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
                '-num_threads', '4'
            ]
            
            subprocess.run(blast_cmd, check=True, capture_output=True)
            
            # Parse BLAST results
            off_target_counts = self._parse_blast_results(tmp_output_path, candidates)
            
            # Add specificity information to candidates
            screened_candidates = []
            for candidate in candidates:
                primer_id = candidate['primer_id']
                off_target_count = off_target_counts.get(primer_id, 0)
                
                # Calculate specificity score
                specificity_score = max(0, 1 - (off_target_count / 100.0))
                
                candidate.update({
                    'off_target_count': off_target_count,
                    'specificity_score': round(specificity_score, 3)
                })
                
                # Only keep candidates with acceptable specificity
                if off_target_count <= max_off_targets:
                    screened_candidates.append(candidate)
            
            self.logger.info(f"Specificity screening: {len(screened_candidates)}/{len(candidates)} candidates passed")
            return screened_candidates
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"BLAST failed: {e}")
            return candidates  # Return unscreened candidates
        
        finally:
            # Cleanup temporary files
            if os.path.exists(tmp_fasta_path):
                os.unlink(tmp_fasta_path)
            if os.path.exists(tmp_output_path):
                os.unlink(tmp_output_path)
    
    def _parse_blast_results(self, blast_output: str, candidates: List[Dict]) -> Dict[str, int]:
        """Parse BLAST output and count off-targets"""
        off_target_counts = defaultdict(int)
        
        # Get target organism for each primer
        primer_targets = {c['primer_id']: c['target_organism'] for c in candidates}
        
        try:
            with open(blast_output, 'r') as f:
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) >= 12:
                        query_id = fields[0]
                        subject_id = fields[1]
                        identity = float(fields[2])
                        
                        # Skip self-hits and low-identity hits
                        if identity < 80:
                            continue
                        
                        # Check if this is a hit to the target organism
                        target_organism = primer_targets.get(query_id, '')
                        if target_organism and target_organism.lower() in subject_id.lower():
                            continue  # Skip hits to target organism
                        
                        # Count as off-target
                        off_target_counts[query_id] += 1
        
        except FileNotFoundError:
            self.logger.warning("BLAST output file not found")
        
        return dict(off_target_counts)
    
    def design_pathogen_panel(self, output_file: str) -> None:
        """Design complete pathogen detection panel"""
        self.logger.info("Starting pathogen panel design")
        
        all_candidates = []
        
        # Design primers for each target
        for target in self.pathogen_targets:
            target_candidates = self.design_primers_for_target(target)
            all_candidates.extend(target_candidates)
        
        self.logger.info(f"Generated {len(all_candidates)} initial candidates")
        
        # Create output DataFrame
        if all_candidates:
            candidates_df = pd.DataFrame(all_candidates)
            
            # Sort by quality score
            candidates_df = candidates_df.sort_values('quality_score', ascending=False)
            
            # Save to file
            candidates_df.to_csv(output_file, sep='\t', index=False)
            self.logger.info(f"Saved {len(candidates_df)} candidates to {output_file}")
        else:
            # Create empty file with headers
            empty_df = pd.DataFrame(columns=[
                'primer_id', 'sequence', 'type', 'length', 'tm_calculated', 
                'gc_content', 'quality_score', 'target_organism', 'target_gene',
                'pathogen_type', 'specificity_level', 'clinical_significance', 'design_mode'
            ])
            empty_df.to_csv(output_file, sep='\t', index=False)
            self.logger.warning("No valid candidates generated")

def main():
    """Main function for Snakemake integration"""
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Get parameters from Snakemake
    config = snakemake.config
    mode = snakemake.wildcards.mode
    
    input_fasta = snakemake.input.pathogen_fasta
    output_candidates = snakemake.output.candidates
    output_fasta = snakemake.output.candidates_fasta
    
    # Initialize designer
    designer = PathogenOligoDesigner(config, mode)
    
    # Load pathogen targets
    designer.load_pathogen_targets(input_fasta)
    
    # Design panel
    designer.design_pathogen_panel(output_candidates)
    
    # Create FASTA output
    if os.path.exists(output_candidates):
        df = pd.read_csv(output_candidates, sep='\t')
        if not df.empty:
            with open(output_fasta, 'w') as f:
                for _, row in df.iterrows():
                    f.write(f">{row['primer_id']}\n{row['sequence']}\n")
        else:
            # Create empty FASTA
            with open(output_fasta, 'w') as f:
                pass

if __name__ == "__main__":
    main()
