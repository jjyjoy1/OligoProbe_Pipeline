#!/usr/bin/env python3
"""
Practical Examples for Oligo Database Queries

This script demonstrates real-world use cases for querying the oligo design
database using PostgreSQL + BioSQL schema with pg_trgm extensions.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from scripts.oligo_db_client import OligoDBClient
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict

class OligoDatabaseAnalyzer:
    """High-level analyzer for oligo database queries"""
    
    def __init__(self, client: OligoDBClient):
        self.client = client
    
    # =================================================================
    # Case 1: Quality Control and Validation
    # =================================================================
    
    def quality_control_report(self, project_name: str) -> Dict:
        """Generate comprehensive QC report for a project"""
        
        print(f"🔍 Quality Control Report for: {project_name}")
        print("=" * 60)
        
        # Basic statistics
        summary = self.client.get_project_summary(project_name)
        print(f"Total oligos: {summary.get('total_oligos', 0)}")
        print(f"Primer pairs: {summary.get('primer_pairs', 0)}")
        print(f"Average Tm: {summary.get('avg_tm', 0):.1f}°C")
        print(f"Average GC: {summary.get('avg_gc', 0):.1f}%")
        
        # Find outliers
        outliers = self.find_quality_outliers(project_name)
        print(f"\n⚠️  Quality outliers: {len(outliers)}")
        
        # Check for problematic sequences
        problematic = self.find_problematic_sequences(project_name)
        print(f"🚨 Problematic sequences: {len(problematic)}")
        
        # Validation status distribution
        validation_dist = self.get_validation_distribution(project_name)
        print(f"\n📊 Validation Status:")
        for status, count in validation_dist.items():
            print(f"  {status}: {count}")
        
        return {
            'summary': summary,
            'outliers': outliers,
            'problematic': problematic,
            'validation_distribution': validation_dist
        }
    
    def find_quality_outliers(self, project_name: str) -> pd.DataFrame:
        """Find oligos with unusual quality metrics"""
        
        query = """
        WITH project_stats AS (
            SELECT 
                AVG(tm_calculated) as mean_tm,
                STDDEV(tm_calculated) as stddev_tm,
                AVG(gc_content) as mean_gc,
                STDDEV(gc_content) as stddev_gc
            FROM v_oligo_complete
            WHERE project_name = %s AND tm_calculated IS NOT NULL
        )
        SELECT 
            oligo_name,
            sequence,
            tm_calculated,
            gc_content,
            ABS(tm_calculated - mean_tm) / NULLIF(stddev_tm, 0) as tm_z_score,
            ABS(gc_content - mean_gc) / NULLIF(stddev_gc, 0) as gc_z_score
        FROM v_oligo_complete, project_stats
        WHERE project_name = %s
          AND (ABS(tm_calculated - mean_tm) / NULLIF(stddev_tm, 0) > 2
               OR ABS(gc_content - mean_gc) / NULLIF(stddev_gc, 0) > 2)
        ORDER BY GREATEST(
            ABS(tm_calculated - mean_tm) / NULLIF(stddev_tm, 0),
            ABS(gc_content - mean_gc) / NULLIF(stddev_gc, 0)
        ) DESC
        """
        
        return self.client.execute_query_df(query, (project_name, project_name))
    
    def find_problematic_sequences(self, project_name: str) -> pd.DataFrame:
        """Find sequences with known problematic patterns"""
        
        query = """
        SELECT 
            oligo_name,
            sequence,
            CASE 
                WHEN sequence ~ 'GGGG+' THEN 'Poly-G run'
                WHEN sequence ~ 'TTTT+' THEN 'Poly-T run'
                WHEN sequence ~ 'AAAA+' THEN 'Poly-A run'
                WHEN sequence ~ 'CCCC+' THEN 'Poly-C run'
                WHEN LENGTH(sequence) - LENGTH(REPLACE(sequence, 'CG', '')) > LENGTH(sequence) * 0.3 THEN 'High CpG'
                WHEN gc_content < 20 OR gc_content > 80 THEN 'Extreme GC'
                ELSE 'Other'
            END as issue_type
        FROM v_oligo_complete
        WHERE project_name = %s
          AND (sequence ~ '[ATGC]{4,}'  -- 4+ consecutive same bases
               OR gc_content < 20 OR gc_content > 80)
        ORDER BY issue_type, oligo_name
        """
        
        return self.client.execute_query_df(query, (project_name,))
    
    def get_validation_distribution(self, project_name: str) -> Dict:
        """Get validation status distribution"""
        
        query = """
        SELECT validation_status, COUNT(*) as count
        FROM v_oligo_complete
        WHERE project_name = %s
        GROUP BY validation_status
        """
        
        df = self.client.execute_query_df(query, (project_name,))
        return dict(zip(df['validation_status'], df['count']))
    
    # =================================================================
    # Case 2: Cross-Reactivity and Specificity Analysis
    # =================================================================
    
    def cross_reactivity_analysis(self, project_name: str, 
                                 similarity_threshold: float = 0.8) -> Dict:
        """Comprehensive cross-reactivity analysis"""
        
        print(f"🔍 Cross-Reactivity Analysis for: {project_name}")
        print("=" * 60)
        
        # Find similar primer sequences
        similar_pairs = self.find_similar_primers(project_name, similarity_threshold)
        print(f"Potentially cross-reactive pairs: {len(similar_pairs)}")
        
        # Check for primer-dimer potential
        dimer_risks = self.check_primer_dimer_risk(project_name)
        print(f"Primer-dimer risks: {len(dimer_risks)}")
        
        # Off-target analysis summary
        off_target_summary = self.summarize_off_targets(project_name)
        print(f"Oligos with high off-targets: {off_target_summary['high_off_target_count']}")
        
        return {
            'similar_pairs': similar_pairs,
            'dimer_risks': dimer_risks,
            'off_target_summary': off_target_summary
        }
    
    def find_similar_primers(self, project_name: str, 
                           similarity_threshold: float = 0.8) -> pd.DataFrame:
        """Find primers with high sequence similarity"""
        
        query = """
        SELECT 
            o1.oligo_name as primer1,
            o2.oligo_name as primer2,
            o1.sequence as seq1,
            o2.sequence as seq2,
            SIMILARITY(o1.sequence, o2.sequence) as similarity,
            o1.oligo_type as type1,
            o2.oligo_type as type2,
            o1.target_region as target1,
            o2.target_region as target2
        FROM v_oligo_complete o1
        JOIN v_oligo_complete o2 ON o1.oligo_design_id < o2.oligo_design_id
        WHERE o1.project_name = %s 
          AND o2.project_name = %s
          AND o1.sequence %% o2.sequence
          AND SIMILARITY(o1.sequence, o2.sequence) >= %s
        ORDER BY similarity DESC
        """
        
        return self.client.execute_query_df(query, (project_name, project_name, similarity_threshold))
    
    def check_primer_dimer_risk(self, project_name: str, 
                               min_overlap: int = 6) -> pd.DataFrame:
        """Check for potential primer-dimer formation"""
        
        query = """
        SELECT 
            o1.oligo_name as primer1,
            o2.oligo_name as primer2,
            o1.sequence as seq1,
            reverse_complement(o2.sequence) as rev_comp2,
            SIMILARITY(o1.sequence, reverse_complement(o2.sequence)) as rc_similarity
        FROM v_oligo_complete o1
        JOIN v_oligo_complete o2 ON o1.oligo_design_id != o2.oligo_design_id
        WHERE o1.project_name = %s 
          AND o2.project_name = %s
          AND o1.sequence %% reverse_complement(o2.sequence)
          AND SIMILARITY(o1.sequence, reverse_complement(o2.sequence)) > 0.6
        ORDER BY rc_similarity DESC
        """
        
        return self.client.execute_query_df(query, (project_name, project_name))
    
    def summarize_off_targets(self, project_name: str) -> Dict:
        """Summarize off-target analysis"""
        
        query = """
        SELECT 
            COUNT(*) as total_oligos,
            COUNT(*) FILTER (WHERE total_off_targets > 10) as high_off_target_count,
            AVG(total_off_targets) as avg_off_targets,
            MAX(total_off_targets) as max_off_targets,
            COUNT(*) FILTER (WHERE high_similarity_hits > 0) as oligos_with_high_sim_hits
        FROM v_oligo_with_off_targets owt
        JOIN v_oligo_complete oc ON owt.oligo_design_id = oc.oligo_design_id
        WHERE oc.project_name = %s
        """
        
        result = self.client.execute_query(query, (project_name,), fetch_all=False)
        columns = ['total_oligos', 'high_off_target_count', 'avg_off_targets', 
                  'max_off_targets', 'oligos_with_high_sim_hits']
        return dict(zip(columns, result))
    
    # =================================================================
    # Case 3: Sequence Pattern Mining and Discovery
    # =================================================================
    
    def sequence_pattern_mining(self, project_name: str) -> Dict:
        """Mine sequences for common patterns and motifs"""
        
        print(f"🔍 Sequence Pattern Mining for: {project_name}")
        print("=" * 60)
        
        # Common motif analysis
        motifs = self.analyze_common_motifs(project_name)
        print(f"Analyzed {len(motifs)} motif patterns")
        
        # Success pattern analysis
        success_patterns = self.find_success_patterns(project_name)
        print(f"Identified patterns in successful oligos")
        
        # Length distribution
        length_dist = self.analyze_length_distribution(project_name)
        print(f"Length distribution analyzed")
        
        return {
            'motifs': motifs,
            'success_patterns': success_patterns,
            'length_distribution': length_dist
        }
    
    def analyze_common_motifs(self, project_name: str) -> pd.DataFrame:
        """Analyze frequency of common sequence motifs"""
        
        motifs = ['CG', 'GC', 'AT', 'TA', 'GGG', 'CCC', 'AAA', 'TTT']
        results = []
        
        for motif in motifs:
            query = f"""
            SELECT 
                '{motif}' as motif,
                COUNT(*) as total_sequences,
                COUNT(*) FILTER (WHERE sequence LIKE '%{motif}%') as contains_motif,
                ROUND(
                    COUNT(*) FILTER (WHERE sequence LIKE '%{motif}%') * 100.0 / COUNT(*), 
                    2
                ) as percentage
            FROM v_oligo_complete
            WHERE project_name = %s
            """
            
            df = self.client.execute_query_df(query, (project_name,))
            results.append(df)
        
        return pd.concat(results, ignore_index=True)
    
    def find_success_patterns(self, project_name: str) -> Dict:
        """Find patterns associated with successful oligos"""
        
        query = """
        SELECT 
            validation_status,
            COUNT(*) as count,
            AVG(LENGTH(sequence)) as avg_length,
            AVG(gc_content) as avg_gc,
            AVG(tm_calculated) as avg_tm,
            -- Starting nucleotides
            COUNT(*) FILTER (WHERE sequence LIKE 'G%') as starts_with_g,
            COUNT(*) FILTER (WHERE sequence LIKE 'C%') as starts_with_c,
            COUNT(*) FILTER (WHERE sequence LIKE 'A%') as starts_with_a,
            COUNT(*) FILTER (WHERE sequence LIKE 'T%') as starts_with_t
        FROM v_oligo_complete
        WHERE project_name = %s
        GROUP BY validation_status
        """
        
        df = self.client.execute_query_df(query, (project_name,))
        return df.set_index('validation_status').to_dict('index')
    
    def analyze_length_distribution(self, project_name: str) -> pd.DataFrame:
        """Analyze sequence length distribution by type"""
        
        query = """
        SELECT 
            oligo_type,
            MIN(length) as min_length,
            MAX(length) as max_length,
            AVG(length) as avg_length,
            PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY length) as median_length,
            STDDEV(length) as stddev_length
        FROM v_oligo_complete
        WHERE project_name = %s
        GROUP BY oligo_type
        """
        
        return self.client.execute_query_df(query, (project_name,))
    
    # =================================================================
    # Case 4: Contamination Detection
    # =================================================================
    
    def contamination_screening(self, project_name: str) -> Dict:
        """Screen for potential contamination sequences"""
        
        print(f"🔍 Contamination Screening for: {project_name}")
        print("=" * 60)
        
        # Common contaminant sequences
        contaminants = {
            'E.coli_16S': 'AGAGTTTGATCCTGGCTCAG',
            'Human_18S': 'GTAACCCGTTGAACCCCATT',
            'Mycoplasma': 'GGGAGCAAACAGGATTAGATACCCT',
            'PhiX174': 'GAGTTTTATCGCTTCCATGACGCAG'
        }
        
        contamination_hits = []
        
        for name, seq in contaminants.items():
            similar = self.client.find_similar_sequences(seq, 0.8, 10)
            project_similar = similar[similar.get('project_name', '') == project_name] if 'project_name' in similar.columns else similar
            
            if not project_similar.empty:
                contamination_hits.append({
                    'contaminant': name,
                    'sequence': seq,
                    'matches': len(project_similar),
                    'best_match': project_similar.iloc[0].to_dict() if not project_similar.empty else None
                })
        
        print(f"Potential contamination sources found: {len(contamination_hits)}")
        
        return {'contamination_hits': contamination_hits}
    
    # =================================================================
    # Case 5: Primer Design Optimization
    # =================================================================
    
    def optimize_primer_pairs(self, project_name: str) -> pd.DataFrame:
        """Find optimal primer pairs based on multiple criteria"""
        
        print(f"🔍 Primer Pair Optimization for: {project_name}")
        print("=" * 60)
        
        query = """
        WITH pair_scores AS (
            SELECT 
                pp.*,
                -- Tm compatibility score (closer Tm is better)
                1.0 - (ABS(forward_tm - reverse_tm) / 10.0) as tm_compatibility,
                -- GC balance score
                1.0 - (ABS(forward_gc - reverse_gc) / 50.0) as gc_balance,
                -- Product size score (prefer 100-300bp)
                CASE 
                    WHEN product_size BETWEEN 100 AND 300 THEN 1.0
                    WHEN product_size BETWEEN 80 AND 400 THEN 0.8
                    ELSE 0.5
                END as size_score,
                -- Overall specificity (average of both primers)
                (COALESCE(f_specificity.specificity_score, 0) + 
                 COALESCE(r_specificity.specificity_score, 0)) / 2.0 as avg_specificity
            FROM v_primer_pairs pp
            LEFT JOIN v_oligo_complete f_specificity ON pp.forward_name = f_specificity.oligo_name
            LEFT JOIN v_oligo_complete r_specificity ON pp.reverse_name = r_specificity.oligo_name
            WHERE f_specificity.project_name = %s
        )
        SELECT 
            *,
            -- Combined optimization score
            (tm_compatibility * 0.3 + gc_balance * 0.2 + size_score * 0.3 + avg_specificity * 0.2) as optimization_score
        FROM pair_scores
        WHERE tm_compatibility > 0.7  -- Basic quality filter
        ORDER BY optimization_score DESC
        LIMIT 20
        """
        
        optimized = self.client.execute_query_df(query, (project_name,))
        print(f"Top optimized primer pairs: {len(optimized)}")
        
        return optimized
    
    # =================================================================
    # Case 6: Performance Analysis and Monitoring
    # =================================================================
    
    def performance_analysis(self) -> Dict:
        """Analyze database performance and query patterns"""
        
        print("🔍 Database Performance Analysis")
        print("=" * 60)
        
        # Index usage statistics
        index_stats = self.get_index_usage_stats()
        print(f"Analyzed {len(index_stats)} indexes")
        
        # Query performance
        query_perf = self.get_query_performance()
        print(f"Query performance data available: {len(query_perf) > 0}")
        
        # Table sizes
        table_sizes = self.get_table_sizes()
        print(f"Analyzed {len(table_sizes)} tables")
        
        return {
            'index_stats': index_stats,
            'query_performance': query_perf,
            'table_sizes': table_sizes
        }
    
    def get_index_usage_stats(self) -> pd.DataFrame:
        """Get index usage statistics"""
        
        query = """
        SELECT 
            schemaname,
            tablename,
            indexname,
            idx_scan,
            idx_tup_read,
            idx_tup_fetch,
            CASE WHEN idx_scan > 0 THEN idx_tup_read::float / idx_scan ELSE 0 END as avg_tuples_per_scan
        FROM pg_stat_user_indexes
        WHERE schemaname = 'public'
        ORDER BY idx_scan DESC
        """
        
        return self.client.execute_query_df(query)
    
    def get_query_performance(self) -> pd.DataFrame:
        """Get query performance statistics if available"""
        
        # Check if query_performance table exists
        check_query = """
        SELECT EXISTS (
            SELECT FROM information_schema.tables 
            WHERE table_schema = 'public' 
            AND table_name = 'query_performance'
        )
        """
        
        exists = self.client.execute_query(check_query, fetch_all=False)[0]
        
        if exists:
            query = """
            SELECT 
                query_type,
                COUNT(*) as execution_count,
                AVG(execution_time_ms) as avg_time_ms,
                MIN(execution_time_ms) as min_time_ms,
                MAX(execution_time_ms) as max_time_ms,
                AVG(rows_returned) as avg_rows
            FROM query_performance
            WHERE executed_at >= CURRENT_DATE - INTERVAL '7 days'
            GROUP BY query_type
            ORDER BY avg_time_ms DESC
            """
            return self.client.execute_query_df(query)
        
        return pd.DataFrame()
    
    def get_table_sizes(self) -> pd.DataFrame:
        """Get table size information"""
        
        query = """
        SELECT 
            schemaname,
            tablename,
            pg_size_pretty(pg_total_relation_size(schemaname||'.'||tablename)) as size,
            pg_total_relation_size(schemaname||'.'||tablename) as size_bytes
        FROM pg_tables
        WHERE schemaname = 'public'
        ORDER BY pg_total_relation_size(schemaname||'.'||tablename) DESC
        """
        
        return self.client.execute_query_df(query)


# =================================================================
# Example Usage and Demonstrations
# =================================================================

def run_comprehensive_analysis(project_name: str):
    """Run a comprehensive analysis of a project"""
    
    # Initialize client and analyzer
    client = OligoDBClient()
    analyzer = OligoDatabaseAnalyzer(client)
    
    print(f"🧬 Comprehensive Oligo Analysis for: {project_name}")
    print("=" * 80)
    
    try:
        # 1. Quality Control
        qc_report = analyzer.quality_control_report(project_name)
        
        # 2. Cross-Reactivity Analysis
        xr_analysis = analyzer.cross_reactivity_analysis(project_name)
        
        # 3. Pattern Mining
        pattern_analysis = analyzer.sequence_pattern_mining(project_name)
        
        # 4. Contamination Screening
        contamination = analyzer.contamination_screening(project_name)
        
        # 5. Primer Optimization
        optimized_pairs = analyzer.optimize_primer_pairs(project_name)
        
        # 6. Performance Analysis
        performance = analyzer.performance_analysis()
        
        print("\n" + "=" * 80)
        print("✅ Analysis Complete!")
        print(f"   - Quality outliers: {len(qc_report['outliers'])}")
        print(f"   - Cross-reactive pairs: {len(xr_analysis['similar_pairs'])}")
        print(f"   - Optimized primer pairs: {len(optimized_pairs)}")
        print(f"   - Contamination hits: {len(contamination['contamination_hits'])}")
        
        return {
            'qc': qc_report,
            'cross_reactivity': xr_analysis,
            'patterns': pattern_analysis,
            'contamination': contamination,
            'optimized_pairs': optimized_pairs,
            'performance': performance
        }
        
    except Exception as e:
        print(f"❌ Analysis failed: {e}")
        return None

def demonstrate_similarity_searches():
    """Demonstrate various similarity search capabilities"""
    
    client = OligoDBClient()
    
    print("🔍 Similarity Search Demonstrations")
    print("=" * 50)
    
    # Example sequence
    test_sequence = "GCATACGTTGTATCCGGGCAT"
    
    # 1. Basic similarity search
    print(f"1. Basic similarity search for: {test_sequence}")
    similar = client.find_similar_sequences(test_sequence, 0.7)
    print(f"   Found {len(similar)} similar sequences")
    
    # 2. Cross-reactivity check
    print(f"\n2. Cross-reactivity analysis")
    cross_reactive = client.find_cross_reactive_primers(0.8)
    print(f"   Found {len(cross_reactive)} potentially cross-reactive pairs")
    
    # 3. Novelty check
    print(f"\n3. Sequence novelty check")
    novelty = client.check_sequence_novelty(test_sequence, 0.9)
    print(f"   Sequence is novel: {novelty['is_novel']}")
    print(f"   Most similar match: {novelty['max_similarity']:.3f}")

def demonstrate_custom_queries():
    """Demonstrate custom query capabilities"""
    
    client = OligoDBClient()
    
    print("🔍 Custom Query Demonstrations")
    print("=" * 50)
    
    # 1. Complex quality analysis
    quality_query = """
    WITH quality_metrics AS (
        SELECT 
            oligo_type,
            COUNT(*) as total,
            AVG(tm_calculated) as avg_tm,
            AVG(gc_content) as avg_gc,
            COUNT(*) FILTER (WHERE off_target_count = 0) as no_off_targets,
            COUNT(*) FILTER (WHERE validation_status = 'validated') as validated
        FROM v_oligo_complete
        GROUP BY oligo_type
    )
    SELECT 
        *,
        ROUND(validated::float / total * 100, 2) as validation_rate,
        ROUND(no_off_targets::float / total * 100, 2) as specificity_rate
    FROM quality_metrics
    ORDER BY total DESC
    """
    
    quality_results = client.execute_query_df(quality_query)
    print("1. Quality metrics by oligo type:")
    print(quality_results.to_string(index=False))
    
    # 2. Advanced similarity with thresholds
    similarity_query = """
    WITH similarity_matrix AS (
        SELECT 
            o1.oligo_name as oligo1,
            o2.oligo_name as oligo2,
            SIMILARITY(o1.sequence, o2.sequence) as similarity,
            o1.oligo_type as type1,
            o2.oligo_type as type2
        FROM v_oligo_complete o1
        JOIN v_oligo_complete o2 ON o1.oligo_design_id < o2.oligo_design_id
        WHERE o1.sequence % o2.sequence
        AND SIMILARITY(o1.sequence, o2.sequence) > 0.8
    )
    SELECT 
        type1,
        type2,
        COUNT(*) as similar_pairs,
        AVG(similarity) as avg_similarity,
        MAX(similarity) as max_similarity
    FROM similarity_matrix
    GROUP BY type1, type2
    ORDER BY similar_pairs DESC
    """
    
    similarity_results = client.execute_query_df(similarity_query)
    print(f"\n2. Similarity analysis by oligo type ({len(similarity_results)} combinations):")
    if not similarity_results.empty:
        print(similarity_results.to_string(index=False))

if __name__ == "__main__":
    # Set up logging
    import logging
    logging.basicConfig(level=logging.INFO)
    
    # Example project name (adjust as needed)
    project_name = "OligoDesign_hg38_2025"
    
    print("🧬 Oligo Database Analysis Examples")
    print("=" * 80)
    
    try:
        # Run comprehensive analysis
        results = run_comprehensive_analysis(project_name)
        
        print("\n" + "=" * 80)
        
        # Demonstrate similarity searches
        demonstrate_similarity_searches()
        
        print("\n" + "=" * 80)
        
        # Demonstrate custom queries
        demonstrate_custom_queries()
        
    except Exception as e:
        print(f"❌ Error running examples: {e}")
        print("Make sure the database is set up and contains data")
