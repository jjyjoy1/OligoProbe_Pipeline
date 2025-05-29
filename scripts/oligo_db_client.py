#!/usr/bin/env python3
"""
OligoDB Client Library

A Python client for interacting with the oligo design PostgreSQL database.
Provides high-level methods for common operations and sequence searches.
"""

import psycopg2
import psycopg2.extras
import pandas as pd
from typing import List, Dict, Optional, Tuple
import logging
from contextlib import contextmanager

class OligoDBClient:
    """Client for interacting with the oligo design database"""
    
    def __init__(self, host="localhost", port=5432, database="oligodb", 
                 user="postgres", password="", **kwargs):
        """Initialize database connection"""
        self.connection_params = {
            'host': host,
            'port': port,
            'database': database,
            'user': user,
            'password': password,
            **kwargs
        }
        self.logger = logging.getLogger(__name__)
        
        # Test connection
        self._test_connection()
    
    def _test_connection(self):
        """Test database connection"""
        try:
            with self.get_connection() as conn:
                with conn.cursor() as cur:
                    cur.execute("SELECT 1")
            self.logger.info("Database connection successful")
        except Exception as e:
            self.logger.error(f"Database connection failed: {e}")
            raise
    
    @contextmanager
    def get_connection(self):
        """Get database connection with context manager"""
        conn = psycopg2.connect(**self.connection_params)
        try:
            yield conn
        finally:
            conn.close()
    
    def execute_query(self, query: str, params: Optional[Tuple] = None, 
                     fetch_all: bool = True) -> List[Tuple]:
        """Execute query and return results"""
        with self.get_connection() as conn:
            with conn.cursor() as cur:
                cur.execute(query, params or ())
                if fetch_all:
                    return cur.fetchall()
                else:
                    return cur.fetchone()
    
    def execute_query_df(self, query: str, params: Optional[Tuple] = None) -> pd.DataFrame:
        """Execute query and return results as DataFrame"""
        with self.get_connection() as conn:
            return pd.read_sql_query(query, conn, params=params)
    
    # =================================================================
    # Project Management
    # =================================================================
    
    def list_projects(self) -> pd.DataFrame:
        """List all projects in the database"""
        query = """
        SELECT 
            biodatabase_id,
            name,
            authority,
            description,
            created_date,
            (SELECT COUNT(*) FROM bioentry WHERE biodatabase_id = bd.biodatabase_id) as oligo_count
        FROM biodatabase bd
        ORDER BY created_date DESC
        """
        return self.execute_query_df(query)
    
    def create_project(self, name: str, description: str = "", 
                      authority: str = "OligoClient") -> int:
        """Create a new project"""
        query = """
        INSERT INTO biodatabase (name, authority, description)
        VALUES (%s, %s, %s)
        RETURNING biodatabase_id
        """
        with self.get_connection() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (name, authority, description))
                project_id = cur.fetchone()[0]
                conn.commit()
        
        self.logger.info(f"Created project '{name}' with ID {project_id}")
        return project_id
    
    def get_project_summary(self, project_name: str) -> Dict:
        """Get summary statistics for a project"""
        query = """
        SELECT 
            bd.name as project_name,
            bd.description,
            COUNT(DISTINCT od.oligo_design_id) as total_oligos,
            COUNT(DISTINCT CASE WHEN ot.name = 'forward_primer' THEN od.oligo_design_id END) as forward_primers,
            COUNT(DISTINCT CASE WHEN ot.name = 'reverse_primer' THEN od.oligo_design_id END) as reverse_primers,
            COUNT(DISTINCT CASE WHEN ot.name = 'capture_probe' THEN od.oligo_design_id END) as capture_probes,
            COUNT(DISTINCT pp.primer_pair_id) as primer_pairs,
            AVG(od.tm_calculated) as avg_tm,
            AVG(od.gc_content) as avg_gc,
            AVG(od.specificity_score) as avg_specificity
        FROM biodatabase bd
        LEFT JOIN bioentry be ON bd.biodatabase_id = be.biodatabase_id
        LEFT JOIN oligo_design od ON be.bioentry_id = od.bioentry_id
        LEFT JOIN oligo_type ot ON od.oligo_type_id = ot.oligo_type_id
        LEFT JOIN primer_pair pp ON (od.oligo_design_id = pp.forward_oligo_id 
                                     OR od.oligo_design_id = pp.reverse_oligo_id)
        WHERE bd.name = %s
        GROUP BY bd.biodatabase_id, bd.name, bd.description
        """
        result = self.execute_query(query, (project_name,), fetch_all=False)
        if result:
            columns = ['project_name', 'description', 'total_oligos', 'forward_primers',
                      'reverse_primers', 'capture_probes', 'primer_pairs', 'avg_tm',
                      'avg_gc', 'avg_specificity']
            return dict(zip(columns, result))
        return {}
    
    # =================================================================
    # Oligo Searches
    # =================================================================
    
    def search_oligos(self, project_name: Optional[str] = None,
                     oligo_type: Optional[str] = None,
                     design_mode: Optional[str] = None,
                     min_tm: Optional[float] = None,
                     max_tm: Optional[float] = None,
                     min_gc: Optional[float] = None,
                     max_gc: Optional[float] = None,
                     min_specificity: Optional[float] = None,
                     validation_status: Optional[str] = None,
                     limit: int = 1000) -> pd.DataFrame:
        """Search oligos with various filters"""
        
        conditions = []
        params = []
        
        if project_name:
            conditions.append("project_name = %s")
            params.append(project_name)
        
        if oligo_type:
            conditions.append("oligo_type = %s")
            params.append(oligo_type)
        
        if design_mode:
            conditions.append("design_mode = %s")
            params.append(design_mode)
        
        if min_tm is not None:
            conditions.append("tm_calculated >= %s")
            params.append(min_tm)
        
        if max_tm is not None:
            conditions.append("tm_calculated <= %s")
            params.append(max_tm)
        
        if min_gc is not None:
            conditions.append("gc_content >= %s")
            params.append(min_gc)
        
        if max_gc is not None:
            conditions.append("gc_content <= %s")
            params.append(max_gc)
        
        if min_specificity is not None:
            conditions.append("specificity_score >= %s")
            params.append(min_specificity)
        
        if validation_status:
            conditions.append("validation_status = %s")
            params.append(validation_status)
        
        where_clause = " AND ".join(conditions) if conditions else "1=1"
        
        query = f"""
        SELECT * FROM v_oligo_complete
        WHERE {where_clause}
        ORDER BY oligo_name
        LIMIT %s
        """
        params.append(limit)
        
        return self.execute_query_df(query, tuple(params))
    
    def search_by_name(self, name_pattern: str, exact_match: bool = False) -> pd.DataFrame:
        """Search oligos by name pattern"""
        if exact_match:
            query = "SELECT * FROM v_oligo_complete WHERE oligo_name = %s"
            params = (name_pattern,)
        else:
            query = "SELECT * FROM v_oligo_complete WHERE oligo_name ILIKE %s ORDER BY oligo_name"
            params = (f"%{name_pattern}%",)
        
        return self.execute_query_df(query, params)
    
    def get_oligo_details(self, oligo_name: str) -> Dict:
        """Get detailed information for a specific oligo"""
        query = """
        SELECT 
            oc.*,
            owt.total_off_targets,
            owt.high_similarity_hits,
            owt.max_off_target_identity
        FROM v_oligo_complete oc
        LEFT JOIN v_oligo_with_off_targets owt ON oc.oligo_design_id = owt.oligo_design_id
        WHERE oc.oligo_name = %s
        """
        df = self.execute_query_df(query, (oligo_name,))
        return df.iloc[0].to_dict() if not df.empty else {}
    
    # =================================================================
    # Sequence Similarity Searches
    # =================================================================
    
    def find_similar_sequences(self, query_sequence: str, 
                              similarity_threshold: float = 0.7,
                              max_results: int = 100) -> pd.DataFrame:
        """Find sequences similar to query using trigram similarity"""
        query = """
        SELECT 
            oligo_name,
            sequence,
            oligo_type,
            design_mode,
            tm_calculated,
            gc_content,
            SIMILARITY(sequence, %s) as similarity_score
        FROM v_oligo_complete 
        WHERE sequence %% %s
          AND SIMILARITY(sequence, %s) >= %s
        ORDER BY similarity_score DESC
        LIMIT %s
        """
        params = (query_sequence, query_sequence, query_sequence, 
                 similarity_threshold, max_results)
        
        return self.execute_query_df(query, params)
    
    def find_cross_reactive_primers(self, similarity_threshold: float = 0.8,
                                  max_results: int = 50) -> pd.DataFrame:
        """Find potentially cross-reactive primer pairs"""
        query = """
        SELECT 
            o1.oligo_name as primer1,
            o2.oligo_name as primer2,
            o1.sequence as seq1,
            o2.sequence as seq2,
            SIMILARITY(o1.sequence, o2.sequence) as similarity,
            o1.oligo_type as type1,
            o2.oligo_type as type2
        FROM v_oligo_complete o1
        JOIN v_oligo_complete o2 ON o1.oligo_design_id < o2.oligo_design_id
        WHERE o1.sequence %% o2.sequence
          AND SIMILARITY(o1.sequence, o2.sequence) >= %s
        ORDER BY similarity DESC
        LIMIT %s
        """
        return self.execute_query_df(query, (similarity_threshold, max_results))
    
    def check_sequence_novelty(self, sequence: str, 
                              similarity_threshold: float = 0.9) -> Dict:
        """Check if a sequence is novel or similar to existing ones"""
        similar = self.find_similar_sequences(sequence, similarity_threshold, 10)
        
        return {
            'sequence': sequence,
            'is_novel': len(similar) == 0,
            'similar_count': len(similar),
            'most_similar': similar.iloc[0].to_dict() if not similar.empty else None,
            'max_similarity': similar['similarity_score'].max() if not similar.empty else 0.0
        }
    
    # =================================================================
    # Primer Pairs
    # =================================================================
    
    def get_primer_pairs(self, project_name: Optional[str] = None,
                        min_product_size: Optional[int] = None,
                        max_product_size: Optional[int] = None,
                        max_tm_difference: Optional[float] = None) -> pd.DataFrame:
        """Get primer pairs with optional filters"""
        
        conditions = []
        params = []
        
        if project_name:
            # Need to join to get project name
            base_query = """
            SELECT 
                pp.*,
                bd.name as project_name
            FROM v_primer_pairs pp
            JOIN v_oligo_complete oc ON pp.forward_name = oc.oligo_name
            JOIN biodatabase bd ON oc.project_name = bd.name
            """
            conditions.append("bd.name = %s")
            params.append(project_name)
        else:
            base_query = "SELECT * FROM v_primer_pairs pp"
        
        if min_product_size is not None:
            conditions.append("pp.product_size >= %s")
            params.append(min_product_size)
        
        if max_product_size is not None:
            conditions.append("pp.product_size <= %s")
            params.append(max_product_size)
        
        if max_tm_difference is not None:
            conditions.append("ABS(pp.forward_tm - pp.reverse_tm) <= %s")
            params.append(max_tm_difference)
        
        where_clause = " AND ".join(conditions) if conditions else "1=1"
        query = f"{base_query} WHERE {where_clause} ORDER BY pp.primer_pair_id"
        
        return self.execute_query_df(query, tuple(params))
    
    def create_primer_pair(self, forward_oligo_name: str, reverse_oligo_name: str,
                          product_size: Optional[int] = None) -> int:
        """Create a new primer pair"""
        
        # Get oligo design IDs
        forward_query = """
        SELECT oligo_design_id FROM v_oligo_complete 
        WHERE oligo_name = %s AND oligo_type = 'forward_primer'
        """
        reverse_query = """
        SELECT oligo_design_id FROM v_oligo_complete 
        WHERE oligo_name = %s AND oligo_type = 'reverse_primer'
        """
        
        forward_result = self.execute_query(forward_query, (forward_oligo_name,), False)
        reverse_result = self.execute_query(reverse_query, (reverse_oligo_name,), False)
        
        if not forward_result or not reverse_result:
            raise ValueError("Could not find specified primers")
        
        forward_id = forward_result[0]
        reverse_id = reverse_result[0]
        
        # Create primer pair
        insert_query = """
        INSERT INTO primer_pair (forward_oligo_id, reverse_oligo_id, product_size)
        VALUES (%s, %s, %s)
        RETURNING primer_pair_id
        """
        
        with self.get_connection() as conn:
            with conn.cursor() as cur:
                cur.execute(insert_query, (forward_id, reverse_id, product_size))
                pair_id = cur.fetchone()[0]
                conn.commit()
        
        return pair_id
    
    # =================================================================
    # Off-target Analysis
    # =================================================================
    
    def get_off_targets(self, oligo_name: str, min_identity: float = 80.0) -> pd.DataFrame:
        """Get off-target hits for a specific oligo"""
        query = """
        SELECT 
            oth.chromosome,
            oth.start_pos,
            oth.end_pos,
            oth.strand,
            oth.identity_percent,
            oth.mismatch_count,
            oth.alignment_length,
            oth.gene_symbol,
            oth.feature_type
        FROM off_target_hit oth
        JOIN oligo_design od ON oth.oligo_design_id = od.oligo_design_id
        JOIN bioentry be ON od.bioentry_id = be.bioentry_id
        WHERE be.name = %s
          AND oth.identity_percent >= %s
        ORDER BY oth.identity_percent DESC
        """
        return self.execute_query_df(query, (oligo_name, min_identity))
    
    def get_high_off_target_oligos(self, max_off_targets: int = 10) -> pd.DataFrame:
        """Get oligos with high off-target counts"""
        query = """
        SELECT 
            oligo_name,
            sequence,
            oligo_type,
            total_off_targets,
            high_similarity_hits,
            max_off_target_identity,
            specificity_score
        FROM v_oligo_with_off_targets
        WHERE total_off_targets > %s
        ORDER BY total_off_targets DESC
        """
        return self.execute_query_df(query, (max_off_targets,))
    
    # =================================================================
    # Statistics and Analysis
    # =================================================================
    
    def get_database_stats(self) -> Dict:
        """Get overall database statistics"""
        query = "SELECT * FROM v_database_stats"
        results = self.execute_query(query)
        return {metric: value for metric, value in results}
    
    def get_quality_distribution(self, project_name: Optional[str] = None) -> Dict:
        """Get quality metric distributions"""
        
        where_clause = "WHERE project_name = %s" if project_name else ""
        params = (project_name,) if project_name else ()
        
        query = f"""
        SELECT 
            oligo_type,
            COUNT(*) as count,
            AVG(tm_calculated) as avg_tm,
            STDDEV(tm_calculated) as stddev_tm,
            AVG(gc_content) as avg_gc,
            STDDEV(gc_content) as stddev_gc,
            AVG(specificity_score) as avg_specificity,
            AVG(off_target_count) as avg_off_targets
        FROM v_oligo_complete
        {where_clause}
        GROUP BY oligo_type
        """
        
        df = self.execute_query_df(query, params)
        return df.set_index('oligo_type').to_dict('index')
    
    def analyze_sequence_patterns(self, pattern: str, 
                                project_name: Optional[str] = None) -> Dict:
        """Analyze frequency of sequence patterns"""
        
        where_clause = "WHERE project_name = %s" if project_name else ""
        params = [pattern]
        if project_name:
            params.append(project_name)
        
        query = f"""
        SELECT 
            oligo_type,
            COUNT(*) as total_oligos,
            COUNT(*) FILTER (WHERE sequence ~ %s) as pattern_matches,
            ROUND(COUNT(*) FILTER (WHERE sequence ~ %s) * 100.0 / COUNT(*), 2) as percentage
        FROM v_oligo_complete
        {where_clause}
        GROUP BY oligo_type
        """
        
        return self.execute_query_df(query, tuple(params)).to_dict('records')
    
    # =================================================================
    # Data Management
    # =================================================================
    
    def update_validation_status(self, oligo_names: List[str], 
                                status: str) -> int:
        """Update validation status for multiple oligos"""
        
        query = """
        UPDATE oligo_design 
        SET validation_status = %s, modified_date = CURRENT_TIMESTAMP
        WHERE bioentry_id IN (
            SELECT bioentry_id FROM bioentry WHERE name = ANY(%s)
        )
        """
        
        with self.get_connection() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (status, oligo_names))
                updated_count = cur.rowcount
                conn.commit()
        
        self.logger.info(f"Updated validation status for {updated_count} oligos")
        return updated_count
    
    def export_project_data(self, project_name: str, format: str = 'csv') -> pd.DataFrame:
        """Export all data for a project"""
        
        query = """
        SELECT 
            oligo_name,
            sequence,
            length,
            oligo_type,
            design_mode,
            tm_calculated,
            gc_content,
            specificity_score,
            off_target_count,
            validation_status,
            target_region,
            created_date
        FROM v_oligo_complete
        WHERE project_name = %s
        ORDER BY oligo_name
        """
        
        df = self.execute_query_df(query, (project_name,))
        
        if format.lower() == 'fasta':
            # Convert to FASTA format
            fasta_lines = []
            for _, row in df.iterrows():
                header = f">{row['oligo_name']} |tm={row['tm_calculated']:.1f}|gc={row['gc_content']:.1f}|type={row['oligo_type']}"
                fasta_lines.append(header)
                fasta_lines.append(row['sequence'])
            return '\n'.join(fasta_lines)
        
        return df
    
    # =================================================================
    # Utility Methods
    # =================================================================
    
    def execute_custom_query(self, query: str, params: Optional[Tuple] = None) -> pd.DataFrame:
        """Execute a custom SQL query and return DataFrame"""
        return self.execute_query_df(query, params)
    
    def close(self):
        """Close database connections (context manager handles this automatically)"""
        pass

# =================================================================
# Convenience Functions
# =================================================================

def create_client_from_config(config_file: str) -> OligoDBClient:
    """Create client from YAML config file"""
    import yaml
    
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    db_config = config.get('database', {})
    return OligoDBClient(**db_config)

# =================================================================
# Example Usage
# =================================================================

if __name__ == "__main__":
    # Example usage
    client = OligoDBClient()
    
    # List projects
    projects = client.list_projects()
    print("Projects:", projects)
    
    # Search for high-quality primers
    quality_primers = client.search_oligos(
        min_tm=58, max_tm=62,
        min_gc=40, max_gc=60,
        min_specificity=0.8,
        validation_status='validated'
    )
    print(f"Found {len(quality_primers)} high-quality primers")
    
    # Find similar sequences
    if not quality_primers.empty:
        sample_sequence = quality_primers.iloc[0]['sequence']
        similar = client.find_similar_sequences(sample_sequence, 0.8)
        print(f"Found {len(similar)} similar sequences")
    
    # Get database statistics
    stats = client.get_database_stats()
    print("Database stats:", stats)


