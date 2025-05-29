# PostgreSQL + BioSQL Database Integration Guide

## 🎯 Overview

The oligo design pipeline integrates with PostgreSQL using the BioSQL schema to provide powerful sequence storage, similarity searching, and analysis capabilities. This integration uses `pg_trgm` extensions for fast trigram-based sequence similarity searches.

## 🔧 Setup

### Database Installation

```bash
# Install PostgreSQL (Ubuntu/Debian)
sudo apt-get install postgresql postgresql-contrib

# Install required extensions
sudo -u postgres psql << EOF
CREATE DATABASE oligodb;
\c oligodb;
CREATE EXTENSION pg_trgm;
CREATE EXTENSION btree_gin;
CREATE EXTENSION btree_gist;
EOF
```

### Pipeline Configuration

Update `config.yaml`:

```yaml
database:
  host: "localhost"
  port: 5432
  database: "oligodb"
  user: "postgres"
  password: ""  # Set via PGPASSWORD environment variable
  project_name: "MyOligoProject_2025"
  project_description: "Custom oligo design project"
```

### Environment Setup

```bash
# Set database password
export PGPASSWORD="your_password"

# Or create .pgpass file
echo "localhost:5432:oligodb:postgres:your_password" > ~/.pgpass
chmod 600 ~/.pgpass
```

## 🚀 Running with Database Integration

### Full Pipeline with Database

```bash
# Run complete pipeline including database storage
snakemake database --cores 8

# Or run everything including database
snakemake all --cores 8
```

### Database-only Operations

```bash
# Setup database schema
snakemake results/database/setup_complete.txt

# Load existing results into database
snakemake results/database/pcr_loaded.txt results/database/panel_loaded.txt

# Generate database summary
snakemake results/database/database_summary.html

# Export data from database
snakemake results/database/exports/all_oligos.fasta
```

## 📊 Database Schema Overview

### Core Tables

| Table | Purpose |
|-------|---------|
| `biodatabase` | Projects/experiments |
| `bioentry` | Individual oligos/primers |
| `biosequence` | Actual sequence data |
| `oligo_design` | Design-specific properties |
| `primer_pair` | PCR primer pairs |
| `off_target_hit` | Off-target analysis results |

### Key Views

| View | Description |
|------|-------------|
| `v_oligo_complete` | Complete oligo information |
| `v_primer_pairs` | Primer pairs with sequences |
| `v_oligo_with_off_targets` | Oligos with off-target summary |

## 🔍 Query Examples

### 1. Basic Sequence Searches

```sql
-- Find all primers for a project
SELECT oligo_name, sequence, oligo_type, tm_calculated
FROM v_oligo_complete 
WHERE project_name = 'MyOligoProject_2025';

-- Search by name pattern
SELECT oligo_name, sequence 
FROM v_oligo_complete 
WHERE oligo_name ILIKE '%GAPDH%';
```

### 2. Quality-Based Filtering

```sql
-- High-quality primers
SELECT oligo_name, sequence, tm_calculated, gc_content, specificity_score
FROM v_oligo_complete 
WHERE tm_calculated BETWEEN 58 AND 62
  AND gc_content BETWEEN 40 AND 60
  AND specificity_score > 0.8
ORDER BY specificity_score DESC;

-- Primers with no off-targets
SELECT oligo_name, sequence, specificity_score
FROM v_oligo_complete 
WHERE off_target_count = 0
  AND validation_status = 'validated';
```

### 3. Sequence Similarity Searches (pg_trgm)

```sql
-- Find sequences similar to query (trigram similarity)
SELECT 
    oligo_name,
    sequence,
    SIMILARITY(sequence, 'GCATACGTTGTATCCGGGCAT') as similarity_score
FROM v_oligo_complete 
WHERE sequence % 'GCATACGTTGTATCCGGGCAT'  -- % operator for trigram matching
  AND SIMILARITY(sequence, 'GCATACGTTGTATCCGGGCAT') > 0.7
ORDER BY similarity_score DESC;

-- Cross-reactivity analysis
SELECT 
    o1.oligo_name as primer1,
    o2.oligo_name as primer2,
    SIMILARITY(o1.sequence, o2.sequence) as similarity
FROM v_oligo_complete o1
JOIN v_oligo_complete o2 ON o1.oligo_design_id < o2.oligo_design_id
WHERE o1.sequence % o2.sequence
  AND SIMILARITY(o1.sequence, o2.sequence) > 0.8
ORDER BY similarity DESC;
```

### 4. Advanced Pattern Matching

```sql
-- Find sequences with specific motifs
SELECT oligo_name, sequence
FROM v_oligo_complete 
WHERE sequence ~ 'G{3,}'  -- 3+ consecutive Gs
   OR sequence LIKE '%CpG%';

-- Reverse complement analysis
SELECT 
    oligo_name,
    sequence,
    reverse_complement(sequence) as rev_comp
FROM v_oligo_complete 
WHERE sequence % reverse_complement(sequence);
```

### 5. Primer Pair Analysis

```sql
-- Well-matched primer pairs
SELECT 
    forward_name,
    reverse_name,
    ABS(forward_tm - reverse_tm) as tm_difference,
    product_size
FROM v_primer_pairs 
WHERE ABS(forward_tm - reverse_tm) <= 2
  AND product_size BETWEEN 100 AND 300
ORDER BY tm_difference;

-- Primer pairs by target region
SELECT COUNT(*) as pair_count, 
       AVG(product_size) as avg_product_size
FROM primer_pair pp
JOIN oligo_design od ON pp.forward_oligo_id = od.oligo_design_id
WHERE od.target_region LIKE 'chr22%'
GROUP BY od.target_region;
```

### 6. Off-Target Analysis

```sql
-- Primers with high off-target potential
SELECT 
    oligo_name,
    total_off_targets,
    high_similarity_hits,
    max_off_target_identity
FROM v_oligo_with_off_targets
WHERE high_similarity_hits > 5
ORDER BY high_similarity_hits DESC;

-- Off-target distribution by chromosome
SELECT 
    chromosome,
    COUNT(*) as hit_count,
    AVG(identity_percent) as avg_identity
FROM off_target_hit
GROUP BY chromosome
ORDER BY hit_count DESC;
```

### 7. Statistical Analysis

```sql
-- Quality metrics by design mode
SELECT 
    design_mode,
    COUNT(*) as total_oligos,
    AVG(tm_calculated) as avg_tm,
    AVG(gc_content) as avg_gc,
    AVG(specificity_score) as avg_specificity
FROM v_oligo_complete
GROUP BY design_mode;

-- Success rate analysis
SELECT 
    validation_status,
    COUNT(*) as count,
    ROUND(COUNT(*) * 100.0 / SUM(COUNT(*)) OVER (), 2) as percentage
FROM v_oligo_complete
GROUP BY validation_status;
```

## 🔧 Advanced Features

### Custom Functions

```sql
-- Use built-in functions
SELECT calculate_gc_content('GCATACGTTGTATCCGGGCAT');
SELECT reverse_complement('ATGCCCGGATACAACGTATGC');

-- Find similar sequences with custom threshold
SELECT * FROM find_similar_sequences('GCATACGTTGTATCCGGGCAT', 0.85);
```

### Performance Optimization

```sql
-- Set similarity threshold for faster searches
SET pg_trgm.similarity_threshold = 0.8;

-- Check index usage
SELECT schemaname, tablename, indexname, idx_scan
FROM pg_stat_user_indexes 
WHERE indexname LIKE '%trgm%';

-- Update table statistics
ANALYZE biosequence;
```

### Batch Operations

```sql
-- Bulk update validation status
UPDATE oligo_design 
SET validation_status = 'auto_approved'
WHERE tm_calculated BETWEEN 58 AND 62
  AND gc_content BETWEEN 40 AND 60
  AND specificity_score > 0.9;

-- Flag problematic sequences
UPDATE oligo_design 
SET validation_status = 'needs_review'
WHERE oligo_design_id IN (
    SELECT od.oligo_design_id
    FROM oligo_design od
    JOIN biosequence bs ON od.bioentry_id = bs.bioentry_id
    WHERE bs.seq ~ '[ATGC]{4,}'  -- 4+ nucleotide repeats
);
```

## 📈 Monitoring and Maintenance

### Database Statistics

```sql
-- Check database size and growth
SELECT 
    pg_size_pretty(pg_database_size('oligodb')) as database_size,
    COUNT(*) as total_oligos
FROM oligo_design;

-- Monitor recent activity
SELECT 
    DATE_TRUNC('day', created_date) as date,
    COUNT(*) as oligos_added
FROM oligo_design
WHERE created_date >= CURRENT_DATE - INTERVAL '7 days'
GROUP BY DATE_TRUNC('day', created_date)
ORDER BY date;
```

### Performance Monitoring

```sql
-- Query performance tracking
SELECT 
    query_type,
    AVG(execution_time_ms) as avg_time_ms,
    COUNT(*) as query_count
FROM query_performance
WHERE executed_at >= CURRENT_DATE - INTERVAL '24 hours'
GROUP BY query_type
ORDER BY avg_time_ms DESC;

-- Index efficiency
SELECT 
    indexname,
    idx_scan,
    idx_tup_read,
    idx_tup_fetch
FROM pg_stat_user_indexes
WHERE schemaname = 'public'
ORDER BY idx_scan DESC;
```

## 🔄 Data Export and Integration

### Export Formats

```bash
# Export all data (automatically done by pipeline)
snakemake results/database/exports/all_oligos.fasta    # FASTA format
snakemake results/database/exports/all_oligos.csv      # CSV format  
snakemake results/database/exports/all_oligos.json     # JSON format
```

### Integration with External Tools

```sql
-- Export for external primer design tools
SELECT 
    '>' || oligo_name || E'\n' || sequence as fasta_format
FROM v_oligo_complete
WHERE validation_status = 'validated';

-- Export primer pairs for PCR setup
SELECT 
    forward_name,
    forward_sequence,
    reverse_name,  
    reverse_sequence,
    product_size
FROM v_primer_pairs
WHERE validation_status = 'validated'
  AND product_size BETWEEN 100 AND 500;
```

## 🛠️ Troubleshooting

### Common Issues

1. **Slow Similarity Searches**
   ```sql
   -- Check trigram threshold
   SHOW pg_trgm.similarity_threshold;
   
   -- Rebuild trigram indexes
   REINDEX INDEX CONCURRENTLY idx_biosequence_seq_gin;
   ```

2. **Connection Errors**
   ```bash
   # Test connection
   psql -h localhost -U postgres -d oligodb -c "SELECT VERSION();"
   
   # Check configuration
   grep -E "host|port" config.yaml
   ```

3. **Missing Data**
   ```sql
   -- Check data loading status
   SELECT * FROM v_database_stats;
   
   -- Verify pipeline completion
   SELECT run_name, status, total_candidates 
   FROM pipeline_run 
   ORDER BY start_time DESC;
   ```

### Performance Tuning

```sql
-- Optimize for large datasets
SET work_mem = '256MB';
SET shared_buffers = '1GB';
SET effective_cache_size = '4GB';

-- For similarity searches
SET pg_trgm.similarity_threshold = 0.8;  -- Higher = faster but less sensitive
```

## 📚 Best Practices

### Query Optimization

1. **Use indexes effectively**
   - Trigram searches with % operator
   - Filter by indexed columns first
   - Use LIMIT for large result sets

2. **Batch operations**
   - Group related INSERTs in transactions
   - Use COPY for bulk data loading
   - Update statistics after large changes

3. **Similarity search tips**
   - Start with higher similarity thresholds
   - Combine with other filters (length, GC content)
   - Use partial indexes for specific use cases

### Data Management

1. **Regular maintenance**
   ```sql
   -- Weekly maintenance
   VACUUM ANALYZE;
   REINDEX INDEX CONCURRENTLY idx_biosequence_seq_gin;
   ```

2. **Backup strategy**
   ```bash
   # Automated backup (in pipeline)
   snakemake results/database/backup_$(date +%Y%m%d).sql
   ```

3. **Monitoring**
   - Track query performance
   - Monitor index usage
   - Watch database growth

This database integration provides powerful sequence analysis capabilities while maintaining the flexibility and performance needed for large-scale oligo design projects.
