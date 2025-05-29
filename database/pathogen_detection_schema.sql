-- =============================================================================
-- Pathogen Detection Database Schema Extensions
-- =============================================================================
-- Extends the base BioSQL schema with pathogen-specific tables and functions
-- for taxonomic classification, host exclusion, and diagnostic applications

-- =============================================================================
-- Taxonomic Information Tables
-- =============================================================================

-- NCBI Taxonomy nodes
CREATE TABLE IF NOT EXISTS taxonomy_node (
    tax_id INTEGER PRIMARY KEY,
    parent_tax_id INTEGER,
    rank VARCHAR(50),
    embl_code VARCHAR(10),
    division_id INTEGER,
    inherited_div_flag BOOLEAN DEFAULT FALSE,
    genetic_code_id INTEGER,
    inherited_gc_flag BOOLEAN DEFAULT FALSE,
    mitochondrial_genetic_code_id INTEGER,
    inherited_mgc_flag BOOLEAN DEFAULT FALSE,
    genbank_hidden_flag BOOLEAN DEFAULT FALSE,
    hidden_subtree_root_flag BOOLEAN DEFAULT FALSE,
    comments TEXT,
    created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Taxonomic names
CREATE TABLE IF NOT EXISTS taxonomy_name (
    tax_id INTEGER REFERENCES taxonomy_node(tax_id),
    name_txt VARCHAR(500) NOT NULL,
    unique_name VARCHAR(500),
    name_class VARCHAR(50) NOT NULL,
    created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (tax_id, name_txt, name_class)
);

-- Pathogen classification
CREATE TABLE IF NOT EXISTS pathogen_classification (
    pathogen_id SERIAL PRIMARY KEY,
    tax_id INTEGER REFERENCES taxonomy_node(tax_id),
    pathogen_type VARCHAR(50) NOT NULL, -- 'bacteria', 'virus', 'fungi', 'parasite'
    pathogenicity_level VARCHAR(20), -- 'BSL-1', 'BSL-2', 'BSL-3', 'BSL-4'
    clinical_significance VARCHAR(100),
    antimicrobial_resistance BOOLEAN DEFAULT FALSE,
    zoonotic BOOLEAN DEFAULT FALSE,
    foodborne BOOLEAN DEFAULT FALSE,
    waterborne BOOLEAN DEFAULT FALSE,
    airborne BOOLEAN DEFAULT FALSE,
    vector_borne BOOLEAN DEFAULT FALSE,
    healthcare_associated BOOLEAN DEFAULT FALSE,
    bioterrorism_agent BOOLEAN DEFAULT FALSE,
    reportable_disease BOOLEAN DEFAULT FALSE,
    created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Host organisms (for exclusion)
CREATE TABLE IF NOT EXISTS host_organism (
    host_id SERIAL PRIMARY KEY,
    tax_id INTEGER REFERENCES taxonomy_node(tax_id),
    common_name VARCHAR(200),
    genome_assembly VARCHAR(100),
    exclusion_priority INTEGER DEFAULT 1, -- Higher priority = more important to exclude
    clinical_relevance VARCHAR(200),
    created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Insert common host organisms
INSERT INTO host_organism (tax_id, common_name, genome_assembly, exclusion_priority, clinical_relevance) VALUES
    (9606, 'Homo sapiens', 'GRCh38', 10, 'Primary clinical host'),
    (10090, 'Mus musculus', 'GRCm39', 8, 'Laboratory model organism'),
    (10116, 'Rattus norvegicus', 'Rnor6.0', 7, 'Laboratory model organism'),
    (9913, 'Bos taurus', 'ARS-UCD1.2', 6, 'Livestock - food safety'),
    (9031, 'Gallus gallus', 'GRCg6a', 6, 'Livestock - food safety'),
    (9823, 'Sus scrofa', 'Sscrofa11.1', 6, 'Livestock - food safety')
ON CONFLICT (tax_id) DO NOTHING;

-- =============================================================================
-- Pathogen-Specific Oligo Design Tables
-- =============================================================================

-- Extend oligo_design table with pathogen-specific fields
ALTER TABLE oligo_design ADD COLUMN IF NOT EXISTS target_tax_id INTEGER REFERENCES taxonomy_node(tax_id);
ALTER TABLE oligo_design ADD COLUMN IF NOT EXISTS target_gene VARCHAR(100);
ALTER TABLE oligo_design ADD COLUMN IF NOT EXISTS target_region_type VARCHAR(50); -- '16S', '23S', 'ITS', 'gyrB', etc.
ALTER TABLE oligo_design ADD COLUMN IF NOT EXISTS specificity_level VARCHAR(50); -- 'strain', 'species', 'genus', 'family'
ALTER TABLE oligo_design ADD COLUMN IF NOT EXISTS clinical_application VARCHAR(100); -- 'diagnostic', 'surveillance', 'research'
ALTER TABLE oligo_design ADD COLUMN IF NOT EXISTS detection_limit_cfu_ml INTEGER; -- CFU/ml for sensitivity
ALTER TABLE oligo_design ADD COLUMN IF NOT EXISTS detection_limit_copies_ul INTEGER; -- copies/μl for qPCR

-- Pathogen target genes/regions
CREATE TABLE IF NOT EXISTS target_gene (
    target_gene_id SERIAL PRIMARY KEY,
    gene_name VARCHAR(100) NOT NULL,
    gene_symbol VARCHAR(20),
    target_type VARCHAR(50), -- 'housekeeping', 'virulence', 'resistance', 'species_specific'
    pathogen_type VARCHAR(50), -- 'bacteria', 'virus', 'fungi', 'parasite'
    conservation_level VARCHAR(50), -- 'highly_conserved', 'moderately_conserved', 'variable'
    copy_number INTEGER, -- typical copy number per genome
    function_description TEXT,
    diagnostic_utility TEXT,
    created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(gene_name, pathogen_type)
);

-- Insert common pathogen target genes
INSERT INTO target_gene (gene_name, gene_symbol, target_type, pathogen_type, conservation_level, copy_number, function_description, diagnostic_utility) VALUES
    ('16S ribosomal RNA', '16S', 'housekeeping', 'bacteria', 'highly_conserved', 3, 'Ribosomal RNA component', 'Universal bacterial detection and identification'),
    ('23S ribosomal RNA', '23S', 'housekeeping', 'bacteria', 'highly_conserved', 3, 'Ribosomal RNA component', 'Bacterial identification with higher resolution than 16S'),
    ('Internal Transcribed Spacer', 'ITS', 'housekeeping', 'fungi', 'moderately_conserved', 1, 'rRNA spacer region', 'Fungal species identification'),
    ('RNA polymerase B', 'rpoB', 'housekeeping', 'bacteria', 'moderately_conserved', 1, 'RNA polymerase subunit', 'Species-level bacterial identification'),
    ('DNA gyrase B', 'gyrB', 'housekeeping', 'bacteria', 'moderately_conserved', 1, 'DNA replication enzyme', 'Species and strain differentiation'),
    ('Heat shock protein 65', 'hsp65', 'housekeeping', 'bacteria', 'moderately_conserved', 1, 'Chaperone protein', 'Mycobacterial identification'),
    ('Elongation factor Tu', 'tuf', 'housekeeping', 'bacteria', 'moderately_conserved', 2, 'Protein synthesis factor', 'Species identification for difficult groups')
ON CONFLICT (gene_name, pathogen_type) DO NOTHING;

-- Cross-reactivity tracking
CREATE TABLE IF NOT EXISTS cross_reactivity (
    cross_reactivity_id SERIAL PRIMARY KEY,
    oligo_design_id INTEGER REFERENCES oligo_design(oligo_design_id),
    cross_reactive_tax_id INTEGER REFERENCES taxonomy_node(tax_id),
    cross_reactive_organism VARCHAR(200),
    similarity_percent DECIMAL(5,2),
    mismatch_count INTEGER,
    clinical_significance VARCHAR(50), -- 'critical', 'moderate', 'low', 'acceptable'
    phylogenetic_distance DECIMAL(8,4), -- evolutionary distance
    notes TEXT,
    created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Host cross-reactivity (critical for clinical applications)
CREATE TABLE IF NOT EXISTS host_cross_reactivity (
    host_cross_reactivity_id SERIAL PRIMARY KEY,
    oligo_design_id INTEGER REFERENCES oligo_design(oligo_design_id),
    host_id INTEGER REFERENCES host_organism(host_id),
    similarity_percent DECIMAL(5,2),
    alignment_length INTEGER,
    hit_location VARCHAR(200), -- genomic location of hit
    clinical_concern_level VARCHAR(20), -- 'critical', 'high', 'medium', 'low'
    mitigation_strategy TEXT,
    created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- =============================================================================
-- Antimicrobial Resistance (AMR) Integration
-- =============================================================================

-- AMR genes database
CREATE TABLE IF NOT EXISTS amr_gene (
    amr_gene_id SERIAL PRIMARY KEY,
    gene_name VARCHAR(100) NOT NULL,
    gene_family VARCHAR(100),
    resistance_mechanism VARCHAR(200),
    antibiotic_class VARCHAR(100),
    specific_antibiotics TEXT[],
    sequence_accession VARCHAR(50),
    reference_sequence TEXT,
    created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Link oligos to AMR genes
CREATE TABLE IF NOT EXISTS oligo_amr_target (
    oligo_design_id INTEGER REFERENCES oligo_design(oligo_design_id),
    amr_gene_id INTEGER REFERENCES amr_gene(amr_gene_id),
    detection_type VARCHAR(50), -- 'presence', 'expression', 'mutation'
    clinical_breakpoint VARCHAR(100),
    PRIMARY KEY (oligo_design_id, amr_gene_id)
);

-- =============================================================================
-- Virulence Factors
-- =============================================================================

-- Virulence factor database
CREATE TABLE IF NOT EXISTS virulence_factor (
    virulence_factor_id SERIAL PRIMARY KEY,
    factor_name VARCHAR(200) NOT NULL,
    factor_type VARCHAR(100), -- 'toxin', 'adhesin', 'invasin', 'immune_evasion'
    pathogen_tax_ids INTEGER[], -- array of tax_ids where this factor is found
    disease_association TEXT,
    sequence_accession VARCHAR(50),
    reference_sequence TEXT,
    created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Link oligos to virulence factors
CREATE TABLE IF NOT EXISTS oligo_virulence_target (
    oligo_design_id INTEGER REFERENCES oligo_design(oligo_design_id),
    virulence_factor_id INTEGER REFERENCES virulence_factor(virulence_factor_id),
    detection_significance VARCHAR(50), -- 'diagnostic', 'epidemiological', 'research'
    PRIMARY KEY (oligo_design_id, virulence_factor_id)
);

-- =============================================================================
-- Clinical Validation and Performance
-- =============================================================================

-- Clinical validation studies
CREATE TABLE IF NOT EXISTS clinical_validation (
    validation_id SERIAL PRIMARY KEY,
    oligo_design_id INTEGER REFERENCES oligo_design(oligo_design_id),
    study_name VARCHAR(200),
    study_type VARCHAR(50), -- 'sensitivity', 'specificity', 'clinical_trial'
    sample_size INTEGER,
    true_positive INTEGER,
    false_positive INTEGER,
    true_negative INTEGER,
    false_negative INTEGER,
    sensitivity_percent DECIMAL(5,2),
    specificity_percent DECIMAL(5,2),
    ppv_percent DECIMAL(5,2), -- positive predictive value
    npv_percent DECIMAL(5,2), -- negative predictive value
    reference_method VARCHAR(200),
    sample_types TEXT[],
    study_population VARCHAR(200),
    validation_date DATE,
    publication_doi VARCHAR(200),
    regulatory_approval VARCHAR(100), -- 'FDA', 'CE-IVD', 'Health Canada', etc.
    notes TEXT,
    created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Quality control for diagnostic use
CREATE TABLE IF NOT EXISTS diagnostic_qc (
    qc_id SERIAL PRIMARY KEY,
    oligo_design_id INTEGER REFERENCES oligo_design(oligo_design_id),
    qc_parameter VARCHAR(100),
    specification_value VARCHAR(100),
    tested_value VARCHAR(100),
    pass_fail VARCHAR(10),
    test_date DATE,
    operator VARCHAR(100),
    instrument VARCHAR(100),
    lot_number VARCHAR(50),
    notes TEXT,
    created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- =============================================================================
-- Pathogen-Specific Functions
-- =============================================================================

-- Function to get taxonomic lineage
CREATE OR REPLACE FUNCTION get_taxonomic_lineage(input_tax_id INTEGER)
RETURNS TABLE(
    tax_id INTEGER,
    rank VARCHAR,
    name VARCHAR
) AS $$
BEGIN
    RETURN QUERY
    WITH RECURSIVE taxonomy_tree AS (
        -- Base case: start with input tax_id
        SELECT 
            tn.tax_id,
            tn.parent_tax_id,
            tn.rank,
            tn2.name_txt as name
        FROM taxonomy_node tn
        JOIN taxonomy_name tn2 ON tn.tax_id = tn2.tax_id 
        WHERE tn.tax_id = input_tax_id AND tn2.name_class = 'scientific name'
        
        UNION ALL
        
        -- Recursive case: get parent nodes
        SELECT 
            tn.tax_id,
            tn.parent_tax_id,
            tn.rank,
            tn2.name_txt as name
        FROM taxonomy_node tn
        JOIN taxonomy_name tn2 ON tn.tax_id = tn2.tax_id
        JOIN taxonomy_tree tt ON tn.tax_id = tt.parent_tax_id
        WHERE tn.tax_id != tn.parent_tax_id AND tn2.name_class = 'scientific name'
    )
    SELECT 
        tt.tax_id,
        tt.rank,
        tt.name
    FROM taxonomy_tree tt
    ORDER BY 
        CASE tt.rank
            WHEN 'superkingdom' THEN 1
            WHEN 'kingdom' THEN 2
            WHEN 'phylum' THEN 3
            WHEN 'class' THEN 4
            WHEN 'order' THEN 5
            WHEN 'family' THEN 6
            WHEN 'genus' THEN 7
            WHEN 'species' THEN 8
            WHEN 'subspecies' THEN 9
            WHEN 'strain' THEN 10
            ELSE 99
        END;
END;
$$ LANGUAGE plpgsql;

-- Function to calculate phylogenetic distance (simplified)
CREATE OR REPLACE FUNCTION calculate_phylogenetic_distance(tax_id1 INTEGER, tax_id2 INTEGER)
RETURNS DECIMAL AS $$
DECLARE
    common_ancestor INTEGER;
    distance1 INTEGER;
    distance2 INTEGER;
BEGIN
    -- Find lowest common ancestor
    WITH lineage1 AS (
        SELECT tax_id FROM get_taxonomic_lineage(tax_id1)
    ),
    lineage2 AS (
        SELECT tax_id FROM get_taxonomic_lineage(tax_id2)
    )
    SELECT l1.tax_id INTO common_ancestor
    FROM lineage1 l1
    JOIN lineage2 l2 ON l1.tax_id = l2.tax_id
    ORDER BY l1.tax_id DESC
    LIMIT 1;
    
    -- Calculate distances to common ancestor
    SELECT COUNT(*) INTO distance1 FROM get_taxonomic_lineage(tax_id1) WHERE tax_id >= common_ancestor;
    SELECT COUNT(*) INTO distance2 FROM get_taxonomic_lineage(tax_id2) WHERE tax_id >= common_ancestor;
    
    RETURN (distance1 + distance2)::DECIMAL;
END;
$$ LANGUAGE plpgsql;

-- Function to check pathogen specificity
CREATE OR REPLACE FUNCTION check_pathogen_specificity(
    oligo_sequence TEXT,
    target_tax_id INTEGER,
    max_cross_reactivity DECIMAL DEFAULT 0.85
)
RETURNS TABLE(
    cross_reactive_organism VARCHAR,
    tax_id INTEGER,
    similarity DECIMAL,
    phylogenetic_distance DECIMAL,
    concern_level VARCHAR
) AS $$
BEGIN
    RETURN QUERY
    SELECT 
        tn.name_txt,
        bs.tax_id,
        SIMILARITY(bs.seq, oligo_sequence) as sim,
        calculate_phylogenetic_distance(target_tax_id, bs.tax_id) as phylo_dist,
        CASE 
            WHEN SIMILARITY(bs.seq, oligo_sequence) > 0.95 THEN 'CRITICAL'
            WHEN SIMILARITY(bs.seq, oligo_sequence) > 0.90 THEN 'HIGH'
            WHEN SIMILARITY(bs.seq, oligo_sequence) > 0.85 THEN 'MEDIUM'
            ELSE 'LOW'
        END as concern
    FROM biosequence bs
    JOIN bioentry be ON bs.bioentry_id = be.bioentry_id
    JOIN taxonomy_name tn ON be.tax_id = tn.tax_id
    WHERE bs.seq % oligo_sequence
      AND SIMILARITY(bs.seq, oligo_sequence) >= max_cross_reactivity
      AND be.tax_id != target_tax_id
      AND tn.name_class = 'scientific name'
    ORDER BY sim DESC;
END;
$$ LANGUAGE plpgsql;

-- =============================================================================
-- Views for Pathogen Detection
-- =============================================================================

-- Complete pathogen oligo view with taxonomic information
CREATE OR REPLACE VIEW v_pathogen_oligo_complete AS
SELECT 
    oc.*,
    tn.name_txt as target_organism,
    tn.rank as taxonomic_rank,
    pc.pathogen_type,
    pc.pathogenicity_level,
    pc.clinical_significance,
    pc.antimicrobial_resistance as has_amr,
    tg.target_type as gene_type,
    tg.conservation_level,
    tg.diagnostic_utility
FROM v_oligo_complete oc
LEFT JOIN taxonomy_name tn ON oc.target_tax_id = tn.tax_id AND tn.name_class = 'scientific name'
LEFT JOIN pathogen_classification pc ON oc.target_tax_id = pc.tax_id
LEFT JOIN target_gene tg ON oc.target_gene = tg.gene_name;

-- Clinical validation summary view
CREATE OR REPLACE VIEW v_clinical_performance AS
SELECT 
    od.oligo_design_id,
    be.name as oligo_name,
    tn.name_txt as target_organism,
    AVG(cv.sensitivity_percent) as avg_sensitivity,
    AVG(cv.specificity_percent) as avg_specificity,
    COUNT(cv.validation_id) as validation_studies,
    MAX(cv.validation_date) as latest_validation,
    COUNT(cv.regulatory_approval) FILTER (WHERE cv.regulatory_approval IS NOT NULL) as regulatory_approvals
FROM oligo_design od
JOIN bioentry be ON od.bioentry_id = be.bioentry_id
LEFT JOIN taxonomy_name tn ON od.target_tax_id = tn.tax_id AND tn.name_class = 'scientific name'
LEFT JOIN clinical_validation cv ON od.oligo_design_id = cv.oligo_design_id
GROUP BY od.oligo_design_id, be.name, tn.name_txt;

-- Cross-reactivity risk assessment view
CREATE OR REPLACE VIEW v_cross_reactivity_risk AS
SELECT 
    od.oligo_design_id,
    be.name as oligo_name,
    COUNT(cr.cross_reactivity_id) as total_cross_reactions,
    COUNT(cr.cross_reactivity_id) FILTER (WHERE cr.clinical_significance = 'critical') as critical_reactions,
    COUNT(hcr.host_cross_reactivity_id) as host_reactions,
    MAX(cr.similarity_percent) as max_pathogen_similarity,
    MAX(hcr.similarity_percent) as max_host_similarity,
    CASE 
        WHEN COUNT(cr.cross_reactivity_id) FILTER (WHERE cr.clinical_significance = 'critical') > 0 THEN 'HIGH RISK'
        WHEN COUNT(hcr.host_cross_reactivity_id) FILTER (WHERE hcr.clinical_concern_level = 'critical') > 0 THEN 'HIGH RISK'
        WHEN COUNT(cr.cross_reactivity_id) > 10 THEN 'MEDIUM RISK'
        ELSE 'LOW RISK'
    END as overall_risk_level
FROM oligo_design od
JOIN bioentry be ON od.bioentry_id = be.bioentry_id
LEFT JOIN cross_reactivity cr ON od.oligo_design_id = cr.oligo_design_id
LEFT JOIN host_cross_reactivity hcr ON od.oligo_design_id = hcr.oligo_design_id
GROUP BY od.oligo_design_id, be.name;

-- =============================================================================
-- Indexes for Pathogen Detection Performance
-- =============================================================================

-- Taxonomic indexes
CREATE INDEX IF NOT EXISTS idx_taxonomy_node_parent ON taxonomy_node(parent_tax_id);
CREATE INDEX IF NOT EXISTS idx_taxonomy_name_tax_id ON taxonomy_name(tax_id);
CREATE INDEX IF NOT EXISTS idx_taxonomy_name_name_class ON taxonomy_name(name_class);
CREATE INDEX IF NOT EXISTS idx_taxonomy_name_name_txt_trgm ON taxonomy_name USING gin(name_txt gin_trgm_ops);

-- Pathogen-specific indexes
CREATE INDEX IF NOT EXISTS idx_oligo_design_target_tax_id ON oligo_design(target_tax_id);
CREATE INDEX IF NOT EXISTS idx_oligo_design_target_gene ON oligo_design(target_gene);
CREATE INDEX IF NOT EXISTS idx_oligo_design_specificity_level ON oligo_design(specificity_level);
CREATE INDEX IF NOT EXISTS idx_pathogen_classification_type ON pathogen_classification(pathogen_type);
CREATE INDEX IF NOT EXISTS idx_cross_reactivity_similarity ON cross_reactivity(similarity_percent);

-- Clinical validation indexes
CREATE INDEX IF NOT EXISTS idx_clinical_validation_sensitivity ON clinical_validation(sensitivity_percent);
CREATE INDEX IF NOT EXISTS idx_clinical_validation_specificity ON clinical_validation(specificity_percent);

-- AMR and virulence indexes
CREATE INDEX IF NOT EXISTS idx_amr_gene_name_trgm ON amr_gene USING gin(gene_name gin_trgm_ops);
CREATE INDEX IF NOT EXISTS idx_virulence_factor_name_trgm ON virulence_factor USING gin(factor_name gin_trgm_ops);

-- =============================================================================
-- Sample Data for Testing
-- =============================================================================

-- Insert some example pathogen classifications
INSERT INTO pathogen_classification (tax_id, pathogen_type, pathogenicity_level, clinical_significance) VALUES
    (1280, 'bacteria', 'BSL-2', 'Major food-borne pathogen - Staphylococcus aureus'),
    (562, 'bacteria', 'BSL-1', 'Indicator organism and opportunistic pathogen - Escherichia coli'),
    (1773, 'bacteria', 'BSL-3', 'Obligate intracellular pathogen - Mycobacterium tuberculosis'),
    (90371, 'bacteria', 'BSL-2', 'Food-borne pathogen - Salmonella enterica'),
    (1352, 'bacteria', 'BSL-2', 'Biofilm-forming opportunistic pathogen - Enterococcus faecalis')
ON CONFLICT DO NOTHING;

-- Performance and maintenance functions
CREATE OR REPLACE FUNCTION update_pathogen_statistics()
RETURNS VOID AS $$
BEGIN
    -- Update table statistics for query optimization
    ANALYZE taxonomy_node;
    ANALYZE taxonomy_name;
    ANALYZE pathogen_classification;
    ANALYZE oligo_design;
    ANALYZE cross_reactivity;
    ANALYZE host_cross_reactivity;
    
    -- Log the update
    INSERT INTO query_performance (query_type, execution_time_ms, rows_returned, query_text)
    VALUES ('maintenance', 0, 0, 'update_pathogen_statistics() - ANALYZE completed');
END;
$$ LANGUAGE plpgsql;

COMMENT ON TABLE taxonomy_node IS 'NCBI taxonomy nodes for pathogen classification';
COMMENT ON TABLE pathogen_classification IS 'Pathogen-specific metadata including BSL levels and clinical significance';
COMMENT ON TABLE host_organism IS 'Host organisms for cross-reactivity exclusion';
COMMENT ON TABLE cross_reactivity IS 'Cross-reactivity tracking between pathogens';
COMMENT ON TABLE host_cross_reactivity IS 'Host cross-reactivity analysis for clinical safety';
COMMENT ON TABLE clinical_validation IS 'Clinical validation studies and performance metrics';
COMMENT ON TABLE amr_gene IS 'Antimicrobial resistance gene database';
COMMENT ON TABLE virulence_factor IS 'Virulence factor database for pathogenesis analysis';


