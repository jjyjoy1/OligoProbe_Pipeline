-- =============================================================================
-- Pathogen Detection Database Query Examples
-- =============================================================================
-- Practical examples for pathogen detection applications using the enhanced
-- PostgreSQL + BioSQL schema with taxonomic and clinical extensions

-- =============================================================================
-- 1. PATHOGEN IDENTIFICATION AND CLASSIFICATION
-- =============================================================================

-- Find all primers targeting specific pathogens
SELECT 
    poc.oligo_name,
    poc.sequence,
    poc.target_organism,
    poc.pathogen_type,
    poc.pathogenicity_level,
    poc.clinical_significance
FROM v_pathogen_oligo_complete poc
WHERE poc.target_organism ILIKE '%salmonella%'
ORDER BY poc.specificity_score DESC;

-- Get taxonomic lineage for a pathogen
SELECT * FROM get_taxonomic_lineage(590); -- Salmonella enterica

-- Find all pathogens in a specific taxonomic family
SELECT DISTINCT
    tn.name_txt as organism_name,
    pc.pathogen_type,
    pc.pathogenicity_level
FROM taxonomy_name tn
JOIN pathogen_classification pc ON tn.tax_id = pc.tax_id
WHERE tn.tax_id IN (
    SELECT DISTINCT tax_id 
    FROM get_taxonomic_lineage(543) -- Enterobacteriaceae family
)
AND tn.name_class = 'scientific name'
ORDER BY organism_name;

-- =============================================================================
-- 2. CROSS-REACTIVITY AND SPECIFICITY ANALYSIS
-- =============================================================================

-- Find primers with potential cross-reactivity issues
SELECT 
    crr.oligo_name,
    crr.total_cross_reactions,
    crr.critical_reactions,
    crr.host_reactions,
    crr.overall_risk_level
FROM v_cross_reactivity_risk crr
WHERE crr.overall_risk_level IN ('HIGH RISK', 'MEDIUM RISK')
ORDER BY crr.critical_reactions DESC, crr.total_cross_reactions DESC;

-- Detailed cross-reactivity analysis for a specific primer
SELECT 
    cr.cross_reactive_organism,
    cr.similarity_percent,
    cr.clinical_significance,
    cr.phylogenetic_distance,
    tn.name_txt as taxonomic_name
FROM cross_reactivity cr
JOIN taxonomy_name tn ON cr.cross_reactive_tax_id = tn.tax_id
WHERE cr.oligo_design_id = (
    SELECT oligo_design_id FROM v_pathogen_oligo_complete 
    WHERE oligo_name = 'SALM_16S_F1'
)
AND tn.name_class = 'scientific name'
ORDER BY cr.similarity_percent DESC;

-- Host cross-reactivity check (critical for clinical use)
SELECT 
    be.name as oligo_name,
    ho.common_name as host_organism,
    hcr.similarity_percent,
    hcr.clinical_concern_level,
    hcr.hit_location
FROM host_cross_reactivity hcr
JOIN oligo_design od ON hcr.oligo_design_id = od.oligo_design_id
JOIN bioentry be ON od.bioentry_id = be.bioentry_id
JOIN host_organism ho ON hcr.host_id = ho.host_id
WHERE hcr.clinical_concern_level IN ('critical', 'high')
ORDER BY hcr.similarity_percent DESC;

-- =============================================================================
-- 3. MULTI-PATHOGEN PANEL DESIGN
-- =============================================================================

-- Design a respiratory pathogen panel
WITH respiratory_pathogens AS (
    SELECT tax_id, name_txt as pathogen_name
    FROM taxonomy_name 
    WHERE name_txt IN (
        'Streptococcus pneumoniae',
        'Haemophilus influenzae', 
        'Mycoplasma pneumoniae',
        'Chlamydophila pneumoniae',
        'Legionella pneumophila'
    )
    AND name_class = 'scientific name'
)
SELECT 
    poc.oligo_name,
    poc.target_organism,
    poc.tm_calculated,
    poc.gc_content,
    poc.specificity_score,
    poc.target_gene
FROM v_pathogen_oligo_complete poc
JOIN respiratory_pathogens rp ON poc.target_tax_id = rp.tax_id
WHERE poc.validation_status = 'validated'
  AND poc.tm_calculated BETWEEN 58 AND 62  -- Uniform Tm for multiplex
  AND poc.specificity_score > 0.9
ORDER BY poc.target_organism, poc.specificity_score DESC;

-- Check primer compatibility for multiplex PCR
WITH panel_primers AS (
    SELECT 
        poc.oligo_name,
        poc.sequence,
        poc.tm_calculated,
        poc.target_organism
    FROM v_pathogen_oligo_complete poc
    WHERE poc.target_organism IN (
        'Escherichia coli',
        'Salmonella enterica', 
        'Listeria monocytogenes'
    )
    AND poc.oligo_type = 'forward_primer'
    AND poc.validation_status = 'validated'
)
SELECT 
    p1.oligo_name as primer1,
    p2.oligo_name as primer2,
    p1.target_organism as target1,
    p2.target_organism as target2,
    ABS(p1.tm_calculated - p2.tm_calculated) as tm_difference,
    SIMILARITY(p1.sequence, p2.sequence) as sequence_similarity
FROM panel_primers p1
CROSS JOIN panel_primers p2
WHERE p1.oligo_name < p2.oligo_name  -- Avoid duplicates
  AND (ABS(p1.tm_calculated - p2.tm_calculated) > 3  -- Tm incompatible
       OR SIMILARITY(p1.sequence, p2.sequence) > 0.7)  -- Too similar
ORDER BY tm_difference DESC, sequence_similarity DESC;

-- =============================================================================
-- 4. ANTIMICROBIAL RESISTANCE (AMR) DETECTION
-- =============================================================================

-- Find primers targeting AMR genes
SELECT 
    be.name as oligo_name,
    ag.gene_name as amr_gene,
    ag.resistance_mechanism,
    ag.antibiotic_class,
    ag.specific_antibiotics,
    od.tm_calculated,
    od.specificity_score
FROM oligo_amr_target oat
JOIN oligo_design od ON oat.oligo_design_id = od.oligo_design_id
JOIN bioentry be ON od.bioentry_id = be.bioentry_id
JOIN amr_gene ag ON oat.amr_gene_id = ag.amr_gene_id
WHERE ag.antibiotic_class = 'beta-lactam'
ORDER BY od.specificity_score DESC;

-- AMR gene prevalence analysis
SELECT 
    ag.gene_name,
    ag.resistance_mechanism,
    COUNT(oat.oligo_design_id) as primer_count,
    AVG(od.specificity_score) as avg_specificity
FROM amr_gene ag
LEFT JOIN oligo_amr_target oat ON ag.amr_gene_id = oat.amr_gene_id
LEFT JOIN oligo_design od ON oat.oligo_design_id = od.oligo_design_id
GROUP BY ag.amr_gene_id, ag.gene_name, ag.resistance_mechanism
HAVING COUNT(oat.oligo_design_id) > 0
ORDER BY primer_count DESC;

-- =============================================================================
-- 5. VIRULENCE FACTOR DETECTION
-- =============================================================================

-- Primers targeting virulence factors
SELECT 
    be.name as oligo_name,
    vf.factor_name,
    vf.factor_type,
    vf.disease_association,
    od.target_organism,
    od.specificity_score
FROM oligo_virulence_target ovt
JOIN oligo_design od ON ovt.oligo_design_id = od.oligo_design_id
JOIN bioentry be ON od.bioentry_id = be.bioentry_id
JOIN virulence_factor vf ON ovt.virulence_factor_id = vf.virulence_factor_id
WHERE vf.factor_type = 'toxin'
ORDER BY od.specificity_score DESC;

-- Comprehensive pathogenicity profile
SELECT 
    poc.target_organism,
    COUNT(DISTINCT ovt.virulence_factor_id) as virulence_factors,
    COUNT(DISTINCT oat.amr_gene_id) as resistance_genes,
    poc.pathogenicity_level,
    poc.clinical_significance
FROM v_pathogen_oligo_complete poc
LEFT JOIN oligo_virulence_target ovt ON poc.oligo_design_id = ovt.oligo_design_id
LEFT JOIN oligo_amr_target oat ON poc.oligo_design_id = oat.oligo_design_id
GROUP BY poc.target_organism, poc.pathogenicity_level, poc.clinical_significance
HAVING COUNT(DISTINCT ovt.virulence_factor_id) > 0 
    OR COUNT(DISTINCT oat.amr_gene_id) > 0
ORDER BY virulence_factors DESC, resistance_genes DESC;

-- =============================================================================
-- 6. CLINICAL VALIDATION AND PERFORMANCE
-- =============================================================================

-- Clinical performance summary
SELECT 
    cp.oligo_name,
    cp.target_organism,
    cp.avg_sensitivity,
    cp.avg_specificity,
    cp.validation_studies,
    cp.regulatory_approvals,
    cp.latest_validation
FROM v_clinical_performance cp
WHERE cp.avg_sensitivity >= 95.0 
  AND cp.avg_specificity >= 98.0
ORDER BY cp.avg_sensitivity DESC, cp.avg_specificity DESC;

-- Failed validation analysis
SELECT 
    be.name as oligo_name,
    poc.target_organism,
    cv.study_name,
    cv.sensitivity_percent,
    cv.specificity_percent,
    cv.false_positive,
    cv.false_negative,
    cv.validation_date
FROM clinical_validation cv
JOIN oligo_design od ON cv.oligo_design_id = od.oligo_design_id
JOIN bioentry be ON od.bioentry_id = be.bioentry_id
JOIN v_pathogen_oligo_complete poc ON od.oligo_design_id = poc.oligo_design_id
WHERE cv.sensitivity_percent < 95.0 OR cv.specificity_percent < 98.0
ORDER BY cv.validation_date DESC;

-- Quality control trends
SELECT 
    be.name as oligo_name,
    qc.qc_parameter,
    qc.specification_value,
    qc.tested_value,
    qc.pass_fail,
    qc.test_date
FROM diagnostic_qc qc
JOIN oligo_design od ON qc.oligo_design_id = od.oligo_design_id
JOIN bioentry be ON od.bioentry_id = be.bioentry_id
WHERE qc.test_date >= CURRENT_DATE - INTERVAL '30 days'
  AND qc.pass_fail = 'FAIL'
ORDER BY qc.test_date DESC;

-- =============================================================================
-- 7. CONTAMINATION AND INTERFERENCE ANALYSIS
-- =============================================================================

-- Common laboratory contaminant screening
WITH lab_contaminants AS (
    SELECT UNNEST(ARRAY[
        'Escherichia coli DH5alpha',
        'Pseudomonas fluorescens', 
        'Bacillus subtilis',
        'Enterococcus faecalis'
    ]) as contaminant_name
)
SELECT 
    poc.oligo_name,
    poc.target_organism,
    lc.contaminant_name,
    SIMILARITY(poc.sequence, bs_cont.seq) as similarity
FROM v_pathogen_oligo_complete poc
CROSS JOIN lab_contaminants lc
JOIN taxonomy_name tn_cont ON lc.contaminant_name = tn_cont.name_txt
JOIN bioentry be_cont ON tn_cont.tax_id = be_cont.tax_id  
JOIN biosequence bs_cont ON be_cont.bioentry_id = bs_cont.bioentry_id
WHERE poc.sequence % bs_cont.seq
  AND SIMILARITY(poc.sequence, bs_cont.seq) > 0.8
  AND tn_cont.name_class = 'scientific name'
ORDER BY similarity DESC;

-- Environmental interference check
SELECT 
    poc.oligo_name,
    poc.target_organism,
    COUNT(DISTINCT oth.chromosome) as chromosomal_hits,
    AVG(oth.identity_percent) as avg_identity,
    MAX(oth.identity_percent) as max_identity
FROM v_pathogen_oligo_complete poc
JOIN off_target_hit oth ON poc.oligo_design_id = oth.oligo_design_id
WHERE oth.gene_symbol IS NOT NULL  -- Annotated genes
  AND oth.identity_percent > 85.0
GROUP BY poc.oligo_name, poc.target_organism
HAVING COUNT(DISTINCT oth.chromosome) > 5  -- Multiple chromosomal locations
ORDER BY max_identity DESC;

-- =============================================================================
-- 8. PHYLOGENETIC AND EVOLUTIONARY ANALYSIS
-- =============================================================================

-- Phylogenetic distance matrix for pathogen panel
WITH pathogen_pairs AS (
    SELECT DISTINCT
        p1.tax_id as tax_id1,
        p1.name_txt as pathogen1,
        p2.tax_id as tax_id2, 
        p2.name_txt as pathogen2
    FROM taxonomy_name p1
    CROSS JOIN taxonomy_name p2
    WHERE p1.tax_id < p2.tax_id
      AND p1.name_class = 'scientific name'
      AND p2.name_class = 'scientific name'
      AND p1.name_txt IN (
          'Escherichia coli',
          'Salmonella enterica',
          'Shigella dysenteriae',
          'Klebsiella pneumoniae'
      )
      AND p2.name_txt IN (
          'Escherichia coli',
          'Salmonella enterica', 
          'Shigella dysenteriae',
          'Klebsiella pneumoniae'
      )
)
SELECT 
    pathogen1,
    pathogen2,
    calculate_phylogenetic_distance(tax_id1, tax_id2) as phylo_distance
FROM pathogen_pairs
ORDER BY phylo_distance;

-- Evolutionary conservation analysis
SELECT 
    tg.gene_name,
    tg.conservation_level,
    COUNT(od.oligo_design_id) as primer_count,
    AVG(od.specificity_score) as avg_specificity,
    STDDEV(od.specificity_score) as specificity_stddev
FROM target_gene tg
JOIN oligo_design od ON tg.gene_name = od.target_gene
WHERE tg.pathogen_type = 'bacteria'
GROUP BY tg.gene_name, tg.conservation_level
ORDER BY conservation_level, avg_specificity DESC;

-- =============================================================================
-- 9. DIAGNOSTIC PANEL OPTIMIZATION
-- =============================================================================

-- Optimal primer selection for foodborne pathogen panel
WITH foodborne_targets AS (
    SELECT UNNEST(ARRAY[
        'Salmonella enterica',
        'Listeria monocytogenes',
        'Escherichia coli O157:H7',
        'Campylobacter jejuni',
        'Staphylococcus aureus'
    ]) as target_pathogen
),
ranked_primers AS (
    SELECT 
        poc.oligo_name,
        poc.target_organism,
        poc.tm_calculated,
        poc.specificity_score,
        poc.off_target_count,
        ROW_NUMBER() OVER (
            PARTITION BY poc.target_organism 
            ORDER BY poc.specificity_score DESC, poc.off_target_count ASC
        ) as primer_rank
    FROM v_pathogen_oligo_complete poc
    JOIN foodborne_targets ft ON poc.target_organism = ft.target_pathogen
    WHERE poc.validation_status = 'validated'
      AND poc.tm_calculated BETWEEN 58 AND 62
      AND poc.oligo_type = 'forward_primer'
)
SELECT * FROM ranked_primers WHERE primer_rank <= 3;

-- Panel performance prediction
SELECT 
    COUNT(DISTINCT poc.target_organism) as pathogens_covered,
    AVG(poc.specificity_score) as avg_panel_specificity,
    MIN(poc.specificity_score) as min_specificity,
    MAX(ABS(poc.tm_calculated - AVG(poc.tm_calculated) OVER())) as max_tm_deviation,
    SUM(CASE WHEN crr.overall_risk_level = 'HIGH RISK' THEN 1 ELSE 0 END) as high_risk_primers
FROM v_pathogen_oligo_complete poc
LEFT JOIN v_cross_reactivity_risk crr ON poc.oligo_design_id = crr.oligo_design_id
WHERE poc.target_organism IN (
    'Salmonella enterica',
    'Listeria monocytogenes', 
    'Escherichia coli O157:H7'
)
AND poc.validation_status = 'validated';

-- =============================================================================
-- 10. REGULATORY AND COMPLIANCE QUERIES
-- =============================================================================

-- FDA-ready primer documentation
SELECT 
    poc.oligo_name,
    poc.target_organism,
    poc.pathogenicity_level,
    cp.avg_sensitivity,
    cp.avg_specificity,
    cp.validation_studies,
    cp.regulatory_approvals,
    STRING_AGG(cv.reference_method, ', ') as validation_methods
FROM v_pathogen_oligo_complete poc
JOIN v_clinical_performance cp ON poc.oligo_design_id = cp.oligo_design_id
LEFT JOIN clinical_validation cv ON poc.oligo_design_id = cv.oligo_design_id
WHERE cp.avg_sensitivity >= 95.0
  AND cp.avg_specificity >= 98.0
  AND cp.validation_studies >= 3
GROUP BY poc.oligo_name, poc.target_organism, poc.pathogenicity_level,
         cp.avg_sensitivity, cp.avg_specificity, cp.validation_studies, 
         cp.regulatory_approvals
ORDER BY cp.avg_sensitivity DESC;

-- Audit trail for validated primers
SELECT 
    be.name as oligo_name,
    od.created_date as design_date,
    od.validation_status,
    od.modified_date as last_modified,
    cv.validation_date,
    cv.study_name,
    cv.regulatory_approval,
    od.created_by
FROM oligo_design od
JOIN bioentry be ON od.bioentry_id = be.bioentry_id
LEFT JOIN clinical_validation cv ON od.oligo_design_id = cv.oligo_design_id
WHERE od.validation_status IN ('validated', 'approved')
ORDER BY od.created_date DESC;

-- =============================================================================
-- 11. PERFORMANCE MONITORING AND MAINTENANCE
-- =============================================================================

-- Database performance metrics
SELECT 
    'Pathogen Oligos' as category,
    COUNT(*) as total_count,
    COUNT(*) FILTER (WHERE validation_status = 'validated') as validated_count,
    AVG(specificity_score) as avg_specificity,
    COUNT(*) FILTER (WHERE off_target_count = 0) as perfect_specificity_count
FROM v_pathogen_oligo_complete

UNION ALL

SELECT 
    'Cross-Reactivity Issues' as category,
    COUNT(*) as total_count,
    COUNT(*) FILTER (WHERE overall_risk_level = 'HIGH RISK') as high_risk_count,
    AVG(total_cross_reactions::numeric) as avg_cross_reactions,
    COUNT(*) FILTER (WHERE critical_reactions > 0) as critical_issues_count
FROM v_cross_reactivity_risk;

-- Taxonomy coverage analysis
SELECT 
    tn.rank,
    COUNT(DISTINCT tn.tax_id) as total_taxa,
    COUNT(DISTINCT poc.target_tax_id) as covered_taxa,
    ROUND(
        COUNT(DISTINCT poc.target_tax_id)::numeric / 
        COUNT(DISTINCT tn.tax_id) * 100, 
        2
    ) as coverage_percentage
FROM taxonomy_node tn
LEFT JOIN v_pathogen_oligo_complete poc ON tn.tax_id = poc.target_tax_id
WHERE tn.rank IN ('species', 'genus', 'family', 'order')
GROUP BY tn.rank
ORDER BY coverage_percentage DESC;

-- Recent validation activity
SELECT 
    DATE_TRUNC('month', cv.validation_date) as validation_month,
    COUNT(*) as validations_completed,
    AVG(cv.sensitivity_percent) as avg_sensitivity,
    AVG(cv.specificity_percent) as avg_specificity,
    COUNT(*) FILTER (WHERE cv.regulatory_approval IS NOT NULL) as regulatory_submissions
FROM clinical_validation cv
WHERE cv.validation_date >= CURRENT_DATE - INTERVAL '12 months'
GROUP BY DATE_TRUNC('month', cv.validation_date)
ORDER BY validation_month DESC;
