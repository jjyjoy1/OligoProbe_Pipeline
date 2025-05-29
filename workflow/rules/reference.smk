# =============================================================================
# Reference Sequence Acquisition & Preparation Rules
# =============================================================================

# Download reference genome
rule download_reference:
    output:
        fasta="resources/{build}.fa.gz"
    params:
        url=config["urls"]["reference_fasta"]
    log:
        "logs/download_reference_{build}.log"
    threads: get_threads("download")
    resources:
        mem_mb=lambda w: get_memory("download").replace("G", "000").replace("M", "")
    shell:
        """
        echo "Downloading reference genome from {params.url}" > {log}
        wget -O {output.fasta} {params.url} 2>> {log}
        echo "Download completed: $(date)" >> {log}
        """

# Uncompress and normalize reference
rule prepare_reference:
    input:
        "resources/{build}.fa.gz"
    output:
        "resources/{build}.fa"
    log:
        "logs/prepare_reference_{build}.log"
    shell:
        """
        module load {config[modules][samtools]}
        
        echo "Uncompressing and normalizing reference..." > {log}
        
        # Uncompress and normalize headers
        gunzip -c {input} | \
        sed 's/^>.*chromosome />chr/; s/^>.*mitochondrion/>chrM/' | \
        sed 's/, GRCh38.*$//' > {output} 2>> {log}
        
        # Index the FASTA
        samtools faidx {output} 2>> {log}
        
        echo "Reference preparation completed: $(date)" >> {log}
        """

# Download VCF files for variant-aware reference
rule download_vcf:
    output:
        vcf="resources/gnomad_chr{chr}.vcf.gz",
        index="resources/gnomad_chr{chr}.vcf.gz.tbi"
    params:
        url=lambda w: config["urls"]["gnomad_vcf"].format(chr=w.chr)
    log:
        "logs/download_vcf_chr{chr}.log"
    threads: get_threads("download")
    shell:
        """
        echo "Downloading VCF for chromosome {wildcards.chr}" > {log}
        
        # Download VCF
        wget -O {output.vcf} {params.url} 2>> {log}
        
        # Download index
        wget -O {output.index} {params.url}.tbi 2>> {log}
        
        echo "VCF download completed: $(date)" >> {log}
        """

# Merge VCF files
rule merge_vcf:
    input:
        vcfs=expand("resources/gnomad_chr{chr}.vcf.gz", chr=config["chromosomes"]),
        indices=expand("resources/gnomad_chr{chr}.vcf.gz.tbi", chr=config["chromosomes"])
    output:
        "resources/gnomad_merged.vcf.gz"
    log:
        "logs/merge_vcf.log"
    threads: get_threads("consensus")
    shell:
        """
        module load {config[modules][bcftools]}
        
        echo "Merging VCF files..." > {log}
        
        # Create file list
        ls resources/gnomad_chr*.vcf.gz > vcf_list.txt
        
        # Merge VCFs
        bcftools concat -f vcf_list.txt -O z -o {output} --threads {threads} 2>> {log}
        
        # Index merged VCF
        bcftools index {output} 2>> {log}
        
        # Cleanup
        rm vcf_list.txt
        
        echo "VCF merge completed: $(date)" >> {log}
        """

# Create variant-aware reference for panel design
rule create_panel_reference:
    input:
        fasta="resources/{build}.fa",
        vcf="resources/gnomad_merged.vcf.gz"
    output:
        "resources/{build}_panel.fa"
    log:
        "logs/create_panel_reference_{build}.log"
    threads: get_threads("consensus")
    resources:
        mem_mb=lambda w: get_memory("consensus").replace("G", "000").replace("M", "")
    shell:
        """
        module load {config[modules][bcftools]}
        module load {config[modules][samtools]}
        
        echo "Creating variant-aware reference..." > {log}
        
        # Create consensus sequence with IUPAC codes for variants
        bcftools consensus -f {input.fasta} {input.vcf} \
            --iupac-coding \
            --threads {threads} \
            -o {output} 2>> {log}
        
        # Index the new reference
        samtools faidx {output} 2>> {log}
        
        echo "Panel reference created: $(date)" >> {log}
        """

# Create test dataset (chr22 only)
rule create_test_data:
    input:
        "resources/{build}.fa"
    output:
        "resources/test_chr22.fa"
    log:
        "logs/create_test_data_{build}.log"
    shell:
        """
        module load {config[modules][samtools]}
        
        echo "Creating test dataset (chr22)..." > {log}
        
        # Extract chromosome 22
        samtools faidx {input} chr22 > {output} 2>> {log}
        
        # Index test file
        samtools faidx {output} 2>> {log}
        
        echo "Test data created: $(date)" >> {log}
        """

# Download annotation (optional)
rule download_annotation:
    output:
        "resources/annotation.gtf.gz"
    params:
        url=config["urls"]["annotation_gtf"]
    log:
        "logs/download_annotation.log"
    shell:
        """
        echo "Downloading annotation from {params.url}" > {log}
        wget -O {output} {params.url} 2>> {log}
        echo "Annotation download completed: $(date)" >> {log}
        """

# Process annotation for exon masking
rule process_annotation:
    input:
        "resources/annotation.gtf.gz"
    output:
        bed="resources/exons.bed",
        introns="resources/introns.bed"
    log:
        "logs/process_annotation.log"
    shell:
        """
        echo "Processing annotation..." > {log}
        
        # Extract exons
        gunzip -c {input} | \
        awk '$3=="exon" {{print $1"\\t"$4-1"\\t"$5"\\t"$10"\\t0\\t"$7}}' | \
        sed 's/[";]//g' | sort -k1,1 -k2,2n > {output.bed} 2>> {log}
        
        # Create intron mask (complement of exons)
        bedtools complement -i {output.bed} -g resources/hg38.fa.fai > {output.introns} 2>> {log}
        
        echo "Annotation processing completed: $(date)" >> {log}
        """

# Quality check for reference files
rule qc_reference:
    input:
        pcr_ref="resources/{build}.fa",
        panel_ref="resources/{build}_panel.fa"
    output:
        "results/qc/reference_qc_{build}.txt"
    log:
        "logs/qc_reference_{build}.log"
    shell:
        """
        echo "Performing reference quality checks..." > {log}
        
        mkdir -p results/qc
        
        # Check file sizes
        PCR_SIZE=$(stat -c%s {input.pcr_ref})
        PANEL_SIZE=$(stat -c%s {input.panel_ref})
        
        echo "Reference QC Report" > {output}
        echo "Generated: $(date)" >> {output}
        echo "===================" >> {output}
        echo "" >> {output}
        echo "File sizes:" >> {output}
        echo "PCR reference: $PCR_SIZE bytes" >> {output}
        echo "Panel reference: $PANEL_SIZE bytes" >> {output}
        echo "" >> {output}
        
        # Check sequence counts
        PCR_SEQS=$(grep -c "^>" {input.pcr_ref})
        PANEL_SEQS=$(grep -c "^>" {input.panel_ref})
        
        echo "Sequence counts:" >> {output}
        echo "PCR reference: $PCR_SEQS sequences" >> {output}
        echo "Panel reference: $PANEL_SEQS sequences" >> {output}
        echo "" >> {output}
        
        # Check for required minimum sizes
        MIN_SIZE={config[qc][min_file_sizes][fasta]}
        if [ $PCR_SIZE -lt $MIN_SIZE ]; then
            echo "ERROR: PCR reference too small ($PCR_SIZE < $MIN_SIZE)" >> {output}
            exit 1
        fi
        
        if [ $PANEL_SIZE -lt $MIN_SIZE ]; then
            echo "ERROR: Panel reference too small ($PANEL_SIZE < $MIN_SIZE)" >> {output}
            exit 1
        fi
        
        echo "QC PASSED: All reference files meet minimum requirements" >> {output}
        echo "Reference QC completed: $(date)" >> {log}
        """

# Cleanup intermediate files
rule cleanup_reference:
    input:
        "resources/{build}.fa",
        "resources/{build}_panel.fa"
    output:
        touch("resources/.cleanup_done_{build}")
    shell:
        """
        # Remove compressed files to save space
        rm -f resources/*.fa.gz
        rm -f resources/gnomad_chr*.vcf.gz*
        echo "Reference cleanup completed"
        """


