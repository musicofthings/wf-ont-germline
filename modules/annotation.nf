/*
========================================================================================
    VARIANT ANNOTATION MODULES
    
    Comprehensive annotation with clinically relevant databases:
    1. VEP (Ensembl Variant Effect Predictor) - Base functional annotation
    2. ClinVar - Clinical significance
    3. gnomAD - Population allele frequencies
    4. dbSNP - Variant identifiers
    5. CADD - Deleteriousness scores
    6. REVEL - Missense pathogenicity
    7. SpliceAI - Splice site predictions
    8. PharmGKB - Pharmacogenomics
    9. COSMIC - Cancer mutations (germline predisposition)
    10. OMIM - Disease associations
========================================================================================
*/

process VEP_ANNOTATE {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/annotation/vep", mode: 'copy'
    
    container 'ensemblorg/ensembl-vep:release_110.1'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi)
    path reference
    path vep_cache
    
    output:
    tuple val(sample_id), path("${sample_id}.vep.vcf.gz"), path("${sample_id}.vep.vcf.gz.tbi"), emit: vcf
    tuple val(sample_id), path("${sample_id}.vep_summary.html"), emit: summary
    path "versions.yml", emit: versions
    
    script:
    def genome_assembly = params.genome_assembly ?: 'GRCh38'
    def species = params.species ?: 'homo_sapiens'
    def extra_args = params.vep_args ?: ""
    """
    echo "[INFO] =============================================="
    echo "[INFO] Starting VEP annotation"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] Assembly: ${genome_assembly}"
    echo "[INFO] Cache: ${vep_cache}"
    echo "[INFO] =============================================="
    
    vep \\
        --input_file ${vcf} \\
        --output_file ${sample_id}.vep.vcf \\
        --vcf \\
        --cache \\
        --dir_cache ${vep_cache} \\
        --species ${species} \\
        --assembly ${genome_assembly} \\
        --fasta ${reference} \\
        --fork ${task.cpus} \\
        --offline \\
        --everything \\
        --force_overwrite \\
        --stats_file ${sample_id}.vep_summary.html \\
        --warning_file ${sample_id}.vep_warnings.txt \\
        --canonical \\
        --hgvs \\
        --symbol \\
        --numbers \\
        --domains \\
        --regulatory \\
        --biotype \\
        --protein \\
        --tsl \\
        --appris \\
        --gene_phenotype \\
        --af \\
        --af_1kg \\
        --af_gnomade \\
        --af_gnomadg \\
        --max_af \\
        --pubmed \\
        --variant_class \\
        --mane \\
        ${extra_args} \\
        2>&1 | tee vep.log
    
    # Compress and index
    bgzip -c ${sample_id}.vep.vcf > ${sample_id}.vep.vcf.gz
    tabix -p vcf ${sample_id}.vep.vcf.gz
    
    echo "[INFO] VEP annotation completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(vep --help 2>&1 | grep "Versions:" -A 1 | tail -1 | sed 's/.*ensembl-vep : //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.vep.vcf.gz
    touch ${sample_id}.vep.vcf.gz.tbi
    touch ${sample_id}.vep_summary.html
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: 110.1
    END_VERSIONS
    """
}

process CLINVAR_ANNOTATE {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/annotation/clinvar", mode: 'copy'
    
    container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi)
    path clinvar_vcf
    path clinvar_tbi
    
    output:
    tuple val(sample_id), path("${sample_id}.clinvar.vcf.gz"), path("${sample_id}.clinvar.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    """
    echo "[INFO] =============================================="
    echo "[INFO] Adding ClinVar annotations"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] ClinVar database: ${clinvar_vcf}"
    echo "[INFO] =============================================="
    
    # Annotate with ClinVar
    bcftools annotate \\
        --threads ${task.cpus} \\
        --annotations ${clinvar_vcf} \\
        --columns INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT,INFO/CLNVC,INFO/CLNVI,INFO/CLNDISDB,INFO/CLNHGVS,INFO/GENEINFO,INFO/MC,INFO/ORIGIN,INFO/RS \\
        --output-type z \\
        --output ${sample_id}.clinvar.vcf.gz \\
        ${vcf}
    
    tabix -p vcf ${sample_id}.clinvar.vcf.gz
    
    echo "[INFO] ClinVar annotation completed"
    echo "[INFO] Pathogenic/Likely pathogenic variants:"
    bcftools query -f '%CHROM\\t%POS\\t%CLNSIG\\n' ${sample_id}.clinvar.vcf.gz | grep -E "Pathogenic|Likely_pathogenic" | wc -l || echo "0"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.clinvar.vcf.gz
    touch ${sample_id}.clinvar.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.18
    END_VERSIONS
    """
}

process GNOMAD_ANNOTATE {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/annotation/gnomad", mode: 'copy'
    
    container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi)
    path gnomad_vcf
    path gnomad_tbi
    
    output:
    tuple val(sample_id), path("${sample_id}.gnomad.vcf.gz"), path("${sample_id}.gnomad.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    """
    echo "[INFO] =============================================="
    echo "[INFO] Adding gnomAD population frequencies"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] gnomAD database: ${gnomad_vcf}"
    echo "[INFO] =============================================="
    
    # Annotate with gnomAD
    bcftools annotate \\
        --threads ${task.cpus} \\
        --annotations ${gnomad_vcf} \\
        --columns INFO/gnomAD_AF,INFO/gnomAD_AF_popmax,INFO/gnomAD_AF_afr,INFO/gnomAD_AF_amr,INFO/gnomAD_AF_asj,INFO/gnomAD_AF_eas,INFO/gnomAD_AF_fin,INFO/gnomAD_AF_nfe,INFO/gnomAD_AF_sas,INFO/gnomAD_nhomalt \\
        --output-type z \\
        --output ${sample_id}.gnomad.vcf.gz \\
        ${vcf}
    
    tabix -p vcf ${sample_id}.gnomad.vcf.gz
    
    echo "[INFO] gnomAD annotation completed"
    echo "[INFO] Rare variants (AF < 0.01):"
    bcftools query -f '%gnomAD_AF\\n' ${sample_id}.gnomad.vcf.gz | awk '\$1 < 0.01 || \$1 == "."' | wc -l || echo "N/A"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.gnomad.vcf.gz
    touch ${sample_id}.gnomad.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.18
    END_VERSIONS
    """
}

process DBSNP_ANNOTATE {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/annotation/dbsnp", mode: 'copy'
    
    container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi)
    path dbsnp_vcf
    path dbsnp_tbi
    
    output:
    tuple val(sample_id), path("${sample_id}.dbsnp.vcf.gz"), path("${sample_id}.dbsnp.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    """
    echo "[INFO] =============================================="
    echo "[INFO] Adding dbSNP rsIDs"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] dbSNP database: ${dbsnp_vcf}"
    echo "[INFO] =============================================="
    
    # Annotate with dbSNP rsIDs
    bcftools annotate \\
        --threads ${task.cpus} \\
        --annotations ${dbsnp_vcf} \\
        --columns ID \\
        --output-type z \\
        --output ${sample_id}.dbsnp.vcf.gz \\
        ${vcf}
    
    tabix -p vcf ${sample_id}.dbsnp.vcf.gz
    
    echo "[INFO] dbSNP annotation completed"
    echo "[INFO] Variants with rsIDs:"
    bcftools query -f '%ID\\n' ${sample_id}.dbsnp.vcf.gz | grep -v "^\\." | wc -l || echo "0"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.dbsnp.vcf.gz
    touch ${sample_id}.dbsnp.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.18
    END_VERSIONS
    """
}

process CADD_ANNOTATE {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/annotation/cadd", mode: 'copy'
    
    container 'quay.io/biocontainers/vcfanno:0.3.3--h9ee0642_0'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi)
    path cadd_snvs
    path cadd_snvs_tbi
    path cadd_indels
    path cadd_indels_tbi
    
    output:
    tuple val(sample_id), path("${sample_id}.cadd.vcf.gz"), path("${sample_id}.cadd.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    """
    echo "[INFO] =============================================="
    echo "[INFO] Adding CADD scores"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] =============================================="
    
    # Create vcfanno config for CADD
    cat > cadd_config.toml << 'EOF'
[[annotation]]
file = "${cadd_snvs}"
fields = ["RawScore", "PHRED"]
names = ["CADD_RAW", "CADD_PHRED"]
ops = ["self", "self"]

[[annotation]]
file = "${cadd_indels}"
fields = ["RawScore", "PHRED"]
names = ["CADD_RAW_INDEL", "CADD_PHRED_INDEL"]
ops = ["self", "self"]
EOF
    
    # Run vcfanno
    vcfanno \\
        -p ${task.cpus} \\
        cadd_config.toml \\
        ${vcf} \\
        | bgzip -c > ${sample_id}.cadd.vcf.gz
    
    tabix -p vcf ${sample_id}.cadd.vcf.gz
    
    echo "[INFO] CADD annotation completed"
    echo "[INFO] Variants with CADD PHRED >= 20 (potentially deleterious):"
    bcftools query -f '%CADD_PHRED\\n' ${sample_id}.cadd.vcf.gz | awk '\$1 >= 20' | wc -l || echo "0"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfanno: \$(vcfanno 2>&1 | head -1 | sed 's/vcfanno version //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.cadd.vcf.gz
    touch ${sample_id}.cadd.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfanno: 0.3.3
    END_VERSIONS
    """
}

process REVEL_ANNOTATE {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/annotation/revel", mode: 'copy'
    
    container 'quay.io/biocontainers/vcfanno:0.3.3--h9ee0642_0'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi)
    path revel_file
    path revel_tbi
    
    output:
    tuple val(sample_id), path("${sample_id}.revel.vcf.gz"), path("${sample_id}.revel.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    """
    echo "[INFO] =============================================="
    echo "[INFO] Adding REVEL missense pathogenicity scores"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] =============================================="
    
    # Create vcfanno config for REVEL
    cat > revel_config.toml << 'EOF'
[[annotation]]
file = "${revel_file}"
fields = ["REVEL"]
names = ["REVEL_score"]
ops = ["self"]
EOF
    
    # Run vcfanno
    vcfanno \\
        -p ${task.cpus} \\
        revel_config.toml \\
        ${vcf} \\
        | bgzip -c > ${sample_id}.revel.vcf.gz
    
    tabix -p vcf ${sample_id}.revel.vcf.gz
    
    echo "[INFO] REVEL annotation completed"
    echo "[INFO] Variants with REVEL >= 0.5 (likely pathogenic missense):"
    bcftools query -f '%REVEL_score\\n' ${sample_id}.revel.vcf.gz | awk '\$1 >= 0.5' | wc -l || echo "0"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfanno: \$(vcfanno 2>&1 | head -1 | sed 's/vcfanno version //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.revel.vcf.gz
    touch ${sample_id}.revel.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfanno: 0.3.3
    END_VERSIONS
    """
}

process SPLICEAI_ANNOTATE {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/annotation/spliceai", mode: 'copy'
    
    container 'quay.io/biocontainers/vcfanno:0.3.3--h9ee0642_0'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi)
    path spliceai_snvs
    path spliceai_snvs_tbi
    path spliceai_indels
    path spliceai_indels_tbi
    
    output:
    tuple val(sample_id), path("${sample_id}.spliceai.vcf.gz"), path("${sample_id}.spliceai.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    """
    echo "[INFO] =============================================="
    echo "[INFO] Adding SpliceAI predictions"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] =============================================="
    
    # Create vcfanno config for SpliceAI
    cat > spliceai_config.toml << 'EOF'
[[annotation]]
file = "${spliceai_snvs}"
fields = ["SpliceAI"]
names = ["SpliceAI_SNV"]
ops = ["self"]

[[annotation]]
file = "${spliceai_indels}"
fields = ["SpliceAI"]
names = ["SpliceAI_INDEL"]
ops = ["self"]
EOF
    
    # Run vcfanno
    vcfanno \\
        -p ${task.cpus} \\
        spliceai_config.toml \\
        ${vcf} \\
        | bgzip -c > ${sample_id}.spliceai.vcf.gz
    
    tabix -p vcf ${sample_id}.spliceai.vcf.gz
    
    echo "[INFO] SpliceAI annotation completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfanno: \$(vcfanno 2>&1 | head -1 | sed 's/vcfanno version //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.spliceai.vcf.gz
    touch ${sample_id}.spliceai.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfanno: 0.3.3
    END_VERSIONS
    """
}

process PHARMGKB_ANNOTATE {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/annotation/pharmgkb", mode: 'copy'
    
    container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi)
    path pharmgkb_vcf
    path pharmgkb_tbi
    
    output:
    tuple val(sample_id), path("${sample_id}.pharmgkb.vcf.gz"), path("${sample_id}.pharmgkb.vcf.gz.tbi"), emit: vcf
    tuple val(sample_id), path("${sample_id}.pgx_variants.tsv"), emit: pgx_report
    path "versions.yml", emit: versions
    
    script:
    """
    echo "[INFO] =============================================="
    echo "[INFO] Adding PharmGKB pharmacogenomics annotations"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] =============================================="
    
    # Annotate with PharmGKB
    bcftools annotate \\
        --threads ${task.cpus} \\
        --annotations ${pharmgkb_vcf} \\
        --columns INFO/PGX_GENE,INFO/PGX_DRUG,INFO/PGX_PHENOTYPE,INFO/PGX_LEVEL \\
        --output-type z \\
        --output ${sample_id}.pharmgkb.vcf.gz \\
        ${vcf}
    
    tabix -p vcf ${sample_id}.pharmgkb.vcf.gz
    
    # Extract PGx variants report
    echo -e "CHROM\\tPOS\\tID\\tREF\\tALT\\tGENE\\tDRUG\\tPHENOTYPE\\tLEVEL" > ${sample_id}.pgx_variants.tsv
    bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%PGX_GENE\\t%PGX_DRUG\\t%PGX_PHENOTYPE\\t%PGX_LEVEL\\n' \\
        -i 'INFO/PGX_GENE!="."' \\
        ${sample_id}.pharmgkb.vcf.gz >> ${sample_id}.pgx_variants.tsv || true
    
    echo "[INFO] PharmGKB annotation completed"
    echo "[INFO] Pharmacogenomic variants found:"
    tail -n +2 ${sample_id}.pgx_variants.tsv | wc -l
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.pharmgkb.vcf.gz
    touch ${sample_id}.pharmgkb.vcf.gz.tbi
    touch ${sample_id}.pgx_variants.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.18
    END_VERSIONS
    """
}

process COSMIC_ANNOTATE {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/annotation/cosmic", mode: 'copy'
    
    container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi)
    path cosmic_vcf
    path cosmic_tbi
    
    output:
    tuple val(sample_id), path("${sample_id}.cosmic.vcf.gz"), path("${sample_id}.cosmic.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    """
    echo "[INFO] =============================================="
    echo "[INFO] Adding COSMIC cancer mutation annotations"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] =============================================="
    
    # Annotate with COSMIC
    bcftools annotate \\
        --threads ${task.cpus} \\
        --annotations ${cosmic_vcf} \\
        --columns INFO/COSMIC_ID,INFO/COSMIC_CNT,INFO/COSMIC_GENE,INFO/COSMIC_TISSUE \\
        --output-type z \\
        --output ${sample_id}.cosmic.vcf.gz \\
        ${vcf}
    
    tabix -p vcf ${sample_id}.cosmic.vcf.gz
    
    echo "[INFO] COSMIC annotation completed"
    echo "[INFO] Variants in COSMIC:"
    bcftools query -f '%COSMIC_ID\\n' ${sample_id}.cosmic.vcf.gz | grep -v "^\\." | wc -l || echo "0"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.cosmic.vcf.gz
    touch ${sample_id}.cosmic.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.18
    END_VERSIONS
    """
}

process OMIM_ANNOTATE {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/annotation/omim", mode: 'copy'
    
    container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi)
    path omim_bed
    
    output:
    tuple val(sample_id), path("${sample_id}.omim.vcf.gz"), path("${sample_id}.omim.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    """
    echo "[INFO] =============================================="
    echo "[INFO] Adding OMIM disease gene annotations"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] =============================================="
    
    # Annotate with OMIM gene regions
    bcftools annotate \\
        --threads ${task.cpus} \\
        --annotations ${omim_bed} \\
        --columns CHROM,FROM,TO,OMIM_GENE,OMIM_MIM,OMIM_PHENOTYPE \\
        --header-lines <(echo '##INFO=<ID=OMIM_GENE,Number=1,Type=String,Description="OMIM gene symbol">') \\
        --output-type z \\
        --output ${sample_id}.omim.vcf.gz \\
        ${vcf}
    
    tabix -p vcf ${sample_id}.omim.vcf.gz
    
    echo "[INFO] OMIM annotation completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.omim.vcf.gz
    touch ${sample_id}.omim.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.18
    END_VERSIONS
    """
}

process MERGE_ANNOTATIONS {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/annotation/final", mode: 'copy'
    
    container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi)
    
    output:
    tuple val(sample_id), path("${sample_id}.fully_annotated.vcf.gz"), path("${sample_id}.fully_annotated.vcf.gz.tbi"), emit: vcf
    tuple val(sample_id), path("${sample_id}.annotation_summary.txt"), emit: summary
    path "versions.yml", emit: versions
    
    script:
    """
    echo "[INFO] =============================================="
    echo "[INFO] Finalizing annotations"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] =============================================="
    
    # Copy final annotated VCF
    cp ${vcf} ${sample_id}.fully_annotated.vcf.gz
    cp ${vcf_tbi} ${sample_id}.fully_annotated.vcf.gz.tbi
    
    # Generate annotation summary
    cat > ${sample_id}.annotation_summary.txt << EOF
================================================================================
VARIANT ANNOTATION SUMMARY
Sample: ${sample_id}
Date: \$(date)
================================================================================

TOTAL VARIANTS:
\$(bcftools view -H ${vcf} | wc -l)

VARIANT TYPES:
\$(bcftools stats ${vcf} | grep "^SN" | cut -f3-)

CLINVAR CLASSIFICATIONS:
\$(bcftools query -f '%CLNSIG\\n' ${vcf} 2>/dev/null | sort | uniq -c | sort -rn | head -10 || echo "N/A")

POPULATION FREQUENCIES (gnomAD):
- Ultra-rare (AF < 0.0001): \$(bcftools query -f '%gnomAD_AF\\n' ${vcf} 2>/dev/null | awk '\$1 < 0.0001 || \$1 == "."' | wc -l || echo "N/A")
- Rare (AF < 0.01): \$(bcftools query -f '%gnomAD_AF\\n' ${vcf} 2>/dev/null | awk '\$1 < 0.01 || \$1 == "."' | wc -l || echo "N/A")
- Common (AF >= 0.01): \$(bcftools query -f '%gnomAD_AF\\n' ${vcf} 2>/dev/null | awk '\$1 >= 0.01' | wc -l || echo "N/A")

CADD SCORES:
- CADD >= 30 (highly deleterious): \$(bcftools query -f '%CADD_PHRED\\n' ${vcf} 2>/dev/null | awk '\$1 >= 30' | wc -l || echo "N/A")
- CADD >= 20 (potentially deleterious): \$(bcftools query -f '%CADD_PHRED\\n' ${vcf} 2>/dev/null | awk '\$1 >= 20' | wc -l || echo "N/A")

PHARMACOGENOMIC VARIANTS:
\$(bcftools query -f '%PGX_GENE\\n' ${vcf} 2>/dev/null | grep -v "^\\." | wc -l || echo "N/A")

================================================================================
EOF
    
    cat ${sample_id}.annotation_summary.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.fully_annotated.vcf.gz
    touch ${sample_id}.fully_annotated.vcf.gz.tbi
    touch ${sample_id}.annotation_summary.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.18
    END_VERSIONS
    """
}

process FILTER_CLINICAL_VARIANTS {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/annotation/clinical", mode: 'copy'
    
    container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi)
    
    output:
    tuple val(sample_id), path("${sample_id}.pathogenic.vcf.gz"), path("${sample_id}.pathogenic.vcf.gz.tbi"), emit: pathogenic_vcf
    tuple val(sample_id), path("${sample_id}.rare_damaging.vcf.gz"), path("${sample_id}.rare_damaging.vcf.gz.tbi"), emit: rare_damaging_vcf
    tuple val(sample_id), path("${sample_id}.clinical_report.tsv"), emit: clinical_report
    path "versions.yml", emit: versions
    
    script:
    def af_threshold = params.rare_af_threshold ?: 0.01
    def cadd_threshold = params.cadd_threshold ?: 20
    """
    echo "[INFO] =============================================="
    echo "[INFO] Filtering clinical variants"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] AF threshold: ${af_threshold}"
    echo "[INFO] CADD threshold: ${cadd_threshold}"
    echo "[INFO] =============================================="
    
    # Filter pathogenic/likely pathogenic variants
    bcftools view \\
        -i 'INFO/CLNSIG~"Pathogenic" || INFO/CLNSIG~"Likely_pathogenic"' \\
        -Oz -o ${sample_id}.pathogenic.vcf.gz \\
        ${vcf} || touch ${sample_id}.pathogenic.vcf.gz
    
    tabix -f -p vcf ${sample_id}.pathogenic.vcf.gz || true
    
    # Filter rare + damaging variants
    bcftools view \\
        -i '(INFO/gnomAD_AF < ${af_threshold} || INFO/gnomAD_AF=".") && INFO/CADD_PHRED >= ${cadd_threshold}' \\
        -Oz -o ${sample_id}.rare_damaging.vcf.gz \\
        ${vcf} || touch ${sample_id}.rare_damaging.vcf.gz
    
    tabix -f -p vcf ${sample_id}.rare_damaging.vcf.gz || true
    
    # Generate clinical report TSV
    echo -e "CHROM\\tPOS\\tID\\tREF\\tALT\\tGENE\\tCONSEQUENCE\\tCLINVAR\\tgnomAD_AF\\tCADD\\tREVEL" > ${sample_id}.clinical_report.tsv
    
    bcftools query \\
        -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%SYMBOL\\t%Consequence\\t%CLNSIG\\t%gnomAD_AF\\t%CADD_PHRED\\t%REVEL_score\\n' \\
        -i 'INFO/CLNSIG~"Pathogenic" || INFO/CLNSIG~"Likely_pathogenic" || (INFO/gnomAD_AF < ${af_threshold} && INFO/CADD_PHRED >= ${cadd_threshold})' \\
        ${vcf} >> ${sample_id}.clinical_report.tsv 2>/dev/null || true
    
    echo "[INFO] Clinical filtering completed"
    echo "[INFO] Pathogenic variants: \$(zcat ${sample_id}.pathogenic.vcf.gz 2>/dev/null | grep -v '^#' | wc -l || echo 0)"
    echo "[INFO] Rare damaging variants: \$(zcat ${sample_id}.rare_damaging.vcf.gz 2>/dev/null | grep -v '^#' | wc -l || echo 0)"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.pathogenic.vcf.gz
    touch ${sample_id}.pathogenic.vcf.gz.tbi
    touch ${sample_id}.rare_damaging.vcf.gz
    touch ${sample_id}.rare_damaging.vcf.gz.tbi
    touch ${sample_id}.clinical_report.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.18
    END_VERSIONS
    """
}
