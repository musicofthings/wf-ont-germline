/*
========================================================================================
    VARIANT CALLING MODULES
========================================================================================
*/

process CLAIR3 {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/variants/small_variants", mode: 'copy'
    
    container 'hkubal/clair3:v1.0.5'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_fai
    val model_path
    
    output:
    tuple val(sample_id), path("${sample_id}.clair3.vcf.gz"), path("${sample_id}.clair3.vcf.gz.tbi"), emit: vcf
    tuple val(sample_id), path("${sample_id}.clair3.g.vcf.gz"), path("${sample_id}.clair3.g.vcf.gz.tbi"), emit: gvcf, optional: true
    path "${sample_id}_clair3/", emit: output_dir
    path "versions.yml", emit: versions
    
    script:
    def gvcf_arg = params.output_gvcf ? "--gvcf" : ""
    def include_all_ctgs = params.include_all_ctgs ? "--include_all_ctgs" : ""
    def haploid_precise = params.haploid_precise ? "--haploid_precise" : ""
    def haploid_sensitive = params.haploid_sensitive ? "--haploid_sensitive" : ""
    def enable_phasing = params.clair3_enable_phasing ? "--enable_phasing" : ""
    def extra_args = params.clair3_args ?: ""
    """
    echo "[INFO] =============================================="
    echo "[INFO] Starting Clair3 variant calling"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] Model: ${model_path}"
    echo "[INFO] Threads: ${task.cpus}"
    echo "[INFO] =============================================="
    
    # Create output directory
    mkdir -p ${sample_id}_clair3
    
    # Run Clair3
    run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${reference} \\
        --threads=${task.cpus} \\
        --platform="ont" \\
        --model_path="${model_path}" \\
        --output=${sample_id}_clair3 \\
        --sample_name=${sample_id} \\
        ${gvcf_arg} \\
        ${include_all_ctgs} \\
        ${haploid_precise} \\
        ${haploid_sensitive} \\
        ${enable_phasing} \\
        ${extra_args} \\
        2>&1 | tee clair3.log
    
    # Copy and rename output files
    cp ${sample_id}_clair3/merge_output.vcf.gz ${sample_id}.clair3.vcf.gz
    cp ${sample_id}_clair3/merge_output.vcf.gz.tbi ${sample_id}.clair3.vcf.gz.tbi
    
    if [ -f "${sample_id}_clair3/merge_output.gvcf.gz" ]; then
        cp ${sample_id}_clair3/merge_output.gvcf.gz ${sample_id}.clair3.g.vcf.gz
        cp ${sample_id}_clair3/merge_output.gvcf.gz.tbi ${sample_id}.clair3.g.vcf.gz.tbi
    fi
    
    # Report variant counts
    echo "[INFO] Clair3 variant calling completed"
    echo "[INFO] Variant counts:"
    bcftools stats ${sample_id}.clair3.vcf.gz | grep "^SN" | head -10
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh --version 2>&1 | head -1 || echo "1.0.5")
    END_VERSIONS
    """
    
    stub:
    """
    mkdir -p ${sample_id}_clair3
    touch ${sample_id}.clair3.vcf.gz
    touch ${sample_id}.clair3.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: 1.0.5
    END_VERSIONS
    """
}

process PEPPER_DEEPVARIANT {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/variants/small_variants", mode: 'copy'
    
    container 'kishwars/pepper_deepvariant:r0.8'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_fai
    
    output:
    tuple val(sample_id), path("${sample_id}.pepper_dv.vcf.gz"), path("${sample_id}.pepper_dv.vcf.gz.tbi"), emit: vcf
    tuple val(sample_id), path("${sample_id}.pepper_dv.g.vcf.gz"), path("${sample_id}.pepper_dv.g.vcf.gz.tbi"), emit: gvcf, optional: true
    tuple val(sample_id), path("${sample_id}.pepper_dv.phased.vcf.gz"), emit: phased_vcf, optional: true
    path "versions.yml", emit: versions
    
    script:
    def gvcf_arg = params.output_gvcf ? "--gvcf" : ""
    def phased_arg = params.pepper_phased_output ? "--phased_output" : ""
    def model_preset = params.pepper_model ?: "ont_r10_q20"
    def extra_args = params.pepper_args ?: ""
    """
    echo "[INFO] =============================================="
    echo "[INFO] Starting PEPPER-DeepVariant variant calling"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] Model preset: ${model_preset}"
    echo "[INFO] Threads: ${task.cpus}"
    echo "[INFO] =============================================="
    
    run_pepper_margin_deepvariant call_variant \\
        -b ${bam} \\
        -f ${reference} \\
        -o . \\
        -p ${sample_id}.pepper_dv \\
        -t ${task.cpus} \\
        --${model_preset} \\
        ${gvcf_arg} \\
        ${phased_arg} \\
        ${extra_args} \\
        2>&1 | tee pepper_dv.log
    
    # Index output VCF
    tabix -p vcf ${sample_id}.pepper_dv.vcf.gz
    
    if [ -f "${sample_id}.pepper_dv.g.vcf.gz" ]; then
        tabix -p vcf ${sample_id}.pepper_dv.g.vcf.gz
    fi
    
    echo "[INFO] PEPPER-DeepVariant completed"
    echo "[INFO] Variant counts:"
    bcftools stats ${sample_id}.pepper_dv.vcf.gz | grep "^SN" | head -10
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pepper_deepvariant: r0.8
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.pepper_dv.vcf.gz
    touch ${sample_id}.pepper_dv.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pepper_deepvariant: r0.8
    END_VERSIONS
    """
}

process SNIFFLES2 {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/variants/structural", mode: 'copy'
    
    container 'quay.io/biocontainers/sniffles:2.2--pyhdfd78af_0'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_fai
    
    output:
    tuple val(sample_id), path("${sample_id}.sniffles.vcf.gz"), path("${sample_id}.sniffles.vcf.gz.tbi"), emit: vcf
    tuple val(sample_id), path("${sample_id}.sniffles.snf"), emit: snf
    path "versions.yml", emit: versions
    
    script:
    def min_sv_len = params.min_sv_length ?: 50
    def min_support = params.min_sv_support ?: 'auto'
    def extra_args = params.sniffles_args ?: ""
    """
    echo "[INFO] =============================================="
    echo "[INFO] Starting Sniffles2 structural variant calling"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] Min SV length: ${min_sv_len}"
    echo "[INFO] Threads: ${task.cpus}"
    echo "[INFO] =============================================="
    
    sniffles \\
        --input ${bam} \\
        --vcf ${sample_id}.sniffles.vcf.gz \\
        --snf ${sample_id}.sniffles.snf \\
        --reference ${reference} \\
        --threads ${task.cpus} \\
        --sample-id ${sample_id} \\
        --minsvlen ${min_sv_len} \\
        --minsupport ${min_support} \\
        --output-rnames \\
        ${extra_args} \\
        2>&1 | tee sniffles.log
    
    # Index VCF
    tabix -p vcf ${sample_id}.sniffles.vcf.gz
    
    echo "[INFO] Sniffles2 completed"
    echo "[INFO] SV counts by type:"
    bcftools query -f '%SVTYPE\\n' ${sample_id}.sniffles.vcf.gz | sort | uniq -c
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles: \$(sniffles --version 2>&1 | sed 's/Sniffles //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.sniffles.vcf.gz
    touch ${sample_id}.sniffles.vcf.gz.tbi
    touch ${sample_id}.sniffles.snf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles: 2.2
    END_VERSIONS
    """
}

process CUTESV {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/variants/structural", mode: 'copy'
    
    container 'quay.io/biocontainers/cutesv:2.1.1--pyhdfd78af_0'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    
    output:
    tuple val(sample_id), path("${sample_id}.cutesv.vcf.gz"), path("${sample_id}.cutesv.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    def min_sv_len = params.min_sv_length ?: 50
    def min_support = params.min_sv_support_cutesv ?: 10
    def max_cluster_bias_INS = params.cutesv_max_cluster_bias_INS ?: 1000
    def diff_ratio_merging_INS = params.cutesv_diff_ratio_merging_INS ?: 0.9
    def max_cluster_bias_DEL = params.cutesv_max_cluster_bias_DEL ?: 1000
    def diff_ratio_merging_DEL = params.cutesv_diff_ratio_merging_DEL ?: 0.5
    """
    echo "[INFO] =============================================="
    echo "[INFO] Starting CuteSV structural variant calling"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] =============================================="
    
    mkdir -p cutesv_workdir
    
    cuteSV \\
        ${bam} \\
        ${reference} \\
        ${sample_id}.cutesv.vcf \\
        cutesv_workdir \\
        --threads ${task.cpus} \\
        --sample ${sample_id} \\
        --min_size ${min_sv_len} \\
        --min_support ${min_support} \\
        --max_cluster_bias_INS ${max_cluster_bias_INS} \\
        --diff_ratio_merging_INS ${diff_ratio_merging_INS} \\
        --max_cluster_bias_DEL ${max_cluster_bias_DEL} \\
        --diff_ratio_merging_DEL ${diff_ratio_merging_DEL} \\
        --genotype \\
        2>&1 | tee cutesv.log
    
    # Compress and index
    bgzip -c ${sample_id}.cutesv.vcf > ${sample_id}.cutesv.vcf.gz
    tabix -p vcf ${sample_id}.cutesv.vcf.gz
    
    echo "[INFO] CuteSV completed"
    echo "[INFO] SV counts by type:"
    bcftools query -f '%SVTYPE\\n' ${sample_id}.cutesv.vcf.gz | sort | uniq -c
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutesv: \$(cuteSV --version 2>&1 | sed 's/cuteSV //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.cutesv.vcf.gz
    touch ${sample_id}.cutesv.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutesv: 2.1.1
    END_VERSIONS
    """
}

process MERGE_SV_CALLS {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/variants/structural", mode: 'copy'
    
    container 'quay.io/biocontainers/survivor:1.0.7--h9f5acd7_2'
    
    input:
    tuple val(sample_id), path(vcfs)
    
    output:
    tuple val(sample_id), path("${sample_id}.sv_merged.vcf.gz"), path("${sample_id}.sv_merged.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    def max_dist = params.sv_merge_max_dist ?: 1000
    def min_callers = params.sv_merge_min_callers ?: 1
    """
    echo "[INFO] Merging SV calls from multiple callers for sample: ${sample_id}"
    
    # Create file list for SURVIVOR
    ls *.vcf.gz > vcf_list.txt
    
    # Decompress VCFs for SURVIVOR
    for vcf in *.vcf.gz; do
        gunzip -c \$vcf > \${vcf%.gz}
    done
    
    ls *.vcf | grep -v merged > vcf_list_unzipped.txt
    
    # Merge with SURVIVOR
    SURVIVOR merge vcf_list_unzipped.txt ${max_dist} ${min_callers} 1 1 0 50 ${sample_id}.sv_merged.vcf
    
    # Compress and index
    bgzip -c ${sample_id}.sv_merged.vcf > ${sample_id}.sv_merged.vcf.gz
    tabix -p vcf ${sample_id}.sv_merged.vcf.gz
    
    echo "[INFO] SV merging completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(SURVIVOR 2>&1 | head -2 | tail -1 | sed 's/Version: //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.sv_merged.vcf.gz
    touch ${sample_id}.sv_merged.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: 1.0.7
    END_VERSIONS
    """
}
