/*
========================================================================================
    HAPLOTYPE PHASING MODULES
========================================================================================
*/

process WHATSHAP_PHASE {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/phasing", mode: 'copy'
    
    container 'quay.io/biocontainers/whatshap:2.2--py39hf95cd2a_0'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi), path(bam), path(bai)
    path reference
    path reference_fai
    
    output:
    tuple val(sample_id), path("${sample_id}.phased.vcf.gz"), path("${sample_id}.phased.vcf.gz.tbi"), emit: vcf
    path "${sample_id}.phasing.log", emit: log
    path "versions.yml", emit: versions
    
    script:
    def indels_arg = params.phase_indels ? "--indels" : ""
    def ignore_read_groups = params.ignore_read_groups ? "--ignore-read-groups" : ""
    def extra_args = params.whatshap_phase_args ?: ""
    """
    echo "[INFO] =============================================="
    echo "[INFO] Starting WhatsHap phasing"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] =============================================="
    
    whatshap phase \\
        --output ${sample_id}.phased.vcf.gz \\
        --reference ${reference} \\
        ${indels_arg} \\
        ${ignore_read_groups} \\
        ${extra_args} \\
        ${vcf} \\
        ${bam} \\
        2>&1 | tee ${sample_id}.phasing.log
    
    # Index phased VCF
    tabix -p vcf ${sample_id}.phased.vcf.gz
    
    echo "[INFO] WhatsHap phasing completed"
    echo "[INFO] Phasing summary:"
    grep -E "^(Found|Phased|Unphased)" ${sample_id}.phasing.log || true
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version 2>&1)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.phased.vcf.gz
    touch ${sample_id}.phased.vcf.gz.tbi
    touch ${sample_id}.phasing.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: 2.2
    END_VERSIONS
    """
}

process WHATSHAP_HAPLOTAG {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/phasing", mode: 'copy'
    
    container 'quay.io/biocontainers/whatshap:2.2--py39hf95cd2a_0'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi), path(bam), path(bai)
    path reference
    path reference_fai
    
    output:
    tuple val(sample_id), path("${sample_id}.haplotagged.bam"), path("${sample_id}.haplotagged.bam.bai"), emit: bam
    path "${sample_id}.haplotag.log", emit: log
    tuple val(sample_id), path("${sample_id}.haplotag_list.tsv.gz"), emit: haplotag_list, optional: true
    path "versions.yml", emit: versions
    
    script:
    def ignore_read_groups = params.ignore_read_groups ? "--ignore-read-groups" : ""
    def output_haplotag_list = params.output_haplotag_list ? "--output-haplotag-list ${sample_id}.haplotag_list.tsv.gz" : ""
    def extra_args = params.whatshap_haplotag_args ?: ""
    """
    echo "[INFO] =============================================="
    echo "[INFO] Starting WhatsHap haplotagging"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] =============================================="
    
    whatshap haplotag \\
        --output ${sample_id}.haplotagged.bam \\
        --reference ${reference} \\
        ${ignore_read_groups} \\
        ${output_haplotag_list} \\
        ${extra_args} \\
        ${vcf} \\
        ${bam} \\
        2>&1 | tee ${sample_id}.haplotag.log
    
    # Index haplotagged BAM
    samtools index -@ ${task.cpus} ${sample_id}.haplotagged.bam
    
    echo "[INFO] WhatsHap haplotagging completed"
    echo "[INFO] Haplotagging summary:"
    grep -E "^(Total|Tagged|Untagged)" ${sample_id}.haplotag.log || true
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version 2>&1)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.haplotagged.bam
    touch ${sample_id}.haplotagged.bam.bai
    touch ${sample_id}.haplotag.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: 2.2
        samtools: 1.18
    END_VERSIONS
    """
}

process WHATSHAP_STATS {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/phasing", mode: 'copy'
    
    container 'quay.io/biocontainers/whatshap:2.2--py39hf95cd2a_0'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi)
    
    output:
    tuple val(sample_id), path("${sample_id}.phasing_stats.tsv"), emit: stats
    tuple val(sample_id), path("${sample_id}.phasing_stats.gtf"), emit: gtf, optional: true
    tuple val(sample_id), path("${sample_id}.block_list.tsv"), emit: block_list, optional: true
    path "versions.yml", emit: versions
    
    script:
    def gtf_output = params.output_phasing_gtf ? "--gtf ${sample_id}.phasing_stats.gtf" : ""
    def block_list = params.output_block_list ? "--block-list ${sample_id}.block_list.tsv" : ""
    """
    echo "[INFO] Generating phasing statistics for sample: ${sample_id}"
    
    whatshap stats \\
        --tsv ${sample_id}.phasing_stats.tsv \\
        ${gtf_output} \\
        ${block_list} \\
        ${vcf}
    
    echo "[INFO] Phasing statistics:"
    cat ${sample_id}.phasing_stats.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version 2>&1)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.phasing_stats.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: 2.2
    END_VERSIONS
    """
}

process WHATSHAP_COMPARE {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/phasing/comparison", mode: 'copy'
    
    container 'quay.io/biocontainers/whatshap:2.2--py39hf95cd2a_0'
    
    input:
    tuple val(sample_id), path(vcf1), path(vcf1_tbi), path(vcf2), path(vcf2_tbi)
    
    output:
    tuple val(sample_id), path("${sample_id}.phasing_comparison.tsv"), emit: comparison
    path "versions.yml", emit: versions
    
    when:
    params.compare_phasing
    
    script:
    """
    echo "[INFO] Comparing phasing between two VCFs for sample: ${sample_id}"
    
    whatshap compare \\
        --tsv-pairwise ${sample_id}.phasing_comparison.tsv \\
        ${vcf1} \\
        ${vcf2}
    
    echo "[INFO] Phasing comparison completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version 2>&1)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.phasing_comparison.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: 2.2
    END_VERSIONS
    """
}

process LONGPHASE {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/phasing", mode: 'copy'
    
    container 'quay.io/biocontainers/longphase:1.5--hf5e1c6e_0'
    
    input:
    tuple val(sample_id), path(snv_vcf), path(snv_vcf_tbi), path(sv_vcf), path(sv_vcf_tbi), path(bam), path(bai)
    path reference
    path reference_fai
    
    output:
    tuple val(sample_id), path("${sample_id}.longphase.snv.vcf.gz"), path("${sample_id}.longphase.snv.vcf.gz.tbi"), emit: snv_vcf
    tuple val(sample_id), path("${sample_id}.longphase.sv.vcf.gz"), path("${sample_id}.longphase.sv.vcf.gz.tbi"), emit: sv_vcf, optional: true
    path "versions.yml", emit: versions
    
    when:
    params.use_longphase
    
    script:
    def sv_arg = sv_vcf ? "-s ${sv_vcf}" : ""
    """
    echo "[INFO] =============================================="
    echo "[INFO] Starting LongPhase phasing"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] =============================================="
    
    longphase phase \\
        -o ${sample_id}.longphase \\
        -r ${reference} \\
        -b ${bam} \\
        -v ${snv_vcf} \\
        ${sv_arg} \\
        -t ${task.cpus} \\
        --ont \\
        2>&1 | tee longphase.log
    
    # Compress and index outputs
    bgzip -c ${sample_id}.longphase.snv.vcf > ${sample_id}.longphase.snv.vcf.gz
    tabix -p vcf ${sample_id}.longphase.snv.vcf.gz
    
    if [ -f "${sample_id}.longphase.sv.vcf" ]; then
        bgzip -c ${sample_id}.longphase.sv.vcf > ${sample_id}.longphase.sv.vcf.gz
        tabix -p vcf ${sample_id}.longphase.sv.vcf.gz
    fi
    
    echo "[INFO] LongPhase completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: \$(longphase --version 2>&1 | head -1)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.longphase.snv.vcf.gz
    touch ${sample_id}.longphase.snv.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: 1.5
    END_VERSIONS
    """
}
