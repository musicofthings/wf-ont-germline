/*
========================================================================================
    ALIGNMENT MODULES
========================================================================================
*/

process MINIMAP2_ALIGN {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: "*.bam*"
    
    container 'quay.io/biocontainers/minimap2:2.26--he4a0461_2'
    
    input:
    tuple val(sample_id), path(fastq)
    path reference
    path reference_fai
    
    output:
    tuple val(sample_id), path("${sample_id}.unsorted.bam"), emit: bam
    path "versions.yml", emit: versions
    
    script:
    def input_files = fastq instanceof List ? fastq.join(' ') : fastq
    def rg_tag = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ONT\\tPU:${sample_id}"
    def preset = params.minimap2_preset ?: 'map-ont'
    def extra_args = params.minimap2_args ?: ''
    """
    echo "[INFO] =============================================="
    echo "[INFO] Starting minimap2 alignment for sample: ${sample_id}"
    echo "[INFO] Reference: ${reference}"
    echo "[INFO] Preset: ${preset}"
    echo "[INFO] Threads: ${task.cpus}"
    echo "[INFO] =============================================="
    
    minimap2 \\
        -ax ${preset} \\
        -t ${task.cpus} \\
        -R '${rg_tag}' \\
        --MD \\
        --eqx \\
        -Y \\
        ${extra_args} \\
        ${reference} \\
        ${input_files} \\
        2> minimap2.log \\
        | samtools view -@ ${task.cpus} -bh -o ${sample_id}.unsorted.bam -
    
    echo "[INFO] minimap2 alignment completed successfully"
    echo "[INFO] Alignment statistics:"
    samtools flagstat ${sample_id}.unsorted.bam | head -5
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.unsorted.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: 2.26
        samtools: 1.18
    END_VERSIONS
    """
}

process SAMTOOLS_SORT_INDEX {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/alignment", mode: 'copy'
    
    container 'quay.io/biocontainers/samtools:1.18--hd87286a_0'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam
    path "${sample_id}.flagstat.txt", emit: flagstat
    path "${sample_id}.idxstats.txt", emit: idxstats
    path "versions.yml", emit: versions
    
    script:
    def mem_per_thread = Math.floor(task.memory.toGiga() / task.cpus).intValue()
    """
    echo "[INFO] =============================================="
    echo "[INFO] Sorting and indexing BAM for sample: ${sample_id}"
    echo "[INFO] Memory per thread: ${mem_per_thread}G"
    echo "[INFO] =============================================="
    
    # Sort BAM
    samtools sort \\
        -@ ${task.cpus} \\
        -m ${mem_per_thread}G \\
        -o ${sample_id}.sorted.bam \\
        ${bam}
    
    echo "[INFO] BAM sorting completed"
    
    # Index BAM
    samtools index \\
        -@ ${task.cpus} \\
        ${sample_id}.sorted.bam
    
    echo "[INFO] BAM indexing completed"
    
    # Generate statistics
    samtools flagstat \\
        -@ ${task.cpus} \\
        ${sample_id}.sorted.bam > ${sample_id}.flagstat.txt
    
    samtools idxstats \\
        ${sample_id}.sorted.bam > ${sample_id}.idxstats.txt
    
    echo "[INFO] Alignment statistics:"
    cat ${sample_id}.flagstat.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.sorted.bam
    touch ${sample_id}.sorted.bam.bai
    touch ${sample_id}.flagstat.txt
    touch ${sample_id}.idxstats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.18
    END_VERSIONS
    """
}

process SAMTOOLS_MERGE {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/alignment", mode: 'copy'
    
    container 'quay.io/biocontainers/samtools:1.18--hd87286a_0'
    
    input:
    tuple val(sample_id), path(bams)
    
    output:
    tuple val(sample_id), path("${sample_id}.merged.bam"), path("${sample_id}.merged.bam.bai"), emit: bam
    path "versions.yml", emit: versions
    
    when:
    bams.size() > 1
    
    script:
    """
    echo "[INFO] Merging ${bams.size()} BAM files for sample: ${sample_id}"
    
    samtools merge \\
        -@ ${task.cpus} \\
        ${sample_id}.merged.bam \\
        ${bams}
    
    samtools index \\
        -@ ${task.cpus} \\
        ${sample_id}.merged.bam
    
    echo "[INFO] BAM merge completed successfully"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.merged.bam
    touch ${sample_id}.merged.bam.bai
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.18
    END_VERSIONS
    """
}

process SAMTOOLS_STATS {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/qc/samtools_stats", mode: 'copy'
    
    container 'quay.io/biocontainers/samtools:1.18--hd87286a_0'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    
    output:
    tuple val(sample_id), path("${sample_id}.stats.txt"), emit: stats
    path "versions.yml", emit: versions
    
    script:
    """
    echo "[INFO] Generating detailed BAM statistics for sample: ${sample_id}"
    
    samtools stats \\
        -@ ${task.cpus} \\
        --reference ${reference} \\
        ${bam} > ${sample_id}.stats.txt
    
    echo "[INFO] Statistics generation completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.stats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.18
    END_VERSIONS
    """
}
