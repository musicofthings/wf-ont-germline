/*
========================================================================================
    QUALITY CONTROL MODULES
========================================================================================
*/

process NANOPLOT {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/qc/nanoplot", mode: 'copy', pattern: "*"
    
    container 'quay.io/biocontainers/nanoplot:1.42.0--pyhdfd78af_0'
    
    input:
    tuple val(sample_id), path(fastq)
    
    output:
    tuple val(sample_id), path("${sample_id}_nanoplot"), emit: report
    path "${sample_id}_nanoplot/NanoStats.txt", emit: stats
    path "versions.yml", emit: versions
    
    script:
    def input_files = fastq instanceof List ? fastq.join(' ') : fastq
    """
    echo "[INFO] Running NanoPlot QC for sample: ${sample_id}"
    echo "[INFO] Input files: ${input_files}"
    
    NanoPlot \\
        --fastq ${input_files} \\
        --outdir ${sample_id}_nanoplot \\
        --threads ${task.cpus} \\
        --plots hex dot \\
        --N50 \\
        --title "${sample_id} - ONT Read QC" \\
        --info_in_report \\
        --store \\
        --raw \\
        2>&1 | tee nanoplot.log
    
    echo "[INFO] NanoPlot completed successfully"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(NanoPlot --version 2>&1 | sed 's/NanoPlot //')
    END_VERSIONS
    """
    
    stub:
    """
    mkdir -p ${sample_id}_nanoplot
    touch ${sample_id}_nanoplot/NanoStats.txt
    touch ${sample_id}_nanoplot/NanoPlot-report.html
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: 1.42.0
    END_VERSIONS
    """
}

process MOSDEPTH {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/qc/mosdepth", mode: 'copy'
    
    container 'quay.io/biocontainers/mosdepth:0.3.6--hd299f12_0'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    
    output:
    tuple val(sample_id), path("${sample_id}.mosdepth.global.dist.txt"), emit: report
    tuple val(sample_id), path("${sample_id}.mosdepth.summary.txt"), emit: summary
    tuple val(sample_id), path("${sample_id}.per-base.bed.gz"), emit: per_base, optional: true
    tuple val(sample_id), path("${sample_id}.regions.bed.gz"), emit: regions, optional: true
    path "versions.yml", emit: versions
    
    script:
    def by_threshold = params.mosdepth_threshold ? "--thresholds ${params.mosdepth_threshold}" : ""
    def window_size = params.mosdepth_window ?: 500
    """
    echo "[INFO] Running mosdepth coverage analysis for sample: ${sample_id}"
    
    mosdepth \\
        --threads ${task.cpus} \\
        --by ${window_size} \\
        --fast-mode \\
        ${by_threshold} \\
        --fasta ${reference} \\
        ${sample_id} \\
        ${bam} \\
        2>&1 | tee mosdepth.log
    
    echo "[INFO] mosdepth completed successfully"
    echo "[INFO] Mean coverage: \$(grep 'total' ${sample_id}.mosdepth.summary.txt | cut -f4)"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.mosdepth.global.dist.txt
    touch ${sample_id}.mosdepth.summary.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: 0.3.6
    END_VERSIONS
    """
}

process FASTQC_LONG {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'
    
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    
    input:
    tuple val(sample_id), path(fastq)
    
    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"), emit: zip
    path "versions.yml", emit: versions
    
    script:
    """
    echo "[INFO] Running FastQC for sample: ${sample_id}"
    
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        --memory ${task.memory.toMega()}M \\
        ${fastq}
    
    echo "[INFO] FastQC completed successfully"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}_fastqc.html
    touch ${sample_id}_fastqc.zip
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: 0.12.1
    END_VERSIONS
    """
}

process PYCOQC {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/qc/pycoqc", mode: 'copy'
    
    container 'quay.io/biocontainers/pycoqc:2.5.2--py_0'
    
    input:
    tuple val(sample_id), path(sequencing_summary)
    
    output:
    tuple val(sample_id), path("${sample_id}_pycoqc.html"), emit: html
    tuple val(sample_id), path("${sample_id}_pycoqc.json"), emit: json
    path "versions.yml", emit: versions
    
    when:
    params.sequencing_summary
    
    script:
    """
    echo "[INFO] Running PycoQC for sample: ${sample_id}"
    
    pycoQC \\
        --summary_file ${sequencing_summary} \\
        --html_outfile ${sample_id}_pycoqc.html \\
        --json_outfile ${sample_id}_pycoqc.json \\
        --min_pass_qual ${params.min_pass_qual ?: 7} \\
        --verbose
    
    echo "[INFO] PycoQC completed successfully"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pycoqc: \$(pycoQC --version 2>&1 | sed 's/pycoQC v//')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}_pycoqc.html
    touch ${sample_id}_pycoqc.json
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pycoqc: 2.5.2
    END_VERSIONS
    """
}

process CRAMINO {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/qc/cramino", mode: 'copy'
    
    container 'quay.io/biocontainers/cramino:0.14.5--h4ac6f70_0'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("${sample_id}_cramino.txt"), emit: report
    path "versions.yml", emit: versions
    
    script:
    """
    echo "[INFO] Running cramino BAM QC for sample: ${sample_id}"
    
    cramino \\
        --threads ${task.cpus} \\
        ${bam} > ${sample_id}_cramino.txt
    
    echo "[INFO] cramino completed successfully"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cramino: \$(cramino --version 2>&1 | sed 's/cramino //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}_cramino.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cramino: 0.14.5
    END_VERSIONS
    """
}
