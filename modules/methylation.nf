/*
========================================================================================
    METHYLATION CALLING MODULES
========================================================================================
*/

process MODKIT_PILEUP {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/methylation", mode: 'copy'
    
    container 'ontresearch/modkit:v0.2.5'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_fai
    
    output:
    tuple val(sample_id), path("${sample_id}.bedmethyl.gz"), emit: bedmethyl
    tuple val(sample_id), path("${sample_id}.bedmethyl.gz.tbi"), emit: bedmethyl_tbi
    tuple val(sample_id), path("${sample_id}_hap1.bedmethyl.gz"), emit: hap1_bedmethyl, optional: true
    tuple val(sample_id), path("${sample_id}_hap2.bedmethyl.gz"), emit: hap2_bedmethyl, optional: true
    path "versions.yml", emit: versions
    
    script:
    def partition_tag = params.partition_by_haplotype ? "--partition-tag HP" : ""
    def motif = params.methylation_motif ?: "CG 0"
    def filter_threshold = params.modkit_filter_threshold ?: 0.75
    def combine_mods = params.combine_mods ? "--combine-mods" : ""
    def extra_args = params.modkit_pileup_args ?: ""
    """
    echo "[INFO] =============================================="
    echo "[INFO] Starting modkit methylation pileup"
    echo "[INFO] Sample: ${sample_id}"
    echo "[INFO] Motif: ${motif}"
    echo "[INFO] Filter threshold: ${filter_threshold}"
    echo "[INFO] Partition by haplotype: ${params.partition_by_haplotype ?: false}"
    echo "[INFO] =============================================="
    
    # Run modkit pileup
    modkit pileup \\
        ${bam} \\
        ${sample_id}.bedmethyl \\
        --ref ${reference} \\
        --threads ${task.cpus} \\
        --motif ${motif} \\
        --filter-threshold ${filter_threshold} \\
        ${partition_tag} \\
        ${combine_mods} \\
        ${extra_args} \\
        --log-filepath modkit_pileup.log \\
        2>&1 | tee modkit.log
    
    # Compress and index main output
    bgzip -c ${sample_id}.bedmethyl > ${sample_id}.bedmethyl.gz
    tabix -p bed ${sample_id}.bedmethyl.gz
    
    # Handle haplotype-partitioned outputs if present
    if [ -f "${sample_id}_1.bedmethyl" ]; then
        bgzip -c ${sample_id}_1.bedmethyl > ${sample_id}_hap1.bedmethyl.gz
        tabix -p bed ${sample_id}_hap1.bedmethyl.gz
    fi
    
    if [ -f "${sample_id}_2.bedmethyl" ]; then
        bgzip -c ${sample_id}_2.bedmethyl > ${sample_id}_hap2.bedmethyl.gz
        tabix -p bed ${sample_id}_hap2.bedmethyl.gz
    fi
    
    echo "[INFO] modkit pileup completed"
    echo "[INFO] Methylation sites:"
    zcat ${sample_id}.bedmethyl.gz | wc -l
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version 2>&1 | sed 's/modkit //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.bedmethyl.gz
    touch ${sample_id}.bedmethyl.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: 0.2.5
    END_VERSIONS
    """
}

process MODKIT_SUMMARY {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/methylation", mode: 'copy'
    
    container 'ontresearch/modkit:v0.2.5'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("${sample_id}.modkit_summary.txt"), emit: summary
    tuple val(sample_id), path("${sample_id}.modkit_summary.tsv"), emit: summary_tsv
    path "versions.yml", emit: versions
    
    script:
    """
    echo "[INFO] Generating modkit summary for sample: ${sample_id}"
    
    # Generate summary
    modkit summary \\
        ${bam} \\
        --threads ${task.cpus} \\
        --log-filepath modkit_summary.log \\
        > ${sample_id}.modkit_summary.txt \\
        2>&1
    
    # Also generate TSV format
    modkit summary \\
        ${bam} \\
        --threads ${task.cpus} \\
        --tsv \\
        > ${sample_id}.modkit_summary.tsv \\
        2>&1
    
    echo "[INFO] modkit summary completed"
    cat ${sample_id}.modkit_summary.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version 2>&1 | sed 's/modkit //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.modkit_summary.txt
    touch ${sample_id}.modkit_summary.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: 0.2.5
    END_VERSIONS
    """
}

process MODKIT_EXTRACT {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/methylation/extract", mode: 'copy'
    
    container 'ontresearch/modkit:v0.2.5'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_fai
    path regions_bed
    
    output:
    tuple val(sample_id), path("${sample_id}.extract.tsv.gz"), emit: extract
    path "versions.yml", emit: versions
    
    when:
    params.extract_methylation
    
    script:
    def regions_arg = regions_bed ? "--include-bed ${regions_bed}" : ""
    """
    echo "[INFO] Extracting per-read methylation for sample: ${sample_id}"
    
    modkit extract \\
        ${bam} \\
        ${sample_id}.extract.tsv \\
        --ref ${reference} \\
        --threads ${task.cpus} \\
        ${regions_arg} \\
        --log-filepath modkit_extract.log \\
        2>&1 | tee extract.log
    
    # Compress output
    gzip ${sample_id}.extract.tsv
    
    echo "[INFO] Extraction completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version 2>&1 | sed 's/modkit //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.extract.tsv.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: 0.2.5
    END_VERSIONS
    """
}

process MODKIT_DMR {
    tag "$comparison"
    label 'process_high'
    publishDir "${params.outdir}/methylation/dmr", mode: 'copy'
    
    container 'ontresearch/modkit:v0.2.5'
    
    input:
    tuple val(comparison), path(bedmethyl_a), path(bedmethyl_b)
    path reference
    path reference_fai
    path regions_bed
    
    output:
    tuple val(comparison), path("${comparison}.dmr.bed"), emit: dmr
    tuple val(comparison), path("${comparison}.dmr_segments.bed"), emit: segments
    path "versions.yml", emit: versions
    
    when:
    params.call_dmr
    
    script:
    def regions_arg = regions_bed ? "--regions ${regions_bed}" : ""
    """
    echo "[INFO] Calling differential methylation regions: ${comparison}"
    
    modkit dmr pair \\
        -a ${bedmethyl_a} \\
        -b ${bedmethyl_b} \\
        -o ${comparison}.dmr.bed \\
        --ref ${reference} \\
        ${regions_arg} \\
        --threads ${task.cpus} \\
        --log-filepath modkit_dmr.log \\
        2>&1 | tee dmr.log
    
    # Segment DMRs
    modkit dmr segment \\
        ${comparison}.dmr.bed \\
        -o ${comparison}.dmr_segments.bed \\
        2>&1
    
    echo "[INFO] DMR calling completed"
    echo "[INFO] DMR count:"
    wc -l ${comparison}.dmr_segments.bed
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version 2>&1 | sed 's/modkit //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${comparison}.dmr.bed
    touch ${comparison}.dmr_segments.bed
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: 0.2.5
    END_VERSIONS
    """
}

process METHYLATION_QC {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/methylation/qc", mode: 'copy'
    
    container 'quay.io/biocontainers/python:3.11'
    
    input:
    tuple val(sample_id), path(bedmethyl)
    
    output:
    tuple val(sample_id), path("${sample_id}_methylation_stats.tsv"), emit: stats
    tuple val(sample_id), path("${sample_id}_methylation_distribution.png"), emit: plot
    path "versions.yml", emit: versions
    
    script:
    """
    echo "[INFO] Generating methylation QC for sample: ${sample_id}"
    
    python3 << 'PYTHON'
import gzip
import pandas as pd
import matplotlib.pyplot as plt

# Read bedmethyl file
data = []
with gzip.open("${bedmethyl}", 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\\t')
        if len(fields) >= 11:
            data.append({
                'chrom': fields[0],
                'start': int(fields[1]),
                'end': int(fields[2]),
                'mod_code': fields[3],
                'score': int(fields[4]),
                'strand': fields[5],
                'coverage': int(fields[9]),
                'percent_modified': float(fields[10])
            })

df = pd.DataFrame(data)

# Calculate statistics
stats = {
    'total_sites': len(df),
    'mean_coverage': df['coverage'].mean(),
    'median_coverage': df['coverage'].median(),
    'mean_methylation': df['percent_modified'].mean(),
    'median_methylation': df['percent_modified'].median(),
    'sites_above_10x': (df['coverage'] >= 10).sum(),
    'highly_methylated': (df['percent_modified'] >= 80).sum(),
    'lowly_methylated': (df['percent_modified'] <= 20).sum()
}

# Save statistics
stats_df = pd.DataFrame([stats])
stats_df.to_csv("${sample_id}_methylation_stats.tsv", sep='\\t', index=False)

# Create distribution plot
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Methylation distribution
axes[0].hist(df['percent_modified'], bins=50, edgecolor='black', alpha=0.7)
axes[0].set_xlabel('Percent Methylation')
axes[0].set_ylabel('Count')
axes[0].set_title('Methylation Distribution')

# Coverage distribution
axes[1].hist(df['coverage'], bins=50, edgecolor='black', alpha=0.7)
axes[1].set_xlabel('Coverage')
axes[1].set_ylabel('Count')
axes[1].set_title('Coverage Distribution')
axes[1].set_xlim(0, df['coverage'].quantile(0.99))

plt.tight_layout()
plt.savefig("${sample_id}_methylation_distribution.png", dpi=150)

print(f"Total sites: {stats['total_sites']}")
print(f"Mean methylation: {stats['mean_methylation']:.2f}%")
print(f"Mean coverage: {stats['mean_coverage']:.1f}x")
PYTHON
    
    echo "[INFO] Methylation QC completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}_methylation_stats.tsv
    touch ${sample_id}_methylation_distribution.png
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}

process BEDMETHYL_TO_BIGWIG {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/methylation/bigwig", mode: 'copy'
    
    container 'quay.io/biocontainers/ucsc-bedgraphtobigwig:447--h954228d_0'
    
    input:
    tuple val(sample_id), path(bedmethyl)
    path chrom_sizes
    
    output:
    tuple val(sample_id), path("${sample_id}.methylation.bw"), emit: bigwig
    path "versions.yml", emit: versions
    
    script:
    """
    echo "[INFO] Converting bedmethyl to BigWig for sample: ${sample_id}"
    
    # Extract methylation percentage as bedGraph
    zcat ${bedmethyl} | \\
        awk 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, \$11}' | \\
        sort -k1,1 -k2,2n > ${sample_id}.bedgraph
    
    # Convert to BigWig
    bedGraphToBigWig ${sample_id}.bedgraph ${chrom_sizes} ${sample_id}.methylation.bw
    
    echo "[INFO] BigWig conversion completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedGraphToBigWig: \$(bedGraphToBigWig 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_id}.methylation.bw
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedGraphToBigWig: 447
    END_VERSIONS
    """
}
