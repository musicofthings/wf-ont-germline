#!/usr/bin/env nextflow

/*
========================================================================================
    wf-ont-germline: ONT Long-Read Germline Variant Calling Pipeline
========================================================================================
    A comprehensive Nextflow pipeline for Oxford Nanopore whole genome sequencing
    germline variant analysis. Compatible with EPI2ME Desktop.
    
    Features:
    - Quality control with NanoPlot/PycoQC
    - Alignment with minimap2
    - SNV/Indel calling with Clair3 or PEPPER-DeepVariant
    - Structural variant calling with Sniffles2
    - Copy number variant calling with Spectre
    - Haplotype phasing with WhatsHap
    - Methylation calling with modkit
    - Comprehensive variant annotation (ClinVar, gnomAD, CADD, etc.)
    
    Author: Biomni Pipeline Generator
    Version: 1.0.0
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// Pipeline version
def pipeline_version = "1.0.0"

/*
========================================================================================
    PRINT PIPELINE INFO
========================================================================================
*/

log.info """
╔═══════════════════════════════════════════════════════════════════════════════════╗
║                                                                                   ║
║   ██╗    ██╗███████╗       ██████╗ ███╗   ██╗████████╗                           ║
║   ██║    ██║██╔════╝      ██╔═══██╗████╗  ██║╚══██╔══╝                           ║
║   ██║ █╗ ██║█████╗  █████╗██║   ██║██╔██╗ ██║   ██║                              ║
║   ██║███╗██║██╔══╝  ╚════╝██║   ██║██║╚██╗██║   ██║                              ║
║   ╚███╔███╔╝██║           ╚██████╔╝██║ ╚████║   ██║                              ║
║    ╚══╝╚══╝ ╚═╝            ╚═════╝ ╚═╝  ╚═══╝   ╚═╝                              ║
║                                                                                   ║
║   GERMLINE VARIANT CALLING PIPELINE v${pipeline_version}                                     ║
║   Oxford Nanopore Long-Read Whole Genome Sequencing Analysis                      ║
║                                                                                   ║
╚═══════════════════════════════════════════════════════════════════════════════════╝

    Pipeline Parameters:
    =====================
    Input/Output:
      --input              : ${params.input}
      --reference          : ${params.reference}
      --outdir             : ${params.outdir}
      --sample_name        : ${params.sample_name}
    
    Analysis Options:
      --variant_caller     : ${params.variant_caller}
      --call_svs           : ${params.call_svs}
      --call_cnvs          : ${params.call_cnvs}
      --phase_variants     : ${params.phase_variants}
      --call_methylation   : ${params.call_methylation}
      --annotate_variants  : ${params.annotate_variants}
    
    Model/Chemistry:
      --basecall_model     : ${params.basecall_model}
      --clair3_model       : ${params.clair3_model}
    
    Resources:
      --max_cpus           : ${params.max_cpus}
      --max_memory         : ${params.max_memory}
    
    =====================
    Starting pipeline...
    =====================
"""

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Check mandatory parameters
def checkPathParamList = [ params.input, params.reference ]
for (param in checkPathParamList) {
    if (param) { file(param, checkIfExists: true) }
}

if (!params.input) {
    exit 1, "ERROR: Input not specified! Please provide --input parameter."
}

if (!params.reference) {
    exit 1, "ERROR: Reference genome not specified! Please provide --reference parameter."
}

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

// Quality Control
include { NANOPLOT } from './modules/qc'
include { MOSDEPTH } from './modules/qc'
include { FASTQC_LONG } from './modules/qc'

// Alignment
include { MINIMAP2_ALIGN } from './modules/alignment'
include { SAMTOOLS_SORT_INDEX } from './modules/alignment'

// Variant Calling
include { CLAIR3 } from './modules/variant_calling'
include { PEPPER_DEEPVARIANT } from './modules/variant_calling'
include { SNIFFLES2 } from './modules/variant_calling'
include { CUTESV } from './modules/variant_calling'

// CNV Calling
include { SPECTRE } from './modules/cnv'
include { QDNASEQ } from './modules/cnv'

// Phasing
include { WHATSHAP_PHASE } from './modules/phasing'
include { WHATSHAP_HAPLOTAG } from './modules/phasing'
include { WHATSHAP_STATS } from './modules/phasing'

// Methylation
include { MODKIT_PILEUP } from './modules/methylation'
include { MODKIT_SUMMARY } from './modules/methylation'

// Annotation
include { VEP_ANNOTATE } from './modules/annotation'
include { CLINVAR_ANNOTATE } from './modules/annotation'
include { GNOMAD_ANNOTATE } from './modules/annotation'
include { CADD_ANNOTATE } from './modules/annotation'
include { SPLICEAI_ANNOTATE } from './modules/annotation'
include { PHARMGKB_ANNOTATE } from './modules/annotation'
include { DBSNP_ANNOTATE } from './modules/annotation'
include { COSMIC_ANNOTATE } from './modules/annotation'
include { MERGE_ANNOTATIONS } from './modules/annotation'

// Reporting
include { MULTIQC } from './modules/reporting'
include { GENERATE_REPORT } from './modules/reporting'

/*
========================================================================================
    HELPER FUNCTIONS
========================================================================================
*/

def logStep(step_name, status) {
    def timestamp = new Date().format("yyyy-MM-dd HH:mm:ss")
    def status_icon = status == "START" ? "🚀" : (status == "DONE" ? "✅" : "⏳")
    log.info "[${timestamp}] ${status_icon} ${step_name} - ${status}"
}

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    
    // Create channels
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    /*
    ================================================================================
        STEP 1: INPUT PREPARATION
    ================================================================================
    */
    logStep("INPUT PREPARATION", "START")
    
    // Handle input - can be FASTQ, BAM, or directory
    if (params.input.endsWith('.bam')) {
        // Input is already aligned BAM
        Channel
            .fromPath(params.input, checkIfExists: true)
            .map { bam -> 
                def sample = params.sample_name ?: bam.simpleName
                tuple(sample, bam, file("${bam}.bai"))
            }
            .set { ch_aligned_bam }
        ch_fastq = Channel.empty()
        skip_alignment = true
    } else if (params.input.endsWith('.fastq') || params.input.endsWith('.fastq.gz') || params.input.endsWith('.fq') || params.input.endsWith('.fq.gz')) {
        // Single FASTQ file
        Channel
            .fromPath(params.input, checkIfExists: true)
            .map { fastq -> 
                def sample = params.sample_name ?: fastq.simpleName.replaceAll(/\.fastq|\.fq|\.gz/, '')
                tuple(sample, fastq)
            }
            .set { ch_fastq }
        ch_aligned_bam = Channel.empty()
        skip_alignment = false
    } else {
        // Directory of FASTQ files
        Channel
            .fromPath("${params.input}/*.{fastq,fastq.gz,fq,fq.gz}", checkIfExists: true)
            .map { fastq -> 
                def sample = params.sample_name ?: fastq.simpleName.replaceAll(/\.fastq|\.fq|\.gz/, '')
                tuple(sample, fastq)
            }
            .groupTuple()
            .map { sample, fastqs -> tuple(sample, fastqs.flatten()) }
            .set { ch_fastq }
        ch_aligned_bam = Channel.empty()
        skip_alignment = false
    }
    
    // Reference genome channel
    ch_reference = Channel.fromPath(params.reference, checkIfExists: true)
    ch_reference_fai = Channel.fromPath("${params.reference}.fai", checkIfExists: true)
    
    logStep("INPUT PREPARATION", "DONE")
    
    /*
    ================================================================================
        STEP 2: QUALITY CONTROL
    ================================================================================
    */
    logStep("QUALITY CONTROL", "START")
    
    if (!skip_alignment) {
        // Run NanoPlot on FASTQ files
        NANOPLOT(ch_fastq)
        ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.report)
        ch_versions = ch_versions.mix(NANOPLOT.out.versions)
    }
    
    logStep("QUALITY CONTROL", "DONE")
    
    /*
    ================================================================================
        STEP 3: ALIGNMENT
    ================================================================================
    */
    if (!skip_alignment) {
        logStep("ALIGNMENT", "START")
        
        // Align with minimap2
        MINIMAP2_ALIGN(
            ch_fastq,
            ch_reference.collect(),
            ch_reference_fai.collect()
        )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
        
        // Sort and index
        SAMTOOLS_SORT_INDEX(MINIMAP2_ALIGN.out.bam)
        ch_aligned_bam = SAMTOOLS_SORT_INDEX.out.bam
        ch_versions = ch_versions.mix(SAMTOOLS_SORT_INDEX.out.versions)
        
        logStep("ALIGNMENT", "DONE")
    }
    
    /*
    ================================================================================
        STEP 4: ALIGNMENT QC
    ================================================================================
    */
    logStep("ALIGNMENT QC", "START")
    
    // Calculate coverage with mosdepth
    MOSDEPTH(
        ch_aligned_bam,
        ch_reference.collect()
    )
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.report)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
    
    logStep("ALIGNMENT QC", "DONE")
    
    /*
    ================================================================================
        STEP 5: SMALL VARIANT CALLING (SNVs/Indels)
    ================================================================================
    */
    logStep("SMALL VARIANT CALLING", "START")
    
    if (params.variant_caller == 'clair3') {
        CLAIR3(
            ch_aligned_bam,
            ch_reference.collect(),
            ch_reference_fai.collect(),
            params.clair3_model
        )
        ch_small_variants = CLAIR3.out.vcf
        ch_versions = ch_versions.mix(CLAIR3.out.versions)
    } else {
        PEPPER_DEEPVARIANT(
            ch_aligned_bam,
            ch_reference.collect(),
            ch_reference_fai.collect()
        )
        ch_small_variants = PEPPER_DEEPVARIANT.out.vcf
        ch_versions = ch_versions.mix(PEPPER_DEEPVARIANT.out.versions)
    }
    
    logStep("SMALL VARIANT CALLING", "DONE")
    
    /*
    ================================================================================
        STEP 6: STRUCTURAL VARIANT CALLING
    ================================================================================
    */
    if (params.call_svs) {
        logStep("STRUCTURAL VARIANT CALLING", "START")
        
        // Sniffles2 for SV calling
        SNIFFLES2(
            ch_aligned_bam,
            ch_reference.collect(),
            ch_reference_fai.collect()
        )
        ch_sv_variants = SNIFFLES2.out.vcf
        ch_versions = ch_versions.mix(SNIFFLES2.out.versions)
        
        // Optional: CuteSV as secondary caller
        if (params.sv_caller_secondary) {
            CUTESV(
                ch_aligned_bam,
                ch_reference.collect()
            )
            ch_versions = ch_versions.mix(CUTESV.out.versions)
        }
        
        logStep("STRUCTURAL VARIANT CALLING", "DONE")
    } else {
        ch_sv_variants = Channel.empty()
    }
    
    /*
    ================================================================================
        STEP 7: COPY NUMBER VARIANT CALLING
    ================================================================================
    */
    if (params.call_cnvs) {
        logStep("CNV CALLING", "START")
        
        SPECTRE(
            ch_aligned_bam,
            ch_reference.collect()
        )
        ch_cnv_variants = SPECTRE.out.vcf
        ch_versions = ch_versions.mix(SPECTRE.out.versions)
        
        logStep("CNV CALLING", "DONE")
    } else {
        ch_cnv_variants = Channel.empty()
    }
    
    /*
    ================================================================================
        STEP 8: HAPLOTYPE PHASING
    ================================================================================
    */
    if (params.phase_variants) {
        logStep("HAPLOTYPE PHASING", "START")
        
        // Phase variants with WhatsHap
        WHATSHAP_PHASE(
            ch_small_variants.join(ch_aligned_bam),
            ch_reference.collect(),
            ch_reference_fai.collect()
        )
        ch_phased_vcf = WHATSHAP_PHASE.out.vcf
        ch_versions = ch_versions.mix(WHATSHAP_PHASE.out.versions)
        
        // Haplotag BAM reads
        WHATSHAP_HAPLOTAG(
            ch_phased_vcf.join(ch_aligned_bam),
            ch_reference.collect(),
            ch_reference_fai.collect()
        )
        ch_haplotagged_bam = WHATSHAP_HAPLOTAG.out.bam
        ch_versions = ch_versions.mix(WHATSHAP_HAPLOTAG.out.versions)
        
        // Generate phasing statistics
        WHATSHAP_STATS(ch_phased_vcf)
        ch_multiqc_files = ch_multiqc_files.mix(WHATSHAP_STATS.out.stats)
        ch_versions = ch_versions.mix(WHATSHAP_STATS.out.versions)
        
        // Use phased VCF for downstream
        ch_variants_for_annotation = ch_phased_vcf
        
        logStep("HAPLOTYPE PHASING", "DONE")
    } else {
        ch_variants_for_annotation = ch_small_variants
        ch_haplotagged_bam = ch_aligned_bam
    }
    
    /*
    ================================================================================
        STEP 9: METHYLATION CALLING
    ================================================================================
    */
    if (params.call_methylation) {
        logStep("METHYLATION CALLING", "START")
        
        // Use haplotagged BAM if available, otherwise aligned BAM
        ch_bam_for_methyl = params.phase_variants ? ch_haplotagged_bam : ch_aligned_bam
        
        // Methylation pileup with modkit
        MODKIT_PILEUP(
            ch_bam_for_methyl,
            ch_reference.collect(),
            ch_reference_fai.collect()
        )
        ch_methylation = MODKIT_PILEUP.out.bedmethyl
        ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions)
        
        // Methylation summary
        MODKIT_SUMMARY(ch_bam_for_methyl)
        ch_multiqc_files = ch_multiqc_files.mix(MODKIT_SUMMARY.out.summary)
        ch_versions = ch_versions.mix(MODKIT_SUMMARY.out.versions)
        
        logStep("METHYLATION CALLING", "DONE")
    }
    
    /*
    ================================================================================
        STEP 10: VARIANT ANNOTATION
    ================================================================================
    */
    if (params.annotate_variants) {
        logStep("VARIANT ANNOTATION", "START")
        
        // Start with unannotated variants
        ch_annotated = ch_variants_for_annotation
        
        // VEP base annotation (only if cache is provided)
        if (params.vep_cache) {
            ch_vep_cache = Channel.fromPath(params.vep_cache, checkIfExists: true).collect()
            VEP_ANNOTATE(
                ch_variants_for_annotation,
                ch_reference.collect(),
                ch_vep_cache
            )
            ch_annotated = VEP_ANNOTATE.out.vcf
            ch_versions = ch_versions.mix(VEP_ANNOTATE.out.versions)
        } else {
            log.warn "[WARN] VEP cache not provided (--vep_cache) - skipping VEP annotation"
        }
        
        // ClinVar annotation
        if (params.clinvar_vcf) {
            CLINVAR_ANNOTATE(
                ch_annotated,
                params.clinvar_vcf
            )
            ch_annotated = CLINVAR_ANNOTATE.out.vcf
            ch_versions = ch_versions.mix(CLINVAR_ANNOTATE.out.versions)
        }
        
        // gnomAD annotation
        if (params.gnomad_vcf) {
            GNOMAD_ANNOTATE(
                ch_annotated,
                params.gnomad_vcf
            )
            ch_annotated = GNOMAD_ANNOTATE.out.vcf
            ch_versions = ch_versions.mix(GNOMAD_ANNOTATE.out.versions)
        }
        
        // dbSNP annotation
        if (params.dbsnp_vcf) {
            DBSNP_ANNOTATE(
                ch_annotated,
                params.dbsnp_vcf
            )
            ch_annotated = DBSNP_ANNOTATE.out.vcf
            ch_versions = ch_versions.mix(DBSNP_ANNOTATE.out.versions)
        }
        
        // CADD annotation
        if (params.cadd_snvs && params.cadd_indels) {
            CADD_ANNOTATE(
                ch_annotated,
                params.cadd_snvs,
                params.cadd_indels
            )
            ch_annotated = CADD_ANNOTATE.out.vcf
            ch_versions = ch_versions.mix(CADD_ANNOTATE.out.versions)
        }
        
        // SpliceAI annotation
        if (params.spliceai_snvs && params.spliceai_indels) {
            SPLICEAI_ANNOTATE(
                ch_annotated,
                params.spliceai_snvs,
                params.spliceai_indels
            )
            ch_annotated = SPLICEAI_ANNOTATE.out.vcf
            ch_versions = ch_versions.mix(SPLICEAI_ANNOTATE.out.versions)
        }
        
        // PharmGKB annotation
        if (params.pharmgkb_vcf) {
            PHARMGKB_ANNOTATE(
                ch_annotated,
                params.pharmgkb_vcf
            )
            ch_annotated = PHARMGKB_ANNOTATE.out.vcf
            ch_versions = ch_versions.mix(PHARMGKB_ANNOTATE.out.versions)
        }
        
        // COSMIC annotation (for cancer predisposition)
        if (params.cosmic_vcf) {
            COSMIC_ANNOTATE(
                ch_annotated,
                params.cosmic_vcf
            )
            ch_annotated = COSMIC_ANNOTATE.out.vcf
            ch_versions = ch_versions.mix(COSMIC_ANNOTATE.out.versions)
        }
        
        ch_final_vcf = ch_annotated
        
        logStep("VARIANT ANNOTATION", "DONE")
    } else {
        ch_final_vcf = ch_variants_for_annotation
    }
    
    /*
    ================================================================================
        STEP 11: REPORTING
    ================================================================================
    */
    logStep("GENERATING REPORTS", "START")
    
    // Collect all versions
    ch_versions
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")
    
    // MultiQC report
    ch_multiqc_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true).collect() : Channel.value([])
    MULTIQC(
        ch_multiqc_files.collect().ifEmpty([]),
        ch_multiqc_config
    )
    
    // Generate final HTML report
    GENERATE_REPORT(
        ch_final_vcf,
        ch_sv_variants.ifEmpty([]),
        ch_cnv_variants.ifEmpty([]),
        MOSDEPTH.out.summary.collect()
    )
    
    logStep("GENERATING REPORTS", "DONE")
    
    /*
    ================================================================================
        COMPLETION
    ================================================================================
    */
    workflow.onComplete {
        def duration = workflow.duration
        def success = workflow.success
        
        log.info """
╔═══════════════════════════════════════════════════════════════════════════════════╗
║                           PIPELINE COMPLETED                                      ║
╠═══════════════════════════════════════════════════════════════════════════════════╣
║  Status      : ${success ? '✅ SUCCESS' : '❌ FAILED'}                                              ║
║  Duration    : ${duration}
║  Output Dir  : ${params.outdir}
║  Work Dir    : ${workflow.workDir}
╚═══════════════════════════════════════════════════════════════════════════════════╝

    Output Files:
    =============
    📁 ${params.outdir}/
    ├── 📁 qc/                    - Quality control reports
    ├── 📁 alignment/             - Aligned BAM files
    ├── 📁 variants/
    │   ├── 📁 small_variants/    - SNVs and Indels (VCF)
    │   ├── 📁 structural/        - Structural variants (VCF)
    │   └── 📁 cnv/               - Copy number variants
    ├── 📁 phasing/               - Phased variants and haplotagged BAMs
    ├── 📁 methylation/           - Methylation calls (BED)
    ├── 📁 annotation/            - Annotated VCF files
    ├── 📁 reports/               - MultiQC and summary reports
    └── 📁 pipeline_info/         - Execution logs and versions

    Thank you for using wf-ont-germline!
"""
    }
    
    workflow.onError {
        log.error """
╔═══════════════════════════════════════════════════════════════════════════════════╗
║                              PIPELINE ERROR                                       ║
╠═══════════════════════════════════════════════════════════════════════════════════╣
║  Error Message: ${workflow.errorMessage}
║  Exit Status  : ${workflow.exitStatus}
║  Work Dir     : ${workflow.workDir}
╚═══════════════════════════════════════════════════════════════════════════════════╝

    Please check the log files in the work directory for more details.
    Common issues:
    - Insufficient memory: Try increasing --max_memory
    - Missing input files: Check --input and --reference paths
    - Container issues: Ensure Docker/Singularity is running
"""
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/
