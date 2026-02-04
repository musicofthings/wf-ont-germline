# wf-ont-germline

**Comprehensive Germline Variant Calling Pipeline for Oxford Nanopore Long-Read WGS**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![EPI2ME Compatible](https://img.shields.io/badge/EPI2ME-Compatible-blue.svg)](https://labs.epi2me.io/)

## Overview

`wf-ont-germline` is a production-ready Nextflow pipeline for comprehensive germline variant analysis from Oxford Nanopore long-read whole genome sequencing data. It is fully compatible with the **EPI2ME Desktop** application for easy execution.

### Features

| Analysis | Tools | Description |
|----------|-------|-------------|
| **Quality Control** | NanoPlot, mosdepth | Read QC and coverage analysis |
| **Alignment** | minimap2 | Long-read alignment to reference |
| **SNV/Indel Calling** | Clair3, PEPPER-DeepVariant | Small variant detection |
| **Structural Variants** | Sniffles2, CuteSV | SV detection (≥50bp) |
| **Copy Number Variants** | Spectre | CNV detection |
| **Haplotype Phasing** | WhatsHap, LongPhase | Variant phasing and read haplotagging |
| **Methylation** | modkit | 5mC methylation calling |
| **Annotation** | VEP + 9 databases | Clinical variant annotation |

### Annotation Databases

The pipeline supports comprehensive variant annotation with clinically relevant databases:

1. **ClinVar** - Clinical significance classifications
2. **gnomAD** - Population allele frequencies
3. **dbSNP** - Variant identifiers (rsIDs)
4. **CADD** - Deleteriousness scores
5. **REVEL** - Missense pathogenicity predictions
6. **SpliceAI** - Splice site impact predictions
7. **PharmGKB** - Pharmacogenomic associations
8. **COSMIC** - Cancer mutation database
9. **OMIM** - Disease gene associations

## Quick Start

### Using EPI2ME Desktop (Recommended)

1. Open EPI2ME Desktop
2. Import this workflow
3. Select your input FASTQ/BAM files
4. Select reference genome
5. Configure analysis options
6. Click "Run"

### Using Command Line

```bash
# Basic usage
nextflow run main.nf \
    --input /path/to/reads.fastq.gz \
    --reference /path/to/GRCh38.fa \
    --outdir ./results \
    -profile docker

# With all features enabled
nextflow run main.nf \
    --input /path/to/reads.fastq.gz \
    --reference /path/to/GRCh38.fa \
    --outdir ./results \
    --variant_caller clair3 \
    --call_svs true \
    --call_cnvs true \
    --phase_variants true \
    --call_methylation true \
    --annotate_variants true \
    -profile docker

# With annotation databases
nextflow run main.nf \
    --input /path/to/reads.fastq.gz \
    --reference /path/to/GRCh38.fa \
    --outdir ./results \
    --vep_cache /path/to/vep_cache \
    --clinvar_vcf /path/to/clinvar.vcf.gz \
    --gnomad_vcf /path/to/gnomad.vcf.gz \
    -profile docker
```

## Input Requirements

### Required Inputs

| Input | Description | Format |
|-------|-------------|--------|
| `--input` | Sequencing reads or aligned BAM | FASTQ, FASTQ.GZ, BAM |
| `--reference` | Reference genome | FASTA (indexed) |

### Optional Inputs

| Input | Description |
|-------|-------------|
| `--sample_name` | Sample identifier |
| `--clair3_model` | Path to Clair3 model |
| `--vep_cache` | Path to VEP cache |
| Annotation databases | See annotation section |

## Output Structure

```
results/
├── qc/
│   ├── nanoplot/           # Read quality reports
│   ├── coverage/           # Coverage statistics
│   └── alignment_stats/    # Alignment metrics
├── alignment/
│   ├── *.sorted.bam        # Aligned reads
│   └── *.sorted.bam.bai    # BAM index
├── variants/
│   ├── small_variants/     # SNVs and indels (VCF)
│   ├── structural/         # Structural variants (VCF)
│   └── cnv/                # Copy number variants
├── phasing/
│   ├── *.phased.vcf.gz     # Phased variants
│   ├── *.haplotagged.bam   # Haplotagged reads
│   └── stats/              # Phasing statistics
├── methylation/
│   ├── *.bedmethyl.gz      # Methylation calls
│   └── summary/            # Methylation statistics
├── annotation/
│   ├── vep/                # VEP annotations
│   ├── final/              # Fully annotated VCF
│   ├── clinical/           # Filtered clinical variants
│   └── pharmacogenomics/   # PGx variants
├── reports/
│   ├── multiqc_report.html # MultiQC summary
│   └── *_final_report.html # Analysis report
└── pipeline_info/
    ├── execution_*.html    # Nextflow reports
    └── software_versions.* # Tool versions
```

## Configuration

### Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run with Docker containers |
| `singularity` | Run with Singularity containers |
| `epi2me` | Optimized for EPI2ME Desktop |
| `test` | Minimal test dataset |
| `slurm` | HPC with SLURM scheduler |

### Key Parameters

#### Analysis Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--variant_caller` | `clair3` | `clair3` or `pepper_deepvariant` |
| `--call_svs` | `true` | Call structural variants |
| `--call_cnvs` | `true` | Call copy number variants |
| `--phase_variants` | `true` | Phase variants |
| `--call_methylation` | `true` | Call methylation |
| `--annotate_variants` | `true` | Annotate variants |

#### Model Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--basecall_model` | `dna_r10.4.1_e8.2_400bps_sup@v4.3.0` | Basecalling model |
| `--clair3_model` | `/opt/models/ont_r10_dorado_sup_4khz` | Clair3 model path |

#### Resource Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | `16` | Maximum CPUs |
| `--max_memory` | `64.GB` | Maximum memory |
| `--max_time` | `240.h` | Maximum runtime |

## Setting Up Annotation Databases

Download annotation databases using the provided script:

```bash
# Download all databases for GRCh38
./bin/download_annotation_databases.sh \
    --output-dir /path/to/annotation_db \
    --genome GRCh38

# Download specific databases
./bin/download_annotation_databases.sh \
    --output-dir /path/to/annotation_db \
    --databases clinvar,gnomad,dbsnp,cadd
```

Then include the generated config:

```bash
nextflow run main.nf \
    --input reads.fastq.gz \
    --reference GRCh38.fa \
    -c /path/to/annotation_db/annotation_db_paths.config \
    -profile docker
```

## Resource Requirements

### Minimum Requirements

- **CPU**: 4 cores
- **Memory**: 16 GB RAM
- **Storage**: 100 GB (excluding databases)

### Recommended Requirements

- **CPU**: 16+ cores
- **Memory**: 64+ GB RAM
- **Storage**: 500+ GB (including databases)

### Runtime Estimates

| Coverage | Time (16 cores) |
|----------|-----------------|
| 10x WGS | 4-8 hours |
| 30x WGS | 12-24 hours |
| 50x WGS | 24-48 hours |

## Troubleshooting

### Common Issues

1. **Out of memory errors**
   - Increase `--max_memory`
   - Reduce parallelism with `-process.maxForks`

2. **Clair3 model not found**
   - Ensure `--clair3_model` points to correct model directory
   - Model must match basecalling chemistry

3. **Missing annotation databases**
   - Run `download_annotation_databases.sh`
   - Check file paths in config

4. **Container issues**
   - Ensure Docker/Singularity is running
   - Check container registry access

### Getting Help

- Check the [FAQ](docs/faq.md)
- Open an issue on GitHub
- Contact support

## Citations

If you use this pipeline, please cite:

- **Clair3**: Zheng et al. (2022) Nature Computational Science. DOI: 10.1038/s43588-022-00387-x
- **Sniffles2**: Smolka et al. (2024) Nature Biotechnology. DOI: 10.1038/s41587-023-02024-y
- **WhatsHap**: Martin et al. (2016) bioRxiv. DOI: 10.1101/085050
- **VEP**: McLaren et al. (2016) Genome Biology. DOI: 10.1186/s13059-016-0974-4

## License

MIT License

## Acknowledgments

- Oxford Nanopore Technologies
- EPI2ME Labs
- Nextflow community
- All tool developers
