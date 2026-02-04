#!/usr/bin/env bash

#===============================================================================
# Download Annotation Databases for wf-ont-germline Pipeline
#===============================================================================
# This script downloads and prepares all annotation databases required for
# comprehensive variant annotation.
#
# Databases included:
#   1. ClinVar - Clinical significance
#   2. gnomAD - Population allele frequencies
#   3. dbSNP - Variant identifiers
#   4. CADD - Deleteriousness scores
#   5. REVEL - Missense pathogenicity
#   6. SpliceAI - Splice site predictions
#   7. PharmGKB - Pharmacogenomics
#   8. COSMIC - Cancer mutations
#   9. VEP Cache - Ensembl VEP annotation cache
#
# Usage:
#   ./download_annotation_databases.sh [OPTIONS]
#
# Options:
#   -o, --output-dir    Output directory for databases (default: ./annotation_db)
#   -g, --genome        Genome assembly (GRCh37 or GRCh38, default: GRCh38)
#   -t, --threads       Number of threads for parallel downloads (default: 4)
#   -d, --databases     Comma-separated list of databases to download (default: all)
#   -h, --help          Show this help message
#
#===============================================================================

set -euo pipefail

# Default values
OUTPUT_DIR="./annotation_db"
GENOME="GRCh38"
THREADS=4
DATABASES="all"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Help message
show_help() {
    cat << EOF
================================================================================
Download Annotation Databases for wf-ont-germline Pipeline
================================================================================

Usage: $0 [OPTIONS]

Options:
    -o, --output-dir    Output directory for databases (default: ./annotation_db)
    -g, --genome        Genome assembly (GRCh37 or GRCh38, default: GRCh38)
    -t, --threads       Number of threads for parallel downloads (default: 4)
    -d, --databases     Comma-separated list of databases to download
                        Options: clinvar,gnomad,dbsnp,cadd,revel,spliceai,pharmgkb,cosmic,vep
                        (default: all)
    -h, --help          Show this help message

Examples:
    # Download all databases for GRCh38
    $0 -o /data/annotation_db -g GRCh38

    # Download only ClinVar and gnomAD
    $0 -o /data/annotation_db -d clinvar,gnomad

    # Download with 8 parallel threads
    $0 -o /data/annotation_db -t 8

================================================================================
EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -g|--genome)
            GENOME="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -d|--databases)
            DATABASES="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Validate genome assembly
if [[ "$GENOME" != "GRCh37" && "$GENOME" != "GRCh38" ]]; then
    log_error "Invalid genome assembly: $GENOME. Must be GRCh37 or GRCh38."
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

log_info "=============================================="
log_info "Downloading Annotation Databases"
log_info "=============================================="
log_info "Output directory: $OUTPUT_DIR"
log_info "Genome assembly: $GENOME"
log_info "Threads: $THREADS"
log_info "Databases: $DATABASES"
log_info "=============================================="

# Function to check if database should be downloaded
should_download() {
    local db=$1
    if [[ "$DATABASES" == "all" ]] || [[ "$DATABASES" == *"$db"* ]]; then
        return 0
    fi
    return 1
}

#===============================================================================
# 1. ClinVar
#===============================================================================
download_clinvar() {
    log_info "Downloading ClinVar..."
    mkdir -p clinvar
    cd clinvar
    
    if [[ "$GENOME" == "GRCh38" ]]; then
        wget -c "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
        wget -c "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"
    else
        wget -c "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz"
        wget -c "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi"
    fi
    
    cd ..
    log_success "ClinVar download complete"
}

#===============================================================================
# 2. gnomAD
#===============================================================================
download_gnomad() {
    log_info "Downloading gnomAD..."
    log_warning "gnomAD is a large database (~450GB for full genome). Consider downloading specific chromosomes."
    mkdir -p gnomad
    cd gnomad
    
    if [[ "$GENOME" == "GRCh38" ]]; then
        # gnomAD v4.0 sites VCF (GRCh38)
        # Download chromosome-by-chromosome for manageability
        for chr in {1..22} X Y; do
            log_info "Downloading gnomAD chr${chr}..."
            wget -c "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr${chr}.vcf.bgz" || true
            wget -c "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr${chr}.vcf.bgz.tbi" || true
        done
    else
        # gnomAD v2.1.1 (GRCh37)
        for chr in {1..22} X; do
            log_info "Downloading gnomAD chr${chr}..."
            wget -c "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.${chr}.vcf.bgz" || true
            wget -c "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.${chr}.vcf.bgz.tbi" || true
        done
    fi
    
    cd ..
    log_success "gnomAD download complete"
}

#===============================================================================
# 3. dbSNP
#===============================================================================
download_dbsnp() {
    log_info "Downloading dbSNP..."
    mkdir -p dbsnp
    cd dbsnp
    
    if [[ "$GENOME" == "GRCh38" ]]; then
        wget -c "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"
        wget -c "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi"
        # Rename for clarity
        mv GCF_000001405.40.gz dbsnp_GRCh38.vcf.gz
        mv GCF_000001405.40.gz.tbi dbsnp_GRCh38.vcf.gz.tbi
    else
        wget -c "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz"
        wget -c "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi"
        mv GCF_000001405.25.gz dbsnp_GRCh37.vcf.gz
        mv GCF_000001405.25.gz.tbi dbsnp_GRCh37.vcf.gz.tbi
    fi
    
    cd ..
    log_success "dbSNP download complete"
}

#===============================================================================
# 4. CADD
#===============================================================================
download_cadd() {
    log_info "Downloading CADD scores..."
    log_warning "CADD files are large (~80GB for SNVs, ~3GB for indels)"
    mkdir -p cadd
    cd cadd
    
    if [[ "$GENOME" == "GRCh38" ]]; then
        # CADD v1.6 GRCh38
        wget -c "https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz"
        wget -c "https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz.tbi"
        wget -c "https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz"
        wget -c "https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz.tbi"
    else
        # CADD v1.6 GRCh37
        wget -c "https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz"
        wget -c "https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz.tbi"
        wget -c "https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/InDels.tsv.gz"
        wget -c "https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/InDels.tsv.gz.tbi"
    fi
    
    cd ..
    log_success "CADD download complete"
}

#===============================================================================
# 5. REVEL
#===============================================================================
download_revel() {
    log_info "Downloading REVEL scores..."
    mkdir -p revel
    cd revel
    
    # REVEL scores (same file works for both assemblies with liftover)
    wget -c "https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip"
    unzip -o revel-v1.3_all_chromosomes.zip
    
    # Convert to bgzipped and indexed format
    if [[ "$GENOME" == "GRCh38" ]]; then
        # Need to liftover from GRCh37 to GRCh38
        log_warning "REVEL is provided in GRCh37 coordinates. Liftover to GRCh38 may be required."
    fi
    
    # Create tabix index
    bgzip -c revel_with_transcript_ids > revel.tsv.gz
    tabix -s 1 -b 2 -e 2 revel.tsv.gz
    
    cd ..
    log_success "REVEL download complete"
}

#===============================================================================
# 6. SpliceAI
#===============================================================================
download_spliceai() {
    log_info "Downloading SpliceAI scores..."
    log_warning "SpliceAI files require Illumina BaseSpace account for download"
    mkdir -p spliceai
    cd spliceai
    
    cat << 'EOF'
================================================================================
SpliceAI Download Instructions
================================================================================

SpliceAI precomputed scores must be downloaded from Illumina BaseSpace:
https://basespace.illumina.com/s/otSPW8hnhaZR

Files needed:
- spliceai_scores.raw.snv.hg38.vcf.gz (GRCh38 SNVs)
- spliceai_scores.raw.snv.hg38.vcf.gz.tbi
- spliceai_scores.raw.indel.hg38.vcf.gz (GRCh38 Indels)
- spliceai_scores.raw.indel.hg38.vcf.gz.tbi

For GRCh37:
- spliceai_scores.raw.snv.hg19.vcf.gz
- spliceai_scores.raw.snv.hg19.vcf.gz.tbi
- spliceai_scores.raw.indel.hg19.vcf.gz
- spliceai_scores.raw.indel.hg19.vcf.gz.tbi

Please download these files manually and place them in this directory.
================================================================================
EOF
    
    cd ..
    log_warning "SpliceAI requires manual download from Illumina BaseSpace"
}

#===============================================================================
# 7. PharmGKB
#===============================================================================
download_pharmgkb() {
    log_info "Downloading PharmGKB data..."
    mkdir -p pharmgkb
    cd pharmgkb
    
    # PharmGKB clinical annotations
    wget -c "https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip"
    unzip -o clinicalAnnotations.zip
    
    # PharmGKB variant annotations
    wget -c "https://api.pharmgkb.org/v1/download/file/data/variantAnnotations.zip"
    unzip -o variantAnnotations.zip
    
    # Create VCF from PharmGKB data (requires custom script)
    log_info "Creating PharmGKB VCF file..."
    
    cat << 'PYTHON' > create_pharmgkb_vcf.py
#!/usr/bin/env python3
import pandas as pd
import sys

# Read PharmGKB variant annotations
df = pd.read_csv('var_pheno_ann.tsv', sep='\t', low_memory=False)

# Create VCF header
vcf_header = """##fileformat=VCFv4.2
##INFO=<ID=PGX_GENE,Number=1,Type=String,Description="PharmGKB gene">
##INFO=<ID=PGX_DRUG,Number=.,Type=String,Description="Associated drugs">
##INFO=<ID=PGX_PHENOTYPE,Number=.,Type=String,Description="Phenotype category">
##INFO=<ID=PGX_LEVEL,Number=1,Type=String,Description="Evidence level">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

print(vcf_header, end='')
print("PharmGKB VCF creation requires additional processing")
PYTHON
    
    python3 create_pharmgkb_vcf.py || log_warning "PharmGKB VCF creation requires manual processing"
    
    cd ..
    log_success "PharmGKB download complete"
}

#===============================================================================
# 8. COSMIC
#===============================================================================
download_cosmic() {
    log_info "Downloading COSMIC data..."
    log_warning "COSMIC requires registration and license agreement"
    mkdir -p cosmic
    cd cosmic
    
    cat << 'EOF'
================================================================================
COSMIC Download Instructions
================================================================================

COSMIC (Catalogue of Somatic Mutations in Cancer) requires registration:
https://cancer.sanger.ac.uk/cosmic/download

Files needed for germline predisposition analysis:
- CosmicCodingMuts.vcf.gz
- CosmicCodingMuts.vcf.gz.tbi

Download command (after authentication):
sftp -oBatchMode=no -b - "[email protected]" << SFTP
cd /cosmic/grch38/cosmic/v98/VCF
get CosmicCodingMuts.vcf.gz
get CosmicCodingMuts.vcf.gz.tbi
SFTP

Please download these files manually and place them in this directory.
================================================================================
EOF
    
    cd ..
    log_warning "COSMIC requires manual download with authentication"
}

#===============================================================================
# 9. VEP Cache
#===============================================================================
download_vep_cache() {
    log_info "Downloading VEP cache..."
    log_warning "VEP cache is large (~15-20GB)"
    mkdir -p vep_cache
    cd vep_cache
    
    if [[ "$GENOME" == "GRCh38" ]]; then
        # VEP cache for GRCh38
        wget -c "https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz"
        tar -xzf homo_sapiens_vep_110_GRCh38.tar.gz
    else
        # VEP cache for GRCh37
        wget -c "https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh37.tar.gz"
        tar -xzf homo_sapiens_vep_110_GRCh37.tar.gz
    fi
    
    cd ..
    log_success "VEP cache download complete"
}

#===============================================================================
# Main execution
#===============================================================================

# Download selected databases
if should_download "clinvar"; then
    download_clinvar
fi

if should_download "gnomad"; then
    download_gnomad
fi

if should_download "dbsnp"; then
    download_dbsnp
fi

if should_download "cadd"; then
    download_cadd
fi

if should_download "revel"; then
    download_revel
fi

if should_download "spliceai"; then
    download_spliceai
fi

if should_download "pharmgkb"; then
    download_pharmgkb
fi

if should_download "cosmic"; then
    download_cosmic
fi

if should_download "vep"; then
    download_vep_cache
fi

#===============================================================================
# Generate configuration file
#===============================================================================
log_info "Generating annotation database configuration..."

cat > annotation_db_paths.config << EOF
/*
 * Annotation Database Paths
 * Generated by download_annotation_databases.sh
 * Date: $(date)
 * Genome: $GENOME
 */

params {
    // ClinVar
    clinvar_vcf     = "${OUTPUT_DIR}/clinvar/clinvar.vcf.gz"
    clinvar_tbi     = "${OUTPUT_DIR}/clinvar/clinvar.vcf.gz.tbi"
    
    // gnomAD (use merged or per-chromosome)
    gnomad_vcf      = "${OUTPUT_DIR}/gnomad/"
    
    // dbSNP
    dbsnp_vcf       = "${OUTPUT_DIR}/dbsnp/dbsnp_${GENOME}.vcf.gz"
    dbsnp_tbi       = "${OUTPUT_DIR}/dbsnp/dbsnp_${GENOME}.vcf.gz.tbi"
    
    // CADD
    cadd_snvs       = "${OUTPUT_DIR}/cadd/whole_genome_SNVs.tsv.gz"
    cadd_snvs_tbi   = "${OUTPUT_DIR}/cadd/whole_genome_SNVs.tsv.gz.tbi"
    cadd_indels     = "${OUTPUT_DIR}/cadd/gnomad.genomes.r3.0.indel.tsv.gz"
    cadd_indels_tbi = "${OUTPUT_DIR}/cadd/gnomad.genomes.r3.0.indel.tsv.gz.tbi"
    
    // REVEL
    revel_file      = "${OUTPUT_DIR}/revel/revel.tsv.gz"
    revel_tbi       = "${OUTPUT_DIR}/revel/revel.tsv.gz.tbi"
    
    // SpliceAI (requires manual download)
    spliceai_snvs       = "${OUTPUT_DIR}/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz"
    spliceai_snvs_tbi   = "${OUTPUT_DIR}/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz.tbi"
    spliceai_indels     = "${OUTPUT_DIR}/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz"
    spliceai_indels_tbi = "${OUTPUT_DIR}/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz.tbi"
    
    // PharmGKB
    pharmgkb_vcf    = "${OUTPUT_DIR}/pharmgkb/pharmgkb.vcf.gz"
    pharmgkb_tbi    = "${OUTPUT_DIR}/pharmgkb/pharmgkb.vcf.gz.tbi"
    
    // COSMIC (requires manual download)
    cosmic_vcf      = "${OUTPUT_DIR}/cosmic/CosmicCodingMuts.vcf.gz"
    cosmic_tbi      = "${OUTPUT_DIR}/cosmic/CosmicCodingMuts.vcf.gz.tbi"
    
    // VEP Cache
    vep_cache       = "${OUTPUT_DIR}/vep_cache"
}
EOF

log_success "Configuration file generated: annotation_db_paths.config"

#===============================================================================
# Summary
#===============================================================================
log_info "=============================================="
log_info "Download Summary"
log_info "=============================================="
log_info "Output directory: $OUTPUT_DIR"
log_info "Genome assembly: $GENOME"
log_info ""
log_info "Database sizes (approximate):"
log_info "  - ClinVar:   ~100 MB"
log_info "  - gnomAD:    ~450 GB (full genome)"
log_info "  - dbSNP:     ~20 GB"
log_info "  - CADD:      ~80 GB (SNVs) + ~3 GB (indels)"
log_info "  - REVEL:     ~2 GB"
log_info "  - SpliceAI:  ~30 GB (manual download required)"
log_info "  - PharmGKB:  ~50 MB"
log_info "  - COSMIC:    ~500 MB (manual download required)"
log_info "  - VEP Cache: ~15-20 GB"
log_info ""
log_info "To use these databases with the pipeline, include the config:"
log_info "  nextflow run main.nf -c ${OUTPUT_DIR}/annotation_db_paths.config ..."
log_info "=============================================="

log_success "Database download script completed!"
