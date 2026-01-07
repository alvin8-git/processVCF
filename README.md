# processVCF Pipeline

A self-contained VCF processing pipeline for TMSP and CEBPA/CEBNX panels. This pipeline automates variant annotation, filtering, and Excel report generation.

## Overview

This pipeline:
- Processes VCF files from the TMSP panel and optionally CEBPA/CEBNX panel
- Runs multiple annotation tools in parallel (ANNOVAR, VEP, snpEff, TransVar, CancerVar)
- Compares variants against reference databases (Clinical, AD, LRM/CAP, NA12878/MOLM14)
- Annotates with population frequencies (SG10K, GenomeAsia)
- Generates filtered Excel reports for each sample

## Directory Structure

```
processVCF/
├── processVCF.sh                    # Main pipeline script (includes vcf_stats & filter_anno)
├── mergeVCFannotation-optimized.sh  # Optimized annotation engine (includes Excel writers)
├── Dockerfile                       # Docker container build file
├── docker-compose.yml               # Docker Compose configuration
└── README.md                        # This file
```

## Prerequisites

### System Tools

The following tools must be installed and in your PATH:

| Tool | Purpose |
|------|---------|
| bcftools | VCF manipulation |
| parallel | GNU Parallel for job distribution |
| rename | Perl rename utility |
| perl | Perl interpreter (with Excel::Writer::XLSX) |
| awk | AWK interpreter |
| java | Java runtime (for snpEff) |
| python3 | Python 3 (for CancerVar) |
| transvar | HGVS annotation |
| vep | Ensembl Variant Effect Predictor |
| bgzip | Block gzip compression |
| tabix | TAB-delimited file indexing |
| vcf-merge | VCF merging (vcftools) |
| vcf-sort | VCF sorting (vcftools) |

### Perl Modules

```bash
cpan install Excel::Writer::XLSX
```

### Software Directories

These should be installed in `$HOME/Software/`:
- `annovar/` - ANNOVAR annotation tool
- `snpEff/` - snpEff annotation tool
- `CancerVar/` - CancerVar annotation tool

### Databases

These should be present in `$HOME/Databases/`:
- `humandb/` - ANNOVAR databases
- `vep/` - VEP cache files
- `WholeGenomeFASTA/` - Reference genome (genome.fa, hg19.sort.fa)
- `Ensembldata/RefSeqSelectTranscript.txt` - RefSeq transcripts
- `SG10K.hg37.vcf/SG10K.TMSPgenes.txt` - SG10K population data
- `genomeAsia/genomeAsia.TMSP54.txt` - GenomeAsia population data
- `TMSPvcf/TSMPclean/` - TMSP reference VCFs (Clinical, AD, LRM, NA12878)
- `TMSPvcf/CEBPA/` - CEBPA reference VCFs (Clinical, AD, CAP, MOLM14)
- `TMSPvcf/TMSPcommon/TMSPcommonlist.txt` - Common variants list

## Installation

### Option 1: Docker (Recommended)

1. Clone the repository:
   ```bash
   git clone https://github.com/alvin8-git/processVCF.git
   cd processVCF
   ```

2. Build the Docker container:
   ```bash
   docker build -t processvcf .
   ```

3. Run with your data and databases:
   ```bash
   docker run -v /path/to/Databases:/home/user/Databases \
              -v /path/to/vcf:/data \
              -it processvcf bash
   ```

### Option 2: Native Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/alvin8-git/processVCF.git
   cd processVCF
   ```

2. Make scripts executable:
   ```bash
   chmod +x *.sh
   ```

3. Verify dependencies:
   ```bash
   ./processVCF.sh --check
   ```

## Usage

### Expected Input Directory Structure

```
/path/to/analysis/
├── vcf/                  <- Run script from here
│   └── *.vcf            <- TMSP VCF files
├── cebpa/
│   └── vcf/
│       └── *.vcf        <- CEBPA VCF files (optional)
├── bam/                  <- BAM files (optional, for IGV)
└── output/               <- Created by script
```

### Running the Pipeline

1. Navigate to the VCF directory:
   ```bash
   cd /path/to/analysis/vcf
   ```

2. Run the pipeline:
   ```bash
   /path/to/processVCF/processVCF.sh
   ```

### Command Line Options

```
Usage: processVCF.sh [options]

Options:
  --check       Check all dependencies and databases
  --help        Show help message
```

## Output

### TMSP Output (`output/annotationTMSP/`)

Per-sample files:
- `*.annotation.txt` - Full annotation table
- `*.Filter.txt` - Filtered variants
- `*.compare.txt` - Comparison with reference databases
- `*.compareVAF.txt` - VAF comparison
- `*.annoCheck.txt` - Annotation concordance check
- `*.VCFstats.txt` - VCF statistics
- `*.xlsx` - Excel report

Merged files:
- `Merge.AnnoAll.txt` - Combined annotations
- `Merge.compare.txt` - Combined comparisons
- `Summary.txt` - Summary of all filtered variants
- `Summary.xlsx` - Summary Excel file

### CEBPA Output (`output/annotationCEBNX/`)

Same structure as TMSP, but:
- No QUAL>=100 filter applied
- Comparison databases: Clinical, AD, CAP, MOLM14
- No Summary.txt/Summary.xlsx generated

## Processing Modes

### TMSP Mode
- Filters: QUAL>=100 AND (PASS OR SB)
- Comparisons: Clinical, AD, LRM, NA12878
- Includes IKZF1/WT1 transcript fix
- Generates Summary files

### CEBPA Mode
- Filters: PASS OR SB (no QUAL filter)
- Comparisons: Clinical, AD, CAP, MOLM14
- No IKZF1/WT1 fix
- No Summary files

## Integrated Functions

The following functions are built into the scripts (no external dependencies):

### vcf_stats
Extracts statistics from VCF files:
- Chromosome, position, ref, alt
- Genotype (GT), allele depth (AD), read depth (DP)
- Quality score (QUAL), variant allele frequency (VAF)

### filter_anno
Filters annotation results based on:
- Remove CUX1 variants
- Remove concordant intronic/synonymous/UTR variants
- Remove variants seen in >10 NA12878 or >15 LRM samples
- Remove variants with DP<100 or VAF<5%
- Keep exception variants (ASXL1 G646W, TP53/JAK2/FLT3/CALR with low VAF)
- Remove common variants from reference list

## Performance

The optimized pipeline runs annotation tools in parallel:
- ANNOVAR, VEP, snpEff, TransVar run simultaneously
- All 4 database comparisons run in parallel
- SG10K and GenomeAsia lookups run in parallel

Typical runtime: ~10-15 minutes for 4-5 samples (vs ~18+ minutes sequential)

## Troubleshooting

### Check Dependencies
```bash
./processVCF.sh --check
```

### Common Issues

1. **bcftools not found**: Install bcftools via conda or package manager
2. **parallel not found**: Install GNU parallel (`apt install parallel`)
3. **VEP errors**: Ensure VEP cache is installed for GRCh37
4. **Excel::Writer::XLSX missing**: Run `cpan install Excel::Writer::XLSX`

## Files Description

| File | Description |
|------|-------------|
| `processVCF.sh` | Main orchestration script with integrated vcf_stats and filter_anno functions |
| `mergeVCFannotation-optimized.sh` | Core annotation engine with parallel execution, integrated Excel writers |
| `Dockerfile` | Docker container build file with all required tools |
| `docker-compose.yml` | Docker Compose configuration for easy deployment |

## Docker Container

The Docker container includes:
- bcftools, tabix, vcftools
- GNU Parallel
- Perl with Excel::Writer::XLSX
- Java (OpenJDK 11)
- Python 3 with TransVar
- ANNOVAR (scripts only)
- snpEff 5.0e
- CancerVar

**Note:** VEP and databases must be mounted at runtime due to their large size.

### Required Volume Mounts

```bash
docker run \
  -v $HOME/Databases:/home/user/Databases:ro \
  -v $HOME/Software/ensembl-vep:/home/user/Software/ensembl-vep:ro \
  -v /path/to/vcf:/data \
  -it processvcf bash
```

## License

MIT License

## Version History

- v2.2 (2025-01): Docker support, integrated Excel writers into main scripts
- v2.1 (2025-01): Integrated VCFstats and filterAnno functions into main scripts
- v2.0 (2025-01): Optimized parallel pipeline, combined TMSP/CEBPA support
- v1.0: Original AWK-based implementation

## Author

Alvin Ng
