# processVCF Pipeline

A self-contained VCF processing pipeline for TMSP and CEBPA/CEBNX panels. This pipeline automates variant annotation, filtering, Excel report generation, and interactive HTML variant reports.

## Overview

This pipeline:
- Processes VCF files from the TMSP panel and optionally CEBPA/CEBNX panel
- Runs multiple annotation tools in parallel (ANNOVAR, VEP, snpEff, TransVar, CancerVar)
- Compares variants against reference databases (Clinical, AD, LRM/CAP, NA12878/MOLM14)
- Annotates with population frequencies (SG10K, GenomeAsia)
- Generates filtered Excel reports for each sample
- Generates interactive HTML variant reports with clinical-grade visualization

## Directory Structure

```
processVCF/
├── processVCF.sh                    # Main pipeline script (includes vcf_stats & filter_anno)
├── mergeVCFannotation-optimized.sh  # Optimized annotation engine (includes Excel writers)
├── excel_to_html_report.py          # HTML variant report generator
├── make_IGV_snapshots.py            # IGV snapshot automator (Python 3)
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
| xvfb-run | Virtual X server for IGV (apt install xvfb) |
| openjdk-8 | Java 8 required for IGV 2.3.81 (apt install openjdk-8-jdk) |

### Perl Modules

```bash
cpan install Excel::Writer::XLSX
```

### Python Modules

```bash
pip install openpyxl
```

### Software Directories

These should be installed in `$HOME/Software/`:
- `annovar/` - ANNOVAR annotation tool
- `snpEff/` - snpEff annotation tool
- `CancerVar/` - CancerVar annotation tool
- `IGV-snapshot-automator/bin/IGV_2.3.81/igv.jar` - IGV for snapshots

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

### IGV Snapshots (`output/SnapShots/`)

Automatically generated IGV snapshots for filtered variants:
- PNG images for each variant position
- Created using xvfb-run (virtual X display)
- Requires BAM files in `../bam/` directory

```
SnapShots/
├── SAMPLE1-TP53-12345.png
├── SAMPLE1-FLT3-54321.png
└── ...
```

BED files used for snapshot generation are saved in `output/IgvBed/`.

### HTML Reports (`output/html_reports/`)

Interactive HTML variant reports for each sample:

```
html_reports/
├── AML-453-ASR-TMSP_S3.html    # Dashboard with variant table
├── variant_1.html              # Individual variant detail page
├── variant_2.html
└── ...
```

**Dashboard Features:**
- Summary statistics (total variants, genes affected, tier counts)
- Clickable stat cards for quick filtering
- Filterable variant table by gene, classification, or search
- Links to individual variant detail pages

**Variant Detail Pages:**
- Gene name and HGVS notation prominently displayed
- Color-coded CancerVar classification badges
- 5 collapsible panels with color-coded headers:
  - **Basic Variant Information** (blue) - Position, genotype, VAF, COSMIC, CancerVar
  - **Sample Comparison** (green) - Comparison with reference databases
  - **Additional Information** (slate) - Consequence, exonic function, cytoBand
  - **Population Databases** (indigo) - Frequencies grouped by database (gnomAD, ExAC, 1000G, etc.)
  - **Computational Predictions** (orange) - ACMG criteria chips, prediction scores, ClinVar

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

## Standalone IGV Snapshot Generation

The IGV snapshot script can be used independently:

```bash
# Generate snapshots for all variants in a BED file
python3 make_IGV_snapshots.py sample.bam -r variants.bed -o SnapShots -bin /path/to/igv.jar

# Use 4th BED column as snapshot filename
python3 make_IGV_snapshots.py sample.bam -r variants.bed -o SnapShots -bin igv.jar -nf4

# Only generate batchscript without running IGV
python3 make_IGV_snapshots.py sample.bam -r variants.bed -o SnapShots -bin igv.jar -nosnap

# Specify custom Java 8 path (required for IGV 2.3.81)
python3 make_IGV_snapshots.py sample.bam -r variants.bed -o SnapShots -bin igv.jar -java /path/to/java8
```

**Note:** IGV 2.3.81 requires Java 8. The script defaults to `/usr/lib/jvm/java-8-openjdk-amd64/bin/java`.

### IGV Command Line Options

```
positional arguments:
  input_files           Input files (BAM, bigwig, etc.)

optional arguments:
  -r REGION_FILE        BED file with regions (default: regions.bed)
  -g GENOME             Reference genome (default: hg19)
  -ht IMAGE_HEIGHT      Track height in pixels (default: 500)
  -o OUTDIR             Output directory (default: IGV_Snapshots)
  -bin IGV_JAR_BIN      Path to IGV jar file
  -java JAVA_PATH       Path to Java 8 executable (default: /usr/lib/jvm/java-8-openjdk-amd64/bin/java)
  -mem IGV_MEM          Memory for IGV in MB (default: 4000)
  -nosnap               Only write batchscript, don't run IGV
  -suffix SUFFIX        Filename suffix for snapshots
  -nf4                  Use 4th BED field as snapshot filename
  -onlysnap ONLYSNAP    Run existing batchscript file
  -s, --group-by-strand Group reads by strand
```

## Standalone HTML Report Generation

The HTML report generator can be used independently on any Excel file from the pipeline:

```bash
# Generate dashboard + variant pages
python3 excel_to_html_report.py sample.xlsx output_dir/

# Generate single-page report (all variants on one page)
python3 excel_to_html_report.py sample.xlsx --single-page

# Generate single-page to specific file
python3 excel_to_html_report.py sample.xlsx --single-page -o report.html
```

### Command Line Options

```
usage: excel_to_html_report.py [-h] [--single-page] [--output OUTPUT] excel_file [output_dir]

positional arguments:
  excel_file            Input Excel file (.xlsx)
  output_dir            Output directory for multi-page reports

optional arguments:
  --single-page, -s     Generate a single HTML page with all variants
  --output, -o OUTPUT   Output file path for single-page report
```

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
| `processVCF.sh` | Main orchestration script with integrated vcf_stats, filter_anno, HTML, and IGV functions |
| `mergeVCFannotation-optimized.sh` | Core annotation engine with parallel execution, integrated Excel writers |
| `excel_to_html_report.py` | Converts Excel reports to interactive HTML variant reports |
| `make_IGV_snapshots.py` | IGV snapshot automator for batch screenshot generation (Python 3) |
| `Dockerfile` | Docker container build file with all required tools |
| `docker-compose.yml` | Docker Compose configuration for easy deployment |

## Docker Container

The Docker container includes:
- bcftools, tabix, vcftools
- GNU Parallel
- Perl with Excel::Writer::XLSX
- Java (OpenJDK 8 for IGV, OpenJDK 11 for other tools)
- Python 3 with TransVar and openpyxl
- ANNOVAR (scripts only)
- snpEff 5.0e
- CancerVar
- IGV 2.3.81 with xvfb for headless snapshots
- HTML report generator

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

- v2.6 (2025-01): IGV snapshot automation
  - Added `make_IGV_snapshots.py` (Python 3) for batch IGV screenshots
  - Integrated IGV snapshot generation into processVCF.sh pipeline
  - Creates BED files from filtered variants automatically
  - Outputs to `output/SnapShots/` and `output/IgvBed/`
  - Uses xvfb-run for headless IGV execution
- v2.5 (2025-01): HTML variant report generation
  - Added `excel_to_html_report.py` for interactive HTML reports
  - Dashboard with filterable variant table
  - Individual variant pages with 5 color-coded panels
  - CancerVar tier badges, ACMG criteria chips, prediction scores
  - Population frequencies grouped by database
- v2.4 (2025-01): Added HTML report integration to processVCF.sh
- v2.3 (2025-01): Added ANNOVAR, snpEff, CancerVar for Docker build
- v2.2 (2025-01): Docker support, integrated Excel writers into main scripts
- v2.1 (2025-01): Integrated VCFstats and filterAnno functions into main scripts
- v2.0 (2025-01): Optimized parallel pipeline, combined TMSP/CEBPA support
- v1.0: Original AWK-based implementation

## Author

Alvin Ng
