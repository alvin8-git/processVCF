#!/bin/bash
# =============================================================================
# processVCF.sh - VCF Processing Pipeline for TMSP and CEBPA
# =============================================================================
# This script automates VCF processing for TMSP and CEBNX panels.
# It uses the optimized parallel annotation pipeline.
#
# Usage: Run from the vcf directory containing TMSP VCF files
#   cd /path/to/analysis/vcf
#   ~/Shared/SCRIPTS/claude/processVCF/processVCF.sh
#
# Expected directory structure:
#   /path/to/analysis/
#   ├── vcf/                  <- Run script from here
#   │   ├── *.vcf            <- TMSP VCF files
#   │   └── annotationTMSP/  <- Created by script
#   ├── cebpa/
#   │   └── vcf/
#   │       ├── *.vcf        <- CEBPA VCF files
#   │       └── annotationCEBNX/  <- Created by script
#   ├── bam/                  <- BAM files for IGV
#   └── output/               <- Final output directory
#
# =============================================================================

set -e  # Exit on error

# =============================================================================
# CONFIGURATION
# =============================================================================

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Database paths
DATABASE_TMSP="$HOME/Databases/TMSPvcf/TSMPclean/Clinical"
DATABASE_CEBPA="$HOME/Databases/TMSPvcf/CEBPA/Clinical"
COMMON_VARIANTS="$HOME/Databases/TMSPvcf/TMSPcommon/TMSPcommonlist.txt"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $1"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" >&2
}

# =============================================================================
# VCFstats FUNCTION (merged from VCFstats.awk)
# =============================================================================
# Extracts useful statistics from a VCF file
# Usage: vcf_stats <vcf_file>
# Output: <sample>.VCFstats.txt

vcf_stats() {
    local vcf="$1"
    local sampleName=$(basename "$vcf" .vcf)

    [ -f "$sampleName.VCFstats.txt" ] && return 0

    paste \
        <(cat \
            <(echo -e "chr\tpos\tRef\tAlt\tGT\tAD\tDP\tQUAL\tVAF") \
            <(bcftools query -f '%CHROM %POS %REF %ALT [%GT %AD %DP %QUAL ]\n' "$vcf" | \
            sed -r 's/ +/\t/g' | cut -f1-8 | sed 's/,/\t/g' | \
            awk '{if($8>0)print $0"\t"$7/$8 ;else print $0"\t-"}' | cut -f1-5,7- )) \
        <(cat \
            <(echo -e "RD,AD") \
            <(bcftools query -f '[%AD]\n' "$vcf")) \
        > "$sampleName.VCFstats.txt"
}

# Export function for use with parallel
export -f vcf_stats

# =============================================================================
# filterAnno FUNCTION (merged from filterAnno.awk)
# =============================================================================
# Filters annotation file based on various criteria
# Usage: filter_anno <annotation_file>
# Output: <sample>.Filter.txt

filter_anno() {
    local FILE="$1"
    local SAMPLE=$(basename "$FILE" .annotation.txt)
    local COMMON="${COMMON_VARIANTS:-$HOME/Databases/TMSPvcf/TMSPcommon/TMSPcommonlist.txt}"

    [ -f "$SAMPLE.Filter.txt" ] && return 0

    # EXCEPTION CASES
    # Extract the ASXL G646Wfs* variant if VAF > 5%
    awk -F"\t" '{if($14>0.05&&$1=="chr20"&&($2==31022441)) print}' "$FILE" \
        > "$SAMPLE.filterASXL1G646W.txt"

    # Extract all TP53 & JAK2 & FLT3 & CALR genes VAF<5%
    awk -F"\t" '{if(($14<0.05&&$6=="TP53")||($14<0.05&&$6=="JAK2")||($14<0.05&&$6=="FLT3")||($14<0.05&&$6=="CALR")) print}' "$FILE" \
        > "$SAMPLE.filterVF5.txt"

    # GENERAL CASES
    # Filter rules:
    # - Remove CUX1 variants
    # - Remove concordant intronic variants (ANNOVAR and VEP agree)
    # - Remove concordant synonymous variants
    # - Remove concordant UTR variants
    # - Remove variants seen in >10 NA12878 samples
    # - Remove variants seen in >15 LRM samples
    # - Remove variants with DP<100
    # - Remove variants with VAF<5%
    # - Keep splice acceptor/donor variants even if intronic

    cat \
        <(awk -F"\t" '{if(NR==1)print}' "$FILE") \
        <(cat \
            <(awk -F"\t" '{if(NR!=1)print}' "$FILE" | \
            awk -F"\t" '{if($6!~/CUX1/)print}' | \
            awk -F"\t" '{if(!($25=="intronic"&&$27=="intron_variant"))print}' | \
            awk -F"\t" '{if(!($25=="intronic"&&$27!~/splice_donor_variant|splice_acceptor_variant/))print}' | \
            awk -F"\t" '{if(!($24=="synonymous SNV"&&$27=="synonymous_variant"))print}' | \
            awk -F"\t" '{if(!($25=="UTR3"&&$27=="3_prime_UTR_variant"))print}' | \
            awk -F"\t" '{if(!($25=="UTR5"&&$27=="5_prime_UTR_variant"))print}' | \
            awk -F"\t" '{if($21<=10)print}' | \
            awk -F"\t" '{if($22<=15)print}' | \
            awk -F"\t" '{if($12>=100)print}' | \
            awk -F"\t" '{if($14>0.05)print}') \
            <(cat "$SAMPLE.filterASXL1G646W.txt") \
            <(cat "$SAMPLE.filterVF5.txt") | \
        sort -k1,1V -k2,2n ) \
        > "$SAMPLE.GFilter.txt"

    # REMOVE VARIANTS IN COMMON LIST
    if [ -f "$COMMON" ]; then
        awk -F"\t" 'NR==FNR{a[$1"_"$2"_"$3"_"$4]=$0;next}!($1"_"$2"_"$3"_"$4 in a){print}' \
            "$COMMON" "$SAMPLE.GFilter.txt" | \
            sort -k1,1V -k2,2n \
            > "$SAMPLE.Filter.txt"
    else
        cp "$SAMPLE.GFilter.txt" "$SAMPLE.Filter.txt"
    fi

    # Cleanup temp files
    rm -f "$SAMPLE.filterASXL1G646W.txt" "$SAMPLE.filterVF5.txt" "$SAMPLE.GFilter.txt"
}

# Export function for use with parallel
export -f filter_anno

# =============================================================================
# TMSP PROCESSING
# =============================================================================

process_tmsp() {
    log_info "######################## PROCESSING TMSP #########################"

    # Check for VCF files
    VCF_COUNT=$(ls *.vcf 2>/dev/null | grep -v "^Merge" | wc -l)
    if [ "$VCF_COUNT" -eq 0 ]; then
        log_error "No VCF files found in current directory"
        return 1
    fi
    log_info "Found $VCF_COUNT VCF files"

    # Rename files if >20 characters
    rename 's/^(.{20}).*.vcf/$1\.vcf/' *.vcf 2>/dev/null || true

    # Clean up VCF by using only QUAL>=100 and PASS/SB variants
    log_info "Filtering VCFs (QUAL>=100 and PASS/SB)..."
    parallel "bcftools view -i \"%QUAL>=100&&(%FILTER='PASS'||%FILTER='SB')\" -o {.}.x.vcf {}; mv {.}.x.vcf {}" ::: $(ls *.vcf | grep -v "^Merge")

    # Create VCFstats files using integrated function
    log_info "Creating VCFstats files..."
    for vcf in $(ls *.vcf | grep -v "^Merge"); do
        vcf_stats "$vcf"
    done

    # Check if sample in Database AND is not empty, else copy to Clinical
    log_info "Checking database for new samples..."
    parallel "[ ! -f $DATABASE_TMSP/{} -a -s {} ] && { echo \"{} absent in Database. Copying...\"; cp {.}.vcf {.}.VCFstats.txt $DATABASE_TMSP/; echo \"Done\"; } || true" ::: $(ls *.vcf | grep -v "^Merge")

    # Start TMSP annotation processing
    if [ ! -f *.xlsx ]; then
        log_info "Running mergeVCFannotation-optimized.sh (TMSP mode)..."
        "$SCRIPT_DIR/mergeVCFannotation-optimized.sh" tmsp
    else
        log_info "Excel files already exist, skipping annotation..."
    fi

    # Move all files to annotation directory
    if [ ! -d ./annotationTMSP ]; then
        log_info "Creating annotationTMSP directory..."
        mkdir ./annotationTMSP
        mv *.txt *.xlsx ./annotationTMSP/ 2>/dev/null || true
        mv transvar/ ./annotationTMSP/ 2>/dev/null || true
        mv *.html ./annotationTMSP/ 2>/dev/null || true
        mv Merge.* ./annotationTMSP/ 2>/dev/null || true
    fi

    log_info "TMSP processing complete"
}

# =============================================================================
# CEBPA/CEBNX PROCESSING
# =============================================================================

process_cebpa() {
    log_info "######################## PROCESSING CEBNX #########################"

    # Change to CEBPA directory
    if [ ! -d ../cebpa/vcf ]; then
        log_info "CEBPA directory not found (../cebpa/vcf), skipping..."
        return 0
    fi

    cd ../cebpa/vcf

    # Check for VCF files
    VCF_COUNT=$(ls *.vcf 2>/dev/null | grep -v "^Merge" | wc -l)
    if [ "$VCF_COUNT" -eq 0 ]; then
        log_info "No VCF files found in CEBPA directory, skipping..."
        cd ../../vcf
        return 0
    fi
    log_info "Found $VCF_COUNT VCF files in CEBPA"

    # Rename files if >20 characters
    rename 's/^(.{20}).*.vcf/$1\.vcf/' *.vcf 2>/dev/null || true

    # Clean up VCF by using only PASS/SB variants (NO QUAL filter for CEBPA)
    log_info "Filtering VCFs (PASS/SB only, no QUAL filter)..."
    parallel "bcftools view -i \"(%FILTER='PASS'||%FILTER='SB')\" -o {.}.x.vcf {}; mv {.}.x.vcf {}" ::: $(ls *.vcf | grep -v "^Merge")

    # Create VCFstats files using integrated function
    log_info "Creating VCFstats files..."
    for vcf in $(ls *.vcf | grep -v "^Merge"); do
        vcf_stats "$vcf"
    done

    # Check if sample in Database AND is NOT empty, else copy to Clinical
    log_info "Checking CEBPA database for new samples..."
    parallel "[ ! -f $DATABASE_CEBPA/{} -a -s {} ] && { echo \"{} absent in Database. Copying...\"; cp {.}.vcf {.}.VCFstats.txt $DATABASE_CEBPA/; echo \"Done\"; } || true" ::: $(ls *.vcf | grep -v "^Merge")

    # Start CEBPA annotation processing
    if [ ! -f *.xlsx ]; then
        log_info "Running mergeVCFannotation-optimized.sh (CEBPA mode)..."
        "$SCRIPT_DIR/mergeVCFannotation-optimized.sh" cebpa
    else
        log_info "Excel files already exist, skipping annotation..."
    fi

    # Move all files to annotation directory
    if [ ! -d ./annotationCEBNX ]; then
        log_info "Creating annotationCEBNX directory..."
        mkdir ./annotationCEBNX
        mv *.txt *.xlsx ./annotationCEBNX/ 2>/dev/null || true
        mv transvar/ ./annotationCEBNX/ 2>/dev/null || true
        mv *.html ./annotationCEBNX/ 2>/dev/null || true
        mv Merge.* ./annotationCEBNX/ 2>/dev/null || true
    fi

    # Return to TMSP directory
    cd ../../vcf

    log_info "CEBPA processing complete"
}

# =============================================================================
# CREATE OUTPUT DIRECTORY
# =============================================================================

create_output() {
    log_info "######################## CREATING OUTPUT #########################"

    # Create output directory if it doesn't exist
    if [ ! -d ../output ]; then
        log_info "Creating output directory..."
        mkdir ../output
    fi

    # Copy TMSP xlsx files and move annotationTMSP
    if [ -d ./annotationTMSP ]; then
        log_info "Copying TMSP files to output..."
        cp ./annotationTMSP/*.xlsx ../output/ 2>/dev/null || true
        mv ./annotationTMSP ../output/ 2>/dev/null || true
    fi

    # Copy CEBPA xlsx files and move annotationCEBNX
    if [ -d ../cebpa/vcf/annotationCEBNX ]; then
        log_info "Copying CEBPA files to output..."
        cp ../cebpa/vcf/annotationCEBNX/*.xlsx ../output/ 2>/dev/null || true
        mv ../cebpa/vcf/annotationCEBNX ../output/ 2>/dev/null || true
    fi

    log_info "Output directory created"
}

# =============================================================================
# GENERATE HTML REPORTS
# =============================================================================

generate_html_reports() {
    log_info "######################## GENERATING HTML REPORTS #########################"

    local output_dir="../output"
    local html_script="$SCRIPT_DIR/excel_to_html_report.py"

    # Check if HTML generation script exists
    if [ ! -f "$html_script" ]; then
        log_error "HTML generation script not found: $html_script"
        log_info "Skipping HTML report generation..."
        return 0
    fi

    # Check if python3 and openpyxl are available
    if ! command -v python3 &> /dev/null; then
        log_error "python3 not found, skipping HTML report generation"
        return 0
    fi

    if ! python3 -c "import openpyxl" 2>/dev/null; then
        log_error "openpyxl module not found. Install with: pip install openpyxl"
        log_info "Skipping HTML report generation..."
        return 0
    fi

    # Generate HTML reports for each TMSP xlsx file (excluding Summary.xlsx and CEBNX files)
    cd "$output_dir"

    # Create html_reports directory
    mkdir -p html_reports

    local xlsx_count=0
    for xlsx in $(ls *-TMSP*.xlsx 2>/dev/null | grep -v "^Summary" | grep -v "^~"); do
        log_info "Generating HTML report for: $xlsx"

        local sample_name=$(basename "$xlsx" .xlsx)

        # Generate single-page report (sample_name.html)
        python3 "$html_script" "$xlsx" --single-page -o "html_reports/${sample_name}.html"

        if [ $? -eq 0 ]; then
            xlsx_count=$((xlsx_count + 1))
            log_info "  -> Generated: html_reports/${sample_name}.html"
        else
            log_error "  -> Failed to generate HTML for $xlsx"
        fi
    done

    # Also generate for CEBNX files if present
    for xlsx in $(ls *-CEBNX*.xlsx 2>/dev/null | grep -v "^Summary" | grep -v "^~"); do
        log_info "Generating HTML report for: $xlsx"

        local sample_name=$(basename "$xlsx" .xlsx)

        # Generate single-page report (sample_name.html)
        python3 "$html_script" "$xlsx" --single-page -o "html_reports/${sample_name}.html"

        if [ $? -eq 0 ]; then
            xlsx_count=$((xlsx_count + 1))
            log_info "  -> Generated: html_reports/${sample_name}.html"
        else
            log_error "  -> Failed to generate HTML for $xlsx"
        fi
    done

    cd - > /dev/null

    if [ $xlsx_count -gt 0 ]; then
        log_info "Generated HTML reports for $xlsx_count sample(s)"
        log_info "HTML reports available at: $output_dir/html_reports/"
    else
        log_info "No Excel files found for HTML report generation"
    fi
}

# =============================================================================
# MAIN FUNCTION
# =============================================================================

main() {
    local start_time=$(date +%s)

    log_info "=========================================="
    log_info "VCF Processing Pipeline (TMSP + CEBPA)"
    log_info "=========================================="
    log_info "Working directory: $(pwd)"
    log_info "Script directory: $SCRIPT_DIR"

    # Export COMMON_VARIANTS for filter_anno function
    export COMMON_VARIANTS

    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --check)
                log_info "Checking dependencies..."
                echo ""
                echo "System tools:"
                for cmd in bcftools parallel rename perl awk java python3 transvar vep bgzip tabix vcf-merge vcf-sort; do
                    if command -v "$cmd" &> /dev/null; then
                        echo "  $cmd: OK"
                    else
                        echo "  $cmd: MISSING"
                    fi
                done
                echo ""
                echo "Local scripts (in $SCRIPT_DIR):"
                for script in writeTMSPtoXLS.pl write1WStoXLS.pl mergeVCFannotation-optimized.sh excel_to_html_report.py; do
                    if [ -f "$SCRIPT_DIR/$script" ]; then
                        echo "  $script: OK"
                    else
                        echo "  $script: MISSING"
                    fi
                done
                echo ""
                echo "Python modules:"
                if python3 -c "import openpyxl" 2>/dev/null; then
                    echo "  openpyxl: OK"
                else
                    echo "  openpyxl: MISSING (pip install openpyxl)"
                fi
                echo ""
                echo "Integrated functions:"
                echo "  vcf_stats: OK (built-in)"
                echo "  filter_anno: OK (built-in)"
                echo ""
                echo "Databases:"
                for db in "$HOME/Databases/TMSPvcf/TSMPclean" "$HOME/Databases/TMSPvcf/CEBPA" "$HOME/Databases/humandb" "$HOME/Databases/vep"; do
                    if [ -d "$db" ]; then
                        echo "  $db: OK"
                    else
                        echo "  $db: MISSING"
                    fi
                done
                exit 0
                ;;
            --help|-h)
                echo "Usage: $0 [options]"
                echo "  Run from the vcf directory containing TMSP VCF files"
                echo ""
                echo "Options:"
                echo "  --check       Check dependencies"
                echo "  --help        Show this help"
                echo ""
                echo "Expected directory structure:"
                echo "  /path/to/analysis/"
                echo "  ├── vcf/                  <- Run script from here"
                echo "  ├── cebpa/vcf/            <- CEBPA VCF files"
                echo "  ├── bam/                  <- BAM files"
                echo "  └── output/               <- Created by script"
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                echo "Use --help for usage information"
                exit 1
                ;;
        esac
    done

    # Process TMSP
    process_tmsp

    # Process CEBPA
    process_cebpa

    # Create output directory and copy files
    create_output

    # Generate HTML reports from Excel files
    generate_html_reports

    local end_time=$(date +%s)
    local duration=$((end_time - start_time))

    log_info "=========================================="
    log_info "######################## DONE #########################"
    log_info "Pipeline completed in ${duration}s"
    log_info "=========================================="
}

# Run main function
main "$@"
