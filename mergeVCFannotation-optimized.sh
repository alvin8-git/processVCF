#!/bin/bash
# =============================================================================
# mergeVCFannotation-optimized.sh - Optimized VCF Annotation Pipeline
# =============================================================================
# This script is an optimized version of mergeVCFannotation.awk that:
# - Runs independent annotations in parallel (ANNOVAR, VEP, snpEff, TransVar)
# - Runs database comparisons in parallel
# - Runs database lookups in parallel (SG10K, genomeAsia)
# - Provides better progress tracking and error handling
# - Supports both TMSP and CEBPA/CEBNX panels
#
# Usage: Run from directory containing VCF files
#   ./mergeVCFannotation-optimized.sh [tmsp|cebpa]
#
# Modes:
#   tmsp  - TMSP panel (default): Clinical, AD, LRM, NA12878 comparisons
#   cebpa - CEBPA/CEBNX panel: Clinical, AD, CAP, MOLM14 comparisons
#
# Author: Optimized from original AWK scripts
# =============================================================================

set -e  # Exit on error

# =============================================================================
# MODE SELECTION
# =============================================================================

MODE="${1:-tmsp}"
MODE=$(echo "$MODE" | tr '[:upper:]' '[:lower:]')

if [[ "$MODE" != "tmsp" && "$MODE" != "cebpa" ]]; then
    echo "Error: Invalid mode '$MODE'. Use 'tmsp' or 'cebpa'"
    echo "Usage: $0 [tmsp|cebpa]"
    exit 1
fi

# =============================================================================
# CONFIGURATION
# =============================================================================

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Paths
HOME_DIR="${HOME:-/home/alvin}"
SCRIPTS_DIR="$SCRIPT_DIR"
SOFTWARE_DIR="$HOME_DIR/Software"
DATABASES_DIR="$HOME_DIR/Databases"

# Annotation tools
ANNOVAR_DIR="$SOFTWARE_DIR/annovar"
SNPEFF_DIR="$SOFTWARE_DIR/snpEff"
CANCERVAR_DIR="$SOFTWARE_DIR/CancerVar"

# Databases
HUMANDB="$DATABASES_DIR/humandb"
VEP_CACHE="$DATABASES_DIR/vep"
HG19_FASTA="$DATABASES_DIR/WholeGenomeFASTA/genome.fa"
HG19_SORTED_FASTA="$DATABASES_DIR/WholeGenomeFASTA/hg19.sort.fa"
REFSEQ_TRANSCRIPTS="$DATABASES_DIR/Ensembldata/RefSeqSelectTranscript.txt"
SG10K_DB="$DATABASES_DIR/SG10K.hg37.vcf/SG10K.TMSPgenes.txt"
GENOMEASIA_DB="$DATABASES_DIR/genomeAsia/genomeAsia.TMSP54.txt"
COMMON_VARIANTS="$DATABASES_DIR/TMSPvcf/TMSPcommon/TMSPcommonlist.txt"

# Mode-specific configuration
if [[ "$MODE" == "tmsp" ]]; then
    VCF_DATABASE="$HOME_DIR/Databases/TMSPvcf/TSMPclean"
    COMPARE_DBS=("Clinical" "AD" "LRM" "NA12878")
    COMPARE_NAMES=("Clin" "AD" "LRM" "NA")
    EXCEL_WRITER="writeTMSPtoXLS.pl"
    DO_IKZF1_WT1_FIX=true
    DO_SUMMARY=true
else
    # CEBPA mode
    VCF_DATABASE="$HOME_DIR/Databases/TMSPvcf/CEBPA"
    COMPARE_DBS=("Clinical" "AD" "CAP" "MOLM14")
    COMPARE_NAMES=("Clin" "AD" "CAP" "MOLM")
    EXCEL_WRITER="writeTMSPtoXLS.pl"
    DO_IKZF1_WT1_FIX=false
    DO_SUMMARY=false
fi

# Processing settings
MAX_PARALLEL_JOBS=4
MEM_SUSPEND="3G"

# =============================================================================
# INTEGRATED FUNCTIONS (merged from VCFstats.awk and filterAnno.awk)
# =============================================================================

# VCFstats function - Extracts useful statistics from a VCF file
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

# Export for parallel
export -f vcf_stats

# filterAnno function - Filters annotation file based on various criteria
# Usage: filter_anno <annotation_file>
# Output: <sample>.Filter.txt
filter_anno() {
    local FILE="$1"
    local SAMPLE=$(basename "$FILE" .annotation.txt)
    local COMMON="${COMMON_VARIANTS:-$HOME/Databases/TMSPvcf/TMSPcommon/TMSPcommonlist.txt}"

    [ -f "$SAMPLE.Filter.txt" ] && return 0

    # EXCEPTION CASES
    awk -F"\t" '{if($14>0.05&&$1=="chr20"&&($2==31022441)) print}' "$FILE" \
        > "$SAMPLE.filterASXL1G646W.txt"

    awk -F"\t" '{if(($14<0.05&&$6=="TP53")||($14<0.05&&$6=="JAK2")||($14<0.05&&$6=="FLT3")||($14<0.05&&$6=="CALR")) print}' "$FILE" \
        > "$SAMPLE.filterVF5.txt"

    # GENERAL FILTERING
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

    # REMOVE COMMON VARIANTS
    if [ -f "$COMMON" ]; then
        awk -F"\t" 'NR==FNR{a[$1"_"$2"_"$3"_"$4]=$0;next}!($1"_"$2"_"$3"_"$4 in a){print}' \
            "$COMMON" "$SAMPLE.GFilter.txt" | \
            sort -k1,1V -k2,2n \
            > "$SAMPLE.Filter.txt"
    else
        cp "$SAMPLE.GFilter.txt" "$SAMPLE.Filter.txt"
    fi

    rm -f "$SAMPLE.filterASXL1G646W.txt" "$SAMPLE.filterVF5.txt" "$SAMPLE.GFilter.txt"
}

# Export for parallel
export -f filter_anno

# =============================================================================
# LOGGING
# =============================================================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${BLUE}[$(date '+%H:%M:%S')]${NC} INFO: $1"
}

log_success() {
    echo -e "${GREEN}[$(date '+%H:%M:%S')]${NC} DONE: $1"
}

log_warn() {
    echo -e "${YELLOW}[$(date '+%H:%M:%S')]${NC} WARN: $1"
}

log_error() {
    echo -e "${RED}[$(date '+%H:%M:%S')]${NC} ERROR: $1" >&2
}

log_step() {
    echo -e "\n${GREEN}======================== $1 ========================${NC}"
}

# =============================================================================
# ANNOTATION FUNCTIONS
# =============================================================================

run_annovar() {
    local vcf="$1"
    local sample=$(basename "$vcf" .vcf)

    [ -f "$sample.annovar.txt" ] && return 0

    log_info "Starting ANNOVAR for $sample"

    # Convert to ANNOVAR input
    perl "$ANNOVAR_DIR/convert2annovar.pl" \
        --format vcf4old "$vcf" \
        --outfile "$sample.avinput" \
        --includeinfo --withzyg 2>/dev/null

    # Run table_annovar
    perl "$ANNOVAR_DIR/table_annovar.pl" "$sample.avinput" "$HUMANDB" \
        --buildver hg19 \
        --outfile "$sample.outfile" \
        --protocol refGene,ensGene,cytoBand,genomicSuperDups,esp6500siv2_all,exac03,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_amr,1000g2015aug_eur,1000g2015aug_sas,intervar_20180118,kaviar_20150923,mcap,revel,gnomad211_exome,gnomad_genome,snp138,ljb26_all,clinvar_20220320,cosmic,dbscsnv11,hrcr1,gme,icgc28,TCGA,civic,tumorportal,gwasCatalog,wgEncodeRegTfbsClustered,wgEncodeRegDnaseClusteredV3,cg69,nci60,hapmapSnpsCEU,hapmapSnpsCHB,hapmapSnpsJPT,hapmapSnpsYRI,wgRna,dgvMerged,phastConsElements46way,phastConsElements100way,tfbsConsSites,targetScanS \
        --operation g,g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,f,f,f,f,f,f,f,r,r,r,f,f,r,r,r,r,r,r,r,r,r,r \
        --argument '-splicing 100',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, \
        --nastring . --otherinfo 2>/dev/null

    mv "$sample.outfile.hg19_multianno.txt" "$sample.annovar.txt"
    rm -f "$sample.outfile."* "$sample.avinput"

    log_success "ANNOVAR completed for $sample"
}

run_vep() {
    local vcf="$1"
    local sample=$(basename "$vcf" .vcf)

    [ -f "$sample.vep.txt" ] && return 0

    if ! command -v vep &> /dev/null; then
        log_warn "VEP not found, skipping"
        return 0
    fi

    log_info "Starting VEP for $sample"

    vep -i "$vcf" \
        --dir_cache "$VEP_CACHE" \
        --cache \
        --assembly GRCh37 \
        --offline \
        --fasta "$HG19_SORTED_FASTA" \
        --everything \
        --force_overwrite \
        --tab \
        --fork 8 \
        --refseq \
        --pick \
        --exclude_predicted \
        --no_escape \
        --use_given_ref \
        -o "$sample.vep.output.txt" 2>/dev/null

    rm -f "${sample}.vep.output.txt_summary.html"
    grep -v "##" "$sample.vep.output.txt" > "$sample.vep.txt"
    rm -f "$sample.vep.output.txt"

    log_success "VEP completed for $sample"
}

run_snpeff() {
    local vcf="$1"
    local sample=$(basename "$vcf" .vcf)

    [ -f "$sample.snpEff.txt" ] && return 0

    log_info "Starting snpEff for $sample"

    local snpEff="$SNPEFF_DIR/snpEff.jar"
    local SnpSift="$SNPEFF_DIR/SnpSift.jar"
    local scriptsDir="$SNPEFF_DIR/scripts"

    # Run snpEff
    java -Xmx2g -jar "$snpEff" GRCh37.p13.RefSeq \
        -canon -onlyProtein -strict -noShiftHgvs "$vcf" \
        > "$sample.snpEff.ref.vcf" 2>/dev/null

    # Split into one effect per line
    cat "$sample.snpEff.ref.vcf" | "$scriptsDir/vcfEffOnePerLine.pl" \
        > "$sample.snpEff.1Eff.vcf"

    # Filter out upstream/downstream/intragenic
    cat "$sample.snpEff.1Eff.vcf" | \
        java -jar "$SnpSift" filter -n \
        "(ANN[*].EFFECT has 'upstream_gene_variant') | \
         (ANN[*].EFFECT has 'downstream_gene_variant') | \
         (ANN[*].EFFECT has 'intragenic_variant')" 2>/dev/null \
        > "$sample.snpEff.filter.vcf"

    # Combine and sort
    cat "$sample.snpEff.filter.vcf" \
        <(awk 'FNR==NR{a[$1"_"$2"_"$4"_"$5]=$0;next}{if(!($1"_"$2"_"$4"_"$5 in a))print}' \
        <(grep -v "^\#" "$sample.snpEff.filter.vcf") \
        <(grep -v "^\#" "$sample.snpEff.ref.vcf")) | vcf-sort -c 2>/dev/null \
        > "$sample.snpEff.vcf"

    # Extract fields
    java -jar "$SnpSift" extractFields \
        -s "," -e "." "$sample.snpEff.vcf" \
        CHROM POS REF ALT FILTER AF AC DP MQ \
        "ANN[*].GENE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" \
        "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].BIOTYPE" "ANN[*].RANK" \
        "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" "ANN[*].CDS_POS" "ANN[*].CDS_POS" \
        "ANN[*].AA_POS" "ANN[*].AA_LEN" \
        "LOF[*].GENE" "LOF[*].NUMTR" "LOF[*].PERC" 2>/dev/null \
        | sed 's/ANN\[\*\]\.//g' \
        > "$sample.presnpEff.txt"

    # Match to VCF order
    cat \
        <(awk -F"\t" '{if(NR==1)print}' "$sample.presnpEff.txt") \
        <(awk -F"\t" 'NR==FNR{a[$1"_"$2"_"$3"_"$4]=$0;next}($1"_"$2"_"$4"_"$5 in a){print a[$1"_"$2"_"$4"_"$5]}' \
        "$sample.presnpEff.txt" <(grep -v "^#" "$vcf")) \
        > "$sample.snpEff.txt"

    rm -f "$sample.snpEff.ref.vcf" "$sample.snpEff.1Eff.vcf" \
          "$sample.snpEff.filter.vcf" "$sample.snpEff.vcf" "$sample.presnpEff.txt"

    log_success "snpEff completed for $sample"
}

run_transvar() {
    local vcf="$1"
    local sample=$(basename "$vcf" .vcf)

    [ -f "$sample.transvar.txt" ] && return 0

    if ! command -v transvar &> /dev/null; then
        log_warn "TransVar not found, skipping"
        return 0
    fi

    log_info "Starting TransVar for $sample"

    mkdir -p transvar

    # Run transvar
    transvar ganno --reference "$HG19_FASTA" --refversion hg19 \
        --vcf "$vcf" --refseq --aa3 --noheader 2>/dev/null | \
        awk '{if($1!~/^\#/)print}' | awk '{if($11!~/^X/)print}' \
        > "./transvar/$sample.transvar.vcf"

    # Extract useful fields
    paste \
        <(awk '{if($1!~/^\##/)print}' "./transvar/$sample.transvar.vcf" | cut -f1-2,4-5) \
        <(awk '{if($1!~/^\##/)print}' "./transvar/$sample.transvar.vcf" | rev | cut -f3-6 | rev) \
        | sed 's/\//\t/g' | awk '{if($5!~/^X/)print}' | sed 's/[[:space:]](protein_coding)//g' \
        > "$sample.noX.txt"

    # Remove version numbers
    awk '{print $5}' "$sample.noX.txt" | cut -d"." -f1 > "$sample.noX.noV.txt"

    # Replace transcripts
    paste \
        <(cut -f1-4 "$sample.noX.txt") \
        <(cut -f1 "$sample.noX.noV.txt") \
        <(cut -f6- "$sample.noX.txt") \
        > "$sample.noX2.txt"

    # Get transcript sizes
    if [ -f "$REFSEQ_TRANSCRIPTS" ]; then
        awk -F"\t" 'NR==FNR{a[$1]=$7;next}{if(($5 in a))print $0"\t"a[$5];else print $0"\t-"}' \
            "$REFSEQ_TRANSCRIPTS" "$sample.noX2.txt" \
            > "$sample.noX.size.txt"
    else
        awk '{print $0"\t-"}' "$sample.noX2.txt" > "$sample.noX.size.txt"
    fi

    # Get VCF positions
    awk '{if($1!~/^\#/)print}' "$vcf" | cut -f1-2,4-5 > "$sample.noheaderVCF.txt"

    # Find longest transcript per variant
    awk '{
        if ( arr[$1"_"$2"_"$3"_"$4] == "" ) {
            arr[$1"_"$2"_"$3"_"$4] = $11
        }
        if ( arr[$1"_"$2"_"$3"_"$4] != "" ) {
            if ( arr[$1"_"$2"_"$3"_"$4] <= $11 ) { arr[$1"_"$2"_"$3"_"$4] = $11 }
        }
    }
    END { for (x in arr) print x"\t"arr[x] }' "$sample.noX.size.txt" \
        | sed 's/\_/\t/g' \
        > "$sample.longTrans.unsort.txt"

    # Order by VCF order
    awk -F"\t" 'NR==FNR{a[$1"\t"$2"\t"$3"\t"$4]=$0;next}{if(($1"\t"$2"\t"$3"\t"$4 in a))print a[$1"\t"$2"\t"$3"\t"$4];else print $1"\t"$2"\t"$3"\t"$4"\t-"}' \
        "$sample.longTrans.unsort.txt" "$sample.noheaderVCF.txt" \
        > "$sample.longTrans.txt"

    # Final output
    echo -e "Chr\tPos\tRef\tAlt\tTranscript\tGene\tStrand\tHGVSg\tHGVSc\tHGVSp\tLength" \
        > "$sample.transvar.txt"
    awk -F"\t" 'NR==FNR{a[$1"\t"$2"\t"$3"\t"$4"\t"$11]=$0;next}{if(($0 in a))print a[$0];else print $1"\t"$2"\t"$3"\t"$4"\t-\t-\t-\t-\t-\t-\t-"}' \
        "$sample.noX.size.txt" "$sample.longTrans.txt" \
        | sed 's/chr[[:alnum:]]*\://g' \
        | sed "s/\(p.*\)X/\1Ter/g" \
        >> "$sample.transvar.txt"

    rm -f "$sample.noX.size.txt" "$sample.noX.txt" "$sample.noX2.txt" \
          "$sample.noX.noV.txt" "$sample.noheaderVCF.txt" \
          "$sample.longTrans.txt" "$sample.longTrans.unsort.txt"

    log_success "TransVar completed for $sample"
}

run_cancervar() {
    local vcf="$1"
    local sample=$(basename "$vcf" .vcf)

    [ -f "$sample.CancerVar.txt" ] && return 0

    if [ ! -f "$CANCERVAR_DIR/CancerVar.py" ]; then
        log_warn "CancerVar not found, skipping"
        return 0
    fi

    log_info "Starting CancerVar for $sample"

    python3 "$CANCERVAR_DIR/CancerVar.py" \
        -b hg19 \
        -i "$vcf" \
        --input_type=VCF \
        -o "$sample.CancerVar" 2>/dev/null

    mv "$sample.CancerVar.hg19_multianno.txt.cancervar" "$sample.CancerVar.txt"
    rm -f "$sample.CancerVar.hg19_multianno.txt"* "$sample.vcf.avinput"

    log_success "CancerVar completed for $sample"
}

run_sg10k() {
    local vcf="$1"
    local sample=$(basename "$vcf" .vcf)

    [ -f "$sample.SG10k.txt" ] && return 0
    [ ! -f "$SG10K_DB" ] && return 0

    log_info "Starting SG10K annotation for $sample"

    awk '{if($1!~/^\#/)print $1"\t"$2"\t"$4"\t"$5}' "$vcf" > "$sample.noheader.vcf"

    awk -F"\t" 'NR==FNR{a["chr"$1"_"$2"_"$3"_"$4]="chr"$0;next}{if(($1"_"$2"_"$3"_"$4 in a))print a[$1"_"$2"_"$3"_"$4];else print $1"\t"$2"\t"$3"\t"$4"\t-\t-\t-\t-\t-\t-\t-\t-"}' \
        "$SG10K_DB" "$sample.noheader.vcf" \
        > "$sample.SG10k.txt"

    sed -i '1s/^/CHR\tBP.B37\tREF\tALT\tAN_All\tAF_All\tAN_CHS\tAF_CHS\tAN_INS\tAF_INS\tAN_MAS\tAF_MAS\n/' "$sample.SG10k.txt"

    rm -f "$sample.noheader.vcf"

    log_success "SG10K completed for $sample"
}

run_genomeasia() {
    local vcf="$1"
    local sample=$(basename "$vcf" .vcf)

    [ -f "$sample.genomeAsia.txt" ] && return 0
    [ ! -f "$GENOMEASIA_DB" ] && return 0

    log_info "Starting GenomeAsia annotation for $sample"

    awk -F"\t" 'NR==FNR{a["chr"$1"_"$2"_"$3"_"$4]="chr"$0;next}{if($1"_"$2"_"$3"_"$4 in a)print a[$1"_"$2"_"$3"_"$4];else print $1"\t"$2"\t"$3"\t"$4"\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-"}' \
        "$GENOMEASIA_DB" \
        <(grep -v "^\#" "$vcf" | cut -f1,2,4,5) | \
        sort -k1,1V -k2,2n | \
        sed '1s/^/CHR\tPOS\tREF\tALT\tSEA_AN\tSEA_AC\tSEA_AF\tSEA_HOM\tNEA_AN\tNEA_AC\tNEA_AF\tNEA_HOM\tSAS_AN\tSAS_AC\tSAS_AF\tSAS_HOM\n/' \
        > "$sample.genomeAsia.txt"

    log_success "GenomeAsia completed for $sample"
}

run_annocheck() {
    local vcf="$1"
    local sample=$(basename "$vcf" .vcf)

    [ -f "$sample.annoCheck.txt" ] && return 0

    # Check required files exist
    for f in "$sample.transvar.txt" "$sample.annovar.txt" "$sample.vep.txt" "$sample.snpEff.txt"; do
        [ ! -f "$f" ] && return 0
    done

    log_info "Running annotation concordance check for $sample"

    paste \
        <(cut -f1-4,5,6,9-10 "$sample.transvar.txt" | sed '1c Chr\tPos\tRef\tAlt\tTrans.Tvar\tGene.Tvar\tHGVSc.Tvar\tHGVSp.Tvar') \
        <(cut -f6,7,9 "$sample.annovar.txt" | sed 's/\./\-/g' | sed '1c Loc.Ann\tGene.Ann\tType.Ann') \
        <(cut -f19,5,43-44,7 "$sample.vep.txt" | awk '{gsub("-",":-",$5)}1' | awk '{gsub(/\..*$/,"",$1)}1' | sed 's/ /\t/g' | sed 's/:/\t/g' | cut -f1-3,5,7 | sed '1c Trans.Vep\tType.Vep\tGene.Vep\tHGVSc.Vep\tHGVSp.Vep') \
        <(cut -f10-13,15 "$sample.snpEff.txt" | awk '{gsub(/\..*$/,"",$2)}1' | sed 's/ /\t/g' | sed '1c Gene.SpEf\tTrans.SpEf\tHGVSc.SpEf\tHGVSp.SpEf\tType.SpEf') \
        > "$sample.PreCheck.txt"

    paste \
        <(cut -f1-4 "$sample.PreCheck.txt") \
        <(awk -F"\t" '{print $6"\t"$10"\t"$14"\t"$17"\t"$5"\t"$12"\t"$18"\t"$7"\t"$15"\t"$19"\t"$8"\t"$16"\t"$20"\t"$9"\t"$11"\t"$13"\t"$21}' "$sample.PreCheck.txt") \
        > "$sample.annoCheck.txt"

    rm -f "$sample.PreCheck.txt"

    log_success "AnnoCheck completed for $sample"
}

# =============================================================================
# COMPARISON FUNCTION
# =============================================================================

run_compare() {
    local sample_vcf="$1"
    local vcf_dir="$2"
    local name="$3"
    local sample=$(basename "$sample_vcf" .vcf)

    [ -f "$sample.compare$name.txt" ] && return 0
    [ ! -d "$vcf_dir" ] && return 0

    log_info "Comparing $sample with $name samples"

    mkdir -p "$sample.compare$name"

    # Compare with each VCF in directory
    for file in "$vcf_dir"/*.vcf; do
        [ ! -f "$file" ] && continue
        local filename=$(basename "$file" .vcf)
        awk 'FNR==NR{a[$1"_"$2"_"$4"_"$5]=$5;next}{if(($1"_"$2"_"$3"_"$4 in a))print a[$1"_"$2"_"$3"_"$4];else print"-"}' \
            "$file" \
            <(grep -v "^#" "$sample_vcf" | cut -f1-2,4-5) \
            > "./$sample.compare$name/$filename.out"
    done

    # Concatenate results
    paste \
        <(grep -v "^#" "$sample_vcf" | cut -f1-2,4-5 | sed '1i #Chr\tPos\tRef\tAlt') \
        <(cat <(echo "$vcf_dir"/*.VCFstats.txt | sed 's/.VCFstats.txt//g' | sed "s|$vcf_dir\/||g" | sed -r "s/\s+/\t/g") \
        <(paste "./$sample.compare$name"/*.out 2>/dev/null)) \
        > "$sample.precompare$name.txt" 2>/dev/null || true

    # Get field count
    local field=$(cut -f5- "$sample.precompare$name.txt" 2>/dev/null | awk '{print NF}' | sort | uniq | head -1)

    # Add count column
    paste \
        <(cut -f1-4 "$sample.precompare$name.txt") \
        <(grep -v "^#" "$sample.precompare$name.txt" | cut -f4- | awk '{print gsub($1,"")-1}' | sed "1i ${name}${field}") \
        <(cut -f5- "$sample.precompare$name.txt") \
        > "$sample.compare$name.txt" 2>/dev/null || true

    # Compare VAF
    for file2 in "$vcf_dir"/*.VCFstats.txt; do
        [ ! -f "$file2" ] && continue
        local filename2=$(basename "$file2" .VCFstats.txt)
        awk -F"\t" 'FNR==NR{a[$1"_"$2"_"$3"_"$4]=$9;next}{if(($1"_"$2"_"$3"_"$4 in a))print a[$1"_"$2"_"$3"_"$4]; else print "-"}' \
            "$file2" \
            <(grep -v "^#" "$sample_vcf" | cut -f1-2,4-5) \
            > "./$sample.compare$name/$filename2.out2"
    done

    paste \
        <(grep -v "^#" "$sample_vcf" | cut -f1-2,4-5 | sed '1i #Chr\tPos\tRef\tAlt') \
        <(cat <(echo "$vcf_dir"/*.VCFstats.txt | sed 's/.VCFstats.txt//g' | sed "s|$vcf_dir\/||g" | sed -r "s/\s+/\t/g") \
        <(paste "./$sample.compare$name"/*.out2 2>/dev/null)) \
        > "$sample.compareVF$name.txt" 2>/dev/null || true

    # Cleanup
    rm -f "$sample.precompare$name.txt"
    rm -rf "$sample.compare$name"

    log_success "Compare $name completed for $sample"
}

# =============================================================================
# MAIN PIPELINE
# =============================================================================

merge_vcfs() {
    local sample="Merge"

    [ -f "$sample.vcf" ] && return 0

    log_step "MERGING VCFs"

    ls *.vcf | grep -v "^Merge" > VCFlist.txt

    # Sort and compress each VCF in parallel
    parallel --memsuspend "$MEM_SUSPEND" "[ -f {.}.sort.gz ] || \
        (bgzip -c <(vcf-sort -c <(bcftools view -a --min-ac=1 --no-update {})) > {.}.sort.gz ; \
        tabix -p vcf {.}.sort.gz)" :::: VCFlist.txt

    # Merge
    vcf-merge -c none -s $(ls *.sort.gz | tr '\n' ' ') | vcf-sort -c > "$sample.vcf"

    # Create combined VCF for CancerVar
    parallel --jobs 1 "egrep -v \"^#\" {} >> $sample.All.vcf" :::: VCFlist.txt
    awk -F"\t" 'NR==FNR{a[$1"_"$2"_"$4"_"$5]=$0;next}($1"_"$2"_"$4"_"$5 in a){print a[$1"_"$2"_"$4"_"$5]}' \
        "$sample.All.vcf" "$sample.vcf" | sort -k1,1V -k2,2n > "$sample.Allsort.vcf"

    rm -f *.sort.gz *.sort.gz.tbi "$sample.All.vcf"

    log_success "VCF merging completed"
}

run_annotations_parallel() {
    local sample="Merge"

    [ -f "$sample.AnnoAll.txt" ] && return 0

    log_step "RUNNING ANNOTATIONS IN PARALLEL"

    # Export functions for parallel execution
    export -f run_annovar run_vep run_snpeff run_transvar log_info log_success log_warn
    export ANNOVAR_DIR SNPEFF_DIR VEP_CACHE HG19_FASTA HG19_SORTED_FASTA REFSEQ_TRANSCRIPTS

    # Run ANNOVAR, VEP, snpEff, TransVar in parallel
    # These are independent and can run simultaneously
    log_info "Starting parallel annotation (ANNOVAR, VEP, snpEff, TransVar)..."

    (
        run_annovar "$sample.vcf" &
        run_vep "$sample.vcf" &
        run_snpeff "$sample.vcf" &
        run_transvar "$sample.vcf" &
        wait
    )

    log_success "Core annotations completed"

    # CancerVar needs to run after (it runs its own ANNOVAR internally)
    log_info "Running CancerVar..."
    run_cancervar "$sample.Allsort.vcf"
    [ -f "$sample.Allsort.CancerVar.txt" ] && mv "$sample.Allsort.CancerVar.txt" "$sample.CancerVar.txt"

    # Run SG10K and genomeAsia in parallel (fast lookups)
    log_info "Running database lookups in parallel (SG10K, GenomeAsia)..."
    (
        run_sg10k "$sample.vcf" &
        run_genomeasia "$sample.vcf" &
        wait
    )

    # Run annotation concordance check
    run_annocheck "$sample.vcf"

    log_success "All annotations completed"
}

run_comparisons_parallel() {
    local sample="Merge"

    [ -f "$sample.compare.txt" ] && return 0

    log_step "RUNNING COMPARISONS IN PARALLEL"
    log_info "Mode: $MODE - Comparing against: ${COMPARE_DBS[*]}"

    # Run all comparisons in parallel using mode-specific databases
    (
        for i in "${!COMPARE_DBS[@]}"; do
            run_compare "$sample.vcf" "$VCF_DATABASE/${COMPARE_DBS[$i]}" "${COMPARE_NAMES[$i]}" &
        done
        wait
    )

    # Combine compare files based on mode
    local name0="${COMPARE_NAMES[0]}"
    local name1="${COMPARE_NAMES[1]}"
    local name2="${COMPARE_NAMES[2]}"
    local name3="${COMPARE_NAMES[3]}"

    if [ -f "$sample.compare${name0}.txt" ]; then
        paste \
            <(cut -f1-5 "$sample.compare${name0}.txt") \
            <(cut -f5 "$sample.compare${name1}.txt" 2>/dev/null || echo "") \
            <(cut -f5 "$sample.compare${name3}.txt" 2>/dev/null || echo "") \
            <(cut -f5 "$sample.compare${name2}.txt" 2>/dev/null || echo "") \
            <(cut -f6- "$sample.compare${name0}.txt") \
            <(cut -f6- "$sample.compare${name1}.txt" 2>/dev/null || echo "") \
            <(cut -f6- "$sample.compare${name2}.txt" 2>/dev/null || echo "") \
            <(cut -f6- "$sample.compare${name3}.txt" 2>/dev/null || echo "") \
            > "$sample.compare.txt" 2>/dev/null || true
        rm -f "$sample.compare${name0}.txt" "$sample.compare${name1}.txt" "$sample.compare${name2}.txt" "$sample.compare${name3}.txt"
    fi

    # Combine VAF compare files
    if [ -f "$sample.compareVF${name0}.txt" ]; then
        paste \
            <(cut -f1- "$sample.compareVF${name0}.txt") \
            <(cut -f5- "$sample.compareVF${name1}.txt" 2>/dev/null || echo "") \
            <(cut -f5- "$sample.compareVF${name2}.txt" 2>/dev/null || echo "") \
            <(cut -f5- "$sample.compareVF${name3}.txt" 2>/dev/null || echo "") \
            > "$sample.compareVAF.txt" 2>/dev/null || true
        rm -f "$sample.compareVF${name0}.txt" "$sample.compareVF${name1}.txt" "$sample.compareVF${name2}.txt" "$sample.compareVF${name3}.txt"
    fi

    log_success "Comparisons completed"
}

combine_annotations() {
    local sample="Merge"

    [ -f "$sample.AnnoAll.txt" ] && return 0

    log_step "COMBINING ANNOTATIONS"

    # Combine annotations
    if [ -f "$sample.annovar.txt" ] && [ -f "$sample.snpEff.txt" ] && [ -f "$sample.transvar.txt" ]; then
        paste \
            <(grep -v "^#" "$sample.vcf" | cut -f1-2,4-5 | sed 's/:\S*//g' | sed '1i Chr\tPos\tRef\tAlt') \
            <(cut -f10 "$sample.snpEff.txt") \
            <(cut -f5,8-10 "$sample.transvar.txt") \
            <(cut -f92,123 "$sample.annovar.txt" 2>/dev/null || echo "") \
            <(cut -f14 "$sample.CancerVar.txt" 2>/dev/null | sed 's/ CancerVar: //g' || echo "") \
            <(cut -f5-8 "$sample.compare.txt" 2>/dev/null || echo "") \
            <(cut -f5 "$sample.snpEff.txt") \
            <(cut -f11,9,16 "$sample.annovar.txt" 2>/dev/null || echo "") \
            <(cut -f7,16,18,39-40 "$sample.vep.txt" 2>/dev/null || echo "") \
            <(cut -f6,8,10,12 "$sample.SG10k.txt" 2>/dev/null || echo "") \
            <(cut -f7,11,15 "$sample.genomeAsia.txt" 2>/dev/null || echo "") \
            <(cut -f18-32,62-64,67,72,74,78,75,77,76,79,73,84-91 "$sample.annovar.txt" 2>/dev/null || echo "") \
            <(cut -f126-139,148-149 "$sample.annovar.txt" 2>/dev/null || echo "") \
            <(cut -f33-61,118-122 "$sample.annovar.txt" 2>/dev/null || echo "") \
            <(cut -f65-66,93-117,156-159 "$sample.annovar.txt" 2>/dev/null || echo "") \
            > "$sample.Anno1.txt" 2>/dev/null || true

        # Fix IKZF1 and WT1 annotations (TMSP mode only)
        if [ -f "$sample.Anno1.txt" ]; then
            if [[ "$DO_IKZF1_WT1_FIX" == true ]]; then
                log_info "Applying IKZF1/WT1 annotation fix (TMSP mode)"
                paste \
                    <(awk '{if($5=="IKZF1" || $5=="WT1")print}' "$sample.Anno1.txt" | cut -f1-5) \
                    <(awk '{if($10=="IKZF1" || $10=="WT1")print}' "$sample.snpEff.txt" | cut -f11) \
                    <(awk '{if($5=="IKZF1" || $5=="WT1")print}' "$sample.Anno1.txt" | cut -f7) \
                    <(awk '{if($10=="IKZF1" || $10=="WT1")print}' "$sample.snpEff.txt" | cut -f12-13) \
                    <(awk '{if($5=="IKZF1" || $5=="WT1")print}' "$sample.Anno1.txt" | cut -f10-) \
                    > "$sample.ReplaceAnno.txt" 2>/dev/null || true

                cat \
                    <(awk '{if(NR==1)print}' "$sample.Anno1.txt") \
                    <(cat <(awk '{if(NR!=1)print}' "$sample.Anno1.txt" | awk '{if($5!="IKZF1"&&$5!="WT1")print}') \
                    <(cat "$sample.ReplaceAnno.txt" 2>/dev/null) | sort -k1,1V -k2,2n) \
                    > "$sample.AnnoAll.txt"

                rm -f "$sample.Anno1.txt" "$sample.ReplaceAnno.txt"
            else
                # CEBPA mode - no IKZF1/WT1 fix needed
                mv "$sample.Anno1.txt" "$sample.AnnoAll.txt"
            fi
        fi
    fi

    log_success "Annotations combined"
}

extract_per_sample() {
    local sample="Merge"

    log_step "EXTRACTING PER-SAMPLE ANNOTATIONS"

    # Generate VCFstats for each sample using integrated function
    export -f vcf_stats
    parallel --memsuspend "$MEM_SUSPEND" "vcf_stats {}" :::: VCFlist.txt

    # Extract per-sample annotations
    if [ -f "$sample.AnnoAll.txt" ]; then
        parallel --memsuspend "$MEM_SUSPEND" "[ -f {.}.annotation.txt ] || \
            (awk '{if(NR==1)print}' $sample.AnnoAll.txt > {.}.Anno.txt; \
            awk -F\"\t\" 'NR==FNR{a[\$1\"_\"\$2\"_\"\$3\"_\"\$4]=\$0;next}(\$1\"_\"\$2\"_\"\$4\"_\"\$5 in a){print a[\$1\"_\"\$2\"_\"\$4\"_\"\$5]}' \
            $sample.AnnoAll.txt {} >> {.}.Anno.txt; \
            paste \
            <(cut -f1-4 {.}.Anno.txt) \
            <(cut -f5 {.}.VCFstats.txt) \
            <(cut -f5-9 {.}.Anno.txt) \
            <(cut -f6-9 {.}.VCFstats.txt) \
            <(cut -f10- {.}.Anno.txt | sed 's/Name=//g') \
            >{.}.Anno2.txt; \
            paste \
            <(cut -f1-15 {.}.Anno2.txt) \
            <(awk -F\"\t\" '{print \$6\":\"\$1\"(GRCh37):\"\$8\"; \"\$7\"(\"\$6\"):\"\$9\";\"\$10\"(\"int(\$14*100+0.5)\"% VAF)\"}' {.}.Anno2.txt) \
            <(cut -f16- {.}.Anno2.txt) \
            >{.}.annotation.txt; \
            rm {.}.Anno.txt {.}.Anno2.txt)" :::: VCFlist.txt
    fi

    # Extract compare files
    if [ -f "$sample.compare.txt" ]; then
        parallel --memsuspend "$MEM_SUSPEND" "[ -f {.}.compare.txt ] || \
            (awk '{if(NR==1)print}' $sample.compare.txt > {.}.compare.txt; \
            awk -F\"\t\" 'NR==FNR{a[\$1\"_\"\$2\"_\"\$3\"_\"\$4]=\$0;next}(\$1\"_\"\$2\"_\"\$4\"_\"\$5 in a){print a[\$1\"_\"\$2\"_\"\$4\"_\"\$5]}' \
            $sample.compare.txt {} >> {.}.compare.txt)" :::: VCFlist.txt
    fi

    # Extract compareVAF files
    if [ -f "$sample.compareVAF.txt" ]; then
        parallel --memsuspend "$MEM_SUSPEND" "[ -f {.}.compareVAF.txt ] || \
            (awk '{if(NR==1)print}' $sample.compareVAF.txt > {.}.compareVAF.txt; \
            awk -F\"\t\" 'NR==FNR{a[\$1\"_\"\$2\"_\"\$3\"_\"\$4]=\$0;next}(\$1\"_\"\$2\"_\"\$4\"_\"\$5 in a){print a[\$1\"_\"\$2\"_\"\$4\"_\"\$5]}' \
            $sample.compareVAF.txt {} >> {.}.compareVAF.txt)" :::: VCFlist.txt
    fi

    # Extract annoCheck files
    if [ -f "$sample.annoCheck.txt" ]; then
        parallel --memsuspend "$MEM_SUSPEND" "[ -f {.}.annoCheck.txt ] || \
            (awk '{if(NR==1)print}' $sample.annoCheck.txt > {.}.annoCheck.txt; \
            awk -F\"\t\" 'NR==FNR{a[\$1\"_\"\$2\"_\"\$3\"_\"\$4]=\$0;next}(\$1\"_\"\$2\"_\"\$4\"_\"\$5 in a){print a[\$1\"_\"\$2\"_\"\$4\"_\"\$5]}' \
            $sample.annoCheck.txt {} >> {.}.annoCheck.txt)" :::: VCFlist.txt
    fi

    log_success "Per-sample extraction completed"
}

filter_and_summarize() {
    log_step "FILTERING AND SUMMARIZING"

    # Filter annotations using integrated function
    export -f filter_anno
    export COMMON_VARIANTS
    parallel --memsuspend "$MEM_SUSPEND" "filter_anno {.}.annotation.txt" :::: VCFlist.txt

    # Generate summary (TMSP mode only)
    if [[ "$DO_SUMMARY" == true ]]; then
        if [ ! -f Summary.txt ]; then
            log_info "Generating Summary.txt (TMSP mode)"
            parallel sed -n \'1p\' {.}.Filter.txt :::: VCFlist.txt 2>/dev/null | uniq | awk '{print "SAMPLE\t"$0}' \
                > Combine.header.Filter.txt
            parallel "awk '{if(NR!=1)print}' {.}.Filter.txt > {.}.body.Filter.txt" :::: VCFlist.txt
            parallel "awk '{print FILENAME,\"\t\"\$0}' {.}.body.Filter.txt" :::: VCFlist.txt | sed 's/.body.Filter.txt//g' \
                > Combine.body.Filter.txt
            cat Combine.header.Filter.txt Combine.body.Filter.txt > Summary.txt
            rm -f *.body.Filter.txt *.header.Filter.txt
        fi
    else
        log_info "Skipping Summary.txt generation (CEBPA mode)"
    fi

    log_success "Filtering and summarizing completed"
}

# =============================================================================
# EXCEL WRITER FUNCTIONS (integrated from writeTMSPtoXLS.pl and write1WStoXLS.pl)
# =============================================================================

# Write multi-worksheet Excel file for a sample
# Usage: write_sample_xlsx <vcf_file>
write_sample_xlsx() {
    local vcf="$1"
    local sample=$(basename "$vcf" .vcf)

    [ -f "$sample.xlsx" ] && return 0

    perl - "$sample" <<'PERL_SCRIPT'
use strict;
use warnings;
use Excel::Writer::XLSX;

my $sample = $ARGV[0];
my $sampleName = $sample;
$sampleName =~ s/^(.{20}).*/$1/;

# Check required files exist
for my $ext (qw(Filter.txt annotation.txt annoCheck.txt compare.txt compareVAF.txt)) {
    die "$sample.$ext not found" unless -f "$sample.$ext";
}

open(FILTER, "$sample.Filter.txt") or die "$sample.Filter.txt: $!";
open(ANNOTATION, "$sample.annotation.txt") or die "$sample.annotation.txt: $!";
open(CHECK, "$sample.annoCheck.txt") or die "$sample.annoCheck.txt: $!";
open(COMPARE, "$sample.compare.txt") or die "$sample.compare.txt: $!";
open(COMPAREVAF, "$sample.compareVAF.txt") or die "$sample.compareVAF.txt: $!";

my $workbook = Excel::Writer::XLSX->new("$sample.xlsx");
my $filter = $workbook->add_worksheet("$sampleName.filter");
my $annotation = $workbook->add_worksheet("$sampleName");
my $check = $workbook->add_worksheet("$sampleName.check");
my $compare = $workbook->add_worksheet("$sampleName.comp");
my $compareVAF = $workbook->add_worksheet("$sampleName.compVF");

my $row = 0;
while (<FILTER>) {
    chomp;
    my @fields = split('\t', $_);
    my $col = 0;
    for my $c (@fields) { $filter->write($row, $col++, $c); }
    $row++;
}

$row = 0;
while (<ANNOTATION>) {
    chomp;
    my @fields = split('\t', $_);
    my $col = 0;
    for my $c (@fields) { $annotation->write($row, $col++, $c); }
    $row++;
}

$row = 0;
while (<CHECK>) {
    chomp;
    my @fields = split('\t', $_);
    my $col = 0;
    for my $c (@fields) { $check->write($row, $col++, $c); }
    $row++;
}

$row = 0;
while (<COMPARE>) {
    chomp;
    my @fields = split('\t', $_);
    my $col = 0;
    for my $c (@fields) { $compare->write($row, $col++, $c); }
    $row++;
}

$row = 0;
while (<COMPAREVAF>) {
    chomp;
    my @fields = split('\t', $_);
    my $col = 0;
    for my $c (@fields) { $compareVAF->write($row, $col++, $c); }
    $row++;
}

$workbook->close();
PERL_SCRIPT
}

# Write single-worksheet Excel file
# Usage: write_summary_xlsx <xlsx_name> <txt_file>
write_summary_xlsx() {
    local xlsname="$1"
    local txtfile="$2"

    [ -f "$xlsname.xlsx" ] && return 0
    [ ! -f "$txtfile" ] && return 1

    perl - "$xlsname" "$txtfile" <<'PERL_SCRIPT'
use strict;
use warnings;
use Excel::Writer::XLSX;

my ($xlsname, $txtfile) = @ARGV;

open(WS1, "$txtfile") or die "$txtfile: $!";

my $workbook = Excel::Writer::XLSX->new("$xlsname.xlsx");
my $sheetname = $txtfile;
$sheetname =~ s/\.txt$//;
my $worksheet = $workbook->add_worksheet("$sheetname");

my $row = 0;
while (<WS1>) {
    chomp;
    my @fields = split('\t', $_);
    my $col = 0;
    for my $c (@fields) { $worksheet->write($row, $col++, $c); }
    $row++;
}

$workbook->close();
PERL_SCRIPT
}

# Export functions for parallel
export -f write_sample_xlsx
export -f write_summary_xlsx

write_excel() {
    log_step "WRITING EXCEL OUTPUT"

    if ! perl -e 'use Excel::Writer::XLSX' 2>/dev/null; then
        log_warn "Excel::Writer::XLSX not found, skipping Excel output"
        return 0
    fi

    log_info "Writing Excel files using integrated functions"

    # Write per-sample Excel files
    for vcf in $(cat VCFlist.txt); do
        write_sample_xlsx "$vcf"
    done

    # Summary Excel only for TMSP mode
    if [[ "$DO_SUMMARY" == true ]] && [ -f Summary.txt ] && [ ! -f Summary.xlsx ]; then
        write_summary_xlsx "Summary" "Summary.txt"
    fi

    log_success "Excel output completed"
}

# =============================================================================
# MAIN
# =============================================================================

main() {
    local start_time=$(date +%s)

    echo ""
    echo "=============================================================="
    echo "  OPTIMIZED VCF ANNOTATION PIPELINE"
    echo "=============================================================="
    echo ""
    echo "  Mode:     $(echo $MODE | tr '[:lower:]' '[:upper:]')"
    echo "  Database: $VCF_DATABASE"
    echo "  Compare:  ${COMPARE_DBS[*]}"
    echo ""
    echo "=============================================================="
    echo ""

    # Check for VCF files
    VCF_COUNT=$(ls *.vcf 2>/dev/null | grep -v "^Merge" | wc -l)
    if [ "$VCF_COUNT" -eq 0 ]; then
        log_error "No VCF files found in current directory"
        exit 1
    fi
    log_info "Found $VCF_COUNT VCF files"

    # Run pipeline
    merge_vcfs
    run_annotations_parallel
    run_comparisons_parallel
    combine_annotations
    extract_per_sample
    filter_and_summarize
    write_excel

    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    local minutes=$((duration / 60))
    local seconds=$((duration % 60))

    echo ""
    echo "=============================================================="
    echo "  PIPELINE COMPLETED in ${minutes}m ${seconds}s"
    echo "=============================================================="
    echo ""
}

main "$@"
