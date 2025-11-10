#!/bin/bash
#
# Mutect2 Longitudinal Joint Analysis Pipeline
# Performs joint filtering on multiple tumor samples from the same patient
# Input: Unfiltered VCF files from individual Mutect2 runs
#

show_usage() {
    echo "Usage: $0 --patient PATIENT_ID [options] TUMOR1_DIR TUMOR2_DIR [TUMOR3_DIR ...]"
    echo ""
    echo "This script performs joint filtering on longitudinal tumor samples."
    echo ""
    echo "Required arguments:"
    echo "  --patient ID          Patient identifier for output files"
    echo "  TUMOR_DIR             Directories containing tumor scatter_gather results"
    echo "                        (must contain: sample.Mut1.vcf, sample.f1r2.tar.gz, etc.)"
    echo ""
    echo "Optional arguments (for joint genotyping):"
    echo "  --normal BAM          Path to normal BAM file (required only for Step 3)"
    echo "  --skip-genotyping     Skip joint genotyping step (use original VCFs)"
    echo "  --ref FASTA           Reference genome (default: from environment)"
    echo "  --bed FILE            Target regions BED file"
    echo "  --output-dir DIR      Output directory (default: PATIENT_ID_longitudinal)"
    echo "  --label LABEL:DIR     Custom label for a tumor directory"
    echo "                        Example: --label primary:tumor1_scatter_gather"
    echo ""
    echo "Examples:"
    echo "  # Basic usage (skip genotyping, use original VCFs)"
    echo "  $0 --patient P001 \\"
    echo "     tumor1_scatter_gather tumor2_scatter_gather tumor3_scatter_gather"
    echo ""
    echo "  # With joint genotyping (requires normal BAM)"
    echo "  $0 --patient P001 --normal normal.bam \\"
    echo "     tumor1_scatter_gather tumor2_scatter_gather tumor3_scatter_gather"
    echo ""
    echo "  # With custom labels"
    echo "  $0 --patient P001 --skip-genotyping \\"
    echo "     --label primary:tumor1_scatter_gather \\"
    echo "     --label relapse:tumor2_scatter_gather \\"
    echo "     --label metastasis:tumor3_scatter_gather"
}

if [ $# -lt 2 ]; then
    show_usage
    exit 1
fi

# Load modules
module reset
module load GCC/12.3.0
module load SAMtools/1.18
module load BCFtools/1.18
module load Java/11.0.20

export RESOURCE=/home/chattopa/projects/common/Shared_data/WGS_resources
export REFDIR=$RESOURCE/Reference
source $RESOURCE/HouseKeeping/source_me

# Default settings
PATIENT_ID=""
NORMAL_BAM=""
REF_FILE="$REF"
TARGET_BED=""
OUTPUT_DIR=""
DICT_FILE="$DICT"
SKIP_GENOTYPING=0

# Arrays for tumor data
declare -A TUMOR_LABELS
declare -A TUMOR_DIRS
declare -A TUMOR_VCFS
TUMOR_ORDER=()

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --patient)
            PATIENT_ID="$2"
            shift 2
            ;;
        --normal)
            NORMAL_BAM="$2"
            shift 2
            ;;
        --ref)
            REF_FILE="$2"
            shift 2
            ;;
        --bed)
            TARGET_BED="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --skip-genotyping)
            SKIP_GENOTYPING=1
            shift
            ;;
        --label)
            IFS=':' read -r label dir <<< "$2"
            TUMOR_LABELS["$dir"]="$label"
            shift 2
            ;;
        --vcf)
            IFS=':' read -r label vcf dir <<< "$2"
            TUMOR_ORDER+=("$label")
            TUMOR_VCFS["$label"]="$vcf"
            TUMOR_DIRS["$label"]="$dir"
            shift 2
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            # Assume it's a tumor directory
            if [ -d "$1" ]; then
                TUMOR_ORDER+=("$1")
                TUMOR_DIRS["$1"]="$1"
            else
                echo "ERROR: Unknown argument or directory not found: $1"
                exit 1
            fi
            shift
            ;;
    esac
done

# Validate inputs
if [ -z "$PATIENT_ID" ]; then
    echo "ERROR: --patient ID is required"
    exit 1
fi

# Check normal BAM only if not skipping genotyping
if [ $SKIP_GENOTYPING -eq 0 ] && [ -n "$NORMAL_BAM" ] && [ ! -f "$NORMAL_BAM" ]; then
    echo "ERROR: Normal BAM file not found: $NORMAL_BAM"
    echo "Use --skip-genotyping if you want to skip joint genotyping step"
    exit 1
fi

if [ ${#TUMOR_ORDER[@]} -lt 2 ]; then
    echo "ERROR: At least 2 tumor samples required for longitudinal analysis"
    exit 1
fi

# Set output directory
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="${PATIENT_ID}_longitudinal"
fi

mkdir -p "$OUTPUT_DIR"

# Create QC log file for all quality control output
QC_LOG="$OUTPUT_DIR/${PATIENT_ID}_QC.log"
exec 3>"$QC_LOG"  # File descriptor 3 for QC logging

# Function to log to both console and QC log
log_qc() {
    echo "$1" >&3      # Write to QC log
    echo "$1"          # Write to console
}

# Function to log only to QC log (suppress console clutter)
log_qc_only() {
    echo "$1" >&3
}

# Get normal sample name if BAM provided
if [ -n "$NORMAL_BAM" ] && [ -f "$NORMAL_BAM" ]; then
    NORMAL_NAME=$(samtools view -H "$NORMAL_BAM" | grep ^@RG | sed -n 's/.*SM:\(\S*\).*/\1/p;' | uniq)
else
    NORMAL_NAME="N/A"
fi

# Initialize QC log
log_qc "================================================================================"
log_qc "MUTECT2 LONGITUDINAL ANALYSIS - QUALITY CONTROL LOG"
log_qc "================================================================================"
log_qc "Patient ID: $PATIENT_ID"
log_qc "Analysis started: $(date)"
log_qc "================================================================================"
log_qc ""

# Initialize statistics file
MAIN_STATS="$OUTPUT_DIR/${PATIENT_ID}_longitudinal_stats.txt"

cat > "$MAIN_STATS" << EOF
================================================================================
MUTECT2 LONGITUDINAL JOINT ANALYSIS STATISTICS
================================================================================
Patient ID: $PATIENT_ID
Normal Sample: $NORMAL_NAME
Analysis Date: $(date)
Output Directory: $OUTPUT_DIR
================================================================================

EOF

echo "================================================================================"
echo "Mutect2 Longitudinal Joint Analysis Pipeline"
echo "================================================================================"
echo "Patient ID: $PATIENT_ID"
echo "Output directory: $OUTPUT_DIR"
echo "QC Log: $QC_LOG"
echo "Statistics: ${PATIENT_ID}_longitudinal_stats.txt"
echo ""
echo "Progress will be logged to: $QC_LOG"
echo "Check that file for detailed quality control information"
echo "================================================================================"
echo ""

# Detect VCF files and validate tumor directories
log_qc "================================================================================"
log_qc "VALIDATION PHASE: Detecting and validating tumor samples"
log_qc "================================================================================"

echo "Validating tumor samples..."

VALID_TUMORS=()
FAILED_TUMORS=()
LOW_QUALITY_TUMORS=()

for tumor_key in "${TUMOR_ORDER[@]}"; do
    tumor_dir="${TUMOR_DIRS[$tumor_key]}"
    
    # Get label
    if [ -n "${TUMOR_LABELS[$tumor_key]}" ]; then
        label="${TUMOR_LABELS[$tumor_key]}"
    elif [ -n "${TUMOR_VCFS[$tumor_key]}" ]; then
        label="$tumor_key"
    else
        label=$(basename "$tumor_dir" | sed 's/_scatter_gather//')
    fi
    
    # Find VCF if not specified
    if [ -z "${TUMOR_VCFS[$label]}" ]; then
        vcf_file=$(find "$tumor_dir" -name "*.Mut1.vcf" | head -1)
        if [ -z "$vcf_file" ]; then
            log_qc "  ✗ $label: Cannot find Mut1.vcf in $tumor_dir"
            FAILED_TUMORS+=("$label:NO_VCF")
            continue
        fi
        TUMOR_VCFS["$label"]="$vcf_file"
    fi
    
    TUMOR_DIRS["$label"]="$tumor_dir"
    vcf_path="${TUMOR_VCFS[$label]}"
    variant_count=$(grep -v "^#" "$vcf_path" 2>/dev/null | wc -l)
    
    log_qc_only "  Checking $label..."
    log_qc_only "    VCF: $vcf_path"
    
    if [ ! -f "$vcf_path" ]; then
        log_qc "  ✗ $label: VCF file not found"
        FAILED_TUMORS+=("$label:VCF_NOT_FOUND")
        continue
    elif [ $variant_count -eq 0 ]; then
        log_qc "  ✗ $label: VCF has ZERO variants"
        FAILED_TUMORS+=("$label:ZERO_VARIANTS")
        continue
    fi
    
    if [ $variant_count -lt 100000 ]; then
        log_qc "  ⚠ $label: LOW QUALITY - Only $variant_count variants"
        LOW_QUALITY_TUMORS+=("$label:$variant_count")
    fi
    
    log_qc "  ✓ $label: Valid ($variant_count variants)"
    VALID_TUMORS+=("$label")
done

log_qc ""
log_qc "Validation Results:"
log_qc "  Total: $((${#VALID_TUMORS[@]} + ${#FAILED_TUMORS[@]}))"
log_qc "  Passed: ${#VALID_TUMORS[@]}"
log_qc "  Failed: ${#FAILED_TUMORS[@]}"

if [ ${#FAILED_TUMORS[@]} -gt 0 ]; then
    log_qc ""
    log_qc "EXCLUDED SAMPLES:"
    for failed in "${FAILED_TUMORS[@]}"; do
        IFS=':' read -r label reason <<< "$failed"
        log_qc "  ✗ $label: $reason"
    done
fi

if [ ${#LOW_QUALITY_TUMORS[@]} -gt 0 ]; then
    log_qc ""
    log_qc "LOW QUALITY WARNINGS:"
    for low_qual in "${LOW_QUALITY_TUMORS[@]}"; do
        IFS=':' read -r label count <<< "$low_qual"
        log_qc "  ⚠ $label: Only $count variants"
    done
fi

log_qc ""
log_qc "================================================================================"
log_qc ""

if [ ${#VALID_TUMORS[@]} -lt 2 ]; then
    log_qc "ERROR: Need at least 2 valid samples"
    echo "ERROR: Insufficient valid samples. Check $QC_LOG"
    exit 1
fi

TUMOR_ORDER=("${VALID_TUMORS[@]}")
echo "Validation complete: ${#VALID_TUMORS[@]} valid samples"

# Update stats with excluded samples
if [ ${#FAILED_TUMORS[@]} -gt 0 ]; then
    cat >> "$MAIN_STATS" << EOF
EXCLUDED SAMPLES (Failed Quality Checks):
------------------------------------------
EOF
    for failed in "${FAILED_TUMORS[@]}"; do
        IFS=':' read -r label reason <<< "$failed"
        echo "  ✗ $label: $reason" >> "$MAIN_STATS"
    done
    echo "" >> "$MAIN_STATS"
fi

echo "" | tee -a "$MAIN_STATS"
log_qc "Tumor samples (valid):" | tee -a "$MAIN_STATS"
for label in "${TUMOR_ORDER[@]}"; do
    echo "  - $label" | tee -a "$MAIN_STATS"
    log_qc "  - $label"
done
echo "" | tee -a "$MAIN_STATS"
log_qc ""

################################################################################
# STEP 1: COUNT RAW VARIANTS
################################################################################

echo "================================================================================"
echo "STEP 1: Analyzing Individual Raw Variant Calls"
echo "================================================================================"

log_qc "================================================================================"
log_qc "STEP 1: RAW VARIANT COUNTS"
log_qc "================================================================================"

cat >> "$MAIN_STATS" << EOF
================================================================================
STEP 1: INDIVIDUAL RAW VARIANT CALLS
================================================================================

EOF

echo "STEP 1: Individual variant counts:" | tee -a "$MAIN_STATS"
log_qc "Individual variant counts:"

for label in "${TUMOR_ORDER[@]}"; do
    vcf="${TUMOR_VCFS[$label]}"
    count=$(grep -v "^#" "$vcf" 2>/dev/null | wc -l)
    echo "  $label: $count variants" | tee -a "$MAIN_STATS"
    log_qc "  $label: $count variants"
done

echo "" | tee -a "$MAIN_STATS"
log_qc ""
log_qc "================================================================================"
log_qc ""

################################################################################
# STEP 2: MERGE STATS AND CREATE UNION
################################################################################

echo "================================================================================"
echo "STEP 2: Merging Stats and Creating Union of Variant Sites"
echo "================================================================================"

cat >> "$MAIN_STATS" << EOF
================================================================================
STEP 2: UNION OF VARIANT SITES
================================================================================

EOF

# Collect stats files
STATS_FILES=()
for label in "${TUMOR_ORDER[@]}"; do
    tumor_dir="${TUMOR_DIRS[$label]}"
    stats_file=$(find "$tumor_dir" -name "*.Mut1.vcf.stats" | head -1)
    if [ -f "$stats_file" ]; then
        STATS_FILES+=("--stats $stats_file")
    fi
done

# Merge stats
echo "Merging Mutect2 statistics..."
MERGED_STATS="$OUTPUT_DIR/${PATIENT_ID}_merged.stats"

$GATK4_PATH --java-options "-Xms8G -Xmx8G" MergeMutectStats \
    ${STATS_FILES[@]} \
    -O "$MERGED_STATS" 2>&1 | grep -v "INFO" | grep -v "WARN" || true

echo "  Merged stats: $MERGED_STATS" | tee -a "$MAIN_STATS"
echo "" | tee -a "$MAIN_STATS"

# Prepare VCFs for merging
echo "Preparing VCF files..."
log_qc "Preparing and indexing VCF files..."

UNION_INPUT_VCFS=()
for label in "${TUMOR_ORDER[@]}"; do
    vcf="${TUMOR_VCFS[$label]}"
    
    log_qc_only "  Processing $label..."
    
    if [[ ! "$vcf" =~ \.gz$ ]]; then
        log_qc_only "    Compressing..."
        bgzip -c "$vcf" > "$OUTPUT_DIR/${label}_input.vcf.gz" 2>&1 | tee -a "$QC_LOG" || true
        
        log_qc_only "    Indexing..."
        tabix -f -p vcf "$OUTPUT_DIR/${label}_input.vcf.gz" 2>&1 | tee -a "$QC_LOG" || true
        
        UNION_INPUT_VCFS+=("$OUTPUT_DIR/${label}_input.vcf.gz")
    else
        if [ ! -f "${vcf}.tbi" ]; then
            tabix -f -p vcf "$vcf" 2>&1 | tee -a "$QC_LOG" || true
        fi
        UNION_INPUT_VCFS+=("$vcf")
    fi
done

echo "Creating union..."
UNION_VCF="$OUTPUT_DIR/${PATIENT_ID}_union_sites.vcf.gz"

bcftools merge \
    ${UNION_INPUT_VCFS[@]} \
    --force-samples \
    -O z -o "$UNION_VCF" 2>&1 | tee -a "$QC_LOG" | grep -v "^#" || true

tabix -f -p vcf "$UNION_VCF" 2>&1 | tee -a "$QC_LOG" || true

UNION_COUNT=$(bcftools view -H "$UNION_VCF" 2>/dev/null | wc -l)
echo "  Union sites: $UNION_COUNT variants" | tee -a "$MAIN_STATS"
echo "" | tee -a "$MAIN_STATS"

################################################################################
# STEP 3: JOINT GENOTYPING (OPTIONAL)
################################################################################

if [ $SKIP_GENOTYPING -eq 1 ]; then
    echo "STEP 3: Joint Genotyping - SKIPPED"
    cat >> "$MAIN_STATS" << EOF
================================================================================
STEP 3: JOINT GENOTYPING - SKIPPED
================================================================================

EOF
fi

################################################################################
# STEP 4: JOINT FILTERING
################################################################################

echo "================================================================================"
echo "STEP 4: Joint Filtering with Shared Artifact Models"
echo "================================================================================"

cat >> "$MAIN_STATS" << EOF
================================================================================
STEP 4: JOINT FILTERING
================================================================================

EOF

# Learn orientation model with BATCHING
echo "Learning read orientation artifacts..."
log_qc "================================================================================"
log_qc "LEARNING ORIENTATION BIAS MODEL"
log_qc "================================================================================"

F1R2_FILES=()
for label in "${TUMOR_ORDER[@]}"; do
    tumor_dir="${TUMOR_DIRS[$label]}"
    while IFS= read -r f1r2_file; do
        F1R2_FILES+=("$f1r2_file")
    done < <(find "$tumor_dir" -name "*.f1r2.tar.gz" 2>/dev/null)
done

F1R2_COUNT=${#F1R2_FILES[@]}
log_qc "Found $F1R2_COUNT f1r2 files"

if [ $F1R2_COUNT -eq 0 ]; then
    log_qc "No f1r2 files found - skipping"
    JOINT_ORIENTATION=""
elif [ $F1R2_COUNT -gt 100 ]; then
    log_qc "Too many files ($F1R2_COUNT) - using batching"
    
    JOINT_ORIENTATION="$OUTPUT_DIR/${PATIENT_ID}_orientation_model.tar.gz"
    BATCH_SIZE=50
    BATCH_MODELS=()
    
    for ((i=0; i<F1R2_COUNT; i+=BATCH_SIZE)); do
        BATCH_NUM=$((i / BATCH_SIZE + 1))
        BATCH_MODEL="$OUTPUT_DIR/orientation_batch_${BATCH_NUM}.tar.gz"
        BATCH_END=$((i + BATCH_SIZE))
        [ $BATCH_END -gt $F1R2_COUNT ] && BATCH_END=$F1R2_COUNT
        BATCH_COUNT=$((BATCH_END - i))
        
        log_qc "  Batch $BATCH_NUM: $BATCH_COUNT files..."
        
        I_ARGS=""
        for ((j=i; j<BATCH_END; j++)); do
            I_ARGS="$I_ARGS -I ${F1R2_FILES[$j]}"
        done
        
        $GATK4_PATH --java-options "-Xms8G -Xmx8G" LearnReadOrientationModel \
            $I_ARGS -O "$BATCH_MODEL" >> "$QC_LOG" 2>&1
        
        if [ -f "$BATCH_MODEL" ]; then
            BATCH_MODELS+=("$BATCH_MODEL")
            log_qc "    ✓ Batch $BATCH_NUM done"
        fi
    done
    
    if [ ${#BATCH_MODELS[@]} -gt 0 ]; then
        cp "${BATCH_MODELS[0]}" "$JOINT_ORIENTATION"
        rm -f "${BATCH_MODELS[@]}"
        log_qc "✓ Orientation model created"
    else
        JOINT_ORIENTATION=""
    fi
else
    log_qc "Processing all $F1R2_COUNT files..."
    
    JOINT_ORIENTATION="$OUTPUT_DIR/${PATIENT_ID}_orientation_model.tar.gz"
    I_ARGS=""
    for f1r2 in "${F1R2_FILES[@]}"; do
        I_ARGS="$I_ARGS -I $f1r2"
    done
    
    $GATK4_PATH --java-options "-Xms16G -Xmx16G" LearnReadOrientationModel \
        $I_ARGS -O "$JOINT_ORIENTATION" >> "$QC_LOG" 2>&1
    
    [ -f "$JOINT_ORIENTATION" ] && log_qc "✓ Orientation model created"
fi

log_qc ""

# Compress and index input VCFs for filtering
echo "Preparing VCFs for filtering..."
for label in "${TUMOR_ORDER[@]}"; do
    input_vcf="${TUMOR_VCFS[$label]}"
    
    if [[ ! "$input_vcf" =~ \.gz$ ]]; then
        compressed_vcf="${input_vcf}.gz"
        bgzip -c "$input_vcf" > "$compressed_vcf" 2>&1 | tee -a "$QC_LOG" || true
        tabix -f -p vcf "$compressed_vcf" 2>&1 | tee -a "$QC_LOG" || true
        TUMOR_VCFS[$label]="$compressed_vcf"
    fi
done

# Filter samples
echo "Filtering samples..."
log_qc "Filtering samples sequentially with high memory allocation..."

FILTER_SUCCESS=()
FILTER_FAILED=()

for label in "${TUMOR_ORDER[@]}"; do
    log_qc "Filtering $label..."
    
    tumor_dir="${TUMOR_DIRS[$label]}"
    input_vcf="${TUMOR_VCFS[$label]}"
    filtered_vcf="$OUTPUT_DIR/${label}_filtered.vcf.gz"
    
    contam_table=$(find "$tumor_dir" -name "*.contamination.table" | head -1)
    seg_file=$(find "$tumor_dir" -name "*.seg" | head -1)
    
    # Use 100GB for large WGS VCFs - filtering is memory intensive
    FILTER_CMD="$GATK4_PATH --java-options \"-Xms100G -Xmx100G -XX:+UseG1GC\" FilterMutectCalls"
    FILTER_CMD="$FILTER_CMD -R $REF_FILE -V $input_vcf --stats $MERGED_STATS"
    
    [ -n "$JOINT_ORIENTATION" ] && [ -f "$JOINT_ORIENTATION" ] && \
        FILTER_CMD="$FILTER_CMD --ob-priors $JOINT_ORIENTATION"
    [ -f "$contam_table" ] && FILTER_CMD="$FILTER_CMD --contamination-table $contam_table"
    [ -f "$seg_file" ] && FILTER_CMD="$FILTER_CMD --tumor-segmentation $seg_file"
    
    FILTER_CMD="$FILTER_CMD -O $filtered_vcf"
    
    log_qc_only "  Command: $FILTER_CMD"
    
    eval "$FILTER_CMD" > "$OUTPUT_DIR/${label}_filter.log" 2>&1
    FILTER_EXIT=$?
    
    if [ $FILTER_EXIT -ne 0 ]; then
        log_qc "  ✗ FilterMutectCalls FAILED (exit code $FILTER_EXIT)"
        # Check for OOM
        if grep -q "OutOfMemoryError\|oom_kill" "$OUTPUT_DIR/${label}_filter.log" 2>/dev/null; then
            log_qc "    CAUSE: Out of Memory - increase Java heap size"
        fi
        tail -20 "$OUTPUT_DIR/${label}_filter.log" >> "$QC_LOG"
        FILTER_FAILED+=("$label")
        continue
    fi
    
    if [ ! -f "$filtered_vcf" ]; then
        log_qc "  ✗ Output not created"
        FILTER_FAILED+=("$label")
        continue
    fi
    
    total=$(bcftools view -H "$filtered_vcf" 2>/dev/null | wc -l)
    passed=$(bcftools view -f PASS -H "$filtered_vcf" 2>/dev/null | wc -l)
    
    if [ $total -eq 0 ]; then
        log_qc "  ✗ Empty output"
        FILTER_FAILED+=("$label")
        continue
    fi
    
    tabix -f -p vcf "$filtered_vcf" 2>&1 | tee -a "$QC_LOG" || true
    log_qc "  ✓ Total=$total, PASS=$passed"
    echo "  $label: Total=$total, PASS=$passed" >> "$MAIN_STATS"
    TUMOR_VCFS["${label}_filtered"]="$filtered_vcf"
    FILTER_SUCCESS+=("$label")
done

if [ ${#FILTER_SUCCESS[@]} -lt 2 ]; then
    echo "ERROR: Too few samples passed filtering"
    exit 1
fi

TUMOR_ORDER=("${FILTER_SUCCESS[@]}")
echo "" >> "$MAIN_STATS"

################################################################################
# STEP 5: EXTRACT PASS VARIANTS
################################################################################

echo "================================================================================"
echo "STEP 5: Creating Multi-Sample VCF"
echo "================================================================================"

cat >> "$MAIN_STATS" << EOF
================================================================================
STEP 5: MULTI-SAMPLE VCF
================================================================================

EOF

PASS_VCFS=()
PASS_SUCCESS=()

for label in "${TUMOR_ORDER[@]}"; do
    filtered_vcf="${TUMOR_VCFS[${label}_filtered]}"
    pass_vcf="$OUTPUT_DIR/${label}_pass.vcf.gz"
    
    bcftools view -f PASS "$filtered_vcf" -O z -o "$pass_vcf" 2>/dev/null
    tabix -p vcf "$pass_vcf" 2>/dev/null
    
    pass_count=$(bcftools view -H "$pass_vcf" 2>/dev/null | wc -l)
    
    if [ $pass_count -gt 0 ]; then
        echo "  $label: $pass_count PASS variants" | tee -a "$MAIN_STATS"
        PASS_VCFS+=("$pass_vcf")
        PASS_SUCCESS+=("$label")
    else
        echo "  $label: 0 PASS (EXCLUDED)" | tee -a "$MAIN_STATS"
    fi
done

if [ ${#PASS_SUCCESS[@]} -lt 2 ]; then
    echo "ERROR: Too few samples with PASS variants"
    exit 1
fi

TUMOR_ORDER=("${PASS_SUCCESS[@]}")
echo "" >> "$MAIN_STATS"

# Merge
MULTISAMPLE_VCF="$OUTPUT_DIR/${PATIENT_ID}_longitudinal.vcf.gz"
bcftools merge ${PASS_VCFS[@]} --force-samples -O z -o "$MULTISAMPLE_VCF"
tabix -p vcf "$MULTISAMPLE_VCF"

multisample_count=$(bcftools view -H "$MULTISAMPLE_VCF" | wc -l)
echo "  Multi-sample VCF: $multisample_count sites" | tee -a "$MAIN_STATS"
echo "" >> "$MAIN_STATS"

################################################################################
# EVOLUTION ANALYSIS
################################################################################

echo "================================================================================"
echo "EVOLUTION ANALYSIS: Clonal Evolution and Temporal Dynamics"
echo "================================================================================"

EVOLUTION_DIR="$OUTPUT_DIR/evolution_analysis"
mkdir -p "$EVOLUTION_DIR"

log_qc "================================================================================"
log_qc "EVOLUTION ANALYSIS"
log_qc "================================================================================"

cat >> "$MAIN_STATS" << EOF
================================================================================
CLONAL EVOLUTION ANALYSIS
================================================================================

EOF

# Trunk mutations (present in ALL samples)
echo "Identifying trunk mutations (present in all samples)..."
log_qc "Trunk mutations (shared by all ${#TUMOR_ORDER[@]} samples)..."

TRUNK_VCF="$EVOLUTION_DIR/${PATIENT_ID}_trunk_mutations.vcf"
bcftools isec -n "=${#TUMOR_ORDER[@]}" "${PASS_VCFS[@]}" -w 1 -O v -o "$TRUNK_VCF" 2>&1 | tee -a "$QC_LOG" || true
trunk_count=$(grep -v "^#" "$TRUNK_VCF" 2>/dev/null | wc -l)

echo "  ✓ Trunk mutations: $trunk_count"
log_qc "  Trunk mutations: $trunk_count variants"

cat >> "$MAIN_STATS" << EOF
Trunk Mutations (shared by all samples):
  Count: $trunk_count variants
  Definition: Variants present in all ${#TUMOR_ORDER[@]} samples
  Interpretation: These represent early clonal events
  Output: ${PATIENT_ID}_trunk_mutations.vcf

EOF

# Private mutations (unique to each sample)
echo "Identifying private mutations (unique to each sample)..."
log_qc ""
log_qc "Private mutations (unique to individual samples)..."

cat >> "$MAIN_STATS" << EOF
Private Mutations (unique to individual samples):
  Definition: Variants present in only one sample
  Interpretation: Sample-specific clonal evolution
EOF

for i in "${!TUMOR_ORDER[@]}"; do
    label="${TUMOR_ORDER[$i]}"
    private_vcf="$EVOLUTION_DIR/${label}_private.vcf"
    
    # Use -n=1 to get variants present in exactly 1 file
    bcftools isec -n=1 "${PASS_VCFS[@]}" -w "$((i+1))" -O v -o "$private_vcf" 2>&1 | tee -a "$QC_LOG" || true
    private_count=$(grep -v "^#" "$private_vcf" 2>/dev/null | wc -l)
    
    echo "  $label: $private_count private variants"
    log_qc "  $label: $private_count private variants"
    echo "  - $label: $private_count variants" >> "$MAIN_STATS"
done

echo "" >> "$MAIN_STATS"

# Pairwise intersection matrix
echo "Calculating pairwise intersection matrix..."
log_qc ""
log_qc "Pairwise intersection matrix..."

MATRIX_FILE="$EVOLUTION_DIR/pairwise_intersection_matrix.txt"
echo "# Pairwise shared variant counts between samples" > "$MATRIX_FILE"
echo "# Diagonal: Total variants per sample" >> "$MATRIX_FILE"
echo "# Off-diagonal: Shared variants between sample pairs" >> "$MATRIX_FILE"
echo "#" >> "$MATRIX_FILE"

# Header row
printf "%-35s" "Sample" >> "$MATRIX_FILE"
for label in "${TUMOR_ORDER[@]}"; do
    short_label=$(echo "$label" | cut -c1-12)
    printf "%-14s" "$short_label" >> "$MATRIX_FILE"
done
printf "\n" >> "$MATRIX_FILE"

# Data rows
for i in "${!TUMOR_ORDER[@]}"; do
    label_i="${TUMOR_ORDER[$i]}"
    vcf_i="${PASS_VCFS[$i]}"
    printf "%-35s" "$label_i" >> "$MATRIX_FILE"
    
    for j in "${!TUMOR_ORDER[@]}"; do
        label_j="${TUMOR_ORDER[$j]}"
        vcf_j="${PASS_VCFS[$j]}"
        
        if [ $i -eq $j ]; then
            # Diagonal: total variants
            count=$(bcftools view -H "$vcf_i" 2>/dev/null | wc -l)
        else
            # Off-diagonal: shared variants (present in both)
            count=$(bcftools isec -n=2 "$vcf_i" "$vcf_j" -w 1 2>/dev/null | bcftools view -H 2>/dev/null | wc -l)
        fi
        
        printf "%-14s" "$count" >> "$MATRIX_FILE"
    done
    printf "\n" >> "$MATRIX_FILE"
done

echo "  ✓ Matrix saved: pairwise_intersection_matrix.txt"
log_qc "  ✓ Pairwise intersection matrix created"

cat >> "$MAIN_STATS" << EOF
Pairwise Intersection Matrix:
  Output: pairwise_intersection_matrix.txt
  Shows shared variants between all sample pairs
  Diagonal = total variants per sample
  Off-diagonal = shared variants between pairs

EOF

# Temporal analysis (first vs last sample)
echo "Analyzing temporal dynamics (first vs last sample)..."
log_qc ""
log_qc "Temporal analysis (first vs last sample)..."

first="${TUMOR_ORDER[0]}"
last="${TUMOR_ORDER[-1]}"

log_qc "  Comparing: $first (first) → $last (last)"

# Gained mutations: present in LAST but not in FIRST
GAINED_VCF="$EVOLUTION_DIR/${first}_to_${last}_GAINED.vcf"
bcftools isec -C "${PASS_VCFS[0]}" "${PASS_VCFS[-1]}" -w 1 -O v -o "$GAINED_VCF" 2>&1 | tee -a "$QC_LOG" || true
gained=$(grep -v "^#" "$GAINED_VCF" 2>/dev/null | wc -l)

# Lost mutations: present in FIRST but not in LAST
LOST_VCF="$EVOLUTION_DIR/${first}_to_${last}_LOST.vcf"
bcftools isec -C "${PASS_VCFS[-1]}" "${PASS_VCFS[0]}" -w 1 -O v -o "$LOST_VCF" 2>&1 | tee -a "$QC_LOG" || true
lost=$(grep -v "^#" "$LOST_VCF" 2>/dev/null | wc -l)

echo "  ✓ Gained: $gained variants"
echo "  ✓ Lost: $lost variants"
log_qc "  Gained mutations: $gained variants"
log_qc "  Lost mutations: $lost variants"

cat >> "$MAIN_STATS" << EOF
Temporal Dynamics ($first → $last):
  Gained mutations: $gained variants
    Definition: Present in last but not in first sample
    Interpretation: Newly acquired during tumor evolution
    Output: ${first}_to_${last}_GAINED.vcf
  
  Lost mutations: $lost variants
    Definition: Present in first but not in last sample
    Interpretation: Lost during evolution (subclonal in first, or technical)
    Output: ${first}_to_${last}_LOST.vcf

EOF

# Summary statistics
total_unique=$(bcftools view -H "$MULTISAMPLE_VCF" 2>/dev/null | wc -l)
shared_ratio=$(awk "BEGIN {printf \"%.1f\", ($trunk_count/$total_unique)*100}")

cat >> "$MAIN_STATS" << EOF
Evolution Summary:
  Total unique variants: $total_unique
  Trunk (all samples): $trunk_count (${shared_ratio}%)
  Temporal changes: +$gained gained, -$lost lost
  
Interpretation Guide:
  - High trunk %: Samples are closely related (shared ancestry)
  - Low trunk %: Samples diverged significantly or have different origins
  - High gained: Active ongoing evolution
  - High lost: Possible subclonal shifts or technical artifacts

EOF

log_qc ""
log_qc "Evolution analysis complete:"
log_qc "  Trunk mutations: $trunk_count"
log_qc "  Private mutations: calculated for all samples"
log_qc "  Temporal changes: +$gained gained, -$lost lost"
log_qc ""

echo ""

################################################################################
# FINAL
################################################################################

cat >> "$MAIN_STATS" << EOF
================================================================================
FINAL SUMMARY
================================================================================
Analysis complete: $(date)
Final analyzed samples: ${#TUMOR_ORDER[@]}

Key Outputs:
  - Multi-sample VCF: ${PATIENT_ID}_longitudinal.vcf.gz
  - Trunk mutations: evolution_analysis/${PATIENT_ID}_trunk_mutations.vcf
  - Private mutations: evolution_analysis/[SAMPLE]_private.vcf
  - Temporal changes: evolution_analysis/[FIRST]_to_[LAST]_GAINED/LOST.vcf
  - Pairwise matrix: evolution_analysis/pairwise_intersection_matrix.txt
  - Statistics: ${PATIENT_ID}_longitudinal_stats.txt
  - QC Log: ${PATIENT_ID}_QC.log

================================================================================
EOF

echo "================================================================================"
echo "Analysis Complete!"
echo "================================================================================"
echo ""
echo "Main Outputs:"
echo "  Multi-sample VCF: $MULTISAMPLE_VCF"
echo "  Statistics: $MAIN_STATS"
echo "  QC Log: $QC_LOG"
echo ""
echo "Evolution Analysis:"
echo "  Trunk mutations: $TRUNK_VCF"
echo "  Private mutations: $EVOLUTION_DIR/[SAMPLE]_private.vcf"
echo "  Temporal analysis: $EVOLUTION_DIR/${first}_to_${last}_GAINED/LOST.vcf"
echo "  Pairwise matrix: $EVOLUTION_DIR/pairwise_intersection_matrix.txt"
echo ""
echo "finito"

exec 3>&-
