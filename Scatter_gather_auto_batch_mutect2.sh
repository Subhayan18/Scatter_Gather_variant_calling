#!/bin/bash
#
# Mutect2 Auto-Batch Submission Wrapper with Node Management and Verification
# Automatically splits tumor samples across multiple SLURM jobs with failure handling
#

show_usage() {
    echo "Usage: $0 [Normal.bam] [Tumor1.bam] [Tumor2.bam] ... [TumorN.bam] [options]"
    echo ""
    echo "This wrapper automatically batches tumor samples across multiple SLURM jobs"
    echo "with intelligent node selection and output verification."
    echo ""
    echo "Required arguments:"
    echo "  Normal.bam            First BAM file (normal sample)"
    echo "  Tumor*.bam            One or more tumor BAM files"
    echo ""
    echo "Optional arguments:"
    echo "  --nodes 'node1,node2' Comma-separated list of nodes to use"
    echo "                        Example: --nodes 'sn02,sn05,sn12'"
    echo "  --tumors-per-job N    Number of tumors to process per job (default: auto-calculate)"
    echo "  --parallel-jobs N     Scatter-gather parallelization per tumor (default: 15)"
    echo "  --time HH:MM:SS       Wall time per job (default: 60:00:00)"
    echo "  --partition NAME      SLURM partition to use"
    echo "  --max-retries N       Maximum retry attempts for failed jobs (default: 2)"
    echo "  --verify-only         Only verify existing outputs, don't submit new jobs"
    echo "  --debug               Show detailed submission commands without submitting"
    echo "  Reference.fasta       Custom reference (must end in .fa or .fasta)"
    echo "  Target.bed            Target regions (must end in .bed)"
    echo "  Reference.dict        Reference dictionary (must end in .dict)"
    echo ""
    echo "Examples:"
    echo "  # Basic usage with node list"
    echo "  $0 normal.bam tumor*.bam --nodes 'sn02,sn05,sn12'"
    echo ""
    echo "  # Force 1 tumor per job with specific nodes"
    echo "  $0 normal.bam tumor*.bam --nodes 'sn05,sn12' --tumors-per-job 1"
    echo ""
    echo "  # Verify existing outputs and resubmit failures"
    echo "  $0 normal.bam tumor*.bam --nodes 'sn02,sn05' --verify-only"
}

if [ $# -lt 2 ]; then
    show_usage
    exit 1
fi

# Default settings
TUMORS_PER_JOB="auto"
PARALLEL_JOBS=15
WALLTIME="60:00:00"
PARTITION=""
NODE_LIST=""
MAX_RETRIES=2
VERIFY_ONLY=0
DEBUG_MODE=0

# Parse arguments
NORMAL_BAM=""
TUMOR_BAMS=()
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case $1 in
        --tumors-per-job)
            TUMORS_PER_JOB="$2"
            shift 2
            ;;
        --parallel-jobs)
            PARALLEL_JOBS="$2"
            shift 2
            ;;
        --time)
            WALLTIME="$2"
            shift 2
            ;;
        --partition)
            PARTITION="$2"
            shift 2
            ;;
        --nodes)
            NODE_LIST="$2"
            shift 2
            ;;
        --max-retries)
            MAX_RETRIES="$2"
            shift 2
            ;;
        --verify-only)
            VERIFY_ONLY=1
            shift
            ;;
        --debug)
            DEBUG_MODE=1
            shift
            ;;
        *.bam)
            if [ -z "$NORMAL_BAM" ]; then
                NORMAL_BAM="$1"
            else
                TUMOR_BAMS+=("$1")
            fi
            shift
            ;;
        *.fa|*.fasta|*.bed|*.dict)
            EXTRA_ARGS+=("$1")
            shift
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

if [ -z "$NORMAL_BAM" ] || [ ${#TUMOR_BAMS[@]} -eq 0 ]; then
    echo "ERROR: Must provide at least one normal and one tumor BAM file"
    exit 1
fi

# Convert node list to array
if [ -n "$NODE_LIST" ]; then
    IFS=',' read -ra NODES <<< "$NODE_LIST"
    # Trim whitespace
    for i in "${!NODES[@]}"; do
        NODES[$i]=$(echo "${NODES[$i]}" | xargs)
    done
else
    echo "WARNING: No node list specified. Jobs may fail if cluster requires --nodelist"
    echo "Specify nodes with: --nodes 'node1,node2,node3'"
    NODES=()
fi

# Calculate optimal tumors per job
if [ "$TUMORS_PER_JOB" == "auto" ]; then
    CPUS_PER_TUMOR=$((PARALLEL_JOBS * 2 + 8))
    TUMORS_PER_JOB=$((48 / CPUS_PER_TUMOR))
    MEM_LIMITED=$((252 / 50))
    
    if [ $MEM_LIMITED -lt $TUMORS_PER_JOB ]; then
        TUMORS_PER_JOB=$MEM_LIMITED
    fi
    
    if [ $TUMORS_PER_JOB -lt 1 ]; then
        TUMORS_PER_JOB=1
    elif [ $TUMORS_PER_JOB -gt 3 ]; then
        TUMORS_PER_JOB=3
    fi
fi

# Function to verify tumor output
verify_tumor_output() {
    local tumor_bam=$1
    local t_name=$(basename $tumor_bam .bam)
    local vcf_file="${t_name}.MUTECT2.m1.vcf"
    local stats_file="${t_name}.filtering_stats.txt"
    
    if [ ! -f "$vcf_file" ]; then
        echo "MISSING"
        return 1
    fi
    
    if [ ! -f "$stats_file" ]; then
        echo "MISSING"
        return 1
    fi
    
    local variant_count=$(grep -v "^#" "$vcf_file" 2>/dev/null | wc -l)
    if [ $variant_count -eq 0 ]; then
        echo "EMPTY"
        return 1
    fi
    
    if ! grep -q "FINAL VARIANT BREAKDOWN" "$stats_file" 2>/dev/null; then
        echo "INCOMPLETE"
        return 1
    fi
    
    echo "OK:$variant_count"
    return 0
}

# Function to get available node (round-robin)
get_next_node() {
    local batch_num=$1
    if [ ${#NODES[@]} -eq 0 ]; then
        echo ""
        return
    fi
    local node_idx=$((batch_num % ${#NODES[@]}))
    echo "${NODES[$node_idx]}"
}

echo "================================================"
echo "Mutect2 Auto-Batch Submission with Verification"
echo "================================================"
echo "Normal BAM: $NORMAL_BAM"
echo "Total tumor samples: ${#TUMOR_BAMS[@]}"
echo "Tumors per job: $TUMORS_PER_JOB"
echo "Parallel jobs per tumor: $PARALLEL_JOBS"

if [ ${#NODES[@]} -gt 0 ]; then
    echo "Available nodes: ${NODES[*]}"
    echo ""
    echo "Checking node availability..."
    for node in "${NODES[@]}"; do
        NODE_INFO=$(sinfo -n "$node" -h -o "%N %T %C" 2>&1)
        if [ $? -eq 0 ] && [ -n "$NODE_INFO" ]; then
            echo "  $node: $NODE_INFO"
        else
            echo "  $node: NOT FOUND or INACCESSIBLE"
        fi
    done
    echo ""
else
    echo "Available nodes: ANY (no restriction)"
    echo "WARNING: No node list specified."
    echo ""
fi

echo "Max retries per job: $MAX_RETRIES"
echo ""

# Verify existing outputs
echo "Verifying existing outputs..."
echo "----------------------------"
FAILED_TUMORS=()
COMPLETE_TUMORS=()

for tumor in "${TUMOR_BAMS[@]}"; do
    t_name=$(basename $tumor .bam)
    status=$(verify_tumor_output "$tumor")
    
    if [[ $status == OK:* ]]; then
        variant_count=${status#OK:}
        echo "  ✓ $t_name: Complete ($variant_count variants)"
        COMPLETE_TUMORS+=("$tumor")
    else
        echo "  ✗ $t_name: $status"
        FAILED_TUMORS+=("$tumor")
    fi
done

echo ""
echo "Summary: ${#COMPLETE_TUMORS[@]} complete, ${#FAILED_TUMORS[@]} need processing"
echo ""

if [ ${#FAILED_TUMORS[@]} -eq 0 ]; then
    echo "All tumor samples have been successfully processed!"
    exit 0
fi

if [ $VERIFY_ONLY -eq 1 ]; then
    echo "Verify-only mode: No new jobs will be submitted"
    echo "Tumors needing processing:"
    for tumor in "${FAILED_TUMORS[@]}"; do
        echo "  - $tumor"
    done
    exit 0
fi

# Update tumor list to only process failed ones
TUMOR_BAMS=("${FAILED_TUMORS[@]}")

# Calculate number of batch jobs needed
TOTAL_TUMORS=${#TUMOR_BAMS[@]}
NUM_BATCH_JOBS=$(( (TOTAL_TUMORS + TUMORS_PER_JOB - 1) / TUMORS_PER_JOB ))

echo "Will submit $NUM_BATCH_JOBS batch job(s) for incomplete samples"
echo "================================================"
echo ""

# Create temporary directory for batch scripts
BATCH_DIR="mutect2_batches_$(date +%Y%m%d_%H%M%S)"
mkdir -p $BATCH_DIR

# Create tracking file
TRACKING_FILE="$BATCH_DIR/job_tracking.txt"
echo "# Batch tracking file created $(date)" > $TRACKING_FILE
echo "# Format: BATCH_NUM|JOB_ID|NODE|STATUS|TUMORS" >> $TRACKING_FILE

# Get script directory for sourcing function
SCRIPT_DIR=$(dirname $(readlink -f $0))

# Process each batch
for ((batch=0; batch<NUM_BATCH_JOBS; batch++)); do
    START_IDX=$((batch * TUMORS_PER_JOB))
    END_IDX=$(( (batch + 1) * TUMORS_PER_JOB ))
    if [ $END_IDX -gt $TOTAL_TUMORS ]; then
        END_IDX=$TOTAL_TUMORS
    fi
    
    # Get tumors for this batch
    BATCH_TUMORS=("${TUMOR_BAMS[@]:$START_IDX:$((END_IDX - START_IDX))}")
    
    # Get node for this batch
    ASSIGNED_NODE=$(get_next_node $batch)
    
    BATCH_NUM=$((batch + 1))
    echo "Batch $BATCH_NUM/$NUM_BATCH_JOBS: ${#BATCH_TUMORS[@]} tumor(s)"
    if [ -n "$ASSIGNED_NODE" ]; then
        echo "  Assigned node: $ASSIGNED_NODE"
    fi
    for tumor in "${BATCH_TUMORS[@]}"; do
        echo "  - $(basename $tumor)"
    done
    
    # Create batch script
    BATCH_SCRIPT="$BATCH_DIR/batch_${BATCH_NUM}.sh"
    
    cat > $BATCH_SCRIPT << EOFSCRIPT
#!/bin/sh
#SBATCH -t $WALLTIME
#SBATCH --cpus-per-task 48
#SBATCH --mem 252G
#SBATCH -J MuTct2_B_${BATCH_NUM}
#SBATCH -o mutect2_batch_${BATCH_NUM}_%j.out
#SBATCH -e mutect2_batch_${BATCH_NUM}_%j.err
EOFSCRIPT

    # Add partition if specified
    if [ -n "$PARTITION" ]; then
        echo "#SBATCH -p $PARTITION" >> $BATCH_SCRIPT
    fi
    
    # Add nodelist if specified
    if [ -n "$ASSIGNED_NODE" ]; then
        echo "#SBATCH --nodelist=$ASSIGNED_NODE" >> $BATCH_SCRIPT
    fi
    
    # Add main script content
    cat >> $BATCH_SCRIPT << 'EOFSCRIPT2'

module reset
module load GCC/12.3.0
module load SAMtools/1.18
module load FastQC/0.11.9-Java-11
module load BCFtools/1.18
module load Java/11.0.20
module load parallel/20230722

export RESOURCE=/home/chattopa/projects/common/Shared_data/WGS_resources
export REFDIR=$RESOURCE/Reference
source $RESOURCE/HouseKeeping/source_me

# Create progress log
PROGRESS_LOG="mutect2_batch_BATCH_NUM_progress_${SLURM_JOB_ID}.log"
exec 3>$PROGRESS_LOG

log_milestone() {
    echo "$1" >&3
    echo "$1"
}

log_milestone "================================================"
log_milestone "Batch BATCH_NUM - Job ID: ${SLURM_JOB_ID}"
log_milestone "Node: $(hostname)"
log_milestone "Started: $(date)"
log_milestone "================================================"

export TMPDIR=${TMPDIR:-/tmp}
mkdir -p $TMPDIR

PARALLEL_JOBS=PARALLEL_JOBS_VAL

REF_FILE="$REF"
B=$REFDIR/Homo_sapiens_assembly38.fasta.bed
DICT_FILE="$DICT"
EOFSCRIPT2

    # Add extra args
    for extra_arg in "${EXTRA_ARGS[@]}"; do
        if [[ "$extra_arg" == *.fa ]] || [[ "$extra_arg" == *.fasta ]]; then
            echo "REF_FILE=\"$extra_arg\"" >> $BATCH_SCRIPT
        elif [[ "$extra_arg" == *.bed ]]; then
            echo "B=\"$extra_arg\"" >> $BATCH_SCRIPT
        elif [[ "$extra_arg" == *.dict ]]; then
            echo "DICT_FILE=\"$extra_arg\"" >> $BATCH_SCRIPT
        fi
    done
    
    # Add memory calculation and tumor BAMs
    cat >> $BATCH_SCRIPT << 'EOFSCRIPT3'

MEMORY_PER_JOB=$((200 / PARALLEL_JOBS))
if [ $MEMORY_PER_JOB -lt 4 ]; then
    MEMORY_PER_JOB=4
fi
JAVA_OPTS="-Xms${MEMORY_PER_JOB}G -Xmx${MEMORY_PER_JOB}G -XX:+UseG1GC -XX:ParallelGCThreads=2 -Djava.io.tmpdir=$TMPDIR"

EOFSCRIPT3

    echo "NORMAL_BAM=\"$NORMAL_BAM\"" >> $BATCH_SCRIPT
    echo "TUMOR_BAMS=(" >> $BATCH_SCRIPT
    for tumor in "${BATCH_TUMORS[@]}"; do
        echo "\"$tumor\"" >> $BATCH_SCRIPT
    done
    echo ")" >> $BATCH_SCRIPT
    
    # Add processing logic
    cat >> $BATCH_SCRIPT << 'EOFSCRIPT4'

n_n=$(samtools view -H $NORMAL_BAM |grep ^@RG | sed -n 's/.*SM:\(\S*\).*/\1/p;' |uniq )

log_milestone "Processing ${#TUMOR_BAMS[@]} tumor(s) in this batch"
for tumor in "${TUMOR_BAMS[@]}"; do
    log_milestone "  - $tumor"
done
log_milestone ""

SCRIPT_DIR="SCRIPT_DIR_VAL"
source $SCRIPT_DIR/mutect2_process_function.sh

export -f process_tumor
export -f log_milestone
export GATK4_PATH REF_FILE B DICT_FILE PARALLEL_JOBS GERMLINE PON BIALLELIC_EXAC DBSNP TMPDIR n_n NORMAL_BAM PROGRESS_LOG JAVA_OPTS

printf "%s\n" "${TUMOR_BAMS[@]}" | \
parallel -j ${#TUMOR_BAMS[@]} --line-buffer "process_tumor {} $NORMAL_BAM" 2>&1 | \
grep -v "INFO" | grep -v "WARN" || true

log_milestone ""
log_milestone "================================================"
log_milestone "Batch processing completed at $(date)"
log_milestone "Verifying outputs..."
log_milestone "================================================"

ALL_SUCCESS=1
for tumor in "${TUMOR_BAMS[@]}"; do
    t_name=$(basename $tumor .bam)
    vcf_file="${t_name}.MUTECT2.m1.vcf"
    stats_file="${t_name}.filtering_stats.txt"
    
    if [ ! -f "$vcf_file" ]; then
        log_milestone "  ✗ $t_name: VCF file missing"
        ALL_SUCCESS=0
    elif [ ! -s "$vcf_file" ]; then
        log_milestone "  ✗ $t_name: VCF file empty"
        ALL_SUCCESS=0
    elif ! grep -q "^#CHROM" "$vcf_file" 2>/dev/null; then
        log_milestone "  ✗ $t_name: VCF file incomplete"
        ALL_SUCCESS=0
    else
        variant_count=$(grep -v "^#" "$vcf_file" 2>/dev/null | wc -l)
        if [ ! -f "$stats_file" ] || ! grep -q "FINAL VARIANT BREAKDOWN" "$stats_file" 2>/dev/null; then
            log_milestone "  ⚠ $t_name: VCF complete ($variant_count variants) but stats incomplete"
        else
            log_milestone "  ✓ $t_name: Complete ($variant_count variants)"
        fi
    fi
done

if [ $ALL_SUCCESS -eq 1 ]; then
    log_milestone ""
    log_milestone "All outputs verified successfully!"
else
    log_milestone ""
    log_milestone "WARNING: Some outputs are missing or incomplete"
    exit 1
fi

exec 3>&-
log_milestone "finito"
EOFSCRIPT4

    # Replace placeholders
    sed -i "s|BATCH_NUM|${BATCH_NUM}|g" $BATCH_SCRIPT
    sed -i "s|PARALLEL_JOBS_VAL|${PARALLEL_JOBS}|g" $BATCH_SCRIPT
    sed -i "s|SCRIPT_DIR_VAL|${SCRIPT_DIR}|g" $BATCH_SCRIPT
    
    chmod +x $BATCH_SCRIPT
    
    # Debug mode
    if [ $DEBUG_MODE -eq 1 ]; then
        echo "  [DEBUG] Would submit: sbatch $BATCH_SCRIPT"
        echo "  [DEBUG] SBATCH directives:"
        grep "^#SBATCH" $BATCH_SCRIPT
        echo ""
        continue
    fi
    
    # Submit job
    SUBMIT_OUTPUT=$(sbatch $BATCH_SCRIPT 2>&1)
    SUBMIT_EXIT=$?
    
    if [ $SUBMIT_EXIT -eq 0 ]; then
        JOB_ID=$(echo "$SUBMIT_OUTPUT" | grep -oE '[0-9]+$')
        if [[ $JOB_ID =~ ^[0-9]+$ ]]; then
            echo "  ✓ Submitted as job $JOB_ID"
            
            TUMOR_NAMES=""
            for tumor in "${BATCH_TUMORS[@]}"; do
                TUMOR_NAMES="${TUMOR_NAMES},$(basename $tumor .bam)"
            done
            TUMOR_NAMES=${TUMOR_NAMES:1}
            
            echo "${BATCH_NUM}|${JOB_ID}|${ASSIGNED_NODE}|SUBMITTED|${TUMOR_NAMES}" >> $TRACKING_FILE
        else
            echo "  ✗ ERROR: Could not parse job ID from: $SUBMIT_OUTPUT"
        fi
    else
        echo "  ✗ ERROR: Submission failed"
        echo "  Output: $SUBMIT_OUTPUT"
        
        if [[ "$SUBMIT_OUTPUT" =~ "available" ]]; then
            echo "  Diagnosis: Node availability issue"
            echo "  Suggestions:"
            echo "    1. Check node names: sinfo -N"
            echo "    2. Try without --nodes flag"
            echo "    3. Use --debug to see full commands"
        fi
    fi
    echo ""
done

echo "================================================"
echo "Submission complete!"
echo "================================================"
echo ""

if [ $DEBUG_MODE -eq 1 ]; then
    echo "DEBUG MODE: No jobs were actually submitted"
    echo "Batch scripts created in: $BATCH_DIR/"
    exit 0
fi

echo "Tracking file: $TRACKING_FILE"
echo "Batch scripts: $BATCH_DIR/"
echo ""

SUCCESSFUL=$(grep -c "SUBMITTED" "$TRACKING_FILE" 2>/dev/null || echo "0")
FAILED=$((NUM_BATCH_JOBS - SUCCESSFUL))

if [ $FAILED -gt 0 ]; then
    echo "WARNING: $FAILED batch job(s) failed to submit"
    echo "Check errors above"
    echo ""
fi

if [ $SUCCESSFUL -gt 0 ]; then
    echo "Successfully submitted: $SUCCESSFUL job(s)"
    echo ""
    echo "Monitor: squeue -u \$USER"
    echo ""
fi
