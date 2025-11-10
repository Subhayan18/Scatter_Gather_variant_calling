# MUTECT2 LONGITUDINAL ANALYSIS PIPELINE - COMPLETE REFERENCE

**Version:** 1.1  
**Date:** November 2025  
**Script:** mutect2_longitudinal.sh  
**Total Lines:** 915  

---

## TABLE OF CONTENTS

1. [Overview](#overview)
2. [Pipeline Architecture](#pipeline-architecture)
3. [Requirements](#requirements)
4. [Usage Guide](#usage-guide)
5. [Command-Line Arguments](#command-line-arguments)
6. [Pipeline Steps (Detailed)](#pipeline-steps-detailed)
7. [Output Files](#output-files)
8. [Internal Functions](#internal-functions)
9. [Quality Control & Validation](#quality-control--validation)
10. [Evolution Analysis](#evolution-analysis)
11. [Memory & Performance](#memory--performance)
12. [Troubleshooting](#troubleshooting)
13. [Examples](#examples)
14. [Technical Notes](#technical-notes)

---

## OVERVIEW

### Purpose
The Mutect2 Longitudinal Analysis Pipeline performs joint variant filtering and evolutionary analysis on multiple tumor samples from the same patient. It identifies trunk mutations (shared across all samples), private mutations (sample-specific), and temporal dynamics (gained/lost variants).

### Key Features
- **Joint filtering** using shared artifact models across all samples
- **Clonal evolution analysis** (trunk, private, temporal mutations)
- **Pairwise intersection matrix** for sample relationships
- **Comprehensive quality control** with detailed logging
- **Optional joint genotyping** for enhanced sensitivity
- **Enhanced statistics** with biological interpretations

### Use Cases
- **Longitudinal tumor analysis**: Track evolution over time
- **Multi-region sampling**: Compare primary and metastatic sites
- **Treatment response**: Analyze pre/post-treatment changes
- **Relapse studies**: Identify variants associated with recurrence

---

## PIPELINE ARCHITECTURE

### Workflow Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                        INPUT VALIDATION                          │
│  • Detect tumor directories                                      │
│  • Find VCF and artifact files                                   │
│  • Validate variant counts                                       │
│  • Flag low-quality samples                                      │
└────────────────────────┬────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────────────┐
│                    STEP 1: RAW VARIANT COUNTS                    │
│  • Count variants per sample                                     │
│  • Report to statistics file                                     │
└────────────────────────┬────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────────────┐
│              STEP 2: MERGE STATS & CREATE UNION                  │
│  • Merge Mutect2 statistics (MergeMutectStats)                  │
│  • Create union of all variant sites                             │
│  • Prepare for joint filtering                                   │
└────────────────────────┬────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────────────┐
│              STEP 3: JOINT GENOTYPING (Optional)                 │
│  • Re-genotype all samples at union sites                        │
│  • Uses normal BAM for improved accuracy                         │
│  • [SKIPPED if --skip-genotyping flag used]                      │
└────────────────────────┬────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────────────┐
│         STEP 4: JOINT FILTERING WITH SHARED MODELS               │
│  • Learn read orientation artifacts (LearnReadOrientationModel) │
│  • Apply joint artifact model to all samples                     │
│  • Filter variants (FilterMutectCalls)                          │
│  • Extract PASS variants                                         │
└────────────────────────┬────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────────────┐
│              STEP 5: CREATE MULTI-SAMPLE VCF                     │
│  • Extract PASS variants per sample                              │
│  • Merge into multi-sample VCF (bcftools merge)                 │
│  • Index final VCF                                               │
└────────────────────────┬────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────────────┐
│                    EVOLUTION ANALYSIS                            │
│  • Trunk mutations (shared by all)                              │
│  • Private mutations (unique to each)                           │
│  • Temporal dynamics (gained/lost)                              │
│  • Pairwise intersection matrix                                 │
│  • Enhanced statistics with interpretations                     │
└────────────────────────┬────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────────────┐
│                       FINAL OUTPUTS                              │
│  • Multi-sample VCF                                              │
│  • Evolution VCF files                                           │
│  • Statistics report                                             │
│  • QC log                                                        │
└─────────────────────────────────────────────────────────────────┘
```

### Data Flow

```
Input Directories → VCF Detection → Validation → Joint Filtering → Evolution Analysis → Outputs
     ↓                   ↓              ↓              ↓                    ↓              ↓
scatter_gather/    *.Mut1.vcf    QC checks     FilterMutect      bcftools isec     Final VCFs
├── VCF            *.f1r2.tar.gz  Counts       Orientation       Private/Trunk      Statistics
├── f1r2           *.stats        Low-quality   Contamination    Gained/Lost        QC logs
├── stats          *.contamination warnings     Segmentation     Pairwise matrix    
└── contamination  *.seg
```

---

## REQUIREMENTS

### Software Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| **GATK4** | 4.3.0.0+ | Mutect2, filtering, genotyping |
| **BCFtools** | 1.18+ | VCF manipulation, merging, intersection |
| **SAMtools** | 1.18+ | BAM file handling, sample name extraction |
| **Java** | 11.0.20+ | GATK4 runtime |
| **Bash** | 4.0+ | Script execution |

### System Requirements

- **Memory:** 100GB+ RAM recommended (for WGS filtering)
- **Storage:** ~10-50GB per sample for intermediate files
- **CPU:** Multi-core recommended for parallel processing

### Input Files Per Sample

Required files in each `scatter_gather` directory:

```
tumor_scatter_gather/
├── sample.Mut1.vcf             # Required: Unfiltered VCF
├── sample.Mut1.vcf.stats       # Required: Mutect2 statistics
├── sample.f1r2.tar.gz          # Required: Orientation artifacts (multiple files)
├── sample.contamination.table  # Required: Contamination estimates
└── sample.seg                  # Required: Segmentation for CNV
```

---

## USAGE GUIDE

### Basic Syntax

```bash
./mutect2_longitudinal.sh --patient PATIENT_ID [OPTIONS] TUMOR_DIR1 TUMOR_DIR2 [TUMOR_DIR3 ...]
```

### Quick Start Examples

#### Example 1: Basic Usage (No Genotyping)
```bash
./mutect2_longitudinal.sh \
  --patient Patient01 \
  --skip-genotyping \
  tumor1_scatter_gather \
  tumor2_scatter_gather \
  tumor3_scatter_gather
```

#### Example 2: With Joint Genotyping
```bash
./mutect2_longitudinal.sh \
  --patient Patient01 \
  --normal /path/to/normal.bam \
  tumor1_scatter_gather \
  tumor2_scatter_gather \
  tumor3_scatter_gather
```

#### Example 3: Custom Labels and Output Directory
```bash
./mutect2_longitudinal.sh \
  --patient Patient01 \
  --skip-genotyping \
  --output-dir Patient01_evolution \
  --label primary:tumor1_scatter_gather \
  --label relapse:tumor2_scatter_gather \
  --label metastasis:tumor3_scatter_gather \
  tumor1_scatter_gather \
  tumor2_scatter_gather \
  tumor3_scatter_gather
```

#### Example 4: With Custom Reference
```bash
./mutect2_longitudinal.sh \
  --patient Patient01 \
  --ref /path/to/reference.fa \
  --skip-genotyping \
  tumor1_scatter_gather \
  tumor2_scatter_gather
```

---

## COMMAND-LINE ARGUMENTS

### Required Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `--patient ID` | Patient identifier for output files | `--patient Patient01` |
| `TUMOR_DIR` | Directories containing scatter_gather results (≥2 required) | `tumor1_scatter_gather` |

### Optional Arguments

| Argument | Description | Default | Example |
|----------|-------------|---------|---------|
| `--normal BAM` | Normal BAM file path (required for Step 3 genotyping) | None | `--normal normal.bam` |
| `--skip-genotyping` | Skip joint genotyping step (use original VCFs) | Off | `--skip-genotyping` |
| `--ref FASTA` | Reference genome FASTA file | From environment (`$REF`) | `--ref hg38.fa` |
| `--bed FILE` | Target regions BED file | None | `--bed targets.bed` |
| `--output-dir DIR` | Output directory | `PATIENT_ID_longitudinal` | `--output-dir results` |
| `--label LABEL:DIR` | Custom label for tumor directory | Auto-generated | `--label primary:tumor1` |
| `--vcf LABEL:VCF:DIR` | Specify VCF directly (advanced) | Auto-detected | See script |
| `-h, --help` | Show usage information | - | `--help` |

### Argument Details

#### `--patient ID`
- **Required:** Yes
- **Purpose:** Patient identifier used in all output filenames
- **Format:** Alphanumeric string (no spaces)
- **Impact:** All output files will be prefixed with this ID

#### `--skip-genotyping`
- **Required:** No
- **Purpose:** Skip Step 3 (joint genotyping)
- **Use when:** 
  - Original VCFs have sufficient quality
  - Normal BAM unavailable
  - Faster processing needed
- **Impact:** Uses original Mut1.vcf files directly for filtering

#### `--label LABEL:DIR`
- **Required:** No
- **Purpose:** Assign custom names to samples
- **Format:** `label:directory_path`
- **Example:** `--label primary:tumor1_scatter_gather`
- **Use cases:**
  - Meaningful biological names (primary, relapse, metastasis)
  - Time points (baseline, month6, month12)
  - Tissue types (liver_met, lung_met)

---

## PIPELINE STEPS (DETAILED)

### VALIDATION PHASE

**Lines:** 215-332  
**Purpose:** Validate input directories and detect required files  

#### Process:
1. **Detect VCF files:**
   - Search for `*.Mut1.vcf` in each directory
   - Extract sample labels from filenames or use custom labels
   
2. **Validate each sample:**
   - Check VCF file exists
   - Count variants (must be > 0)
   - Flag samples with < 100,000 variants as low quality
   
3. **Classification:**
   - **Valid:** VCF exists, variants > 0
   - **Failed:** VCF missing, zero variants
   - **Low Quality:** < 100,000 variants (warning, not excluded)

#### Output:
```
Validation Results:
  Total: 8
  Passed: 8
  Failed: 0

LOW QUALITY WARNINGS:
  ⚠ Sample1: Only 25000 variants
```

#### Quality Thresholds:
- **Minimum variants:** > 0 (hard requirement)
- **Warning threshold:** < 100,000 variants
- **Minimum samples:** ≥ 2 valid samples required

---

### STEP 1: RAW VARIANT COUNTS

**Lines:** 334-367  
**Purpose:** Count and report raw variant calls per sample  
**Tools:** grep, wc  

#### Process:
```bash
for label in "${TUMOR_ORDER[@]}"; do
    vcf="${TUMOR_VCFS[$label]}"
    count=$(grep -v "^#" "$vcf" 2>/dev/null | wc -l)
    echo "  $label: $count variants"
done
```

#### Output Example:
```
STEP 1: Individual variant counts:
  Sample1: 25000 variants
  Sample2: 750000 variants
  Sample3: 700000 variants
```

#### Purpose:
- Document raw variant calls before filtering
- Identify potential quality issues
- Baseline for filtering effectiveness calculation

---

### STEP 2: MERGE STATS & CREATE UNION

**Lines:** 368-450  
**Purpose:** Merge Mutect2 statistics and create union of variant sites  
**Tools:** GATK MergeMutectStats, bcftools  
**Memory:** 8GB  

#### Process:

**2a. Merge Statistics:**
```bash
GATK MergeMutectStats \
  --stats sample1.stats \
  --stats sample2.stats \
  --stats sample3.stats \
  -O merged.stats
```

**Purpose:** Combine filtering statistics across all samples for joint filtering

**2b. Create Union VCF:**
```bash
# Compress and index input VCFs
bgzip -c sample1.vcf > sample1.vcf.gz
tabix -p vcf sample1.vcf.gz

# Create union of all sites
bcftools merge sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz \
  --force-samples -O v -o union.vcf
```

**Purpose:** Identify all unique variant sites across samples for joint genotyping

#### Output:
```
Merged stats: PatientX_merged.stats
Union sites: 3000000 variants
```

#### Technical Details:
- **Compression:** bgzip for bgzf format (required for tabix)
- **Indexing:** tabix creates .tbi index files
- **Force samples:** Handles duplicate sample names
- **Output format:** Uncompressed VCF for downstream processing

---

### STEP 3: JOINT GENOTYPING (Optional)

**Lines:** 452-514  
**Purpose:** Re-genotype all samples at union sites using normal BAM  
**Tools:** GATK Mutect2 (genotyping mode)  
**Memory:** 30GB per sample  
**Status:** SKIPPED if `--skip-genotyping` flag used  

#### When to Use:
- ✅ Normal BAM available
- ✅ Want maximum sensitivity
- ✅ Samples have sufficient coverage
- ❌ Skip if: fast processing needed, normal BAM unavailable

#### Process:
```bash
for each sample:
  GATK Mutect2 \
    -I tumor.bam \
    -I normal.bam \
    -normal NORMAL_NAME \
    -R reference.fa \
    --alleles union.vcf \
    --genotype-germline-sites true \
    --genotype-pon-sites true \
    -O sample_genotyped.vcf
```

#### Parameters:
- `--alleles`: Force genotyping at union sites
- `--genotype-germline-sites`: Genotype germline variants
- `--genotype-pon-sites`: Genotype panel of normals sites
- `-normal`: Specify normal sample name

#### Output:
```
STEP 3: Joint Genotyping - SKIPPED
(or)
STEP 3: Genotyping samples at union sites...
  Sample1: 3000000 sites genotyped
  Sample2: 3000000 sites genotyped
  Sample3: 3000000 sites genotyped
```

---

### STEP 4: JOINT FILTERING

**Lines:** 515-632  
**Purpose:** Apply joint artifact models and filter variants  
**Tools:** GATK LearnReadOrientationModel, FilterMutectCalls  
**Memory:** 16GB for orientation, 100GB for filtering  

#### Process:

**4a. Learn Read Orientation Artifacts:**

```bash
# Collect all f1r2 files
find scatter_gather_dirs -name "*.f1r2.tar.gz"

# If > 50 files, use batching (50 files per batch)
GATK LearnReadOrientationModel \
  -I batch1_file1.f1r2.tar.gz \
  -I batch1_file2.f1r2.tar.gz \
  ... (50 files)
  -O batch1_model.tar.gz

# Merge batch models
cp batch1_model.tar.gz final_orientation_model.tar.gz
```

**Why batching?** GATK has issues with >50 input files

**4b. Filter Each Sample:**

```bash
GATK FilterMutectCalls \
  -R reference.fa \
  -V sample.vcf.gz \
  --stats merged.stats \
  --ob-priors orientation_model.tar.gz \
  --contamination-table sample.contamination.table \
  --tumor-segmentation sample.seg \
  -O sample_filtered.vcf.gz
```

#### Key Parameters:
- `--stats`: Merged statistics from all samples (joint filtering)
- `--ob-priors`: Read orientation artifact model
- `--contamination-table`: Sample-specific contamination
- `--tumor-segmentation`: CNV segmentation

#### Memory Requirements:
- **Orientation model:** 16GB
- **Filtering (WGS):** 100GB (increased from default due to large VCF size)
- **Filtering (WES):** 30-50GB typically sufficient

#### Output:
```
Filtering samples sequentially with high memory allocation...
Filtering Sample1...
  ✓ Total=25000, PASS=150
Filtering Sample2...
  ✓ Total=750000, PASS=5500
Filtering Sample3...
  ✓ Total=700000, PASS=8500
```

#### Error Handling:
- Exit code checking
- OOM (Out Of Memory) detection
- Log tail capture for failures
- Automatic sample exclusion if filtering fails

---

### STEP 5: CREATE MULTI-SAMPLE VCF

**Lines:** 634-686  
**Purpose:** Extract PASS variants and create merged multi-sample VCF  
**Tools:** bcftools view, bcftools merge  

#### Process:

**5a. Extract PASS Variants:**
```bash
for each sample:
  bcftools view -f PASS filtered.vcf.gz -O z -o pass.vcf.gz
  tabix -p vcf pass.vcf.gz
```

**5b. Merge PASS VCFs:**
```bash
bcftools merge \
  sample1_pass.vcf.gz \
  sample2_pass.vcf.gz \
  sample3_pass.vcf.gz \
  --force-samples \
  -O z \
  -o PATIENT_longitudinal.vcf.gz

tabix -p vcf PATIENT_longitudinal.vcf.gz
```

#### Sample Exclusion:
Samples with 0 PASS variants are automatically excluded:
```
Sample1: 150 PASS variants
Sample2: 0 PASS (EXCLUDED)
Sample3: 5500 PASS variants
```

#### Output:
```
Multi-sample VCF: 14500 sites
(Sites = unique positions across all samples)
```

#### Quality Control:
- Minimum 2 samples with PASS variants required
- Pipeline exits if insufficient PASS samples
- All exclusions logged to QC log

---

### EVOLUTION ANALYSIS

**Lines:** 687-885  
**Purpose:** Analyze clonal evolution and temporal dynamics  
**Tools:** bcftools isec  

#### Components:

##### 1. Trunk Mutations (Shared by All)

**Command:**
```bash
bcftools isec \
  -n "=${#SAMPLES[@]}" \  # Exactly N samples
  sample1.vcf.gz \
  sample2.vcf.gz \
  sample3.vcf.gz \
  -w 1 \                  # Write first file
  -O v \
  -o trunk_mutations.vcf
```

**Logic:** `-n =N` means "present in exactly N files"

**Output:**
```
Trunk (all samples): 45 variants
```

**Interpretation:**
- **High count:** Samples are closely related (shared ancestry)
- **Low count:** Significant divergence or different origins
- **0:** No variants shared by all (highly diverged or technical issues)

##### 2. Private Mutations (Unique to Each Sample)

**Command:**
```bash
for i in "${!SAMPLES[@]}"; do
  bcftools isec \
    -n=1 \                    # Present in exactly 1 file
    sample1.vcf.gz \
    sample2.vcf.gz \
    sample3.vcf.gz \
    -w $((i+1)) \             # Write the i-th sample
    -O v \
    -o sample${i}_private.vcf
done
```

**Logic:** `-n=1` identifies variants present in exactly one file

**Output:**
```
Private Mutations:
  - Sample1: 12 variants
  - Sample2: 8 variants
  - Sample3: 200 variants
  - Sample4: 1500 variants
```

**Interpretation:**
- **High count:** Sample-specific evolution, subclonal expansion
- **Low count:** Sample closely related to others
- **Pattern matters:** Increasing private variants → ongoing evolution

##### 3. Temporal Dynamics (Gained/Lost)

**Gained Mutations:**
```bash
bcftools isec \
  -C \                      # Complement mode
  first_sample.vcf.gz \     # Reference sample
  last_sample.vcf.gz \      # Query sample
  -w 1 \                    # Write complement (must be 1)
  -O v \
  -o GAINED.vcf
```

**Logic:** Complement gives variants in last NOT in first = gained

**Lost Mutations:**
```bash
bcftools isec \
  -C \
  last_sample.vcf.gz \      # Reverse order
  first_sample.vcf.gz \
  -w 1 \
  -O v \
  -o LOST.vcf
```

**Logic:** Variants in first NOT in last = lost

**Output:**
```
Temporal Dynamics (Sample1 → Sample8):
  Gained mutations: 8200 variants
  Lost mutations: 18 variants
```

**Interpretation:**
- **Gained:** Newly acquired variants (active evolution)
- **Lost:** Subclonal shifts or technical artifacts
- **Ratio:** High gained/lost ratio = progressive evolution

##### 4. Pairwise Intersection Matrix

**Purpose:** Show variant overlap between all sample pairs

**Command:**
```bash
for i in samples:
  for j in samples:
    if i == j:
      # Diagonal: total variants
      count = total_variants(sample_i)
    else:
      # Off-diagonal: shared variants
      count = bcftools isec -n=2 sample_i sample_j -w 1 | count
```

**Output Example:**
```
Sample               Sample1    Sample2    Sample3    Sample4
Sample1              150        130        110        115
Sample2              130        5500       4300       4200
Sample3              110        4300       8500       7900
Sample4              115        4200       7900       6100
```

**Interpretation:**
- **Diagonal:** Total PASS variants per sample
- **Off-diagonal:** Shared variants between pairs
- **Symmetric matrix:** Overlap(i,j) = Overlap(j,i)
- **Pattern:** Cluster of high values = closely related samples

##### 5. Evolution Summary Statistics

**Output:**
```
Evolution Summary:
  Total unique variants: 14500
  Trunk (all samples): 45 (0.3%)
  Temporal changes: +8200 gained, -18 lost

Interpretation Guide:
  - High trunk %: Samples closely related (shared ancestry)
  - Low trunk %: Significant divergence or different origins
  - High gained: Active ongoing evolution
  - High lost: Possible subclonal shifts or technical artifacts
```

---

## OUTPUT FILES

### Directory Structure

```
PATIENT_longitudinal/
│
├── PATIENT_longitudinal.vcf.gz         # Multi-sample VCF (all PASS variants)
├── PATIENT_longitudinal.vcf.gz.tbi     # VCF index
├── PATIENT_longitudinal_stats.txt      # Main statistics report
├── PATIENT_QC.log                      # Quality control log
├── PATIENT_merged.stats                # Merged Mutect2 statistics
│
├── sample1_filtered.vcf.gz             # Per-sample filtered VCFs
├── sample2_filtered.vcf.gz
├── sample3_filtered.vcf.gz
│
├── sample1_pass.vcf.gz                 # Per-sample PASS-only VCFs
├── sample2_pass.vcf.gz
├── sample3_pass.vcf.gz
│
├── sample1_filter.log                  # Per-sample filter logs
├── sample2_filter.log
├── sample3_filter.log
│
└── evolution_analysis/
    ├── PATIENT_trunk_mutations.vcf                    # Trunk (all samples)
    ├── sample1_private.vcf                            # Private per sample
    ├── sample2_private.vcf
    ├── sample3_private.vcf
    ├── sample1_to_sample3_GAINED.vcf                  # Temporal gained
    ├── sample1_to_sample3_LOST.vcf                    # Temporal lost
    └── pairwise_intersection_matrix.txt               # Sample relationships
```

### Key Output Files

#### 1. Multi-sample VCF
**File:** `PATIENT_longitudinal.vcf.gz`  
**Format:** Compressed VCF (bgzip)  
**Content:** All PASS variants across all samples  
**Usage:** Primary output for downstream analysis  

**Features:**
- Contains genotype information for all samples
- Only PASS-filtered variants
- Suitable for: annotation, functional analysis, visualization

#### 2. Statistics Report
**File:** `PATIENT_longitudinal_stats.txt`  
**Format:** Plain text  
**Content:** Comprehensive analysis summary  

**Sections:**
- Input validation results
- Raw variant counts
- Filtering statistics
- Evolution analysis results
- Interpretation guidelines

**Example:**
```
================================================================================
CLONAL EVOLUTION ANALYSIS
================================================================================

Trunk Mutations (shared by all samples):
  Count: 45 variants
  Definition: Variants present in all 8 samples
  Interpretation: These represent early clonal events
  Output: PatientX_trunk_mutations.vcf

Private Mutations (unique to individual samples):
  Definition: Variants present in only one sample
  Interpretation: Sample-specific clonal evolution
  - Sample1: 12 variants
  - Sample2: 8 variants
  - Sample3: 200 variants
  ...
```

#### 3. QC Log
**File:** `PATIENT_QC.log`  
**Format:** Plain text  
**Content:** Detailed execution log  

**Contains:**
- Validation checks
- Command outputs
- Error messages
- Progress tracking
- Quality warnings

**Use for:**
- Troubleshooting failures
- Quality assessment
- Performance monitoring

#### 4. Evolution VCFs

| File | Description | Use Case |
|------|-------------|----------|
| `trunk_mutations.vcf` | Variants in ALL samples | Identify early clonal events |
| `*_private.vcf` | Variants unique to sample | Find sample-specific mutations |
| `*_GAINED.vcf` | New variants in later sample | Track acquired mutations |
| `*_LOST.vcf` | Variants lost over time | Identify subclonal shifts |

#### 5. Pairwise Matrix
**File:** `pairwise_intersection_matrix.txt`  
**Format:** Tab-delimited text  
**Content:** Variant overlap between all pairs  

**Usage:**
```R
# Load in R for visualization
matrix <- read.table("pairwise_intersection_matrix.txt", header=TRUE, row.names=1)
heatmap(as.matrix(matrix))
```

---

## INTERNAL FUNCTIONS

### Logging Functions

#### `log_qc()`
**Purpose:** Write to both console and QC log  
**Usage:** Important messages visible to user  
```bash
log_qc "Processing sample: $SAMPLE_NAME"
```

#### `log_qc_only()`
**Purpose:** Write only to QC log (suppress console)  
**Usage:** Detailed debug information  
```bash
log_qc_only "  Command: $FULL_COMMAND"
```

### File Descriptor Management
```bash
exec 3>"$QC_LOG"        # Open FD 3 for QC log
echo "text" >&3         # Write to FD 3 (QC log)
exec 3>&-               # Close FD 3 at end
```

---

## QUALITY CONTROL & VALIDATION

### Validation Checks

| Check | Threshold | Action |
|-------|-----------|--------|
| VCF exists | Required | Exclude if missing |
| Variant count > 0 | Required | Exclude if zero |
| Variant count < 100K | Warning | Flag but include |
| Min samples | ≥ 2 | Exit if insufficient |
| PASS variants > 0 | Per-sample | Exclude from evolution |

### Quality Flags

**✓ Valid:** Meets all requirements  
**⚠ Low Quality:** < 100,000 variants (included with warning)  
**✗ Failed:** VCF missing or zero variants (excluded)  

### Error Handling

#### Filtering Failures
```bash
if [ $FILTER_EXIT -ne 0 ]; then
    log_qc "  ✗ FilterMutectCalls FAILED (exit code $FILTER_EXIT)"
    if grep -q "OutOfMemoryError" "$LOG_FILE"; then
        log_qc "    CAUSE: Out of Memory - increase Java heap size"
    fi
    FILTER_FAILED+=("$label")
    continue
fi
```

#### Automatic Recovery
- Failed samples excluded from evolution analysis
- Pipeline continues with remaining samples
- All failures logged to QC log

---

## MEMORY & PERFORMANCE

### Memory Allocation by Step

| Step | Tool | Memory | Duration (WGS) |
|------|------|--------|----------------|
| Validation | bash | 1GB | 1 min |
| Merge Stats | GATK | 8GB | 1 min |
| Union VCF | bcftools | 4GB | 5-10 min |
| Joint Genotyping | GATK | 30GB | 1-3 hours per sample |
| Orientation Model | GATK | 16GB | 10-30 min |
| Filtering | GATK | **100GB** | 30-60 min per sample |
| Evolution | bcftools | 4GB | 5-10 min |

### Performance Optimization

#### Memory
```bash
# For WGS filtering (line 578)
--java-options "-Xms100G -Xmx100G -XX:+UseG1GC"
```

**Why 100GB?**
- WGS VCFs can be very large (700K+ variants)
- FilterMutectCalls loads entire VCF into memory
- G1GC provides better garbage collection for large heaps

#### Parallelization
- Filtering is **sequential** (one sample at a time)
- Can parallelize by running multiple patients simultaneously
- Joint genotyping could be parallelized (not implemented)

#### Batching
```bash
# f1r2 files batched at 50 files (line 517)
if [ $F1R2_COUNT -gt 50 ]; then
    # Process in batches of 50
fi
```

**Reason:** GATK has issues with >50 input files

---

## TROUBLESHOOTING

### Common Issues

#### Issue 1: Out of Memory During Filtering
**Error:** `OutOfMemoryError` or `oom_kill` in logs  
**Solution:**
```bash
# Increase memory (line 578)
--java-options "-Xms150G -Xmx150G"
```

Or request more memory in SLURM:
```bash
#SBATCH --mem=200G
```

#### Issue 2: No PASS Variants
**Symptom:** All samples excluded in Step 5  
**Causes:**
- Over-stringent filtering
- Very low quality samples
- Contamination issues

**Solutions:**
1. Check raw variant counts (Step 1)
2. Review filter logs
3. Adjust contamination thresholds
4. Check sequencing quality metrics

#### Issue 3: VCF Not Found
**Error:** `Cannot find Mut1.vcf in directory`  
**Causes:**
- Wrong directory path
- Different filename pattern
- Missing scatter_gather step

**Solutions:**
```bash
# Check directory contents
ls -la tumor_scatter_gather/

# Look for VCF files
find tumor_scatter_gather/ -name "*.vcf"

# Use --vcf flag to specify directly
--vcf sample1:path/to/sample.vcf:scatter_gather_dir
```

#### Issue 4: Duplicate Sample Names
**Symptom:** bcftools merge warnings  
**Solution:** Use `--force-samples` flag (already in script)

#### Issue 5: Low Trunk Mutations
**Not necessarily an error!**  
**Possible reasons:**
- Samples are from different patients (contamination/swap)
- Very divergent samples
- Different tumor origins (multi-focal)
- Very stringent filtering

**Actions:**
1. Check pairwise matrix for sample relationships
2. Review patient metadata
3. Check sample identity (SNP fingerprinting)
4. Consider less stringent filtering

#### Issue 6: BCFtools Version Issues
**Error:** Unexpected bcftools behavior  
**Solution:**
```bash
# Check version
bcftools --version

# Requires BCFtools 1.18+
# Update if necessary
```

---

## EXAMPLES

### Example 1: Standard Workflow (Skip Genotyping)

```bash
#!/bin/bash
#SBATCH --job-name=mutect2_long
#SBATCH --mem=120G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00

./mutect2_longitudinal.sh \
  --patient Patient01 \
  --skip-genotyping \
  --output-dir Patient01_evolution \
  /path/to/tumor1_scatter_gather \
  /path/to/tumor2_scatter_gather \
  /path/to/tumor3_scatter_gather
```

**Expected Runtime:** 3-5 hours for 3 WGS samples

### Example 2: With Joint Genotyping

```bash
./mutect2_longitudinal.sh \
  --patient Patient01 \
  --normal /path/to/normal.bam \
  --ref /path/to/hg38.fa \
  --output-dir Patient01_genotyped \
  tumor1_scatter_gather \
  tumor2_scatter_gather
```

**Expected Runtime:** 6-12 hours for 2 WGS samples

### Example 3: Time Series with Labels

```bash
./mutect2_longitudinal.sh \
  --patient Patient01 \
  --skip-genotyping \
  --label baseline:tumor_day0 \
  --label month3:tumor_day90 \
  --label month6:tumor_day180 \
  --label month12:tumor_day365 \
  tumor_day0 \
  tumor_day90 \
  tumor_day180 \
  tumor_day365
```

### Example 4: Multiple Metastases

```bash
./mutect2_longitudinal.sh \
  --patient Patient01 \
  --skip-genotyping \
  --label primary:tumor_primary \
  --label liver_met1:tumor_liver1 \
  --label liver_met2:tumor_liver2 \
  --label lung_met:tumor_lung \
  --label brain_met:tumor_brain \
  tumor_primary \
  tumor_liver1 \
  tumor_liver2 \
  tumor_lung \
  tumor_brain
```

### Example 5: Treatment Response Analysis

```bash
./mutect2_longitudinal.sh \
  --patient Patient01 \
  --skip-genotyping \
  --label pre_treatment:tumor_pre \
  --label on_treatment_week4:tumor_week4 \
  --label on_treatment_week8:tumor_week8 \
  --label post_treatment:tumor_post \
  tumor_pre \
  tumor_week4 \
  tumor_week8 \
  tumor_post
```

---

## TECHNICAL NOTES

### BCFtools isec Reference

#### `-n` Flag (Number of files)
```bash
-n=1        # Exactly 1 file
-n=2        # At least 2 files
-n=3        # At least 3 files
-n=N        # Exactly N files (use with =)
```

#### `-C` Flag (Complement)
```bash
# Complement: positions in file2 but NOT in file1
bcftools isec -C file1.vcf file2.vcf -w 1

# CRITICAL: -C only works with -w 1
# Why? Complement is defined relative to first file
```

#### `-w` Flag (Write file)
```bash
-w 1        # Write first input file
-w 2        # Write second input file
-w N        # Write Nth input file
```

### GATK FilterMutectCalls

#### Key Filters
- **weak_evidence:** Low evidence for variant
- **strand_bias:** Strand bias detected
- **orientation:** Read orientation artifact
- **contamination:** Sample contamination
- **germline:** Likely germline variant
- **panel_of_normals:** In panel of normals

#### Adjusting Stringency
Default is balanced. To adjust:
```bash
# More sensitive (keep more variants)
--f-score-beta 0.5

# More specific (stricter filtering)
--f-score-beta 2.0
```

### VCF Format Notes

#### Genotype Fields (FORMAT)
- **GT:** Genotype (0/0, 0/1, 1/1)
- **AD:** Allelic depths (ref,alt)
- **AF:** Allele frequency
- **DP:** Read depth
- **F1R2, F2R1:** Strand counts

#### Filter Field
- **PASS:** Passed all filters
- **filter_name:** Failed specific filter
- **.** : Not filtered

---

## BIOLOGICAL INTERPRETATION GUIDE

### Trunk Mutations
**Definition:** Variants present in all samples

**High Percentage (>10%):**
- Early clonal events
- Strong shared ancestry
- Samples from same tumor lineage

**Low Percentage (<1%):**
- Significant divergence
- Possible multi-focal origins
- Long evolutionary distance

### Private Mutations
**Definition:** Variants unique to one sample

**High Count:**
- Active ongoing evolution
- Subclonal expansion
- Sample-specific selective pressure

**Increasing Pattern:**
- Progressive evolution
- Accumulating mutations over time
- May indicate treatment resistance

### Gained/Lost Mutations

**High Gained Count:**
- Active tumor evolution
- New subclonal populations
- Possible resistance mechanisms

**High Lost Count:**
- Subclonal shifts
- Sampling artifacts
- Population bottlenecks

---

## LICENSE

This pipeline is provided as-is for research use. Please check individual tool licenses:
- GATK: BSD 3-Clause License
- BCFtools: MIT/GPL
- SAMtools: MIT/GPL

---

## SUPPORT & CONTRIBUTIONS

For issues, questions, or contributions:
- Check the troubleshooting section
- Review QC logs for detailed error messages
- Verify software versions meet requirements
- Ensure adequate system resources

---

## VERSION HISTORY

### v1.1 (Current)
- Enhanced statistics with biological interpretations
- Pairwise intersection matrix
- Improved error handling and logging
- Optimized memory management for WGS
- Added percentage calculations
- Comprehensive documentation

### v1.0
- Initial release
- Core pipeline functionality
- Basic evolution analysis
- QC logging

---

## APPENDIX: File Format Examples

### Statistics File Format

```
================================================================================
MUTECT2 LONGITUDINAL JOINT ANALYSIS STATISTICS
================================================================================
Patient ID: PatientX
Normal Sample: NormalX
Analysis Date: Mon Nov 2025
Output Directory: PatientX_longitudinal
================================================================================

Tumor samples (valid):
  - Sample1
  - Sample2
  - Sample3

================================================================================
STEP 1: INDIVIDUAL RAW VARIANT CALLS
================================================================================

Individual variant counts:
  Sample1: 25000 variants
  Sample2: 750000 variants
  Sample3: 700000 variants

... [additional sections]
```

### Pairwise Matrix Format

```
# Pairwise shared variant counts between samples
# Diagonal: Total variants per sample
# Off-diagonal: Shared variants between sample pairs
#
Sample                              Sample1        Sample2        Sample3
Sample1                             150            130            110
Sample2                             130            5500           4300
Sample3                             110            4300           8500
```

---

## REFERENCES

### Related Protocols
- GATK Best Practices: Somatic SNVs + Indels
- Longitudinal tumor analysis guidelines
- Clonal evolution analysis methods

---

**Document Version:** 1.1  
**Last Updated:** November 2025  
**Script Version:** 1.1  

---

END OF REFERENCE DOCUMENT
