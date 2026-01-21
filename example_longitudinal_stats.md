# MUTECT2 LONGITUDINAL JOINT ANALYSIS STATISTICS

**Patient ID:** P1  
**Normal Sample:** N/A  
**Analysis Date:** Mon Nov 10 02:36:54 PM CET 2025  
**Output Directory:** P1_longitudinal

---

## Tumor Samples (Valid)

- CMD659A1268_24-257902_S1
- CMD659A1269_24-259314_S1
- CMD659A1270_24-261422_S1
- CMD659A1271_24-261422_S2
- CMD659A1272_24-262943_S1
- CMD659A1273_24-262943_S2

---

## STEP 1: Individual Raw Variant Calls

| Sample | Variant Count |
|--------|---------------|
| CMD659A1268_24-257902_S1 | 736,350 |
| CMD659A1269_24-259314_S1 | 706,050 |
| CMD659A1270_24-261422_S1 | 813,399 |
| CMD659A1271_24-261422_S2 | 744,362 |
| CMD659A1272_24-262943_S1 | 731,282 |
| CMD659A1273_24-262943_S2 | 784,218 |

---

## STEP 2: Union of Variant Sites

**Merged stats:** P1_longitudinal/P1_merged.stats  
**Union sites:** 2,937,449 variants

---

## STEP 3: Joint Genotyping

**Status:** SKIPPED

---

## STEP 4: Joint Filtering

| Sample | Total | PASS |
|--------|-------|------|
| CMD659A1268_24-257902_S1 | 736,350 | 5,492 |
| CMD659A1269_24-259314_S1 | 706,050 | 8,459 |
| CMD659A1270_24-261422_S1 | 813,399 | 6,009 |
| CMD659A1271_24-261422_S2 | 744,362 | 2,938 |
| CMD659A1272_24-262943_S1 | 731,282 | 5,478 |
| CMD659A1273_24-262943_S2 | 784,218 | 8,206 |

---

## STEP 5: Multi-Sample VCF

| Sample | PASS Variants |
|--------|---------------|
| CMD659A1268_24-257902_S1 | 5,492 |
| CMD659A1269_24-259314_S1 | 8,459 |
| CMD659A1270_24-261422_S1 | 6,009 |
| CMD659A1271_24-261422_S2 | 2,938 |
| CMD659A1272_24-262943_S1 | 5,478 |
| CMD659A1273_24-262943_S2 | 8,206 |

**Multi-sample VCF:** 14,573 sites

---

## Clonal Evolution Analysis

### Trunk Mutations (Shared by All Samples)

- **Count:** 719 variants
- **Definition:** Variants present in all 6 samples
- **Interpretation:** These represent early clonal events
- **Output:** P1_trunk_mutations.vcf

### Private Mutations (Unique to Individual Samples)

**Definition:** Variants present in only one sample  
**Interpretation:** Sample-specific clonal evolution

| Sample | Private Variants |
|--------|------------------|
| CMD659A1268_24-257902_S1 | 837 |
| CMD659A1269_24-259314_S1 | 1,511 |
| CMD659A1270_24-261422_S1 | 1,054 |
| CMD659A1271_24-261422_S2 | 1,080 |
| CMD659A1272_24-262943_S1 | 1,021 |
| CMD659A1273_24-262943_S2 | 1,153 |

### Pairwise Intersection Matrix

- **Output:** pairwise_intersection_matrix.txt
- Shows shared variants between all sample pairs
- Diagonal = total variants per sample
- Off-diagonal = shared variants between pairs

### Temporal Dynamics

**From:** CMD659A1268_24-257902_S1 → CMD659A1273_24-262943_S2

#### Gained Mutations: 1,210 variants
- **Definition:** Present in last but not in first sample
- **Interpretation:** Newly acquired during tumor evolution
- **Output:** CMD659A1268_24-257902_S1_to_CMD659A1273_24-262943_S2_GAINED.vcf

#### Lost Mutations: 3,924 variants
- **Definition:** Present in first but not in last sample
- **Interpretation:** Lost during evolution (subclonal in first, or technical)
- **Output:** CMD659A1268_24-257902_S1_to_CMD659A1273_24-262943_S2_LOST.vcf

### Evolution Summary

- **Total unique variants:** 14,573
- **Trunk (all samples):** 719 (4.9%)
- **Temporal changes:** +1,210 gained, -3,924 lost

### Interpretation Guide

- **High trunk %:** Samples are closely related (shared ancestry)
- **Low trunk %:** Samples diverged significantly or have different origins
- **High gained:** Active ongoing evolution
- **High lost:** Possible subclonal shifts or technical artifacts

---

## Final Summary

**Analysis complete:** Mon Nov 10 04:02:02 PM CET 2025  
**Final analyzed samples:** 6

### Key Outputs

- **Multi-sample VCF:** P1_longitudinal.vcf.gz
- **Trunk mutations:** evolution_analysis/P1_trunk_mutations.vcf
- **Private mutations:** evolution_analysis/[SAMPLE]_private.vcf
- **Temporal changes:** evolution_analysis/[FIRST]_to_[LAST]_GAINED/LOST.vcf
- **Pairwise matrix:** evolution_analysis/pairwise_intersection_matrix.txt
- **Statistics:** P1_longitudinal_stats.txt
- **QC Log:** P1_QC.log
