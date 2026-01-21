================================================================================
MUTECT2 LONGITUDINAL JOINT ANALYSIS STATISTICS
================================================================================
Patient ID: P1
Normal Sample: N/A
Analysis Date: Mon Nov 10 02:36:54 PM CET 2025
Output Directory: P1_longitudinal
================================================================================


Tumor samples (valid):
  - CMD659A1268_24-257902_S1
  - CMD659A1269_24-259314_S1
  - CMD659A1270_24-261422_S1
  - CMD659A1271_24-261422_S2
  - CMD659A1272_24-262943_S1
  - CMD659A1273_24-262943_S2

================================================================================
STEP 1: INDIVIDUAL RAW VARIANT CALLS
================================================================================

STEP 1: Individual variant counts:
  CMD659A1268_24-257902_S1: 736350 variants
  CMD659A1269_24-259314_S1: 706050 variants
  CMD659A1270_24-261422_S1: 813399 variants
  CMD659A1271_24-261422_S2: 744362 variants
  CMD659A1272_24-262943_S1: 731282 variants
  CMD659A1273_24-262943_S2: 784218 variants

================================================================================
STEP 2: UNION OF VARIANT SITES
================================================================================

  Merged stats: P1_longitudinal/P1_merged.stats

  Union sites: 2937449 variants

================================================================================
STEP 3: JOINT GENOTYPING - SKIPPED
================================================================================

================================================================================
STEP 4: JOINT FILTERING
================================================================================

  CMD659A1268_24-257902_S1: Total=736350, PASS=5492
  CMD659A1269_24-259314_S1: Total=706050, PASS=8459
  CMD659A1270_24-261422_S1: Total=813399, PASS=6009
  CMD659A1271_24-261422_S2: Total=744362, PASS=2938
  CMD659A1272_24-262943_S1: Total=731282, PASS=5478
  CMD659A1273_24-262943_S2: Total=784218, PASS=8206

================================================================================
STEP 5: MULTI-SAMPLE VCF
================================================================================

  CMD659A1268_24-257902_S1: 5492 PASS variants
  CMD659A1269_24-259314_S1: 8459 PASS variants
  CMD659A1270_24-261422_S1: 6009 PASS variants
  CMD659A1271_24-261422_S2: 2938 PASS variants
  CMD659A1272_24-262943_S1: 5478 PASS variants
  CMD659A1273_24-262943_S2: 8206 PASS variants

  Multi-sample VCF: 14573 sites

================================================================================
CLONAL EVOLUTION ANALYSIS
================================================================================

Trunk Mutations (shared by all samples):
  Count: 719 variants
  Definition: Variants present in all 6 samples
  Interpretation: These represent early clonal events
  Output: P1_trunk_mutations.vcf

Private Mutations (unique to individual samples):
  Definition: Variants present in only one sample
  Interpretation: Sample-specific clonal evolution
  - CMD659A1268_24-257902_S1: 837 variants
  - CMD659A1269_24-259314_S1: 1511 variants
  - CMD659A1270_24-261422_S1: 1054 variants
  - CMD659A1271_24-261422_S2: 1080 variants
  - CMD659A1272_24-262943_S1: 1021 variants
  - CMD659A1273_24-262943_S2: 1153 variants

Pairwise Intersection Matrix:
  Output: pairwise_intersection_matrix.txt
  Shows shared variants between all sample pairs
  Diagonal = total variants per sample
  Off-diagonal = shared variants between pairs

Temporal Dynamics (CMD659A1268_24-257902_S1 → CMD659A1273_24-262943_S2):
  Gained mutations: 1210 variants
    Definition: Present in last but not in first sample
    Interpretation: Newly acquired during tumor evolution
    Output: CMD659A1268_24-257902_S1_to_CMD659A1273_24-262943_S2_GAINED.vcf
  
  Lost mutations: 3924 variants
    Definition: Present in first but not in last sample
    Interpretation: Lost during evolution (subclonal in first, or technical)
    Output: CMD659A1268_24-257902_S1_to_CMD659A1273_24-262943_S2_LOST.vcf

Evolution Summary:
  Total unique variants: 14573
  Trunk (all samples): 719 (4.9%)
  Temporal changes: +1210 gained, -3924 lost
  
Interpretation Guide:
  - High trunk %: Samples are closely related (shared ancestry)
  - Low trunk %: Samples diverged significantly or have different origins
  - High gained: Active ongoing evolution
  - High lost: Possible subclonal shifts or technical artifacts

================================================================================
FINAL SUMMARY
================================================================================
Analysis complete: Mon Nov 10 04:02:02 PM CET 2025
Final analyzed samples: 6

Key Outputs:
  - Multi-sample VCF: P1_longitudinal.vcf.gz
  - Trunk mutations: evolution_analysis/P1_trunk_mutations.vcf
  - Private mutations: evolution_analysis/[SAMPLE]_private.vcf
  - Temporal changes: evolution_analysis/[FIRST]_to_[LAST]_GAINED/LOST.vcf
  - Pairwise matrix: evolution_analysis/pairwise_intersection_matrix.txt
  - Statistics: P1_longitudinal_stats.txt
  - QC Log: P1_QC.log

================================================================================
