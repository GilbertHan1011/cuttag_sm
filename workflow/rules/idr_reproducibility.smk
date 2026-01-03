# ============================================================================
# IDR-based Reproducibility Analysis Workflow
# ============================================================================
# This module implements comprehensive IDR analysis for CUT&Tag replicates:
#   1. BAM to tagAlign conversion
#   2. Pseudo-replicate generation and peak calling
#   3. IDR analysis (self-consistency, pooled pseudo-replicates, true replicates)
#   4. Reproducibility QC and optimal peak selection
# ============================================================================

import os
import math

DATA_DIR = config["output_base_dir"].rstrip("/")


def get_all_true_replicate_pairs():
    """
    Generate all pairwise combinations of replicates for IDR analysis.
    Returns a list of tuples: (group, rep1, rep2) for all pairs.
    Used to expand outputs in rule all.
    """
    pairs = []
    for group in get_sample_reps():
        samples = get_reproducibility_sample(group)
        for i in range(len(samples)):
            for j in range(i+1, len(samples)):
                pairs.append((group, samples[i], samples[j]))
    return pairs


def is_broad_mark_for_sample(sample):
    """
    Determine if a sample uses broad peak calling based on BROAD_MARKS config.
    Returns True if the sample's mark is in BROAD_MARKS, False otherwise.
    """
    row = st[st['sample'] == sample]
    if row.empty:
        return False
    mark = str(row['mark'].iloc[0])
    mark_lower = mark.lower()
    broad_marks = set(m.lower() for m in config.get('BROAD_MARKS', []))
    return mark_lower in broad_marks


def get_broad_params_str():
    """
    Get broad peak parameters as a string for passing to the script.
    """
    broad_params = config.get('BROAD_PARAMS', {})
    return f"slocal:{broad_params.get('slocal', 10000)},max-gap:{broad_params.get('max-gap', 1000)},broad-cutoff:{broad_params.get('broad-cutoff', 0.01)}"


def is_broad_mark_for_group(group):
    """
    Determine if a replicate group uses broad peak calling based on BROAD_MARKS config.
    Checks the first sample in the group to determine the mark.
    """
    samples = get_reproducibility_sample(group)
    if not samples:
        return False
    # Check the first sample in the group
    return is_broad_mark_for_sample(samples[0])

# ============================================================================
# Configuration & Constants
# ============================================================================

# Genome references
GENOME_SIZE = config.get("reference", {}).get("genome_size_bp", "hs")
CHROM_SIZES = config.get("reference", {}).get("chrom_sizes", "")

# Replicate settings (with defaults)
REP_OPTS = config.get("pipeline_steps", {}).get("replicates", {})
PVAL_THRESH = REP_OPTS.get("pval_thresh", 1e-3)
SMOOTH_WIN = REP_OPTS.get("smooth_win", 150)
CAP_NUM_PEAK = REP_OPTS.get("cap_num_peak", 300000)
PSEUDO_RANDOM_SEED = REP_OPTS.get("pseudoreplication_random_seed", 0)

# IDR settings
IDR_THRESH = REP_OPTS.get("idr_thresh", 0.05)
IDR_RANK = REP_OPTS.get("idr_rank", "signal.value")

# Calculate shiftsize for MACS2 (negative half of smooth window)
SHIFT_SIZE = -1 * (int(SMOOTH_WIN) // 2)

# Calculate IDR rank column (narrowPeak format)
# Column 7 = signalValue, Column 8 = pValue, Column 9 = qValue
IDR_RANK_COL = {
    "signal.value": "7",
    "p.value": "8",
    "q.value": "9"
}.get(IDR_RANK, "7")

# Calculate negative log10 threshold for IDR filtering
NEG_LOG10_THRESH = f"{-math.log10(float(IDR_THRESH)):.6f}"

# Resource defaults
MEM_MB_DEFAULT = config.get("workflow_resources", {}).get("compute", {}).get("memory_gb_default", 32) * 1024
THREADS_DEFAULT = config.get("workflow_resources", {}).get("compute", {}).get("threads_default", 8)

# ============================================================================
# Rule 1: Convert BAM to tagAlign format
# ============================================================================

rule bam_to_tagalign:
    """
    Convert filtered BAM to tagAlign format to capture individual Tn5 insertion sites.
    Uses tagAlign conversion approach with proper paired-end handling and TN5 shift.
    """
    input:
        bam = os.path.join(DATA_DIR, "Important_processed", "Bam", "{sample}.sorted.markd.bam"),
        bai = os.path.join(DATA_DIR, "Important_processed", "Bam", "{sample}.sorted.markd.bam.bai"),
    output:
        bed = os.path.join(DATA_DIR, "middle_files", "bed", "{sample}.tagAlign.gz"),
    params:
        script = os.path.join(workflow.basedir, "src", "bam_to_tagalign.sh"),
        bed_dir = os.path.join(DATA_DIR, "middle_files", "bed"),
        is_paired = lambda w: st.loc[st['sample'] == w.sample, 'read_type'].iloc[0] == "paired-end" if 'read_type' in st.columns else True,
        disable_tn5_shift = REP_OPTS.get("disable_tn5_shift", False),
    resources:
        mem_mb = MEM_MB_DEFAULT,
        runtime = 300,
    threads: THREADS_DEFAULT
    conda:
        "../envs/pybedtools.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "bam_to_tagalign_{sample}.log")
    shell:
        """
        mkdir -p {params.bed_dir}
        {params.script} \\
            {input.bam} \\
            {output.bed} \\
            {wildcards.sample} \\
            {params.is_paired} \\
            {params.disable_tn5_shift} \\
            {threads} \\
            {params.bed_dir} \\
            {log}
        """

# ============================================================================
# PHASE 1: Individual Sample Pseudo-Replicate Analysis
# ============================================================================

# Rule 2: Split each individual sample's tagAlign into pseudo-replicates
rule split_pseudoreplicates_replicate:
    """
    Split each individual sample's tagAlign into two pseudo-replicates (pr1, pr2).
    For paired-end data, keeps pairs together during splitting.
    Uses external script split_reproducibility.sh.
    """
    input:
        tagalign = lambda w: os.path.join(DATA_DIR, "middle_files", "bed", f"{w.replicate}.tagAlign.gz"),
    output:
        pr1 = os.path.join(DATA_DIR, "middle_files", "replicates", "{replicate}", "{replicate}.pr1.tagAlign.gz"),
        pr2 = os.path.join(DATA_DIR, "middle_files", "replicates", "{replicate}", "{replicate}.pr2.tagAlign.gz"),
    params:
        script = os.path.join(workflow.basedir, "src", "split_reproducibility.sh"),
        out_dir = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", w.replicate),
        is_paired = lambda w: st.loc[st['sample'] == w.replicate, 'read_type'].iloc[0] == "paired-end" if 'read_type' in st.columns else True,
        random_seed = PSEUDO_RANDOM_SEED,
    resources:
        mem_mb = 4 * MEM_MB_DEFAULT,
        runtime = 80,
    threads: 1
    conda:
        "../envs/pybedtools.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "replicates", "split_pseudoreplicates_{replicate}.log")
    shell:
        """
        {params.script} \\
            {input.tagalign} \\
            {output.pr1} \\
            {output.pr2} \\
            {wildcards.replicate} \\
            {params.is_paired} \\
            {params.random_seed} \\
            {log}
        """

# Rule 3: Call peaks on individual sample pseudo-replicates (pr1, pr2)
rule call_peaks_pseudoreplicate:
    """
    Call peaks on individual sample pseudo-replicates using MACS2.
    Uses broad or narrow peak calling based on the sample's marker.
    Note: Outputs narrowPeak format for IDR compatibility, but uses broad parameters for broad marks.
    """
    input:
        tagalign = os.path.join(DATA_DIR, "middle_files", "replicates", "{replicate}", "{replicate}.{pr}.tagAlign.gz"),
    output:
        narrowpeak = os.path.join(DATA_DIR, "middle_files", "replicates", "{replicate}", "{replicate}.{pr}.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", w.replicate),
        prefix = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", w.replicate, f"{w.replicate}.{w.pr}"),
        genome_size = GENOME_SIZE,
        pval_thresh = PVAL_THRESH,
        smooth_win = SMOOTH_WIN,
        cap_num_peak = CAP_NUM_PEAK,
        chrom_sizes = CHROM_SIZES,
        shiftsize = SHIFT_SIZE,
        is_broad = lambda w: is_broad_mark_for_sample(w.replicate),
        broad_params = get_broad_params_str(),
    resources:
        mem_mb = MEM_MB_DEFAULT,
        runtime = 600,
    threads: THREADS_DEFAULT
    wildcard_constraints:
        pr = "pr1|pr2"
    conda:
        "../envs/reproducibility.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "replicates", "call_peaks_{replicate}_{pr}.log")
    script:
        "../src/call_peaks_macs2.sh"

# Rule 4: Call peaks on full individual sample (for pooled peak reference)
rule call_peaks_individual_sample:
    """
    Call peaks on full individual sample tagAlign (for use as pooled peak reference).
    Uses broad or narrow peak calling based on the sample's marker.
    Note: Outputs narrowPeak format for IDR compatibility, but uses broad parameters for broad marks.
    """
    input:
        tagalign = lambda w: os.path.join(DATA_DIR, "middle_files", "bed", f"{w.sample}.tagAlign.gz"),
    output:
        narrowpeak = os.path.join(DATA_DIR, "middle_files", "replicates", "{sample}", "{sample}.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", w.sample),
        prefix = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", w.sample, f"{w.sample}"),
        genome_size = GENOME_SIZE,
        pval_thresh = PVAL_THRESH,
        smooth_win = SMOOTH_WIN,
        cap_num_peak = CAP_NUM_PEAK,
        chrom_sizes = CHROM_SIZES,
        shiftsize = SHIFT_SIZE,
        is_broad = lambda w: is_broad_mark_for_sample(w.sample),
        broad_params = get_broad_params_str(),
    resources:
        mem_mb = MEM_MB_DEFAULT,
        runtime = 80,
    threads: THREADS_DEFAULT
    conda:
        "../envs/reproducibility.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "replicates", "call_peaks_individual_{sample}.log")
    script:
        "../src/call_peaks_macs2.sh"

# Rule 5: Pool peaks from individual samples in a group for IDR reference
rule pool_sample_peaks_for_group:
    """
    Pool peaks from individual samples in a group for IDR reference peak list.
    """
    input:
        peaks = lambda w: [
            os.path.join(DATA_DIR, "middle_files", "replicates", sample, f"{sample}.narrowPeak.gz")
            for sample in get_reproducibility_sample(w.group)
        ],
    output:
        pooled_peaks = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-samples.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", "groups", w.group),
    resources:
        mem_mb = MEM_MB_DEFAULT // 2,
        runtime = 60,
    threads: 1
    conda:
        "../envs/reproducibility.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "replicates", "pool_sample_peaks_{group}.log")
    shell:
        """
        mkdir -p {params.out_dir}
        # Concatenate all peaks, sort, and merge overlapping peaks
        zcat {input.peaks} | sort -k1,1 -k2,2n | gzip -nc > {output.pooled_peaks}
        """

# Rule 6: IDR between pr1 and pr2 for each individual sample (self-consistency)
rule idr_self_consistency:
    """
    Calculate IDR between pr1 and pr2 for each individual sample (self-consistency check).
    """
    input:
        pr1_peaks = os.path.join(DATA_DIR, "middle_files", "replicates", "{sample}", "{sample}.pr1.narrowPeak.gz"),
        pr2_peaks = os.path.join(DATA_DIR, "middle_files", "replicates", "{sample}", "{sample}.pr2.narrowPeak.gz"),
        pooled_peaks = os.path.join(DATA_DIR, "middle_files", "replicates", "{sample}", "{sample}.narrowPeak.gz"),
    output:
        idr_peak = os.path.join(DATA_DIR, "middle_files", "replicates", "{sample}", "{sample}-pr1_vs_{sample}-pr2.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", w.sample),
        prefix = lambda w: f"{w.sample}-pr1_vs_{w.sample}-pr2",
        idr_thresh = IDR_THRESH,
        idr_rank = IDR_RANK,
        chrom_sizes = CHROM_SIZES,
        idr_rank_col = IDR_RANK_COL,
        neg_log10_thresh = NEG_LOG10_THRESH,
        method = REP_OPTS.get("peak_combination_method", "idr"),
        nonamecheck = "-nonamecheck" if REP_OPTS.get("nonamecheck", False) else "",
    resources:
        mem_mb = 2 * MEM_MB_DEFAULT,
        runtime = 80,
    threads: THREADS_DEFAULT
    conda:
        "../envs/reproducibility.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "replicates", "idr_self_consistency_{sample}.log")
    script:
        "../src/idr_analysis.sh"

# ============================================================================
# PHASE 2: Pooled Pseudo-Replicate Analysis (for replicate groups with >1 sample)
# ============================================================================

# Rule 7: Pool all pr1 tagAligns across all samples in a group
rule pool_group_pr1_tagaligns:
    """
    Pool all pr1 tagAligns across all samples in a replicate group.
    """
    input:
        pr1_tagaligns = lambda w: get_group_pr1_tagaligns(w),
    output:
        pooled_pr1 = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr1.tagAlign.gz"),
    params:
        out_dir = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", "groups", w.group),
    resources:
        mem_mb = MEM_MB_DEFAULT,
        runtime = 120,
    threads: 1
    conda:
        "../envs/pybedtools.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "replicates", "pool_group_pr1_{group}.log")
    shell:
        """
        mkdir -p {params.out_dir}
        zcat {input.pr1_tagaligns} | gzip -nc > {output.pooled_pr1}
        """

# Rule 8: Pool all pr2 tagAligns across all samples in a group
rule pool_group_pr2_tagaligns:
    """
    Pool all pr2 tagAligns across all samples in a replicate group.
    """
    input:
        pr2_tagaligns = lambda w: get_group_pr2_tagaligns(w),
    output:
        pooled_pr2 = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr2.tagAlign.gz"),
    params:
        out_dir = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", "groups", w.group),
    resources:
        mem_mb = MEM_MB_DEFAULT,
        runtime = 120,
    threads: 1
    conda:
        "../envs/pybedtools.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "replicates", "pool_group_pr2_{group}.log")
    shell:
        """
        mkdir -p {params.out_dir}
        zcat {input.pr2_tagaligns} | gzip -nc > {output.pooled_pr2}
        """

# Rule 9: Pool all replicate tagAligns in a group (for peak calling on all pooled data)
rule pool_group_all_tagaligns:
    """
    Pool all replicate tagAligns in a group for peak calling on all pooled data.
    """
    input:
        tagaligns = lambda w: get_replicate_tagaligns(w),
    output:
        pooled_all = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-all.tagAlign.gz"),
    params:
        out_dir = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", "groups", w.group),
    resources:
        mem_mb = MEM_MB_DEFAULT,
        runtime = 120,
    threads: 1
    conda:
        "../envs/pybedtools.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "replicates", "pool_group_all_{group}.log")
    shell:
        """
        mkdir -p {params.out_dir}
        zcat {input.tagaligns} | gzip -nc > {output.pooled_all}
        """

# Rule 10: Call peaks on pooled-pr1
rule call_peaks_pooled_pr1:
    """
    Call peaks on pooled-pr1 tagAlign using MACS2.
    Uses broad or narrow peak calling based on the group's marker.
    Note: Outputs narrowPeak format for IDR compatibility, but uses broad parameters for broad marks.
    """
    input:
        tagalign = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr1.tagAlign.gz"),
    output:
        narrowpeak = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr1.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", "groups", w.group),
        prefix = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", "groups", w.group, f"{w.group}.pooled-pr1"),
        genome_size = GENOME_SIZE,
        pval_thresh = PVAL_THRESH,
        smooth_win = SMOOTH_WIN,
        cap_num_peak = CAP_NUM_PEAK,
        chrom_sizes = CHROM_SIZES,
        shiftsize = SHIFT_SIZE,
        is_broad = lambda w: is_broad_mark_for_group(w.group),
        broad_params = get_broad_params_str(),
    resources:
        mem_mb = MEM_MB_DEFAULT,
        runtime = 80,
    threads: THREADS_DEFAULT
    conda:
        "../envs/reproducibility.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "replicates", "call_peaks_pooled_pr1_{group}.log")
    script:
        "../src/call_peaks_macs2.sh"

# Rule 11: Call peaks on pooled-pr2
rule call_peaks_pooled_pr2:
    """
    Call peaks on pooled-pr2 tagAlign using MACS2.
    Uses broad or narrow peak calling based on the group's marker.
    Note: Outputs narrowPeak format for IDR compatibility, but uses broad parameters for broad marks.
    """
    input:
        tagalign = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr2.tagAlign.gz"),
    output:
        narrowpeak = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr2.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", "groups", w.group),
        prefix = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", "groups", w.group, f"{w.group}.pooled-pr2"),
        genome_size = GENOME_SIZE,
        pval_thresh = PVAL_THRESH,
        smooth_win = SMOOTH_WIN,
        cap_num_peak = CAP_NUM_PEAK,
        chrom_sizes = CHROM_SIZES,
        shiftsize = SHIFT_SIZE,
        is_broad = lambda w: is_broad_mark_for_group(w.group),
        broad_params = get_broad_params_str(),
    resources:
        mem_mb = MEM_MB_DEFAULT,
        runtime = 80,
    threads: THREADS_DEFAULT
    conda:
        "../envs/reproducibility.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "replicates", "call_peaks_pooled_pr2_{group}.log")
    script:
        "../src/call_peaks_macs2.sh"

# Rule 12: Call peaks on all pooled tagAligns
rule call_peaks_pooled_all:
    """
    Call peaks on all pooled tagAligns using MACS2.
    Uses broad or narrow peak calling based on the group's marker.
    Note: Outputs narrowPeak format for IDR compatibility, but uses broad parameters for broad marks.
    """
    input:
        tagalign = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-all.tagAlign.gz"),
    output:
        narrowpeak = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-all.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", "groups", w.group),
        prefix = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", "groups", w.group, f"{w.group}.pooled-all"),
        genome_size = GENOME_SIZE,
        pval_thresh = PVAL_THRESH,
        smooth_win = SMOOTH_WIN,
        cap_num_peak = CAP_NUM_PEAK,
        chrom_sizes = CHROM_SIZES,
        shiftsize = SHIFT_SIZE,
        is_broad = lambda w: is_broad_mark_for_group(w.group),
        broad_params = get_broad_params_str(),
    resources:
        mem_mb = 4 * MEM_MB_DEFAULT,
        runtime = 80,
    threads: THREADS_DEFAULT
    conda:
        "../envs/reproducibility.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "replicates", "call_peaks_pooled_all_{group}.log")
    script:
        "../src/call_peaks_macs2.sh"

# Rule 13: IDR between pooled-pr1 vs pooled-pr2
rule idr_pooled_pseudoreplicates:
    """
    Calculate IDR between pooled-pr1 vs pooled-pr2 pseudo-replicates.
    """
    input:
        pr1_peaks = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr1.narrowPeak.gz"),
        pr2_peaks = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr2.narrowPeak.gz"),
        pooled_peaks = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-all.narrowPeak.gz"),
    output:
        idr_peak = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "{group}.pooled-pr1_vs_pooled-pr2.idr.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", "groups", w.group),
        prefix = lambda w: f"{w.group}.pooled-pr1_vs_pooled-pr2",
        idr_thresh = IDR_THRESH,
        idr_rank = IDR_RANK,
        chrom_sizes = CHROM_SIZES,
        idr_rank_col = IDR_RANK_COL,
        neg_log10_thresh = NEG_LOG10_THRESH,
        method = REP_OPTS.get("peak_combination_method", "idr"),
        nonamecheck = "-nonamecheck" if REP_OPTS.get("nonamecheck", False) else "",
    resources:
        mem_mb = 2 * MEM_MB_DEFAULT,
        runtime = 80,
    threads: THREADS_DEFAULT
    conda:
        "../envs/reproducibility.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "replicates", "idr_pooled_pr_{group}.log")
    script:
        "../src/idr_analysis.sh"

# ============================================================================
# PHASE 3: True Replicate Comparisons (pairwise between biological replicates)
# ============================================================================

# Rule 14: IDR between pairs of true replicates
rule idr_true_replicates:
    """
    Calculate IDR between pairs of true biological replicates.
    """
    input:
        pr1_peaks = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", w.rep1, f"{w.rep1}.narrowPeak.gz"),
        pr2_peaks = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", w.rep2, f"{w.rep2}.narrowPeak.gz"),
        pooled_peaks = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", "groups", w.group, f"{w.group}.pooled-all.narrowPeak.gz"),
    output:
        idr_peak = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "true_replicates", "{rep1}_vs_{rep2}.idr.narrowPeak.gz"),
    params:
        out_dir = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", "groups", w.group, "true_replicates"),
        prefix = lambda w: f"{w.rep1}_vs_{w.rep2}",
        idr_thresh = IDR_THRESH,
        idr_rank = IDR_RANK,
        chrom_sizes = CHROM_SIZES,
        idr_rank_col = IDR_RANK_COL,
        neg_log10_thresh = NEG_LOG10_THRESH,
        method = REP_OPTS.get("peak_combination_method", "idr"),
        nonamecheck = "-nonamecheck" if REP_OPTS.get("nonamecheck", False) else "",
    resources:
        mem_mb = 2 * MEM_MB_DEFAULT,
        runtime = 80,
    threads: THREADS_DEFAULT
    conda:
        "../envs/reproducibility.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "replicates", "idr_true_replicates_{group}_{rep1}_vs_{rep2}.log")
    script:
        "../src/idr_analysis.sh"

# ============================================================================
# PHASE 4: Correlation Analysis (Replicate Correlation Metrics)
# ============================================================================

# Rule 15: Merge peaks from all replicates in a group for correlation analysis
rule merge_group_peaks_for_correlation:
    """
    Merge peaks from all replicates in a group to create a unified peak set
    for correlation analysis. This ensures we compare read counts in the same
    genomic regions across replicates.
    """
    input:
        peaks = lambda w: [
            os.path.join(DATA_DIR, "middle_files", "replicates", sample, f"{sample}.narrowPeak.gz")
            for sample in get_reproducibility_sample(w.group)
        ],
    output:
        merged_peaks = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "{group}.merged_peaks.bed"),
    params:
        out_dir = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", "groups", w.group),
    resources:
        mem_mb = MEM_MB_DEFAULT // 2,
        runtime = 60,
    threads: 1
    conda:
        "../envs/reproducibility.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "replicates", "merge_group_peaks_correlation_{group}.log")
    shell:
        """
        mkdir -p {params.out_dir}
        # Extract peaks from narrowPeak files (first 3 columns: chrom, start, end)
        # Sort and merge overlapping peaks
        zcat {input.peaks} | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i stdin > {output.merged_peaks}
        """

# Rule 16: Calculate correlation between replicate pairs
rule calculate_replicate_correlation:
    """
    Calculate Pearson and Spearman correlation between biological replicate pairs
    based on read counts in merged peaks. Uses CPM normalization and log2 transform.
    """
    input:
        merged_peaks = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "{group}.merged_peaks.bed"),
        rep1_tagalign = lambda w: os.path.join(DATA_DIR, "middle_files", "bed", f"{w.rep1}.tagAlign.gz"),
        rep2_tagalign = lambda w: os.path.join(DATA_DIR, "middle_files", "bed", f"{w.rep2}.tagAlign.gz"),
    output:
        counts = temp(os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "true_replicates", "{rep1}_vs_{rep2}.counts.tsv")),
        tsv = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "true_replicates", "{rep1}_vs_{rep2}.correlation.tsv"),
        png = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "true_replicates", "{rep1}_vs_{rep2}.correlation.png"),
    params:
        out_dir = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", "groups", w.group, "true_replicates"),
        script = os.path.join(workflow.basedir, "src", "calculate_replicate_correlation.py"),
    resources:
        mem_mb = MEM_MB_DEFAULT,
        runtime = 120,
    threads: THREADS_DEFAULT
    conda:
        "../envs/reproducibility.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "replicates", "calculate_correlation_{group}_{rep1}_vs_{rep2}.log")
    shell:
        """
        mkdir -p {params.out_dir}
        
        # Use bedtools coverage to count reads in merged peaks for each replicate
        # Add a name column (chrom:start-end) for joining
        bedtools coverage -counts -a {input.merged_peaks} -b {input.rep1_tagalign} | \\
            awk '{{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$4}}' | \\
            sort -k4,4 > {output.counts}.rep1
        
        bedtools coverage -counts -a {input.merged_peaks} -b {input.rep2_tagalign} | \\
            awk '{{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$4}}' | \\
            sort -k4,4 > {output.counts}.rep2
        
        # Join the two count files on the name column (column 4)
        # Output: chrom, start, end, name, rep1_count, rep2_count
        join -t $'\t' -1 4 -2 4 -o 1.1,1.2,1.3,1.4,1.5,2.5 \\
            {output.counts}.rep1 \\
            {output.counts}.rep2 > {output.counts}
        
        # Calculate correlation using Python script
        python {params.script} \\
            --input {output.counts} \\
            --out-tsv {output.tsv} \\
            --out-png {output.png} 2>> {log}
        """

# Rule 17: Aggregate correlation results for all pairs in a group
rule aggregate_replicate_correlations:
    """
    Aggregate correlation statistics from all replicate pairs in a group
    into a single summary file.
    """
    input:
        correlation_files = lambda w: get_correlation_pair_files(w),
    output:
        summary_tsv = os.path.join(DATA_DIR, "Report", "reproducibility", "{group}.replicate_correlations.tsv"),
    params:
        pairs = lambda w: [
            (rep1, rep2, os.path.join(
                DATA_DIR, "middle_files", "replicates", "groups", 
                w.group, "true_replicates", 
                f"{rep1}_vs_{rep2}.correlation.tsv"
            ))
            for rep1, rep2 in get_unique_replicate_pairs(w.group)
        ],
    resources:
        mem_mb = MEM_MB_DEFAULT // 4,
        runtime = 30,
    threads: 1
    conda:
        "../envs/reproducibility.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "replicates", "aggregate_correlations_{group}.log")
    script:
        "../src/aggregate_correlations.py"

# ============================================================================
# PHASE 5: Reproducibility QC and Optimal Peak Selection
# ============================================================================

# Rule 18: Reproducibility QC - Evaluate consistency and select optimal peaks
rule reproducibility_qc:
    """
    Evaluate reproducibility consistency and select optimal peaks based on IDR results.
    """
    input:
        peaks_pr = lambda w: [
            os.path.join(DATA_DIR, "middle_files", "replicates", rep, f"{rep}-pr1_vs_{rep}-pr2.narrowPeak.gz")
            for rep in get_reproducibility_sample(w.group)
        ],
        peak_ppr = lambda w: os.path.join(DATA_DIR, "middle_files", "replicates", "groups", w.group, f"{w.group}.pooled-pr1_vs_pooled-pr2.idr.narrowPeak.gz"),
        # Ensure all true replicate pairs exist
        true_rep_pairs = lambda w: get_true_replicate_pair_files(w),
    output:
        qc_json = os.path.join(DATA_DIR, "Report", "reproducibility", "{group}.reproducibility.qc.json"),
        optimal_peak = os.path.join(DATA_DIR, "middle_files", "replicates", "groups", "{group}", "{group}.optimal_peak.narrowPeak.gz"),
    params:
        prefix = lambda w: w.group,
    resources:
        mem_mb = 4000,
        runtime = 60,
    threads: 1
    conda:
        "../envs/reproducibility.yml"
    log:
        os.path.join(DATA_DIR, "logs", "rules", "replicates", "reproducibility_qc_{group}.log")
    script:
        "../src/reproducibility_qc.py"

