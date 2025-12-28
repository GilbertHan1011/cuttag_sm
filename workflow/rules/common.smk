# ============================================================================
# Common helper functions for Snakemake workflow
# ============================================================================
# Functions are organized by category:
#   1. Core utilities (data directory, config helpers)
#   2. Sample table queries (samples, marks, runs)
#   3. Input file path builders (FASTQ, BAM, peaks)
#   4. Peak calling helpers (MACS2, broad/narrow detection)
#   5. Reproducibility helpers
#   6. Internal helpers (prefixed with _)
# ============================================================================
import pandas as pd
import os

# ============================================================================
# 1. CORE UTILITIES
# ============================================================================

def get_data_dir():
    """Return the base output directory."""
    return config["output_base_dir"].rstrip("/")


def _get_broad_marks_set():
    """
    Helper to get broad marks set from config (cached for reuse).
    Returns a set of lowercase mark names.
    """
    return set(m.lower() for m in config.get('BROAD_MARKS', []))


# ============================================================================
# 2. SAMPLE TABLE QUERIES
# ============================================================================

def get_samples():
    """Return sorted list of unique sample identifiers."""
    return sorted(list(pd.Index(st.index).unique()))


def get_marks():
    """Return sorted list of unique marks from sample table."""
    return sorted(list(st['mark'].astype(str).unique()))


def get_mark_conditions():
    """Return list of unique mark_condition combinations."""
    st['mark_condition'] = st['mark'].astype(str) + "_" + st['condition']
    return st['mark_condition'].unique().tolist()


def get_runs_for_sample(wildcards):
    """
    Return list of run identifiers for a biological sample.
    Accepts wildcards object with 'sample' attribute or sample string.
    """
    sample = wildcards.sample if hasattr(wildcards, 'sample') else wildcards
    runs = st[st['sample'] == sample]['run'].astype(str).tolist()
    return runs


def get_sample_reps():
    """
    Return list of sample_base values that have more than 1 sample (replicates).
    These are the sample_rep values used in reproducibility rules.
    """
    if "sample_base" not in st.columns:
        return []
    # Get sample_base values that have more than 1 sample
    sample_base_counts = st.groupby("sample_base")["sample"].nunique()
    sample_reps = sample_base_counts[sample_base_counts > 1].index.tolist()
    return sorted(sample_reps)


def get_reads():
    """Return list of all read identifiers (sample.run_R1, sample.run_R2)."""
    reads = []
    for _, row in st.iterrows():
        sample = row['sample']
        run = str(row['run']) if 'run' in row else '1'
        reads.append(f"{sample}.{run}_R1")
        reads.append(f"{sample}.{run}_R2")
    return reads


# ============================================================================
# 3. INPUT FILE PATH BUILDERS
# ============================================================================

def get_bowtie2_input(wildcards):
    """
    Return R1 and R2 FASTQ paths for a sample (first run if multiple).
    """
    try:
        r1 = st.loc[wildcards.sample]['R1']
        r2 = st.loc[wildcards.sample]['R2']
        return r1, r2
    except Exception:
        row = st[st['sample'] == wildcards.sample].iloc[0]
        return row['R1'], row['R2']


def get_bowtie2_input_by_run(wildcards):
    """
    Return R1 and R2 FASTQ paths for a specific sample/run combination.
    Requires 'sample' and 'run' wildcards.
    """
    row = st[(st['sample'] == wildcards.sample) & (st['run'] == wildcards.run)]
    if row.empty:
        raise ValueError(f"No samplesheet row for sample={wildcards.sample}, run={wildcards.run}")
    row = row.iloc[0]
    return row['R1'], row['R2']


def _all_trimmed_fastqs(sample):
    """
    Return lists of all trimmed FASTQ files (R1 and R2) for a sample across all runs.
    Used internally by alignment rules.
    Format matches preprocess.smk output: {sample}.{run}_R1.fastq.gz
    """
    runs = get_runs_for_sample(sample)
    data_dir = get_data_dir()
    r1 = [f"{data_dir}/middle_file/Trimmed_fastq/fastp/{sample}.{r}_R1.fastq.gz" for r in runs]
    r2 = [f"{data_dir}/middle_file/Trimmed_fastq/fastp/{sample}.{r}_R2.fastq.gz" for r in runs]
    return r1, r2


def _all_raw_fastqs(sample):
    """
    Return lists of all raw FASTQ files (R1 and R2) for a sample across all runs.
    Used internally by alignment rules.
    """
    rows = st[st['sample'] == sample]
    r1 = rows["R1"].astype(str).tolist()
    r2 = rows["R2"].astype(str).tolist()
    return r1, r2


def get_sorted_bams_for_sample(wildcards):
    """Build paths to per-run sorted BAMs for a sample."""
    runs = get_runs_for_sample(wildcards)
    data_dir = get_data_dir()
    return [
        f"{data_dir}/middle_file/aligned/{wildcards.sample}.{r}.sort.bam"
        for r in runs
    ]


def get_callpeaks(wildcards):
    """Return BAM file path for peak calling."""
    data_dir = get_data_dir()
    bam = f"{data_dir}/Important_processed/Bam/{wildcards.sample}.sorted.markd.bam"
    return [bam]


def get_igg(wildcards):
    """
    Return IgG control BAM path as MACS2 -c parameter string.
    Returns empty string if USEIGG is False or sample is IgG itself.
    """
    if not config.get('USEIGG', False):
        return ""
    
    row = st[st["sample"] == wildcards.sample]
    if row.empty:
        return ""
    
    igg = str(row["igg"].iloc[0])
    if not igg or config["IGG"] in wildcards.sample:
        return ""
    
    data_dir = get_data_dir()
    iggbam = f"{data_dir}/Important_processed/Bam/{igg}.sorted.markd.bam"
    return f'-c {iggbam}'


def get_tracks_by_mark_condition(wildcards):
    """Return list of BigWig track paths for a mark_condition."""
    st['mark_condition'] = st['mark'].astype(str) + "_" + st['condition']
    samps = st.groupby(["mark_condition"])["sample"].apply(list)[wildcards.mark_condition]
    data_dir = get_data_dir()
    return [f"{data_dir}/Important_processed/Track/tracks/{s}.bw" for s in samps]


def _get_peaks_by_mark_condition_base(wildcards, use_blacklist=False):
    """
    Base function to get peak files by mark_condition.
    Used by get_peaks_by_mark_condition and get_peaks_by_mark_condition_blacklist.
    """
    st["mark_condition"] = st["mark"].astype(str) + "_" + st["condition"]
    samps = st.groupby(["mark_condition"])["sample"].apply(list)[wildcards.mark_condition]
    data_dir = get_data_dir()
    suffix = "_peaks_noBlacklist.bed" if use_blacklist else "_peaks.bed"
    return [
        f"{data_dir}/Important_processed/Peaks/callpeaks/{s}{suffix}"
        for s in samps
    ]


def get_peaks_by_mark_condition(wildcards):
    """Return list of peak BED files for a mark_condition."""
    return _get_peaks_by_mark_condition_base(wildcards, use_blacklist=False)


def get_peaks_by_mark_condition_blacklist(wildcards):
    """Return list of blacklist-filtered peak BED files for a mark_condition."""
    return _get_peaks_by_mark_condition_base(wildcards, use_blacklist=True)


# ============================================================================
# 4. PEAK CALLING HELPERS
# ============================================================================

def is_broad_mark(wildcards):
    """
    Determine if a sample uses broad peak calling based on BROAD_MARKS config.
    """
    row = st[st['sample'] == wildcards.sample]
    if row.empty:
        return False
    mark = str(row['mark'].iloc[0])
    mark_lower = mark.lower()
    broad_marks = _get_broad_marks_set()
    return mark_lower in broad_marks


def get_peak_file_for_sample(sample):
    """
    Get the appropriate MACS2 peak file (broadPeak or narrowPeak) for a given sample.
    """
    row = st[st['sample'] == sample]
    if row.empty:
        return None
    mark = str(row['mark'].iloc[0])
    mark_lower = mark.lower()
    broad_marks = _get_broad_marks_set()
    data_dir = get_data_dir()
    
    if mark_lower in broad_marks:
        return os.path.join(
            data_dir,
            "Important_processed",
            "Peaks",
            "callpeaks",
            f"macs2_broad_{sample}_peaks.broadPeak"
        )
    else:
        return os.path.join(
            data_dir,
            "Important_processed",
            "Peaks",
            "callpeaks",
            f"macs2_narrow_{sample}_peaks.narrowPeak"
        )


def get_qvalue_for_mark(wildcards):
    """
    Return the q-value threshold for a given sample based on its marker.
    Looks up the marker from the sample table and matches it with config['PEAK_QVAL'].
    Returns a default of 0.01 if not specified.
    """
    row = st[st['sample'] == wildcards.sample]
    if row.empty:
        return 0.01
    
    mark = str(row['mark'].iloc[0])
    peak_qval = config.get('PEAK_QVAL', {})
    
    # Try exact match first
    if mark in peak_qval:
        return peak_qval[mark]
    
    # Try case-insensitive match
    mark_lower = mark.lower()
    for key, value in peak_qval.items():
        if key.lower() == mark_lower:
            return value
    
    # Default q-value
    return 0.01


def macs2_extra_params(wildcards):
    """
    Return MACS2 extra parameters based on marker type from BROAD_MARKS config.
    """
    row = st[st['sample'] == wildcards.sample]
    if row.empty:
        return "--nomodel --keep-dup all --format BAMPE"
    
    mark = str(row['mark'].iloc[0])
    mark_lower = mark.lower()
    base = "--nomodel --keep-dup all --format BAMPE"
    broad_marks = _get_broad_marks_set()
    
    if mark_lower in broad_marks:
        return base + " --broad"
    return base


def get_macs2_outputs(data_dir, igg_name="IgG"):
    """
    Get MACS2 output files based on marker type.
    Returns broad peak outputs for markers in BROAD_MARKS config, narrow peak outputs for others.
    Excludes IgG control samples.
    """
    outputs = []
    peak_dir = os.path.join(data_dir, "Important_processed", "Peaks", "callpeaks")
    broad_marks = _get_broad_marks_set()
    igg_lower = igg_name.lower()
    
    for sample in st['sample'].unique():
        # Skip IgG samples
        if igg_lower in str(sample).lower():
            continue
        
        row = st[st['sample'] == sample]
        if row.empty:
            continue
        
        mark = str(row['mark'].iloc[0])
        mark_lower = mark.lower()
        
        # Skip if mark is IgG
        if igg_lower in mark_lower:
            continue
            
        # Check if it's a broad peak marker based on config
        if mark_lower in broad_marks:
            # Broad peak outputs
            outputs.extend([
                f"{peak_dir}/macs2_broad_{sample}_peaks.xls",
                f"{peak_dir}/macs2_broad_{sample}_peaks.broadPeak",
                f"{peak_dir}/macs2_broad_{sample}_peaks.gappedPeak",
            ])
        else:
            # Narrow peak outputs
            outputs.extend([
                f"{peak_dir}/macs2_narrow_{sample}_peaks.xls",
                f"{peak_dir}/macs2_narrow_{sample}_peaks.narrowPeak",
                f"{peak_dir}/macs2_narrow_{sample}_summits.bed",
            ])
    return outputs


# ============================================================================
# 5. REPRODUCIBILITY HELPERS
# ============================================================================
# (Currently handled by get_sample_reps() in section 2)


# ============================================================================
# 6. INTERNAL HELPERS (for config and data processing)
# ============================================================================

def defect_mode(wildcards, attempt):
    """Return defect mode flag for preseq based on retry attempt."""
    if attempt == 1:
        return ""
    elif attempt > 1:
        return "-D"
    return ""


def _resolve_config_paths(node, paths_dict):
    """
    Recursively replace placeholders such as {raw_data} using values defined in
    the PEP `paths` block.
    """
    if not paths_dict:
        return node
    if isinstance(node, dict):
        for key, value in node.items():
            node[key] = _resolve_config_paths(value, paths_dict)
        return node
    if isinstance(node, list):
        return [_resolve_config_paths(item, paths_dict) for item in node]
    if isinstance(node, str) and "{" in node and "}" in node:
        resolved = node
        for placeholder, replacement in paths_dict.items():
            if replacement is None:
                continue
            resolved = resolved.replace(f"{{{placeholder}}}", str(replacement))
        return resolved
    return node


def _explode_list_columns(df):
    """
    Expand list-valued columns that can be produced by PEP modifiers.
    """
    list_cols = [
        col
        for col in df.columns
        if df[col].apply(lambda x: isinstance(x, (list, tuple))).any()
    ]
    if not list_cols:
        return df
    expanded = df.copy()
    for col in list_cols:
        expanded = expanded.explode(col, ignore_index=True)
    return expanded


def _ensure_sample_columns(df):
    """
    Normalize column names expected by downstream logic.
    """
    normalized = df.copy()
    if "sample" not in normalized.columns:
        if "sample_name" in normalized.columns:
            normalized["sample"] = normalized["sample_name"].astype(str)
        elif normalized.index.name:
            normalized["sample"] = normalized.index.astype(str)
        else:
            raise ValueError(
                "PEP sample table is missing a 'sample' or 'sample_name' column."
            )
    normalized["sample"] = normalized["sample"].astype(str)
    if "run" not in normalized.columns:
        normalized["run"] = 1
    normalized["run"] = pd.to_numeric(normalized["run"], errors="raise").astype(int)
    return normalized
