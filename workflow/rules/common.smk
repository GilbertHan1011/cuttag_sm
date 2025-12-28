# common holds pyhton function to be used in the snakefile
import pandas as pd
import os

def get_data_dir():
    return config["output_base_dir"].rstrip("/")

# map samples to fastqs
def get_samples():
    """
    return list of samples from samplesheet.tsv
    """
    return sorted(list(pd.Index(st.index).unique()))

def get_marks():
    """
    return list of marks from samplesheet.tsv
    """
    return sorted(list(st['mark'].astype(str).unique()))

def get_mark_conditions():
    """
    return list of samples by condition
    """
    st['mark_condition']=st['mark'].astype(str)+"_"+st['condition']
    return st['mark_condition'].unique().tolist()

def get_tracks_by_mark_condition(wildcards):
    """
    return list of tracks by mark_condition
    """
    st['mark_condition']=st['mark'].astype(str)+"_"+st['condition']
    samps = st.groupby(["mark_condition"])["sample"].apply(list)[wildcards.mark_condition]
    return [f"{get_data_dir()}/Important_processed/Track/tracks/{s}.bw" for s in samps]
    
def get_peaks_by_mark_condition(wildcards):
    """
    return list of peaks by mark_condition
    """
    st["mark_condition"] = st["mark"].astype(str) + "_" + st["condition"]
    samps = st.groupby(["mark_condition"])["sample"].apply(list)[wildcards.mark_condition]
    return [
        f"{get_data_dir()}/Important_processed/Peaks/callpeaks/{s}_peaks.bed"
        for s in samps
    ]

def get_peaks_by_mark_condition_blacklist(wildcards):
    """
    return list of peaks by mark_condition
    """
    st["mark_condition"] = st["mark"].astype(str) + "_" + st["condition"]
    samps = st.groupby(["mark_condition"])["sample"].apply(list)[wildcards.mark_condition]
    return [
        f"{get_data_dir()}/Important_processed/Peaks/callpeaks/{s}_peaks_noBlacklist.bed"
        for s in samps
    ]


def get_bowtie2_input(wildcards):
    """
    returns reads associated with a sample
    """
    try:
        r1=st.loc[wildcards.sample]['R1']
        r2=st.loc[wildcards.sample]['R2']
        return r1,r2
    except Exception:
        row = st[st['sample'] == wildcards.sample].iloc[0]
        return row['R1'], row['R2']

def get_bowtie2_input_by_run(wildcards):
    """
    returns reads associated with a specific sample/run row
    requires 'sample' and 'run' wildcards
    """
    row = st[(st['sample'] == wildcards.sample) & (st['run'] == wildcards.run)]
    if row.empty:
        raise ValueError(f"No samplesheet row for sample={wildcards.sample}, run={wildcards.run}")
    row = row.iloc[0]
    return row['R1'], row['R2']

def get_runs_for_sample(wildcards):
    """
    return list of run identifiers for a biological sample
    """
    runs = st[st['sample'] == wildcards.sample]['run'].astype(str).tolist()
    return runs

def get_sorted_bams_for_sample(wildcards):
    """
    Build paths to per-run sorted BAMs for a sample
    """
    runs = get_runs_for_sample(wildcards)
    return [
        f"{get_data_dir()}/middle_file/aligned/{wildcards.sample}.{r}.sort.bam"
        for r in runs
    ]

def get_reads():
    """
    get list of all reads
    """
    reads = []
    for _, row in st.iterrows():
        sample = row['sample']
        run = str(row['run']) if 'run' in row else '1'
        reads.append(f"{sample}.{run}_R1")
        reads.append(f"{sample}.{run}_R2")
    return reads

def get_igg(wildcards):
    """
    Returns the igg file for the sample unless
    the sample is IgG then no control file is used.
    """ 
    if config['USEIGG']:
        row = st[st["sample"] == wildcards.sample]
        igg = str(row["igg"].iloc[0]) if not row.empty else ""
        iggbam = (
            f"{get_data_dir()}/Important_processed/Bam/{igg}.sorted.markd.bam"
        )
        isigg = config["IGG"] in wildcards.sample
        if not isigg:
            return f'-c {iggbam}'
        else:
            return ""
    else:
        return ""

def get_callpeaks(wildcards):
    """
    Returns the callpeaks input files
    """
    bam = f"{get_data_dir()}/Important_processed/Bam/{wildcards.sample}.sorted.markd.bam"
    # Only the BAM is needed by gopeaks; index will be created by previous rule
    return [bam]


def macs2_extra_params(wildcards):
    """
    Returns MACS2 parameters based on marker type from BROAD_MARKS config.
    """
    row = st[st['sample'] == wildcards.sample]
    mark = str(row['mark'].iloc[0]) if not row.empty else ""
    mark_lower = mark.lower()
    base = "--nomodel --keep-dup all --format BAMPE"
    # Get broad marks from config, convert to lowercase for comparison
    broad_marks = set(m.lower() for m in config.get('BROAD_MARKS', []))
    if mark_lower in broad_marks:
        return base + " --broad"
    return base

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

def get_macs2_outputs(data_dir, igg_name="IgG"):
    """
    Get MACS2 output files based on marker type.
    Returns broad peak outputs for markers in BROAD_MARKS config, narrow peak outputs for others.
    Excludes IgG control samples.
    """
    outputs = []
    peak_dir = os.path.join(data_dir, "Important_processed", "Peaks", "callpeaks")
    # Get broad marks from config, convert to lowercase for comparison
    broad_marks = set(m.lower() for m in config.get('BROAD_MARKS', []))
    
    for sample in st['sample'].unique():
        # Skip IgG samples
        if igg_name.lower() in str(sample).lower():
            continue
        
        row = st[st['sample'] == sample]
        if row.empty:
            continue
        mark = str(row['mark'].iloc[0])
        mark_lower = mark.lower()
        
        # Skip if mark is IgG
        if igg_name.lower() in mark_lower:
            continue
            
        # Check if it's a broad peak marker based on config
        if mark_lower in broad_marks:
            # Broad peak outputs
            outputs.extend(
                [
                    f"{peak_dir}/macs2_broad_{sample}_peaks.xls",
                    f"{peak_dir}/macs2_broad_{sample}_peaks.broadPeak",
                    f"{peak_dir}/macs2_broad_{sample}_peaks.gappedPeak",
                ]
            )
        else:
            # Narrow peak outputs
            outputs.extend(
                [
                    f"{peak_dir}/macs2_narrow_{sample}_peaks.xls",
                    f"{peak_dir}/macs2_narrow_{sample}_peaks.narrowPeak",
                    f"{peak_dir}/macs2_narrow_{sample}_summits.bed",
                ]
            )
    return outputs

def is_broad_mark(wildcards):
    """
    Determines if a sample uses broad peak calling based on BROAD_MARKS config.
    """
    row = st[st['sample'] == wildcards.sample]
    if row.empty:
        return False
    mark = str(row['mark'].iloc[0])
    mark_lower = mark.lower()
    # Get broad marks from config, convert to lowercase for comparison
    broad_marks = set(m.lower() for m in config.get('BROAD_MARKS', []))
    return mark_lower in broad_marks

def defect_mode(wildcards, attempt):
    if attempt == 1:
        return ""
    elif attempt > 1:
        return "-D"

def get_qvalue_for_mark(wildcards):
    """
    Returns the q-value threshold for a given sample based on its marker.
    Looks up the marker from the sample table and matches it with config['PEAK_QVAL'].
    Returns a default of 0.01 if not specified.
    """
    row = st[st['sample'] == wildcards.sample]
    if row.empty:
        return 0.01
    
    mark = str(row['mark'].iloc[0])
    # Try exact match first
    if mark in config.get('PEAK_QVAL', {}):
        return config['PEAK_QVAL'][mark]
    
    # Try case-insensitive match
    mark_lower = mark.lower()
    for key, value in config.get('PEAK_QVAL', {}).items():
        if key.lower() == mark_lower:
            return value
    
    # Default q-value
    return 0.01

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



def get_peak_file_for_sample(sample):
    """
    Get the appropriate MACS2 peak file (broadPeak or narrowPeak) for a given sample.
    """
    row = st[st['sample'] == sample]
    if row.empty:
        return None
    mark = str(row['mark'].iloc[0])
    mark_lower = mark.lower()
    # Get broad marks from config, convert to lowercase for comparison
    broad_marks = set(m.lower() for m in config.get('BROAD_MARKS', []))
    
    if mark_lower in broad_marks:
        return os.path.join(
            DATA_DIR,
            "Important_processed",
            "Peaks",
            "callpeaks",
            f"macs2_broad_{sample}_peaks.broadPeak"
        )
    else:
        return os.path.join(
            DATA_DIR,
            "Important_processed",
            "Peaks",
            "callpeaks",
            f"macs2_narrow_{sample}_peaks.narrowPeak"
        )
