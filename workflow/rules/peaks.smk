
DATA_DIR = config["output_base_dir"].rstrip("/")

rule callpeaks_macs2_broad:
    input:
        treatment=lambda wc: get_callpeaks(wc)[0]
    output:
        f"{DATA_DIR}/Important_processed/Peaks/callpeaks/macs2_broad_{{sample}}_peaks.xls",
        f"{DATA_DIR}/Important_processed/Peaks/callpeaks/macs2_broad_{{sample}}_peaks.broadPeak",
        f"{DATA_DIR}/Important_processed/Peaks/callpeaks/macs2_broad_{{sample}}_peaks.gappedPeak"
    log:
        f"{DATA_DIR}/logs/{{sample}}_macs2peaks_broad.json"
    params:
        extra=(lambda wc: (
            f"-f BAMPE -g " + (
                "hs" if "hg" in str(config.get("GENES", "")).lower() else (
                    "mm" if "mm" in str(config.get("GENES", "")).lower() else "hs"
                )
            ) + 
            f" --nomodel --keep-dup all --broad" +
            (f" --slocal {config['BROAD_PARAMS']['slocal']}" if 'slocal' in config.get('BROAD_PARAMS', {}) else "") +
            (f" --max-gap {config['BROAD_PARAMS']['max-gap']}" if 'max-gap' in config.get('BROAD_PARAMS', {}) else "") +
            (f" --broad-cutoff {config['BROAD_PARAMS']['broad-cutoff']}" if 'broad-cutoff' in config.get('BROAD_PARAMS', {}) else " --broad-cutoff 0.1")
        ))
    threads: 4
    resources:
        mem_mb=32000,
        runtime = 60,
    wrapper:
        "v2.9.1/bio/macs2/callpeak"

rule callpeaks_macs2_narrow:
    input:
        treatment=lambda wc: get_callpeaks(wc)[0]
    output:
        f"{DATA_DIR}/Important_processed/Peaks/callpeaks/macs2_narrow_{{sample}}_peaks.xls",
        f"{DATA_DIR}/Important_processed/Peaks/callpeaks/macs2_narrow_{{sample}}_peaks.narrowPeak",
        f"{DATA_DIR}/Important_processed/Peaks/callpeaks/macs2_narrow_{{sample}}_summits.bed"
    log:
        f"{DATA_DIR}/logs/{{sample}}_macs2peaks_narrow.json"
    params:
        extra=(lambda wc: (
            f"-f BAMPE -g " + (
                "hs" if "hg" in str(config.get("GENES", "")).lower() else (
                    "mm" if "mm" in str(config.get("GENES", "")).lower() else "hs"
                )
            ) + f" --nomodel --keep-dup all -q {get_qvalue_for_mark(wc)}"
        ))
    threads: 4
    resources:
        mem_mb=32000,
        runtime = 60,
    wrapper:
        "v2.9.1/bio/macs2/callpeak"

rule process_peaks:
    """
    Process MACS2 peak files: remove blacklist if present, otherwise just convert to BED format.
    Uses get_peak_file_for_sample to get the appropriate MACS2 output (broadPeak or narrowPeak).
    """
    input:
        peak_file = lambda wc: get_peak_file_for_sample(wc.sample)
    output:
        bed_file = f"{DATA_DIR}/Important_processed/Peaks/callpeaks/{{sample}}_peaks.bed"
    params:
        blacklist = blacklist_file if os.path.isfile(blacklist_file) else None
    conda:
        "../envs/bedtools.yml"
    threads: 1
    shell:
        """
        if [ -n "{params.blacklist}" ] && [ -f "{params.blacklist}" ]; then
            # Remove blacklist regions and convert to BED (first 3 columns)
            bedtools intersect -v -a {input.peak_file} -b {params.blacklist} | cut -f 1-3 > {output.bed_file}
        else
            # Just convert to BED format (first 3 columns)
            cut -f 1-3 {input.peak_file} > {output.bed_file}
        fi
        """


# get consensus
rule consensus:
    input:
        expand(f"{DATA_DIR}/Important_processed/Peaks/callpeaks/{{sample}}_peaks.bed", sample=sample_noigg)
    output:
        consensus_counts = f"{DATA_DIR}/Important_processed/Peaks/counts/{{mark}}_counts.tsv",
        consensus_bed = f"{DATA_DIR}/Important_processed/Peaks/counts/{{mark}}_consensus.bed"
    params:
        blacklist_flag = "-b" if os.path.isfile(blacklist_file) else ""
    conda:
        "../envs/bedtools.yml"
    shell:
        f"OUTPUT_BASE_DIR=\"{DATA_DIR}\" bash workflow/src/consensus_peaks.sh -m {{wildcards.mark}} -n {{config[N_INTERSECTS]}} -o {{output.consensus_counts}} {{params.blacklist_flag}}"
