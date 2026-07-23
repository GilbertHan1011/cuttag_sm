
DATA_DIR = config["output_base_dir"].rstrip("/")
rule plotFinger:
    input:
        lambda wc: canonical_bam(wc.sample),
        lambda wc: canonical_bai(wc.sample)
    output:
        f"{DATA_DIR}/Report/dtools/fingerprint_{{sample}}.tsv"
    conda:
        "../envs/dtools.yml"
    log:
        f"{DATA_DIR}/logs/fingerprint_{{sample}}.log"
    shell:
        "plotFingerprint -b {input[0]} --smartLabels --outRawCounts {output} > {log} 2>&1"

rule frip_plot:
    input:
        expand(f"{DATA_DIR}/Report/plotEnrichment/frip_{{sample}}.tsv", sample = sample_noigg)
    output:
        f"{DATA_DIR}/Report/plotEnrichment/frip.html"
    conda:
        "../envs/plot_report.yml"
    script:
        "../src/frip_plot.py"

rule preseq:
    input:
        lambda wc: canonical_bam(wc.sample)
    output:
        f"{DATA_DIR}/preseq/estimates_{{sample}}.txt"
    resources:
        defect_mode = defect_mode
    conda:
        "../envs/preseq.yml"
    log:
        f"{DATA_DIR}/logs/preseq_{{sample}}.log"
    shell:
        "preseq c_curve -B {resources.defect_mode} -l 1000000000 -o {output} {input} > {log} 2>&1"

rule preseq_lcextrap:
    input:
        lambda wc: canonical_bam(wc.sample)
    output:
        f"{DATA_DIR}/Report/preseq/lcextrap_{{sample}}.txt"
    resources:
        defect_mode = defect_mode
    conda:
        "../envs/preseq.yml"
    log:
        f"{DATA_DIR}/logs/preseq_{{sample}}.log"
    shell:
        "preseq lc_extrap -B {resources.defect_mode} -l 1000000000 -e 1000000000 -o {output} {input} > {log} 2>&1 || (echo 'preseq lc_extrap failed; creating empty output' >> {log} 2>&1; : > {output})"


rule frip:
    input:
        peaks=f"{DATA_DIR}/Important_processed/Peaks/callpeaks/{{sample}}_peaks.bed",
        bam=lambda wc: canonical_bam(wc.sample),
        bai=lambda wc: canonical_bai(wc.sample)
    output:
        f"{DATA_DIR}/Report/plotEnrichment/frip_{{sample}}.png", f"{DATA_DIR}/Report/plotEnrichment/frip_{{sample}}.tsv"
    conda:
        "../envs/dtools.yml"
    resources:
        mem_mb=16000,
        runtime = 120,
    log:
        f"{DATA_DIR}/logs/plotEnrichment_{{sample}}.log"
    shell:
        "bash workflow/src/skip_frip.sh {input.peaks} {input.bam} {output[0]} {output[1]} {log}"

rule genomic_coverage:
    input:
        peaks=lambda wc: get_peak_file_for_sample(wc.sample),
        chrom_sizes=config["CSIZES"]
    output:
        f"{DATA_DIR}/Report/peak_stat/coverage/{{sample}}_coverage.tsv"
    conda:
        "../envs/bedtools.yml"
    resources:
        mem_mb=32000,
        runtime = 120,
    log:
        f"{DATA_DIR}/logs/coverage_{{sample}}.log"
    shell:
        """
        set -euo pipefail
        exec >{log} 2>&1
        # Calculate total genome size
        TOTAL_GENOME_SIZE=$(awk -F'\t' '{{sum += $2}} END {{print sum}}' {input.chrom_sizes})
        
        # Calculate covered bases (merge overlapping intervals first)
        COVERED_BASES=$(sort -k1,1 -k2,2n {input.peaks} | \
                        bedtools merge -i stdin | \
                        awk -F'\t' '{{sum += $3 - $2}} END {{print sum}}')
        
        # Calculate percentage coverage
        COVERAGE_PERCENT=$(echo "scale=6; ($COVERED_BASES / $TOTAL_GENOME_SIZE) * 100" | bc)
        
        # Write output with header
        echo -e "Sample\tCoverage_Percent\tCovered_Bases\tTotal_Genome_Size" > {output}
        echo -e "{wildcards.sample}\t$COVERAGE_PERCENT\t$COVERED_BASES\t$TOTAL_GENOME_SIZE" >> {output}
        """
        
rule count_peaks_per_sample:
    """
    Count peaks for each sample separately.
    """
    input:
        peak_file = lambda wc: get_peak_file_for_sample(wc.sample)
    output:
        count_file = f"{DATA_DIR}/Report/peak_stat/peakcount/{{sample}}_peakcount.txt"
    log:
        f"{DATA_DIR}/logs/count_peaks_{{sample}}.log"
    shell:
        """
        if [ -f {input.peak_file} ] && [ -s {input.peak_file} ]; then
            count=$(wc -l < {input.peak_file})
            echo "{wildcards.sample} $count" > {output.count_file}
            echo "Sample {wildcards.sample}: $count peaks" > {log}
        else
            echo "{wildcards.sample} 0" > {output.count_file}
            echo "Sample {wildcards.sample}: 0 peaks (file not found or empty)" > {log}
        fi
        """


rule bam_correlation_to_multiqc:
    """
    Convert BAM correlation matrices to MultiQC custom content format.
    Processes one correlation matrix per sample_rep, creating separate MultiQC files.
    """
    input:
        matrix_file = os.path.join(DATA_DIR, "Report", "bamReproducibility", "{sample_rep}_global_rep_cor.txt")
    output:
        multiqc_yaml = os.path.join(DATA_DIR, "Report", "bamReproducibility", "{sample_rep}_bam_correlation_mqc.yaml"),
        stats_tsv = os.path.join(DATA_DIR, "Report", "bamReproducibility", "{sample_rep}_bam_correlation_stats_mqc.tsv"),
    conda:
        "../envs/reproducibility.yml"
    log:
        f"{DATA_DIR}/logs/bam_correlation_to_multiqc_{{sample_rep}}.log"
    shell:
        """
        set -euo pipefail
        exec >{log} 2>&1
        
        mkdir -p $(dirname {output.multiqc_yaml})
        python workflow/src/bam_correlation_to_multiqc.py \
            {input.matrix_file} \
            --output-yaml {output.multiqc_yaml} \
            --output-stats {output.stats_tsv} \
            2>> {log}
        """

rule multiqc:
    input:
        expand(f"{DATA_DIR}/Report/plotEnrichment/frip_{{sample}}.tsv", sample=sample_noigg),
        expand(f"{DATA_DIR}/Report/preseq/lcextrap_{{sample}}.txt", sample=samps),
        expand(f"{DATA_DIR}/Report/peak_stat/peakcount/{{sample}}_peakcount.txt", sample=sample_noigg),
        expand(f"{DATA_DIR}/Report/peak_stat/coverage/{{sample}}_coverage.tsv", sample=sample_noigg),
    output:
        f"{DATA_DIR}/Report/multiqc/multiqc_report.html",
        f"{DATA_DIR}/Report/multiqc/multiqc_data/multiqc_data.json"
    conda:
        "../envs/multiqc.yml"
    log:
        f"{DATA_DIR}/logs/multiqc.log"
    shell:
        """
        export LC_ALL=C.UTF-8; export LANG=C.UTF-8;
        
        multiqc {DATA_DIR}/Report/ \
            --ignore {DATA_DIR}/Report/multiqc \
            -o {DATA_DIR}/Report/multiqc \
            -f \
            >> {log} 2>&1
        """
