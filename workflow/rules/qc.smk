
DATA_DIR = config["output_base_dir"].rstrip("/")
rule plotFinger:
    input:
        rules.markdup.output.bam,
        rules.index_bam.output.bai
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
        rules.markdup.output
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
        rules.markdup.output
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

rule multiqc:
    input:
        expand(f"{DATA_DIR}/Report/plotEnrichment/frip_{{sample}}.tsv", sample=sample_noigg),
        expand(f"{DATA_DIR}/Report/preseq/lcextrap_{{sample}}.txt", sample=samps),
        expand(f"{DATA_DIR}/Report/peak_stat/peakcount/{{sample}}_peakcount.txt", sample=sample_noigg),
        f"{DATA_DIR}/Report/peak_stat/coverage_report.tsv"
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
        
        multiqc {DATA_DIR}/Report/ {DATA_DIR}/logs/ \
            --ignore {DATA_DIR}/Report/multiqc \
            -o {DATA_DIR}/Report/multiqc \
            -f \
            -c workflow/src/multiqc_conf.yml \
            --cl-config "annotation: {SAMPLE_SHEET}" \
            >> {log} 2>&1
        """

