
DATA_DIR = config["output_base_dir"].rstrip("/")
ALIGNMENT_BACKEND = config["alignment_backend"]

rule tracks:
    input:
        bam = lambda wc: canonical_bam(wc.sample),
        bai = lambda wc: canonical_bai(wc.sample),
    output:
        f"{DATA_DIR}/Important_processed/Track/tracks/{{sample}}.bw"
    conda:
        "../envs/dtools.yml"
    resources:
        mem_mb=32000,
        runtime = 200,
    threads: 8
    shell:
        "bamCoverage -b {input[0]} -o {output} --binSize 10 --smoothLength 50 --normalizeUsing CPM -p {threads} "

if ALIGNMENT_BACKEND == "legacy":
    rule bam_to_ucf:
        input:
            bam = lambda wc: canonical_bam(wc.sample),
            chrom_sizes = config["CSIZES"]
        output:
            ucf = f"{DATA_DIR}/Important_processed/ucf/{{sample}}.sorted.markd.ucf"
        params:
            ucf_cli = config.get("UCF_CLI", "ucf-cli"),
            outdir = lambda wc, output: os.path.dirname(output.ucf)
        resources:
            mem_mb=8000,
            runtime=240,
        threads: 1
        shell:
            "mkdir -p {params.outdir:q} && "
            "{params.ucf_cli:q} build event-1d "
            "--input {input.bam:q} "
            "--output {output.ucf:q} "
            "--chrom-sizes {input.chrom_sizes:q} "
            "--assay cut-and-tag"

rule merge_bw:
    input:
        get_tracks_by_mark_condition
    output:
        f"{DATA_DIR}/mergebw/{{mark_condition}}.bw"
    conda:
        "../envs/mergebw.yml"
    resources:
        mem_mb=32000,
        runtime = 300,
    shell:
        "bash workflow/src/mergebw.sh -c {config[CSIZES]} -o {output} {input}"

rule fraglength:
    input:
        lambda wc: canonical_bam(wc.sample)
    output:
        f"{DATA_DIR}/Important_processed/Bam/{{sample}}.sorted.markd.fraglen.tsv"
    conda:
        "../envs/align.yml"
    threads: 1
    shell:
        "workflow/src/fraglen-dist.sh {input} {output}"

# This should plot in multiqc
rule fraglength_plot:
    input:
        expand(f"{DATA_DIR}/Important_processed/Bam/{{sample}}.sorted.markd.fraglen.tsv", sample = samps)
    output:
        f"{DATA_DIR}/Report/fraglen.html"
    container: None
    threads: 1
    run:
        pd.options.plotting.backend = "plotly"
        dfs = []
        for i in input:
            cond_marker = [os.path.basename(i).split(".")[0]]
            temp_df = pd.read_csv(i, sep = "\t", index_col = 0, names = cond_marker)
            dfs.append(temp_df)
        df = pd.concat(dfs, axis = 1)
        fraglen = df.plot()
        fraglen.update_layout( 
            title='Fragment Length Distribution', 
            xaxis_title='Fragment Length (bp)', 
            yaxis_title='Counts', 
            legend_title_text='Samples')
        fraglen.write_html(str(output))
