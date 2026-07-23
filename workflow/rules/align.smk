from shlex import quote

DATA_DIR = config["output_base_dir"].rstrip("/")
ALIGNMENT_BACKEND = config["alignment_backend"]
TACHYON_CONFIG = config["pipeline_steps"]["alignment"].get("tachyon_upstream", {})

if ALIGNMENT_BACKEND == "legacy":
    if config.get("aligner", "bowtie2") == "bowtie2":
        rule align:
            input:
                r1 = lambda wc: _all_trimmed_fastqs(wc.sample)[0],
                r2 = lambda wc: _all_trimmed_fastqs(wc.sample)[1]
            output:
                bam = temp(f"{DATA_DIR}/middle_file/aligned/{{sample}}.sort.bam")
            log:
                err=f"{DATA_DIR}/Report/align/bowtie2_{{sample}}.err"
            conda:
                "../envs/align.yml"
            threads: 8
            resources:
                mem_mb=align_mem_mb,
                runtime = 1000,
            params:
                bowtie2_input = lambda w, input: f"-1 {','.join(input.r1)} -2 {','.join(input.r2)}"
            shell:
                (
                "bowtie2 --local --very-sensitive-local "
                "--no-unal --no-mixed --threads {threads} "
                "--no-discordant --phred33 "
                "-I 10 -X 700 -x {config[GENOME]} "
                    "{params.bowtie2_input} 2>{log.err} | "
                    "samtools view -@ {threads} -Sbh - | "
                    "samtools sort -@ {threads} -T {DATA_DIR}/middle_file/aligned/{wildcards.sample}.tmp -o {output.bam} -"
                )
    else:
        rule align:
            input:
                r1 = lambda wc: _all_trimmed_fastqs(wc.sample)[0],
                r2 = lambda wc: _all_trimmed_fastqs(wc.sample)[1]
            output:
                bam = temp(f"{DATA_DIR}/middle_file/aligned/{{sample}}.sort.bam")
            log:
                err=f"{DATA_DIR}/logs/bwa_mem2_{{sample}}.err"
            params:
                rg = lambda wc: f"@RG\\tID:{wc.sample}\\tSM:{wc.sample}\\tPL:{config.get('sequencing_platform', 'ILLUMINA')}",
                ref = lambda wc: config["FASTA"],
                bwa_args = config.get("bwa_args", ""),
                bwa_input = lambda w, input: " ".join(input.r1 + input.r2)
            resources:
                mem_mb=align_mem_mb,
                runtime = 800,
            conda:
                "../envs/bwa.yml"
            threads: 8
            shell:
                (
                "bwa-mem2 mem {params.bwa_args} -t {threads} -R \"{params.rg}\" {params.ref} "
                    "{params.bwa_input} 2>{log.err} | "
                    "samtools view -@ {threads} -Sbh - | "
                    "samtools sort -@ {threads} -T {DATA_DIR}/middle_file/aligned/{wildcards.sample}.tmp -o {output.bam} -"
                )

    rule markdup:
        input:
            bam = rules.align.output.bam
        output:
            bam = f"{DATA_DIR}/Important_processed/Bam/{{sample}}.sorted.markd.bam",
        conda:
            "../envs/sambamba.yml"
        threads: 4
        log:
            f"{DATA_DIR}/Report/sambamba/sambamba_markdup_{{sample}}.log"
        shell:
            f"sambamba markdup --tmpdir={DATA_DIR}/Important_processed/Bam -t {{threads}} {{input.bam}} {{output.bam}} > {{log}} 2>&1"

    rule index_bam:
        input:
            rules.markdup.output.bam
        output:
            bai = f"{DATA_DIR}/Important_processed/Bam/{{sample}}.sorted.markd.bam.bai"
        conda:
            "../envs/samtools.yml"
        threads: 8
        resources:
            runtime = 360,
            mem_mb = 4*config.get("mem", 8000)
        shell:
            "samtools index -@ {threads} {input} {output}"
else:
    rule align_tachyon_upstream:
        input:
            r1 = lambda wc: _paired_raw_fastqs(wc.sample)[0],
            r2 = lambda wc: _paired_raw_fastqs(wc.sample)[1],
            index = lambda wc: [f"{TACHYON_CONFIG['index']}{suffix}" for suffix in (".amb", ".ann", ".pac", ".bwt.2bit.64", ".0123")]
        output:
            bam = f"{DATA_DIR}/Important_processed/Bam/{{sample}}.sorted.markd.bam",
            bai = f"{DATA_DIR}/Important_processed/Bam/{{sample}}.sorted.markd.bam.bai",
            ucf = f"{DATA_DIR}/Important_processed/ucf/{{sample}}.sorted.markd.ucf",
            stats = f"{DATA_DIR}/Report/tachyon_upstream/{{sample}}.stats.json"
        log:
            f"{DATA_DIR}/logs/tachyon_upstream_{{sample}}.log"
        threads: 8
        resources:
            mem_mb=align_mem_mb,
            runtime=1000,
        params:
            executable = TACHYON_CONFIG["executable"],
            adapter = TACHYON_CONFIG.get("adapter", "illumina"),
            index = TACHYON_CONFIG["index"],
            r1_args = lambda wc, input: " ".join(f"--r1 {quote(path)}" for path in input.r1),
            r2_args = lambda wc, input: " ".join(f"--r2 {quote(path)}" for path in input.r2),
            outdir = lambda wc, output: os.path.dirname(output.bam),
            ucf_dir = lambda wc, output: os.path.dirname(output.ucf),
            report_dir = lambda wc, output: os.path.dirname(output.stats),
            log_dir = f"{DATA_DIR}/logs"
        shell:
            "mkdir -p {params.outdir:q} {params.ucf_dir:q} {params.report_dir:q} {params.log_dir:q} && "
            "{params.executable:q} {params.r1_args} {params.r2_args} "
            "--assay cut-and-tag --adapter {params.adapter:q} --index {params.index:q} "
            "--bam {output.bam:q} --ucf {output.ucf:q} --stats {output.stats:q} "
            "--threads {threads} > {log:q} 2>&1"
