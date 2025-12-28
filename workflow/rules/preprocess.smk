

DATA_DIR = config["output_base_dir"].rstrip("/")


rule fastp:
    input:
        r1 = lambda wc: get_bowtie2_input_by_run(wc)[0],
        r2 = lambda wc: get_bowtie2_input_by_run(wc)[1]
    output:
        r1 = temp(f"{DATA_DIR}/middle_file/Trimmed_fastq/fastp/{{sample}}.{{run}}_R1.fastq.gz"),
        r2 = temp(f"{DATA_DIR}/middle_file/Trimmed_fastq/fastp/{{sample}}.{{run}}_R2.fastq.gz"),
        report_html = os.path.join(DATA_DIR,"Report","fastp","{sample}_run{run}_fastp.html"),
        report_json = os.path.join(DATA_DIR,"Report","fastp","{sample}_run{run}_fastp.json"),
    params:
        adapter_fasta_file = config["ADAPTER_FASTA"]
    conda: 
        "../envs/fastp.yml"
    threads: 4
    resources:
        runtime = 120,
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --detect_adapter_for_pe --trim_poly_g --adapter_fasta {params.adapter_fasta_file} --thread {threads} -j {output.report_json} -h {output.report_html}"
