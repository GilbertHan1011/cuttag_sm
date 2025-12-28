
DATA_DIR = config["output_base_dir"].rstrip("/")


def get_reproducibility_sample(sample_rep):
    """
    Return all samples that belong to a given reproducibility group.
    Uses the 'sample_base' column from the annotation table to group replicates.
    Only returns a list if there are >1 samples in the group.
    """
    if "sample_base" not in annot.columns:
        return []
    samples = (
        annot.loc[annot["sample_base"] == sample_rep, "sample"]
        .astype(str)
        .unique()
        .tolist()
    )
    return samples if len(samples) > 1 else []

rule bamReproducibility:
    """
    Global reproducibility / correlation heatmap using deeptools:
    - multiBamSummary bins
    - plotCorrelation heatmap (Spearman)

    This is a Snakemake translation of the provided SLURM script.
    """
    input:
        bams = lambda w: [
            os.path.join(
                DATA_DIR,
                "Important_processed",
                "Bam",
                f"{s}.sorted.markd.bam",
            )
            for s in get_reproducibility_sample(w.sample_rep)
        ]
    output:
        npz = temp(
            os.path.join(
                DATA_DIR,
                "Report",
                "bamReproducibility",
                "{sample_rep}_global_rep.npz",
            )
        ),
        heatmap = os.path.join(
            DATA_DIR,
            "Report",
            "bamReproducibility",
            "{sample_rep}_global_rep_heatmap.pdf",
        ),
        matrix = os.path.join(
            DATA_DIR,
            "Report",
            "bamReproducibility",
            "{sample_rep}_global_rep_cor.txt",
        ),
    params:
        bin_size = 1000,
        prefix = lambda w: w.sample_rep,
    log:
        os.path.join(
            DATA_DIR,
            "logs",
            "bamReproducibility",
            "{sample_rep}_global_rep.log",
        )
    conda:
        "../envs/dtools.yml"
    resources:
        mem_mb = 28000,
        runtime = 480,
    threads: 20
    shell:
        """
        mkdir -p $(dirname {output.npz})
        mkdir -p $(dirname {log})

        echo "Starting multiBamSummary..." >&2
        multiBamSummary bins \
          --bamfiles {input.bams} \
          --binSize {params.bin_size} \
          --numberOfProcessors {threads} \
          --outFileName {output.npz} 2> {log}

        echo "Starting plotCorrelation..." >&2
        plotCorrelation \
            --corData {output.npz} \
            --corMethod spearman \
            --whatToPlot heatmap \
            --plotTitle "Spearman Correlation of {params.prefix}" \
            --plotFile {output.heatmap} \
            --outFileCorMatrix {output.matrix} 2>> {log}
        """

rule bedtools_jaccard:
    input:
        peaks = lambda w: [
            get_peak_file_for_sample(s)
            for s in get_reproducibility_sample(w.sample_rep)
        ]
    output:
        jaccard = os.path.join(
            DATA_DIR,
            "Report",
            "bedtools_jaccard",
            "{sample_rep}_jaccard.txt",
        )
    log:
        os.path.join(
            DATA_DIR,
            "logs",
            "bedtools_jaccard",
            "{sample_rep}_jaccard.log",
        )
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools jaccard -a {input.peaks[0]} -b {input.peaks[1]} > {output.jaccard}
        """