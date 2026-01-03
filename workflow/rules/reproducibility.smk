
DATA_DIR = config["output_base_dir"].rstrip("/")


def get_reproducibility_sample(sample_rep):
    """
    Return all samples that belong to a given reproducibility group.
    Uses the 'replicate_sample_name' column from the annotation table to group replicates.
    Only returns a list if there are >1 samples in the group.
    """
    if "replicate_sample_name" not in st.columns:
        return []
    samples = (
        st.loc[st["replicate_sample_name"] == sample_rep, "sample"]
        .astype(str)
        .unique()
        .tolist()
    )
    return samples if len(samples) > 1 else []

rule bedtools_jaccard:
    """
    Calculate pairwise Jaccard similarity for all pairs of samples in a replicate group.
    Uses processed peaks (with blacklist removed if applicable) from process_peaks rule.
    Outputs a table with all pairwise comparisons.
    """
    input:
        peaks = lambda w: [
            os.path.join(
                DATA_DIR,
                "Important_processed",
                "Peaks",
                "callpeaks",
                f"{s}_peaks.bed"
            )
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
        set -euo pipefail
        exec >{log} 2>&1
        
        # Extract sample names from processed peak file paths
        # Peak files are in format: SAMPLE_peaks.bed
        PEAK_FILES=({input.peaks})
        NUM_SAMPLES=${{#PEAK_FILES[@]}}
        
        # Extract sample names from file paths
        declare -a SAMPLES
        for peak_file in "${{PEAK_FILES[@]}}"; do
            # Extract sample name: remove path and _peaks.bed suffix
            SAMPLE=$(basename "$peak_file" | sed 's/_peaks\\.bed$//')
            SAMPLES+=("$SAMPLE")
        done
        
        # Write header
        echo -e "Sample1\tSample2\tJaccard_Similarity\tIntersection\tUnion_Intersection" > {output.jaccard}
        
        # Calculate pairwise Jaccard for all pairs
        for ((i=0; i<NUM_SAMPLES; i++)); do
            for ((j=i+1; j<NUM_SAMPLES; j++)); do
                SAMPLE1="${{SAMPLES[i]}}"
                SAMPLE2="${{SAMPLES[j]}}"
                PEAK1="${{PEAK_FILES[i]}}"
                PEAK2="${{PEAK_FILES[j]}}"
                
                echo "Comparing $SAMPLE1 vs $SAMPLE2..." >&2
                
                # Run bedtools jaccard and extract values
                # bedtools jaccard output format: 
                #   intersection    union-intersection    jaccard    n_intersections
                JACCARD_OUT=$(bedtools jaccard -a "$PEAK1" -b "$PEAK2" 2>&1 | tail -n 1)
                
                if [ -n "$JACCARD_OUT" ] && ! echo "$JACCARD_OUT" | grep -q "ERROR"; then
                    INTERSECTION=$(echo "$JACCARD_OUT" | awk '{{print $1}}')
                    UNION=$(echo "$JACCARD_OUT" | awk '{{print $2}}')
                    JACCARD=$(echo "$JACCARD_OUT" | awk '{{print $3}}')
                    
                    echo -e "$SAMPLE1\t$SAMPLE2\t$JACCARD\t$INTERSECTION\t$UNION" >> {output.jaccard}
                else
                    echo "Warning: Failed to calculate Jaccard for $SAMPLE1 vs $SAMPLE2" >&2
                    echo -e "$SAMPLE1\t$SAMPLE2\tNA\tNA\tNA" >> {output.jaccard}
                fi
            done
        done
        
        echo "Pairwise Jaccard calculation complete. Results written to {output.jaccard}" >&2
        """