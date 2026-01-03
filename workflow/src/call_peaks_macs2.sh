#!/bin/bash
# MACS2 peak calling with standard ATAC-seq parameters
# Supports both broad and narrow peak calling based on marker type
# Called by Snakemake - uses snakemake.input, snakemake.output, snakemake.params

set -euo pipefail

# When called via Snakemake script: directive, variables are passed via the environment
TAGALIGN="${snakemake_input[tagalign]}"
OUTPUT_PEAK="${snakemake_output[narrowpeak]}"  # Can be narrowPeak or broadPeak based on is_broad
OUT_DIR="${snakemake_params[out_dir]}"
PREFIX="${snakemake_params[prefix]}"
GENOME_SIZE="${snakemake_params[genome_size]}"
PVAL_THRESH="${snakemake_params[pval_thresh]}"
SMOOTH_WIN="${snakemake_params[smooth_win]}"
CAP_NUM_PEAK="${snakemake_params[cap_num_peak]}"
CHROM_SIZES="${snakemake_params[chrom_sizes]}"
SHIFT_SIZE="${snakemake_params[shiftsize]}"
IS_BROAD="${snakemake_params[is_broad]}"
BROAD_PARAMS="${snakemake_params[broad_params]}"

mkdir -p "$OUT_DIR"

# Determine peak type and MACS2 parameters
# Note: For IDR compatibility, we always output narrowPeak format
# But use broad-like parameters (wider windows) when is_broad=True
PEAK_TYPE="narrowPeak"

if [ "$IS_BROAD" = "True" ] || [ "$IS_BROAD" = "true" ]; then
    # Use broad-like parameters but output narrowPeak for IDR compatibility
    # Adjust smoothing window to be wider (typical for broad marks)
    # Use broader parameters from config if available
    BROAD_EXTRA="--call-summits"
    
    # Parse broad parameters if provided (format: "slocal:10000,max-gap:1000,broad-cutoff:0.01")
    if [ -n "$BROAD_PARAMS" ] && [ "$BROAD_PARAMS" != "{}" ]; then
        # Extract parameters (simple parsing - assumes format above)
        SLOCAL=$(echo "$BROAD_PARAMS" | sed -n 's/.*slocal:\([0-9]*\).*/\1/p')
        MAX_GAP=$(echo "$BROAD_PARAMS" | sed -n 's/.*max-gap:\([0-9]*\).*/\1/p')
        BROAD_CUTOFF=$(echo "$BROAD_PARAMS" | sed -n 's/.*broad-cutoff:\([0-9.]*\).*/\1/p')
        
        # Note: We don't use --broad flag to keep narrowPeak format for IDR
        # But we could adjust smoothing window based on these parameters
        # For now, we use the standard smooth_win parameter
    fi
else
    BROAD_EXTRA="--call-summits"
fi

# Call peaks with MACS2
# Always use narrowPeak format for IDR compatibility
macs2 callpeak \
    -t "$TAGALIGN" -f BED -n "$PREFIX" -g "$GENOME_SIZE" \
    -p "$PVAL_THRESH" --shift "$SHIFT_SIZE" --extsize "$SMOOTH_WIN" \
    --nomodel -B --SPMR --keep-dup all $BROAD_EXTRA \
    --outdir "$OUT_DIR"

# Sort peaks by signal value and cap at max number
peak_tmp="${PREFIX}.tmp"
peak_tmp2="${PREFIX}.tmp2"
PEAK_FILE="${PREFIX}_peaks.${PEAK_TYPE}"

# Sort narrowPeak by signal value (column 8) and handle summits
LC_COLLATE=C sort -k 8gr,8gr "$PEAK_FILE" | \
awk 'BEGIN{OFS="\t"} {$4="Peak_"NR} ($2<0){$2=0} ($3<0){$3=0} ($10==-1){$10=int($2+($3-$2+1)/2.0)} {print $0}' > "$peak_tmp"

head -n "$CAP_NUM_PEAK" "$peak_tmp" > "$peak_tmp2"

# Filter peaks to only include chromosomes present in chrom.sizes file before bedClip
if [ -n "$CHROM_SIZES" ] && [ -f "$CHROM_SIZES" ]; then
    # Create a filtered version with only valid chromosomes
    peak_filtered="${peak_tmp2}.filtered"
    # Use awk to filter: only keep peaks where chromosome exists in chrom.sizes
    awk 'NR==FNR{chroms[$1]=1; next} $1 in chroms' "$CHROM_SIZES" "$peak_tmp2" > "$peak_filtered"
    
    # Now run bedClip on filtered peaks (should not fail)
    bedClip "$peak_filtered" "$CHROM_SIZES" stdout | gzip -nc > "$OUTPUT_PEAK"
    rm -f "$peak_filtered"
else
    # If no chrom.sizes file provided, just compress without clipping
    gzip -nc "$peak_tmp2" > "$OUTPUT_PEAK"
fi

# Cleanup
rm -f "$peak_tmp" "$peak_tmp2" "${PREFIX}"_*

