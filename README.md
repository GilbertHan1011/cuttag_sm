# CUT&Tag Snakemake Pipeline

A comprehensive Snakemake workflow for processing CUT&Tag sequencing data, from raw FASTQ files to peak calling, quality control, and reproducibility analysis.

## Overview

This pipeline processes CUT&Tag (Cleavage Under Targets and Tagmentation) sequencing data through the following stages:

1. **Preprocessing**: Adapter trimming and quality control (optional)
2. **Alignment**: Read alignment to reference genome (Bowtie2 or BWA-MEM2)
3. **Post-alignment**: Duplicate marking, indexing, fragment length analysis
4. **Peak Calling**: MACS2 peak calling with marker-specific parameters
5. **Quality Control**: FRiP, library complexity, reproducibility metrics
6. **Reporting**: MultiQC integration with custom statistics

## Features

- **Flexible alignment**: Supports Bowtie2 and BWA-MEM2 aligners
- **Marker-specific peak calling**: Configurable q-values and broad/narrow peak parameters per histone mark
- **Comprehensive QC**: FRiP, library complexity (Preseq), genomic coverage, reproducibility analysis
- **MultiQC integration**: Custom statistics including peak counts, genomic coverage, and Preseq metrics
- **Reproducibility assessment**: Jaccard similarity and correlation analysis for replicates
- **Blacklist filtering**: Optional blacklist region filtering
- **Conda/Singularity support**: Reproducible environments

## Requirements

### Software Dependencies

- **Snakemake** ≥ 8.20.1
- **Python** ≥ 3.7
- **Conda** or **Singularity** (for environment management)
- **SLURM** (optional, for cluster execution)

### Reference Data

- Reference genome FASTA file
- Genome index (Bowtie2 or BWA-MEM2)
- Chromosome sizes file
- Gene annotation BED file (optional)
- Adapter sequences FASTA (if trimming)
- Blacklist regions BED (optional)

## Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd cuttag_sm
```

2. Install Snakemake (if not already installed):
```bash
conda install -c bioconda -c conda-forge snakemake
```

3. Configure your environment:
   - Edit `config/config.yml` with your project-specific settings
   - Prepare your sample sheet (see Configuration section)

## Configuration

### Main Configuration File (`config/config.yml`)

The main configuration file contains project information, pipeline settings, and reference paths:

```yaml
project_info:
  name: "your_project_name"
  description: "Project description"
  genome: "hg38"  # or "mm10", etc.
  samples:
    table: "config/your_samplesheet.csv"
    pep_config: "config/pep/project_config.yaml"
  outputs:
    root: "{process_dir}/{project_name}/"

pipeline_steps:
  alignment:
    aligner: "bowtie2"  # or "bwa-mem2"
  peaks:
    use_igg: false
    igg_token: "IgG"
    n_intersects: 2
  peak_qval:
    CTCF: 0.0001
    H3K4me3: 0.001
    H3K27ac: 0.001
    # ... marker-specific q-values
  broad_marks:
    - H3K27me3
    - H3K9me3
  broad_params:
    slocal: 10000
    max-gap: 1000
    broad-cutoff: 0.01

reference:
  genome_index: "{database_dir}/hg38/indices_for_Bowtie2/hg38"
  genome_fasta: "{database_dir}/hg38/hg38.fa"
  chrom_sizes: "genes/hg38.chrom.sizes"
  adapter_fasta: "genes/adapter_seqs.fa"
  blacklist: ""  # Optional: path to blacklist BED file
```

### Sample Sheet Format

The sample sheet should be a CSV file with the following columns:

- `sample`: Sample identifier
- `run`: Run/replicate number
- `R1`: Path to R1 FASTQ file
- `R2`: Path to R2 FASTQ file
- `mark`: Histone mark or protein (e.g., "H3K27ac", "CTCF")
- `condition`: Experimental condition
- `igg`: IgG control sample (if applicable)
- `sample_base`: Base sample name for grouping replicates

Example:
```csv
sample,run,R1,R2,mark,condition,igg,gopeaks,sample_base
Sample1,1,/path/to/Sample1_R1.fastq.gz,/path/to/Sample1_R2.fastq.gz,H3K27ac,Condition1,Sample1_IgG,-,Sample1
Sample1,2,/path/to/Sample1_rep2_R1.fastq.gz,/path/to/Sample1_rep2_R2.fastq.gz,H3K27ac,Condition1,Sample1_IgG,-,Sample1
```

## Usage

### Basic Execution

Run the pipeline with Snakemake:

```bash
snakemake --configfile config/config.yml -c <number_of_cores>
```

### Cluster Execution (SLURM)

For cluster environments, use a profile:

```bash
snakemake --configfile config/config.yml --profile slurm
```

### Dry Run

Test the workflow without executing:

```bash
snakemake --configfile config/config.yml -n
```

### Specific Targets

Run specific outputs:

```bash
snakemake --configfile config/config.yml multiqc_report.html
```

## Pipeline Steps

### 1. Preprocessing (`rules/preprocess.smk`)

- **Adapter trimming**: Optional trimming using fastp
- **Quality control**: FastQC-style reports

### 2. Alignment (`rules/align.smk`)

- **Read alignment**: Bowtie2 or BWA-MEM2
- **Sorting**: Coordinate-sorted BAM files
- **Duplicate marking**: Sambamba markdup
- **Indexing**: BAM index generation

### 3. Post-alignment (`rules/postalign.smk`)

- **Fragment length analysis**: Calculate fragment size distributions
- **BigWig track generation**: Normalized coverage tracks

### 4. Peak Calling (`rules/peaks.smk`)

- **MACS2 peak calling**:
  - **Broad peaks**: For marks like H3K27me3, H3K9me3
  - **Narrow peaks**: For other marks (CTCF, H3K4me3, etc.)
- **Marker-specific parameters**:
  - Q-values configured per mark
  - Broad peak parameters (slocal, max-gap, broad-cutoff)
- **Blacklist filtering**: Optional removal of blacklisted regions
- **Peak processing**: Conversion to BED format

### 5. Quality Control (`rules/qc.smk`)

- **FRiP (Fraction of Reads in Peaks)**: DeepTools plotEnrichment
- **Library complexity**: Preseq lc_extrap
- **Fingerprinting**: DeepTools plotFingerprint
- **MultiQC report**: Aggregated QC metrics

### 6. Reproducibility (`rules/reproducibility.smk`)

- **BAM correlation**: DeepTools multiBamSummary and plotCorrelation
- **Peak similarity**: bedtools jaccard for replicate comparison

### 7. Statistics

- **Peak counting**: Per-sample peak counts
- **Genomic coverage**: Percentage of genome covered by peaks

## Output Structure

```
{output_base_dir}/
├── Important_processed/
│   ├── Bam/
│   │   ├── {sample}.sorted.markd.bam
│   │   └── {sample}.sorted.markd.bam.bai
│   ├── Track/
│   │   └── tracks/{sample}.bw
│   └── Peaks/
│       └── callpeaks/
│           ├── macs2_broad_{sample}_peaks.broadPeak
│           ├── macs2_narrow_{sample}_peaks.narrowPeak
│           └── {sample}_peaks.bed
├── Report/
│   ├── multiqc/
│   │   └── multiqc_report.html
│   ├── plotEnrichment/
│   │   ├── frip_{sample}.tsv
│   │   └── frip.html
│   ├── preseq/
│   │   └── lcextrap_{sample}.txt
│   ├── peak_stat/
│   │   ├── peakcount/
│   │   │   └── {sample}_peakcount.txt
│   │   └── coverage_report.tsv
│   ├── bamReproducibility/
│   │   ├── {sample_rep}_global_rep_cor.txt
│   │   └── {sample_rep}_global_rep_heatmap.pdf
│   └── bedtools_jaccard/
│       └── {sample_rep}_jaccard.txt
└── logs/
    └── {rule}_{sample}.log
```

## Quality Control Metrics

The pipeline generates comprehensive QC metrics:

1. **FRiP (Fraction of Reads in Peaks)**: Percentage of reads overlapping called peaks
2. **Library Complexity**: Expected distinct reads at various sequencing depths (Preseq)
3. **Genomic Coverage**: Percentage of genome covered by peaks
4. **Peak Counts**: Number of peaks per sample
5. **Reproducibility**:
   - Spearman correlation between replicates (BAM level)
   - Jaccard similarity between peak sets

All metrics are integrated into the MultiQC report for easy visualization.

## Customization

### Marker-Specific Peak Calling

Configure q-values and broad peak parameters in `config/config.yml`:

```yaml
peak_qval:
  CTCF: 0.0001
  H3K4me3: 0.001
  H3K27ac: 0.001

broad_marks:
  - H3K27me3
  - H3K9me3

broad_params:
  slocal: 10000
  max-gap: 1000
  broad-cutoff: 0.01
```

### Aligner Selection

Switch between aligners in the config:

```yaml
pipeline_steps:
  alignment:
    aligner: "bowtie2"  # or "bwa-mem2"
```

### Adapter Trimming

Enable/disable adapter trimming:

```yaml
TRIM_ADAPTERS: true  # or false
```

## Troubleshooting

### Common Issues

1. **Missing input files**: Ensure sample sheet paths are correct and files exist
2. **Memory errors**: Adjust `mem_mb` in rule resources or use cluster execution
3. **Wildcard errors**: Check that sample names are consistent across config and sample sheet
4. **Conda environment issues**: Rebuild environments with `--use-conda`

### Log Files

Check rule-specific logs in `{output_base_dir}/logs/` for detailed error messages.

## Citation

If you use this pipeline, please cite:

- CUT&Tag method: [Kaya-Okur et al., 2019](https://www.nature.com/articles/s41467-019-09982-5)
- Snakemake: [Mölder et al., 2021](https://doi.org/10.12688/f1000research.29032.2)
- MACS2: [Zhang et al., 2008](https://doi.org/10.1186/gb-2008-9-9-r137)

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Support

For issues, questions, or contributions, please open an issue on the repository.
