# Optional Tachyon-Upstream Backend Plan

## Goal

Add an opt-in `pipeline_steps.alignment.backend: tachyon_upstream` to replace
only paired FASTQ preprocessing and alignment. Preserve `legacy` as default
and retain the canonical `Important_processed/Bam/{sample}.sorted.markd.bam`
and BAI contract for all downstream work.

## Phases

1. Preserve the pre-existing dirty UCF work and write static backend-contract
   tests before changing workflow rules.
2. Completed: added one raw-FASTQ-to-canonical-BAM Tachyon producer with native
   UCF/stats side outputs. Downstream QC implementations remain shared.
3. Completed locally: legacy and Tachyon DAGs select the intended producer;
   the integrated Tachyon fixture producer passed structural checks.
4. Completed isolated paired pilot and Slurm comparison. Repeat on the expanded
   cohort before promotion because peak-set gates did not pass on this library.

## Acceptance gates

- Legacy remains the default and its producer DAG is unchanged.
- Tachyon schedules no fastp, Bowtie2/BWA legacy align, Sambamba markdup, or
  separate index job.
- Both branches create readable, coordinate-sorted, indexed canonical BAMs.
- Peak, track, FRiP, fingerprint, Preseq, reproducibility, and IDR continue to
  consume canonical paths without backend-specific logic.
- Real-data promotion requires frozen paired inputs and the documented CUT&Tag
  mapping, fragment, peak, signal, and QC agreement gates.

## Blockers

- Existing dirty `workflow/Snakefile`, `workflow/rules/postalign.smk`, and
  `test/` UCF changes need separate review and must not be overwritten.
- Full-data pilot selection awaits an immutable manifest with accessible paired
  raw FASTQs, BWA index, and a pinned HPC Tachyon executable.

## Local evidence (2026-07-23)

- `pytest -q test`: 8 passed.
- Legacy dry-run reached `fastp`; Tachyon dry-run reached only
  `align_tachyon_upstream`. Both stop locally at unavailable archived inputs.
- The Tachyon fixture producer created a coordinate-sorted, indexed, mapped-only
  BAM (4,971 records), valid 13-section UCF, and stats identical to the pinned
  native fixture result.
