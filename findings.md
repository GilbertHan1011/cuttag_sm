# Findings

- The current legacy chain is fastp -> Bowtie2/BWA-MEM2 -> Sambamba markdup ->
  samtools index, with canonical BAM/BAI under `Important_processed/Bam/`.
- The existing downstream rules consume those canonical paths; this is the
  backend boundary and should remain unchanged.
- Tachyon Upstream accepts raw paired FASTQs and emits mapped-only sorted BAM,
  BAI, native UCF, and stats JSON, but does not emit trimmed FASTQs or fastp
  reports and uses BWA-MEM2 rather than Bowtie2.
- The CUT&Tag repository contains unrelated uncommitted UCF work in the
  Snakefile/postalign/test paths. Preserve it while integrating the backend.
- The backend is config-selected after Snakefile hydration. Legacy rules are
  absent in the Tachyon DAG; shared downstream rules consume canonical paths.
- The pinned Tachyon binary is source commit `d30dc24` with SHA256
  `19a2011d0ecd1a38ccbb172692eece50c99ad85222cfbb4ca2386581fa8d3046`.
- Integrated fixture outputs match native fixture stats exactly: 2,500 input
  pairs, 4,971 mapped BAM records, 2,477 UCF rows, and 38 duplicate pairs.
