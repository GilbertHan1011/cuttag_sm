# Progress

## 2026-07-23

- Recorded the Tachyon-Upstream integration plan.
- Confirmed the intended seam is raw paired FASTQs to the existing canonical
  BAM/BAI outputs, leaving downstream QC and analysis rules shared.
- Implementation starts with static backend contract tests and local dry-runs;
  real-data HPC validation awaits immutable, accessible input manifests.
- Added the opt-in backend with `legacy` as default; all downstream consumers
  now use backend-neutral canonical BAM/BAI paths.
- Local tests: 8 passed. Legacy/Tachyon dry-runs resolve their respective
  producers. The Tachyon fixture run passed BAM/BAI/UCF checks and reproduced
  pinned native stats exactly.
- Isolated HPC pilot completed for `B51-1_CUT22`: both canonical BAMs pass
  integrity/index checks and downstream targets completed. BAM records differ
  by 2.44%, fragment Spearman is 0.99107, FRiP differs by 0.59 pp, and BigWig
  Pearson is 0.97953. Peak overlap/Jaccard misses the predeclared gate, so the
  backend remains opt-in pending the expanded cohort. Details:
  `test/hpc_validation/RESULTS.md`.
