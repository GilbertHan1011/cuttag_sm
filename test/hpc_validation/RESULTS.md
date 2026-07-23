# Tachyon backend pilot: B51-1_CUT22

Validation root: `/storage/zhangkaiLab/hanlitian/macrophage/validation/20260723_tachyon_backend`.
Tachyon source was `d30dc24`; the glibc-2.34-compatible binary SHA256 was
`896b606e092c2239fac23de4a2061917390da47aa743a11dba73d5664c116e3e`.

The isolated legacy Bowtie2 and Tachyon BWA arms used the same declared paired
FASTQs. Both canonical BAMs pass `samtools quickcheck`, are coordinate sorted,
and support indexed `chr1:1-100000` queries.

| Metric | Legacy | Tachyon | Result |
|---|---:|---:|---|
| BAM records | 51,420,036 | 52,676,684 | +2.44% |
| chrM fraction | 0.4400% | 0.4380% | -0.0021 pp |
| duplicate rate | 29.37% | 29.51% | +0.14 pp |
| median fragment | 171 bp | 168 bp | pass |
| fragment histogram Spearman (≤1 kb) | — | 0.99107 | pass |
| FRiP | 4.29% | 3.70% | -0.59 pp |
| 10 kb BigWig Pearson | — | 0.97953 | pass |
| peak count | 9,786 | 7,788 | -20.42% |
| legacy/Tachyon peak overlap | 64.99% / 81.95% | — | below gate |
| peak bp Jaccard | — | 0.5414 | below gate |

The structural, mapping-proxy, fragment, FRiP, and signal checks are close.
Peak-set agreement narrowly misses the predeclared CUT&Tag gates, consistent
with the Bowtie2-versus-BWA geometry difference. Do not relax the gate: repeat
on the planned expanded cohort before considering promotion. The backend remains
opt-in.

Machine-readable results: `$validation_root/comparison/cuttag_pilot.json`.
