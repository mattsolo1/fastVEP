# Benchmark run versions

Each row of the `output_v*/` directories represents one end-to-end run
of the ClinVar 2-star+ benchmark on the same 673,660-variant truth set.
Successive runs differ only in (a) which supplementary annotation
databases were loaded into `--sa-dir` and (b) which classifier / output
bugs had been fixed.

## SA stack per run

|                                     |  v1  |  v2  |  v4  |  v5  |  v6  |
|-------------------------------------|:----:|:----:|:----:|:----:|:----:|
| **Variant-level (.osa)**            |      |      |      |      |      |
| ClinVar                             |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |
| gnomAD v4.1 exomes (per-chrom)      |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |
| REVEL v1.3 (per-chrom)              |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |
| **PhyloP** (per-chrom)              |  ❌  |  ✅  |  ✅  |  ✅  |  ✅  |
| **SpliceAI** (per-chrom)            |  ❌  |  ✅  |  ✅  |  ✅  |  ✅  |
| **Gene-level (.oga)**               |      |      |      |      |      |
| ClinVar protein                     |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |
| gnomAD gene constraints             |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |
| **ClinGen Gene-Disease Validity**   |  ❌  |  ✅  |  ✅  |  ✅  |  ✅  |

(PhyloP and SpliceAI are distilled from gnomAD v4 INFO fields
`phylop` and `spliceai_ds_max` rather than re-downloaded; the gnomAD
v4 sites VCF already includes them. ClinGen GDV substitutes for
OMIM `genemap2.txt` per ClinGen SVI / Abou Tayoun 2018 — same `.oga`
schema, `omim` json_key, but with a multi-curator scored rubric and
explicit Definitive/Strong/Moderate filtering.)

## Code fixes per run

|                                                              |  v1  |  v2  |  v4  |  v5  |  v6  |
|--------------------------------------------------------------|:----:|:----:|:----:|:----:|:----:|
| SpliceAI `spliceAI` json_key recognised in classifier        |      |      |  ✅  |  ✅  |  ✅  |
| PhyloP read from `allele_supplementary` (CLI's actual route) |      |      |  ✅  |  ✅  |  ✅  |
| VCF + bgzip output (vs 25 GB pretty JSON)                    |      |      |  ✅  |  ✅  |  ✅  |
| `vep_allele(ref, alt)` indel matching in concordance script  |      |      |      |  ✅  |  ✅  |
| **PM2 fires when variant absent from gnomAD** (`pm2_absent_when_no_record`) |      |      |      |      |  ✅  |
| **BP4-splice gated to non-PVS1 consequences** (Walker 2023)  |      |      |      |      |  ✅  |

v3 was a partial run (PhyloP+SpliceAI loaded but bugs still latent);
its results are functionally indistinguishable from v2 and were
overwritten before being preserved.

## Headline metrics per run

|                            |     v1     |     v5     |     v6     |  Δ v1→v6  |
|----------------------------|-----------:|-----------:|-----------:|----------:|
| Same-direction concordance |   54.7 %   |   65.1 %   | **70.3 %** |**+15.6 pp**|
| Exact match                |   52.7 %   |   56.0 %   | **56.8 %** | +4.1 pp   |
| Opposite direction         |   0.005 %  |   0.06 %   |   0.05 %   | (≈0)      |
| NoCall                     |   0.0 %    |   0.0 %    |   0.0 %    | —         |
| **Pathogenic recall**      |   **15.7 %** | 20.6 %   | **63.8 %** |**+48 pp** |
| **Likely_pathogenic recall** | **20.9 %** | 26.7 %   | **51.8 %** |**+31 pp** |
| VUS recall                 |   96.6 %   |   92.6 %   |   91.5 %   | -5 pp     |
| **Likely_benign recall**   |   **3.2 %**|   42.4 %   |   42.4 %   |**+39 pp** |
| Benign recall              |   33.2 %   |   58.0 %   |   58.0 %   | +25 pp    |

## Driver of each lift

- **+39 pp LB recall, +25 pp B recall** (v1 → v5): BP7 went from **0**
  → **81,706 fires** once PhyloP+SpliceAI were loaded *and* both
  wiring bugs were fixed. (Walker 2023: BP7 needs synonymous + low
  SpliceAI + low PhyloP.)
- **+48 pp Pathogenic recall, +31 pp LP recall** (v5 → v6): the
  classifier's PM2 evaluator previously refused to fire when
  `input.gnomad` was `None` (no gnomAD record at the variant). For
  ~78 % of the truth's pathogenic class — most rare disease variants
  aren't in gnomAD — that meant PM2 NotEvaluated, so PVS1 had no
  partner criterion and the PVS+≥1 PP (SVI Sept 2020) → LP rule
  couldn't trigger. Default config flag
  `pm2_absent_when_no_record = true` interprets a missing record as
  "absent from gnomAD" per ClinGen SVI v1.0. PM2_Supporting fires
  jumped from ~12K → 340K, unlocking the LP rule for ~50K Pathogenic
  variants. PVS1 also nearly doubled (27K → 50K P+LP) because
  BP4-splice is no longer (incorrectly) firing on frameshift / null
  variants — Walker 2023 explicitly scopes BP4-splice to
  splice-territory consequences.
- **VUS recall slight drop (-5 pp)**: by design — when more benign
  evidence fires, some variants previously called VUS now correctly
  drop to LB / B (which doesn't match a VUS truth). Same-direction
  rate still rises because the P/LP/B/LB gains far outweigh the VUS
  loss.

## Where to find each version

- v1 baseline: `output_v1/concordance_matrix.csv` +
  `output_v1/README.md` (raw outputs were overwritten; matrix
  reconstructed from documentation)
- v6 current: `output_v6/` (full outputs + figures + raw VCF.gz)
