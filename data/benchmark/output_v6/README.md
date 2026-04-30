# v6 — Current run (full SA stack + 5 fixes)

This is the **current run** with all SA sources loaded and five
classifier / script fixes layered on top of v1.

## What was loaded

| SA source | Loaded? | Notes |
|-----------|:-------:|-------|
| ClinVar (.osa)            | ✅ | 4,402,501 records |
| ClinVar protein (.oga)    | ✅ | 4,554 genes |
| gnomAD v4.1 exomes (.osa, per-chrom)         | ✅ | 25 chromosomes |
| gnomAD v4.1 gene constraints (.oga)         | ✅ | 18,173 genes |
| REVEL v1.3 (.osa, per-chrom)        | ✅ | 24 chromosomes |
| **PhyloP** (.osa, per-chrom)        | ✅ | distilled from gnomAD v4 INFO `phylop` (Zoonomia 241-mammal score) |
| **SpliceAI** (.osa, per-chrom)      | ✅ | distilled from gnomAD v4 INFO `spliceai_ds_max` |
| **ClinGen Gene-Disease Validity (.oga)**     | ✅ | 2,419 Definitive/Strong/Moderate genes — preferred over OMIM per ClinGen SVI / Abou Tayoun 2018 |

## Code fixes vs v1

1. **SpliceAI camelCase mismatch** in `sa_extract.rs`: writer wrote
   `spliceAI` but classifier matched `spliceai|spliceAi|splice_ai`.
   Added `spliceAI` to the match arms. *(introduced in v4)*
2. **PhyloP routing** in `sa_extract.rs`: PhyloP was attached to
   `aa.supplementary` (allele-level) but the classifier read it from
   `variant_supplementary` (variant-level). Read from both.
   *(introduced in v4)*
3. **Indel allele matching** in
   `analysis/acmg_benchmark/real_data/03_evaluate_concordance.py`:
   added `vep_allele(ref, alt)` to convert VCF (REF, ALT) → VEP CSQ
   Allele convention (`-` for deletion; insertion = right portion
   only) before grouping CSQ entries by allele. Without this all
   48,539 indels in the truth fell into NoCall. *(introduced in v5)*
4. **PM2 fires when variant is absent from gnomAD** (new config flag
   `pm2_absent_when_no_record = true`, default). Per ClinGen SVI v1.0
   "absent or extremely rare in population databases" — a missing
   gnomAD record is interpreted as "gnomAD never observed it, i.e.
   absent". For ~78 % of pathogenic ClinVar truth, the variant has
   no gnomAD record (most rare disease variants aren't in gnomAD);
   the prior strict-coverage stance left PM2 NotEvaluated, so PVS1
   had no partner and the PVS+PP SVI rule couldn't fire. PM2
   _Supporting fires went from ~12K → 340K. *(introduced in v6)*
5. **BP4-splice gated to non-PVS1 consequences** in `benign_supporting.rs`.
   Walker 2023 BP4-splice (SpliceAI ≤ 0.1 → Supporting benign)
   applies to splice-territory consequences, not to LOF variants
   whose pathogenic mechanism is the truncation itself. Frameshift /
   stop_gained / start_lost / transcript_ablation / canonical splice
   ±1/±2 no longer get BP4 just because their SpliceAI happens to
   be low. Cleared the BP4+PVS1 conflict on ~5K pathogenic
   frameshifts. *(introduced in v6)*

## Headline metrics

| Metric | Value | Δ vs v1 |
|--------|------:|--------:|
| Same-direction concordance | **70.3 %** | **+15.6 pp** |
| Exact match | **56.8 %** | +4.1 pp |
| Opposite direction | 0.05 % | (unchanged ~0) |
| **Pathogenic recall** | **63.8 %** | **+48 pp** |
| **Likely_pathogenic recall** | **51.8 %** | **+31 pp** |
| Likely_benign recall | **42.4 %** | **+39 pp** |
| Benign recall | **58.0 %** | +25 pp |

## Files

- `clinvar_2star.fastvep.vcf.gz` — bgzipped VCF with ACMG in CSQ INFO
  (~70 MB; vs ~25 GB for the prior pretty-printed JSON output)
- `concordance_matrix.csv` — 5-class truth × predicted matrix
- `concordance_summary.txt` — text rollup
- `concordance_by_chrom.csv` — per-chromosome breakdown
- `concordance_by_consequence.csv` — top consequences × class
- `criterion_firing_rates.csv` — per-criterion fire counts by truth class
- `rule_distribution.csv` — top criteria-set signatures
- `discrepancies.tsv` — opposite-direction calls (top 10k)
- `figures/` — 6 PNG + PDF figures (incl. v1 vs v6 comparison panels)

For the v1 baseline (no PhyloP / no SpliceAI / no ClinGen GDV, no
fixes), see `../output_v1/`. For the version-by-version SA stack and
code-fix diff, see `../RUN_VERSIONS.md`.
