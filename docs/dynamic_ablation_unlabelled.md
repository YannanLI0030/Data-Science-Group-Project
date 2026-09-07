# Dynamic unlabelled ablation: panel and output guide

## Purpose and boundaries

This stage separates two questions that require different gene panels:

1. `ablation_100genes_unfiltered_random_*` measures behaviour under naturally
   sparse multi-omics coverage (missing-data stress test).
2. `ablation_100genes_coverage_stratified_*` makes protein, confidence, and
   adaptive-trust interventions identifiable across pre-declared coverage
   strata.

Neither panel currently has independently verified relevance labels. This run
therefore measures **ranking sensitivity and robustness**, not recommendation
accuracy. It does not compute NDCG, Recall, MRR, or select a winning weight
profile.

The V4 candidate pool and gold remain development evidence. After configuration
definitions are frozen, a separate 12-20-gene V5 subset can be pooled, blinded,
judged, and evaluated once.

## Inputs

- Panel lists and long tables: `benchmarks/panels/`
- Matching sample universe: local `9_DepMap_sample_info.csv`
- RNA-protein correlation table
- Latest dynamic production scorer, used only for A0 numeric-parity checks

Duplicate `(gene, source, DepMap_ID)` measurements are aggregated by arithmetic
mean of `value_raw` and `value_std`, matching the production dynamic adapter.
Rows whose DepMap IDs are absent from the matching sample-info snapshot are
removed and counted in the manifest.

## Controlled configurations

The runner includes the original B0/A0/A1a/A1b/A2/A3/A4/A5/A6 configurations
and adds:

- `A1c_no_protein_confidence`: direct protein retained, protein-derived
  confidence removed;
- `A2a_no_conf_completeness`;
- `A2b_no_conf_source_support`;
- `A2c_no_conf_consistency`.

For each confidence-component ablation, retained internal weights are
renormalised to one. A4 uses the frozen breakpoints: correlation below 0.25 is
low trust, 0.25-0.50 is the transition, and at least 0.50 is full RNA trust.

## Local command

```bash
cd "/Users/liyannan/Desktop/Data-Science-Group-Project"

/Users/liyannan/miniconda3/envs/cellline-v3/bin/python \
  scripts/dynamic_ablation_runner.py \
  --panels-dir "benchmarks/panels" \
  --panel both \
  --sample-info "/Users/liyannan/Desktop/cellline_selector_v3/data_s3/nomenclature/9_DepMap_sample_info.csv" \
  --corr "/Users/liyannan/Desktop/Data Science MSc Graduation Project/gene_rna_protein_correlations.csv" \
  --production-script "/Users/liyannan/Desktop/cellline_selector_v3/dynamic_cellline_selector_gene_protein.py" \
  --out-dir "results/dynamic_ablation_unlabelled" \
  --top-k 10
```

Use `--force` only for a deliberate regeneration. Run `--self-test` after code
changes.

## Outputs

- `dynamic_ablation_unlabelled_summary.csv`: aggregate sensitivity by panel,
  stratum, and configuration, plus `ALL` rows.
- `dynamic_ablation_unlabelled_by_gene_config.csv`: one row per
  panel/gene/configuration with coverage, tie, score-variance, and A0-comparison
  diagnostics.
- `dynamic_ablation_unlabelled_gene_diagnostics.csv`: configuration-collapse
  count for each gene.
- `dynamic_ablation_unlabelled_top10.csv`: Top-10 with all exact-score ties
  retained; this can exceed ten rows per gene/configuration.
- `dynamic_ablation_unlabelled_manifest.json`: code/data hashes, formulae,
  aggregation, parity, and interpretation restrictions.

## Current run: structural findings

All 200 panel genes passed numeric A0 parity against the latest dynamic scorer.

In the unfiltered random panel:

- only 19/100 genes had protein evidence;
- 43/100 genes produced one identical full ranking across all 14 configs;
- A1a and A0 had identical full rankings for 81% of genes and identical Top-10
  sets for 85%;
- A0 confidence was constant within 67% of genes.

This is the expected sparse-data stress-test result: many interventions cannot
act because their required modality is absent.

In the coverage-stratified panel:

- 50/100 genes had protein evidence by design;
- 20 of the 25 DepMap-RNA-only genes, and no genes in the other strata,
  collapsed to one full ranking across every configuration;
- A1a and A0 had identical full rankings for 50% and identical Top-10 sets for
  63%, showing that the protein-rich half makes direct-protein removal visible;
- A4 was identical to A0 for all high-correlation genes, as defined by the
  intervention;
- A4 changed the Top-10 for 70% of mid-correlation genes and all low-correlation
  genes; mean Top-10 Jaccard versus A0 was 0.741 and 0.556 respectively.

For confidence subcomponents over the full stratified panel, Top-10 identity
versus A0 was 66% after removing completeness, 62% after removing source
support, and 78% after removing consistency. These figures describe structural
influence only; without gold labels they do not say which component improves
biological correctness.

Exact four-decimal score rounding can split or merge ties even when a
transformation is monotonic. Use Top-10-with-ties, Jaccard, and rank correlation
together rather than interpreting full-order identity alone.

## Next decision gate

Use these outputs to verify that every intended intervention is active in the
appropriate stratum and benign when its modality is absent. Then freeze the
configuration definitions. Only after freezing should a separate
`candidate_pool_v5_holdout_review.csv` be generated and independently labelled.
