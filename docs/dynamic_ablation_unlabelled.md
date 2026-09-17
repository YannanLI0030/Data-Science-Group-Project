# Dynamic unlabelled ablation

## What this run tests

The dynamic ablation uses two 100-gene panels for different purposes. The
unfiltered random panel retains the sparsity of the integrated data and is used
as a missing-data stress test. The coverage-stratified panel deliberately
includes Protein, Confidence, and RNA--Protein correlation regimes in which the
corresponding interventions can affect a ranking.

Neither panel has independently verified relevance labels. The reported
quantities describe ranking sensitivity and robustness, not recommendation
accuracy. This run does not calculate NDCG, Recall, or MRR and does not select
a winning weight profile. V4 remains development evidence. V5 is reviewed only
after its configurations and 12--20 genes have been frozen.

## Inputs and data handling

The runner uses:

- panel lists and long tables from `benchmarks/panels/`;
- the matching local `9_DepMap_sample_info.csv`;
- the RNA--Protein correlation table;
- the current dynamic production scorer, used only to check numeric A0 parity.

Repeated `(gene, source, DepMap_ID)` measurements are reduced to the arithmetic
mean of `value_raw` and `value_std`, matching the production adapter. DepMap IDs
that do not appear in the matching sample-information snapshot are removed and
counted in the manifest.

## Configurations

The 14-configuration runner retains B0, A0, A1a, A1b, A2, A3, A3b, A4, A5,
and A6 and adds four controlled ablations:

- `A1c_no_protein_confidence`: retains direct Protein scoring but removes
  Protein-derived Confidence;
- `A2a_no_conf_completeness`;
- `A2b_no_conf_source_support`;
- `A2c_no_conf_consistency`.

When one Confidence component is removed, the remaining internal weights are
renormalised to one. A4 uses the frozen correlation thresholds: below 0.25 is
low RNA trust, 0.25--0.50 is the transition range, and 0.50 or above is full
RNA trust.

## Running the analysis

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

Run `--self-test` after changing the code. Use `--force` only when deliberately
regenerating an existing result set.

## Output files

- `dynamic_ablation_unlabelled_summary.csv`: panel-, stratum-, and
  configuration-level summaries, including `ALL` rows;
- `dynamic_ablation_unlabelled_by_gene_config.csv`: coverage, ties, score
  variance, and A0 comparisons for each panel/gene/configuration combination;
- `dynamic_ablation_unlabelled_gene_diagnostics.csv`: the number of collapsed
  configurations for each gene;
- `dynamic_ablation_unlabelled_top10.csv`: Top-10 candidates with every
  exact-score tie at the cutoff retained, so some groups contain more than ten
  rows;
- `dynamic_ablation_unlabelled_manifest.json`: hashes, formulae, aggregation
  rules, parity checks, and interpretation limits.

## Results

All 200 panel genes matched the current dynamic scorer numerically under A0.

The unfiltered random panel contained Protein evidence for only 19/100 genes.
Across the fourteen configurations, 43/100 genes produced one unchanged full
ranking. A1a and A0 were identical for 81% of full rankings and 85% of Top-10
sets, while A0 Confidence was constant within 67% of genes. These results show
how often an ablation becomes inactive because its input is absent or does not
vary within the gene.

The stratified panel increased Protein availability to 50/100 genes. In the
DepMap-RNA-only stratum, 20 of 25 genes collapsed to one full ranking under every
configuration; no gene in the other strata did so. A1a and A0 were identical
for 50% of full rankings and 63% of Top-10 sets, making direct-Protein removal
observable in the Protein-rich half of the panel.

A4 matched A0 for every high-correlation gene, as required by its definition.
It changed the Top-10 for 70% of intermediate-correlation genes and for all
low-correlation genes. Mean Top-10 Jaccard similarity to A0 was 0.741 and 0.556
in those two strata. Across the full stratified panel, Top-10 identity to A0
was 66% after removing completeness, 62% after removing source support, and
78% after removing consistency.

These figures show where each intervention changes the ranking. Without Gold
labels, the direction of that change cannot be classified as helpful or
harmful. Four-decimal score rounding can also split or merge ties under a
monotonic transformation. Top-10-with-ties, Jaccard similarity, and rank
correlation should be read together instead of relying on full-order identity
alone.

## Next step

The run is complete once each intervention has been shown to act in its intended
stratum and remain benign when its required modality is absent. Configuration
definitions can then be frozen. Only after that freeze should
`candidate_pool_v5_holdout_review.csv` be generated and independently labelled.
