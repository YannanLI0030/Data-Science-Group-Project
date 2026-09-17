# V4 development-stage weight search

## Status and purpose

This package records the constrained weight-search experiment performed on the
V4 development benchmark. It is retained as an auditable experiment and as
evidence of the model-selection work carried out during development. It is not
the team's later label-free sensitivity analysis, and it did not establish a
final or universally optimal set of weights.

The experiment varied only the outer RNA, direct-Protein, and Confidence
weights. Confidence subweights (completeness 0.40, source support 0.35, and
RNA--Protein consistency 0.25), the exclusion-penalty cap (0.30), and adaptive
trust remained fixed. The declared 0.05 grid contained 32 configurations:

- RNA weight: 0.50--0.90;
- direct-Protein weight: 0.00--0.30;
- Confidence weight: 0.10--0.30;
- the three outer weights always summed to one.

A zero direct-Protein weight used the `direct_off` interpretation: Protein did
not enter the biological score directly but could still contribute to
Protein-derived Confidence features.

## Data and selection procedure

The search used `benchmarks/gold_standard_v4_development.csv`: 50 verified
labels across ten genes, comprising 48 positives and two negatives. Because
most labels came from model-generated candidate pools and the negative sample
was very small, this is a development benchmark rather than an untouched test
set.

Mean NDCG@5 was the primary metric. Each configuration was evaluated for all
ten genes; descriptive paired intervals used 5,000 bootstrap resamples with
seed 42. Leave-one-gene-out (LOGO) folds applied a one-standard-error rule. The
rule selected, among configurations within one standard error of the best
training NDCG@5, the configuration closest to the predeclared A1a anchor
(RNA 0.85, direct Protein 0.00, Confidence 0.15).

## Main result and evidence boundary

The unconstrained highest mean NDCG@5 was obtained by RNA 0.70, direct Protein
0.00, and Confidence 0.30 (NDCG@5 0.5385), but its Recall@10 was 0.7656. The
one-standard-error procedure returned the A1a anchor (NDCG@5 0.5027 and
Recall@10 0.8481). The production A0 reference (RNA 0.55, Protein 0.30,
Confidence 0.15) achieved NDCG@5 0.4595 and Recall@10 0.8772.

All 32 configurations lay inside the one-standard-error band. The returned
candidate was therefore determined by the conservative anchor rule, not by
clear empirical separation between weight settings. The candidate JSON is
explicitly marked `candidate_not_yet_frozen` and is not a production
configuration file.

The later frozen V5 evaluation did not reproduce a benefit from removing
direct Protein and retained A0 as the production configuration. See
`docs/v5_holdout_evaluation.md`. Consequently, this V4 package should be
described as incomplete development-stage weight-search evidence rather than
as final weight optimisation.

## Files

- `scripts/weight_search_v4_logo.py`: declared grid, bootstrap summaries,
  one-standard-error selection, and LOGO evaluation;
- `results/weight_search_v4_development_summary.csv`: 32 configuration-level
  summaries;
- `results/weight_search_v4_development_by_gene.csv`: 320 configuration--gene
  records;
- `results/weight_search_v4_development_logo.csv`: ten held-out-gene folds;
- `config/scoring_selected_v4_development_candidate.json`: archival candidate
  metadata, not a production-ready configuration;
- `results/weight_search_v4_development_manifest.json`: file hashes, inputs,
  parameters, and integrity counts.

## Reproduction

Run the self-test first:

```bash
python scripts/weight_search_v4_logo.py --self-test
```

To avoid replacing the archived outputs, write a reproduction run to a new
directory:

```bash
python scripts/weight_search_v4_logo.py \
  --data-dir "/path/to/merged" \
  --gold benchmarks/gold_standard_v4_development.csv \
  --corr "/path/to/gene_rna_protein_correlations.csv" \
  --summary-out results/weight_search_v4_reproduction/summary.csv \
  --detail-out results/weight_search_v4_reproduction/by_gene.csv \
  --logo-out results/weight_search_v4_reproduction/logo.csv \
  --config-out results/weight_search_v4_reproduction/candidate_config.json \
  --bootstrap 5000 \
  --seed 42
```

The integrated data and correlation table are external project inputs and are
not duplicated in this repository. Their hashes for the archived run are
recorded in the manifest.
