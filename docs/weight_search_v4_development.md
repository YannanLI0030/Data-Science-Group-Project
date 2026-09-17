# V4 development weight search

## Scope

This package preserves the constrained weight search run on the V4 development
benchmark. It records the model-selection work carried out during development;
it is separate from the later label-free sensitivity analysis and does not
define the final production weights.

Only the outer RNA, direct-Protein, and Confidence weights were varied. The
following settings remained fixed:

- Confidence composition: completeness 0.40, source support 0.35, and
  RNA--Protein consistency 0.25;
- maximum exclusion penalty: 0.30;
- adaptive trust: disabled.

The predeclared 0.05 grid contained 32 configurations. RNA ranged from
0.50--0.90, direct Protein from 0.00--0.30, and Confidence from 0.10--0.30;
the three weights always summed to one. A direct-Protein weight of zero used
`direct_off`: Protein was removed from the biological score but could still
contribute to Protein-derived Confidence.

## Benchmark and selection rule

`benchmarks/gold_standard_v4_development.csv` contains 50 verified labels
across ten genes: 48 positive and two negative. Most labels came from
model-generated candidate pools, and the negative sample is small. The data
support development comparisons but not an untouched test claim.

Mean NDCG@5 was the primary metric. The summaries use 5,000 paired bootstrap
resamples with seed 42. In each leave-one-gene-out (LOGO) fold, the
one-standard-error rule first identified configurations within one standard
error of the best training NDCG@5. It then chose the configuration closest to
the predeclared A1a anchor: RNA 0.85, direct Protein 0.00, and Confidence 0.15.

## Result

The unconstrained highest mean NDCG@5 came from RNA 0.70, direct Protein 0.00,
and Confidence 0.30. Its NDCG@5 was 0.5385 and its Recall@10 was 0.7656. The
one-standard-error rule returned the A1a anchor, which reached NDCG@5 0.5027
and Recall@10 0.8481. For comparison, the production A0 setting (RNA 0.55,
Protein 0.30, Confidence 0.15) reached NDCG@5 0.4595 and Recall@10 0.8772.

All 32 configurations were inside the one-standard-error band. The conservative
anchor rule, rather than a clear separation in NDCG@5, selected A1a. Its JSON
status is
`candidate_not_yet_frozen`; the file is an experiment record, not a
production-ready configuration.

The later frozen V5 comparison did not reproduce an advantage from removing
direct Protein, so A0 remained the production configuration. The V5 record is
in `docs/v5_holdout_evaluation.md`. The V4 search should be reported as
incomplete development-stage model selection, not as final weight
optimisation.

## Files

- `scripts/weight_search_v4_logo.py`: grid generation, bootstrap summaries,
  one-standard-error selection, and LOGO evaluation;
- `results/weight_search_v4_development_summary.csv`: summaries for 32
  configurations;
- `results/weight_search_v4_development_by_gene.csv`: 320
  configuration--gene records;
- `results/weight_search_v4_development_logo.csv`: ten held-out-gene folds;
- `config/scoring_selected_v4_development_candidate.json`: archival candidate
  metadata;
- `results/weight_search_v4_development_manifest.json`: inputs, parameters,
  hashes, and integrity counts.

## Reproduction

Run the self-test first:

```bash
python scripts/weight_search_v4_logo.py --self-test
```

Run `python scripts/weight_search_v4_logo.py --help` for the full input and
output options. Reproduction runs should use a new output directory so that
the archived files remain unchanged.

The integrated data and correlation table are external inputs and are not
duplicated in this repository. The manifest records their hashes for the
archived run.
