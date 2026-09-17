# Frozen V5 holdout evaluation

V5 is the one-time labelled comparison run after the blinded reviewer file was
completed and frozen.

## Frozen inputs

- `benchmarks/candidate_pool_v5_holdout_review_completed.csv` contains 189
  reviewed gene--DepMap-ID pairs: 99 positive and 90 unknown, with no verified
  negative.
- `benchmarks/gold_standard_v5_holdout.csv` contains the 99 verified positives
  used in the evaluation.
- `benchmarks/gold_standard_v5_holdout_manifest.json` records the label-freeze
  checks, hashes, per-gene counts, and interpretation limits.
- `results/v5_holdout_internal/candidate_pool_v5_holdout_audit.csv` contains
  the frozen rankings and configuration provenance. It was opened only after
  reviewer labels had been frozen.

The untouched reviewer template remains at
`benchmarks/candidate_pool_v5_holdout_review.csv`. It is the pre-review record
and should not be replaced by the completed file.

## Results and provenance

`results/v5_holdout_evaluation/` contains the macro summary, per-gene metrics,
controlled contrasts, exact-score tie audit, evaluation manifest, and a
post-evaluation provenance clarification. The manifest is the reference for
input and output hashes, bootstrap settings, metric definitions, and permitted
claims.

The two scripts under `scripts/v5/` are the exact copies used to freeze labels
and run the evaluation. Their original local paths are retained so that their
hashes continue to match the provenance record. They are archival evidence,
not portable command-line tools.

## How to read V5

V5 measures early ranking within the candidate union produced by ten frozen
configurations. Unknown candidates remain unlabelled rather than being counted
as negative. Precision is a verified lower bound, and negative-sink behaviour
cannot be evaluated.

The comparison is independent of the V4 model-selection labels but remains
targeted to the frozen union. If V5 is later used to redesign a configuration,
it becomes development evidence for that change. Confirmation would then
require another untouched benchmark.
