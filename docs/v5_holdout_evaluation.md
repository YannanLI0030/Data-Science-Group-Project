# Frozen V5 holdout evaluation

This directory set records the one-time V5 component comparison performed
after the blinded reviewer labels were frozen.

## Frozen inputs

- `benchmarks/candidate_pool_v5_holdout_review_completed.csv`: 189 reviewed
  gene--DepMap-ID candidates (99 positive, 90 unknown, no verified negative).
- `benchmarks/gold_standard_v5_holdout.csv`: the 99 verified positive labels
  used by the evaluation.
- `benchmarks/gold_standard_v5_holdout_manifest.json`: label-freeze checks,
  hashes, per-gene counts, and interpretation limits.
- `results/v5_holdout_internal/candidate_pool_v5_holdout_audit.csv`: the
  frozen rankings and configuration provenance opened only after labels were
  frozen.

The original blank reviewer pool remains at
`benchmarks/candidate_pool_v5_holdout_review.csv`; it is retained as the
pre-review baseline and should not be replaced by the completed reviewer file.

## Evaluation outputs

`results/v5_holdout_evaluation/` contains the macro summary, per-gene results,
controlled contrasts, exact-score tie audit, evaluation manifest, and the
post-evaluation provenance clarification. The manifest is the source of truth
for input and output hashes, metric definitions, bootstrap settings, and
claim boundaries.

The two files in `scripts/v5/` are archival copies of the exact scripts used
to freeze the labels and calculate the one-time evaluation. They retain the
original local paths so that their hashes continue to match the provenance
record; they are included for auditability rather than as portable command-line
entry points.

## Interpretation boundary

V5 compares early ranking within the frozen ten-configuration candidate union.
Unknown candidates are unlabelled rather than negative, so precision is a
verified lower bound and negative-sink behaviour is unavailable. If these
results motivate a later configuration change, V5 must be treated as
development evidence and the change confirmed on a new untouched benchmark.
