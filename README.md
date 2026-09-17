# CellLineSelector

MSc Data Science group project: an auditable multi-omics recommendation system
for selecting human cell lines by target gene, optional exclusion gene and
optional disease/tissue context.

## Current scoring deliverable

The production scoring entry point is `src/merged_cellline_selector.py`. It
reads `master_table.csv` and `cellline_annotations.csv`, ranks candidates and
writes:

- a complete ranked CSV for audit and evaluation;
- a structured JSON payload for the UI/RAG integration;
- score components, confidence components, evidence flags, mutation/fusion
  annotations, data gaps and alternative cell lines.

The retained production configuration, A0, assigns 0.55 to RNA, 0.30 to
Protein and 0.15 to Confidence. Confidence combines completeness (0.40),
source support (0.35) and RNA--Protein consistency (0.25). The configuration
file and internal name still contain `provisional` for compatibility with
earlier outputs. The frozen V5 comparison did not support replacing A0 with
A1a, but it also did not establish A0 as a universally optimal weight set.

## Setup

Python 3.10 or newer is recommended.

```bash
python -m pip install -r requirements.txt
```

The committed tables and reports can be inspected without downloading the raw
data. Re-running the full experiments also requires the local integrated data,
the RNA--Protein correlation table and, for the Penalty experiment, the latest
dynamic production scorer and caches. AWS is not required.

The merged data is intentionally not committed. Either pass its directory on
each run or set an environment variable:

```bash
export CELLLINESELECTOR_DATA_DIR="/absolute/path/to/merged"
```

## Run a recommendation

Pan-cancer, non-interactive example:

```bash
python src/merged_cellline_selector.py \
  --target_gene EGFR \
  --no_prompt
```

Disease-filtered example with an exclusion gene:

```bash
python src/merged_cellline_selector.py \
  --target_gene EGFR \
  --exclusion_gene ABCB1 \
  --disease lung \
  --data_dir "/absolute/path/to/merged" \
  --no_prompt
```

To test a future selected weight set without editing Python code, copy and edit
`config/scoring_a0_provisional.json`, then add:

```bash
--scoring_config config/my_selected_weights.json
```

Relative output directories are resolved from the current working directory.
The default `results/` output contains a full ranking CSV and a structured JSON
file whose names include the target and query context.

## Verify scoring and run ablations

```bash
python scripts/ablation_runner.py --self-test
```

This test verifies, among other invariants, that production A0 and ablation A0
produce the same scores and ordering.

## Evaluation record

The experiments answer different questions and should not be merged into one
performance claim:

| Stage | Purpose | Main boundary |
|---|---|---|
| V4 development | Generate component and weight hypotheses from 50 labels across 10 genes | Model-informed pool; only 2 verified negatives |
| Dynamic ablation | Check whether 14 interventions alter rankings across two 100-gene panels | No relevance labels; structural sensitivity only |
| Frozen V5 | Compare 10 configurations on 12 reviewed genes | 99 positives and 90 unknowns; no verified negatives |
| Penalty V1 | Test P00/P15/P30/P45 mechanics on 5 fixed queries | No pair-specific Gold; no optimal-cap claim |

V5 retained A0 because the V4 advantage of A1a did not reproduce. Removing the
full Confidence block or source support caused clearer NDCG@5 reductions, while
other variants showed endpoint- or gene-specific trade-offs. The Penalty run
passed all 155 implementation checks and showed that P30 changes rankings, but
it did not show that 0.30 is biologically optimal.

Detailed records are in:

- `docs/dynamic_ablation_unlabelled.md`;
- `docs/dynamic_ablation_v5_freeze.md`;
- `docs/v5_holdout_evaluation.md`;
- `docs/exclusion_penalty_structural_v1.md`;
- `docs/weight_search_v4_development.md`;
- `docs/experiment_provenance.md`.

## Scope and limitations

- The current merged master table contains a selected gene panel, not all human
  genes.
- Protein coverage is sparse; missing protein changes completeness and
  confidence and triggers row-wise biological-weight redistribution.
- Mutation and fusion are auditable annotations, not universally positive or
  negative ranking signals.
- The current alternative-line similarity is a target RNA/protein component
  fallback, not full-transcriptome or full-multi-omics similarity.
- V5 evaluates early ranking within a frozen candidate union; it is not an
  open-world estimate across all genes and cell lines.
- V5 contains no verified negatives, so it cannot support a negative-sink
  conclusion.
- The archived V4 weight search is development-stage model selection, not a
  completed optimisation of production weights.

See `docs/scoring_module_handoff.md` for the scoring formula and output
contract. Later evaluation decisions are recorded in the V5 and Penalty
documents listed above.
