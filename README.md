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

The default `A0_team_baseline_provisional` configuration is RNA 0.55, protein
0.30 and confidence 0.15. Confidence is completeness 0.40, source support 0.35
and RNA-protein consistency 0.25. These values remain **provisional** until the
gold standard is expanded and the planned weight search is complete.

## Setup

Python 3.10 or newer is recommended.

```bash
python -m pip install -r requirements.txt
```

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
produce the same scores and ordering. A real-data ablation run is documented in
`docs/scoring_module_handoff.md`.

## Scope and limitations

- The current merged master table contains a selected gene panel, not all human
  genes.
- Protein coverage is sparse; missing protein changes completeness and
  confidence and triggers row-wise biological-weight redistribution.
- Mutation and fusion are auditable annotations, not universally positive or
  negative ranking signals.
- The current alternative-line similarity is a target RNA/protein component
  fallback, not full-transcriptome or full-multi-omics similarity.
- The present gold standard is too small to claim that A1a or any other
  configuration is the final optimum.

See `docs/scoring_module_handoff.md` for the scoring formula, output contract,
handoff boundary and the remaining priority-two work.
