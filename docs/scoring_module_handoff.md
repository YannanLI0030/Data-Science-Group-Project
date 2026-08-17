# Scoring module handoff

Status: interface frozen for team integration; default weight values are
provisional pending benchmark expansion and weight search.

## 1. Supported query

Required input:

- `target_gene`: a gene present in the merged gene panel.

Optional input:

- `exclusion_gene`: applies an RNA-expression penalty when the gene is in the
  panel; if absent from the panel, the program warns and applies no penalty.
- `disease`: a hard disease/lineage filter; when omitted, the query is
  pan-cancer.
- `scoring_config`: JSON file matching `ScoringConfig`.

The program requires a directory containing `master_table.csv` and
`cellline_annotations.csv`, supplied by `--data_dir` or
`CELLLINESELECTOR_DATA_DIR`.

## 2. Provisional A0 formula

For candidate cell line `i`:

```text
RNA_i = min-max scaled mean of available standardized DepMap/HPA/GEO RNA
Protein_i = min-max scaled standardized CCLE-Gygi protein, when present

Confidence_i =
    0.40 * completeness_i
  + 0.35 * source_support_i
  + 0.25 * RNA_protein_consistency_i

Biology_i = weighted mean of available RNA/protein layers

Final_i = clip(
    0.85 * Biology_i
  + 0.15 * Confidence_i
  - exclusion_penalty_i,
  0, 1
)
```

Within `Biology`, RNA/protein relative weights are 0.55/0.30. If protein is
missing, the available biological weight is redistributed to RNA; confidence
still records the missing layer through completeness, source support and the
single-layer consistency convention. The maximum exclusion penalty is 0.30.

All seven numeric values are validated as non-negative; outer weights and
confidence-component weights must each sum to 1.0.

## 3. Output contract (schema 1.0)

Each query writes two complementary artifacts:

1. `*_ranked_recommendations.csv`: the full ordered candidate table.
2. `*_recommendations.json`: the stable handoff payload for API/UI/RAG.

Top-level JSON keys:

| Key | Meaning |
| --- | --- |
| `schemaVersion` | Output-contract version. |
| `query` | Target, exclusion, disease filter and pan-cancer/filtered scope. |
| `scoringConfig` | Exact weights and configuration name used for the run. |
| `recommendations` | Top-N ranked rows with score breakdown and evidence flags. |
| `alternatives` | Similar alternatives to the top recommendation. |
| `topRecommendationEvidenceTrace` | Dataset-level evidence for the top row. |
| `artifacts` | Name of the complete ranking CSV. |
| `knownLimitations` | Machine-readable caveats for downstream reporting. |

The ranked rows expose the exact scoring config name and weights, stable
identifiers (`DepMap_ID`), display name, lineage, disease, RNA/protein and
confidence components, final score, recommendation level, source-presence
flags, raw source values, mutation/fusion annotations, assay-readiness score
and risk flags. This keeps the CSV auditable even if it is separated from its
JSON payload.

## 4. Reproducible checks

From the repository root:

```bash
python scripts/ablation_runner.py --self-test
```

The self-test covers production-A0/ablation-A0 parity, adaptive-trust branches,
complete protein removal, missing gold positives and exact-score ties.

Real-data ablation command:

```bash
python scripts/ablation_runner.py \
  --data-dir "/absolute/path/to/merged" \
  --gold benchmarks/gold_standard_v2.csv \
  --corr "/absolute/path/to/gene_rna_protein_correlations.csv" \
  --out results/ablation_results_v2.csv \
  --detail-out results/ablation_results_v2_by_gene.csv \
  --bootstrap 5000 \
  --seed 42
```

## 5. Handoff boundary

Ready for teammates now:

- CLI and callable `score_candidates(rows, config=...)` interface;
- validated external JSON weight configuration;
- deterministic ranking and production/ablation parity guard;
- pan-cancer or disease-filtered queries;
- full CSV plus structured JSON output;
- explicit evidence, confidence and limitation fields.

Not yet a defensible final model claim:

- final numerical weight selection;
- superiority of A1a, A1b or A4 on the current 6-gene/18-positive benchmark;
- full-omics cell-line similarity;
- universal interpretation of mutation/fusion directionality.

## 6. Next priority before the report results section

1. Build the union of Top-10 candidates from the main configurations for each
   benchmark gene.
2. Label candidates as positive, negative or unknown; do not treat all
   unlisted candidates as negatives.
3. Verify/add missing known positives, beginning with MET/MKN45 and then the
   NUGC4/SNU5 review list.
4. Add more genes, particularly low RNA-protein-correlation genes needed to
   test adaptive trust.
5. Run a constrained weight search with held-out or leave-one-gene-out
   evaluation, then save the selected values as a new JSON config without
   changing the production interface.
