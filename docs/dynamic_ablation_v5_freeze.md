# Dynamic ablation analysis and V5 design freeze

## 1. What the four supplementary configurations mean

The dynamic runner contains the previous ten configurations plus four
supplementary controlled ablations, for fourteen configurations in total.  The
new configurations are not a weight search and are not four proposed production
models.

| Configuration | Change relative to A0 | Controlled question |
|---|---|---|
| `A1c_no_protein_confidence` | Keep direct Protein scoring; remove Protein-derived completeness, source-support denominator, and RNA-Protein consistency | Does Protein affect ranking through Confidence even when direct Protein scoring is retained? |
| `A2a_no_conf_completeness` | Set completeness weight to zero; renormalise source support and consistency | Does the composition of Confidence depend on completeness? |
| `A2b_no_conf_source_support` | Set source-support weight to zero; renormalise completeness and consistency | Does the composition of Confidence depend on source support? |
| `A2c_no_conf_consistency` | Set RNA-Protein consistency weight to zero; renormalise completeness and source support | Does the composition of Confidence depend on consistency? |

Together, `A0`, `A1a`, `A1b`, and `A1c` distinguish direct Protein scoring
from Protein-derived Confidence.  `A2`, `A2a`, `A2b`, and `A2c` distinguish the
whole Confidence block from its internal components.  Removing an internal
component necessarily changes the relative weights of the retained components;
the retained weights are therefore renormalised to sum to one.

## 2. What the unlabelled dynamic run establishes

All 200 gene runs (100 old random-panel genes plus 100 coverage-stratified
genes) passed numerical parity between `A0_team_baseline` and the current
production scorer.

### Old unfiltered random panel: missing-data stress test

- Only 19/100 genes had Protein evidence.
- 43/100 genes produced the same full ranking under all fourteen configurations.
- `A1a` and `A0` had identical full rankings for 81% of genes and identical
  Top-10 sets for 85%.
- A0 Confidence was constant within 67% of genes.

This is not evidence that the removed components are unimportant.  It shows
that an intervention cannot act when the required modality is absent or when a
Confidence component is constant within a gene.

### Coverage-stratified panel: component identifiability

- Protein evidence was present for 50/100 genes by design.
- `A1a` and `A0` had identical full rankings for 50% of genes and identical
  Top-10 sets for 63%; direct Protein removal is therefore observable in the
  protein-rich half of the panel.
- In the high-correlation stratum, `A4` and `A0` had identical Top-10 sets for
  100% of genes, as required by the adaptive-trust definition.
- In the mid-correlation and low-correlation strata, the corresponding identity
  rates were 30% and 0%; the adaptive intervention is active where intended.
- Removing completeness, source support, or consistency produced A0-identical
  Top-10 sets for 66%, 62%, and 78% of all stratified genes, respectively.

These figures establish structural influence, not biological correctness.
There were no independently verified V5 labels in this run, so NDCG, Recall,
MRR, or a winning model cannot be inferred from rank movement alone.

## 3. Frozen configuration scope for V5

The V5 candidate union uses ten primary controlled configurations:

1. `B0_rna_mean`
2. `A0_team_baseline`
3. `A1a_no_direct_protein`
4. `A1b_no_protein_evidence`
5. `A1c_no_protein_confidence`
6. `A2_no_confidence`
7. `A2a_no_conf_completeness`
8. `A2b_no_conf_source_support`
9. `A2c_no_conf_consistency`
10. `A4_adaptive_trust`

Four configurations are also frozen, but remain secondary unlabelled
sensitivity analyses: `A3_equal_bio`, `A3b_reversed_bio`, `A5_v3_structure`,
and `A6_v3_full`.  They change biological weights or score architecture rather
than removing one component.  Excluding their candidates from the V5 union
keeps the labelled experiment focused and limits manual-review growth.  They
must not receive V5 accuracy claims because their complete Top-10 unions are
not included in the reviewer pool.

The machine-readable freeze is
`config/ablation_v5_frozen_design.json`.  A configuration change after V5
labels are inspected converts the data from holdout evidence into development
evidence and requires a new version.

## 4. Preselected V5 genes

Twelve genes were chosen without V5 candidate labels.  Selection used coverage
stratum balance, human interpretability, and unlabelled structural
discriminability.  This is a targeted component holdout, not a random sample of
all human genes.

| Stratum | Genes | Candidate rows in frozen union |
|---|---|---:|
| Low reliable RNA-Protein correlation | MFN2, CALML3 | 37 |
| Mid reliable RNA-Protein correlation | OGDH, PSME4 | 37 |
| High reliable RNA-Protein correlation | ERBB3, CHKA, ARID3A | 44 |
| Protein present, no reliable correlation | GALNT7, SIGLEC9, CALCB | 48 |
| Multiple RNA sources, no Protein | POU2F3, KRT76 | 23 |
| **Total** | **12 genes** | **189 rows** |

The DepMap-RNA-only stratum remains part of the missing-data stress test and is
not sent to literature review: configurations mostly collapse there, and the
selected non-coding or pseudogene-like entries would make external cell-line
verification disproportionately difficult.

## 5. Candidate-pool separation and blinding

`benchmarks/candidate_pool_v5_holdout_review.csv` is independent of the V4
development pool.  It is the union of Top-10 candidates, including exact-score
ties, from all ten frozen primary configurations for the twelve selected genes.

The reviewer-facing file contains 189 unique gene-DepMap pairs.  It deliberately
omits configuration names, retrieval counts, ranks, scores, and coverage strata.
Rows are grouped by frozen gene order but deterministically shuffled within each
gene.  Every row starts with `judgement=unknown` and `verified=no`.

Retrieval provenance is stored separately in
`results/v5_holdout_internal/candidate_pool_v5_holdout_audit.csv`.  This file
must not be consulted during labelling.  It should be opened only after labels
are frozen, for metric computation and error analysis.

The candidate query is gene-only, across all cell lines in the matching
sample-info snapshot, with no disease filter and no exclusion gene.  Therefore,
V5 conclusions apply to this ablation query mode and must not be presented as a
validation of disease-context, exclusion-gene, protein-only, or combined-query
workflows.

## 6. Review protocol and next gate

For every row, fill the same seven fields used previously:

- `judgement`: `positive`, `negative`, or `unknown`;
- `benchmark_task`;
- `evidence_type`;
- `source_url`;
- `evidence_summary`;
- `verified`: `yes` only when the external evidence supports the judgement;
- `review_notes`.

Unknown is not a negative.  Dataset values and model ranks are not external
ground truth.  Identity evidence alone can verify the cell-line identity, but
does not establish positive or negative target expression.

After all 189 rows are reviewed, the next gate is to validate fields and
duplicates, freeze the review-file hash, export a separate
`gold_standard_v5_holdout.csv`, and evaluate the ten primary configurations
once.  If the results are used to alter configurations, V5 must be reported as
development rather than holdout evidence.

## 7. Reproduction and verification

```bash
cd "/Users/liyannan/Desktop/Data-Science-Group-Project"

/Users/liyannan/miniconda3/envs/cellline-v3/bin/python \
  scripts/build_candidate_pool_v5_holdout.py \
  --sample-info "/Users/liyannan/Desktop/cellline_selector_v3/data_s3/nomenclature/9_DepMap_sample_info.csv"
```

The builder refuses to overwrite existing outputs unless `--force` is used,
and refuses even with `--force` once manual review content is detected.  Verify
the frozen hashes without regenerating files using:

```bash
/Users/liyannan/miniconda3/envs/cellline-v3/bin/python \
  scripts/build_candidate_pool_v5_holdout.py \
  --sample-info "/Users/liyannan/Desktop/cellline_selector_v3/data_s3/nomenclature/9_DepMap_sample_info.csv" \
  --verify-only
```
