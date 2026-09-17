# Dynamic ablation and the V5 design freeze

## Supplementary controlled ablations

The dynamic runner evaluates fourteen configurations: the previous ten plus
four component-level ablations. The four additions are diagnostic comparisons,
not a weight search or proposed production models.

| Configuration | Change from A0 | Question |
|---|---|---|
| `A1c_no_protein_confidence` | Retain direct Protein scoring; remove Protein-derived completeness, the source-support denominator, and RNA--Protein consistency | Can Protein change the ranking through Confidence when direct Protein scoring remains active? |
| `A2a_no_conf_completeness` | Remove completeness; renormalise source support and consistency | How much does completeness contribute to the Confidence ordering? |
| `A2b_no_conf_source_support` | Remove source support; renormalise completeness and consistency | How much does source support contribute to the Confidence ordering? |
| `A2c_no_conf_consistency` | Remove RNA--Protein consistency; renormalise completeness and source support | How much does consistency contribute to the Confidence ordering? |

Together, A0, A1a, A1b, and A1c separate direct Protein scoring from
Protein-derived Confidence. A2, A2a, A2b, and A2c compare the complete
Confidence block with its three components. Removing one component changes the
relative contribution of the other two, so their retained weights are
renormalised to sum to one.

## What was learned before the freeze

All 200 runs, comprising 100 genes from the unfiltered panel and 100 from the
coverage-stratified panel, reproduced the current production A0 scores.

The unfiltered panel behaved as a missing-data stress test. Only 19/100 genes
had Protein evidence, and 43/100 produced the same full ranking under all
fourteen configurations. A1a and A0 matched for 81% of full rankings and 85%
of Top-10 sets. A0 Confidence was constant within 67% of genes. Lack of
movement here mainly reflects unavailable or non-discriminative inputs, not
evidence that the removed component is unnecessary.

The stratified panel made the interventions easier to observe. Protein was
available for 50/100 genes. A1a and A0 matched for 50% of full rankings and
63% of Top-10 sets. A4 matched A0 for 100% of high-correlation genes, which is
the expected behaviour of the adaptive rule. Top-10 identity fell to 30% in the
intermediate-correlation stratum and 0% in the low-correlation stratum. Removing
completeness, source support, or consistency left 66%, 62%, and 78% of Top-10
sets identical to A0 across the full panel.

These are structural results. No independently verified V5 labels were used,
so ranking movement cannot be translated into NDCG, Recall, MRR, or a claim
that one configuration is more accurate.

## Configurations frozen for V5

The labelled comparison uses ten primary configurations:

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

A3, A3b, A5, and A6 are also frozen but remain secondary unlabelled
sensitivity analyses. They change biological weights or score architecture
rather than isolating one component. Their complete Top-10 candidate sets are
not part of the V5 reviewer pool, so V5 accuracy results do not apply to them.

The machine-readable record is `config/ablation_v5_frozen_design.json`. Once
reviewers have seen V5 labels, changing a configuration turns V5 into
development evidence and requires a new version for confirmation.

## Gene selection

Twelve genes were selected before any V5 candidate was labelled. Selection
balanced coverage strata, biological interpretability, and structural
discriminability in the unlabelled run. The result is a targeted component
holdout rather than a random sample of all human genes.

| Stratum | Genes | Candidate rows in frozen union |
|---|---|---:|
| Low reliable RNA--Protein correlation | MFN2, CALML3 | 37 |
| Intermediate reliable RNA--Protein correlation | OGDH, PSME4 | 37 |
| High reliable RNA--Protein correlation | ERBB3, CHKA, ARID3A | 44 |
| Protein present, no reliable correlation | GALNT7, SIGLEC9, CALCB | 48 |
| Multiple RNA sources, no Protein | POU2F3, KRT76 | 23 |
| **Total** | **12 genes** | **189 rows** |

The DepMap-RNA-only stratum remains in the missing-data stress test but was not
sent for literature review. Configurations mostly collapse in this stratum,
and the selected non-coding or pseudogene-like entries would be unusually hard
to verify at cell-line level.

## Candidate pool and blinding

`benchmarks/candidate_pool_v5_holdout_review.csv` is separate from the V4
development pool. It contains the union of Top-10 candidates, including exact
ties at the cutoff, returned by the ten primary configurations for the twelve
genes.

The reviewer file has 189 unique gene--DepMap pairs. It contains no
configuration names, retrieval counts, scores, ranks, or coverage strata. Rows
follow the frozen gene order and are deterministically shuffled within each
gene. Every entry begins with `judgement=unknown` and `verified=no`.

Configuration provenance is stored in
`results/v5_holdout_internal/candidate_pool_v5_holdout_audit.csv`. Reviewers do
not consult this file while assigning labels. It is opened only after the
labels are frozen, when metrics and error analyses are calculated.

The candidate query is gene-only and covers all cell lines in the matching
sample-information snapshot. It has no disease filter or exclusion gene. V5
applies to this ablation query mode, not disease-context,
exclusion-gene, Protein-only, or combined-query workflows.

## Review and freeze procedure

Reviewers complete seven fields for every row:

- `judgement`: `positive`, `negative`, or `unknown`;
- `benchmark_task`;
- `evidence_type`;
- `source_url`;
- `evidence_summary`;
- `verified`: `yes` only when external evidence supports the judgement;
- `review_notes`.

An unknown entry is not treated as negative. Local dataset values and model
ranks are not external Gold evidence. A source that verifies cell-line
identity alone does not verify target expression.

After all 189 rows are reviewed, the fields and duplicate keys are checked,
the reviewer-file hash is frozen, and a separate
`gold_standard_v5_holdout.csv` is exported. The ten primary configurations are
then evaluated once. If that result is used to revise a configuration, V5 must
be described as development rather than holdout evidence.

## Reproduction and verification

From the repository root, run
`python scripts/build_candidate_pool_v5_holdout.py --help` and supply the local
DepMap sample-information file.

The builder will not replace existing outputs unless `--force` is supplied.
Once it detects manual review content, replacement is refused even with
`--force`. Use `--verify-only` with the same sample-information input to check
the frozen hashes without regenerating any file.
