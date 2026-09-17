# Exclusion-penalty structural ablation V1

## Purpose

This supplementary experiment checks the exclusion-penalty implementation and
measures how different caps change a fixed ranking. There is no pair-specific
Gold standard, so the analysis does not estimate NDCG, Recall, MRR, or
biological accuracy. In particular, it cannot identify 0.30 as an optimal cap.

The experiment was designed after the V5 target-gene holdout had been
completed. V5 labels, rankings, and internal audit fields were not used to
choose the queries or calculate their scores. The runner reads only the V5
design-file hash as an integrity check.

## Frozen design

The design and query definitions are stored in
`config/exclusion_penalty_structural_v1_frozen_design.json` and
`benchmarks/penalty/penalty_structural_v1_queries.csv`.

Every profile uses the production target-scoring weights: RNA 0.55, Protein
0.30, and Confidence 0.15. Confidence retains completeness 0.40, source support
0.35, and RNA--Protein consistency 0.25. Only the maximum penalty changes.

| Profile | Maximum penalty | Role |
|---|---:|---|
| `P00_no_penalty` | 0.00 | Negative control |
| `P15_weak` | 0.15 | Weaker penalty |
| `P30_team_baseline` | 0.30 | Production setting |
| `P45_stress` | 0.45 | Stress condition |

Five disease-context queries were selected from genes already present in the
local cache. They were frozen before the complete P00/P15/P30/P45 output was
generated. Selection considered cache availability, candidate count, exclusion
RNA coverage, and variation in disease context; it did not use labels or
profile-comparison results.

| Query | Target | Exclusion | Disease | Candidates | Exclusion RNA observed |
|---|---|---|---|---:|---:|
| Q01 | FGFR2 | MET | Gastric | 43 | 43 |
| Q02 | ERBB2 | ESR1 | Breast | 73 | 73 |
| Q03 | CD3E | CD86 | Blood | 239 | 239 |
| Q04 | ASGR1 | AFP | Liver | 25 | 25 |
| Q05 | MSLN | MUC1 | Lung | 236 | 236 |

This is a purposive mechanics panel, not a random or representative sample.
The pairs are not claims that MET, ESR1, CD86, AFP, or MUC1 is undesirable in
every experiment involving the corresponding target. Familiar genes can be
reused here because the endpoint is score behaviour rather than held-out
accuracy.

## Scoring behaviour

For an exclusion gene, the production path averages the available standardised
DepMap, HPA, and GEO RNA values and applies Min--Max scaling within the
disease-filtered candidate set. The score is

\[
\operatorname{Final}_\alpha =
\operatorname{clip}\left(
0.85\,\operatorname{Biology}
+0.15\,\operatorname{Confidence}
-\alpha\,\operatorname{ExclusionScaled},
0,1
\right).
\]

`ExclusionScaled` is relative to the current query. A value of 0.8 in one
disease context is not an absolute threshold and cannot be compared directly
with 0.8 from another context.

When exclusion RNA is missing, the production implementation applies zero
penalty. The runner keeps that rule but records the value as missing rather
than treating it as confirmed low expression.

Each query is fetched once. The four profiles then score copies of the same
candidate rows. The cap is changed sequentially in memory and restored after
every scorer call; no production source or cache file is edited.

## Implementation audit

The run produced 155 substantive checks, all marked `PASS`. They cover:

- identical candidate IDs and unchanged RNA, Protein, Biology, and Confidence
  components across the four profiles;
- zero Penalty for every P00 row;
- Penalty equal to the rounded cap multiplied by scaled exclusion;
- P30 equal to twice P15 and P45 equal to three times P15 within the required
  four-decimal rounding tolerance;
- non-decreasing candidate Penalty and non-increasing candidate final score as
  the cap rises;
- score loss equal to the reported Penalty unless clipping at zero limits the
  realised change;
- exact P30 parity with the unmodified production scorer for all 5 queries;
- deterministic repeated runs;
- unchanged hashes for the production scorer, team data modules, sample
  information, frozen designs, query table, and ten gene caches.

All 5 real queries have exclusion RNA for every candidate, so the real-data
runs do not exercise missing exclusion evidence. Synthetic checks cover all-missing
and partly missing values, lower-bound clipping, constant positive, zero, and
negative values, and the `GENE_MULTIOMICS`, `PROTEIN_ONLY`, and `COMBINED`
subtraction paths. The query-mode checks inspect both the reported Penalty and
the subtraction from the final score. A separate cutoff case verifies that
exact final-score ties across positions 10 and 11 are both retained.

Assertions that do not apply to a real query, such as checking for missing
values when none are missing, are not counted as automatic passes. One
synthetic case also records an existing Min--Max edge condition: if all
observed exclusion values are the same positive number, the production scaler
assigns the full cap to every observed row. The experiment documents this
behaviour but does not change it.

## Structural results

The table reports macro means over the 5 frozen queries. Top-10 sets retain
every exact four-decimal final-score tie at the cutoff.

| Profile | Mean Top-10 Jaccard vs P00 | Queries with Top-10 change | Queries with Top-1 change | Mean Spearman vs P00 | Mean absolute rank change | Zero-score rows | Mean Top-10 scaled exclusion |
|---|---:|---:|---:|---:|---:|---:|---:|
| P00 | 1.0000 | 0/5 | 0/5 | 1.0000 | 0.00 | 0 | 0.4021 |
| P15 | 0.7463 | 4/5 | 0/5 | 0.9475 | 9.46 | 7 | 0.3160 |
| P30 | 0.7160 | 4/5 | 3/5 | 0.8778 | 15.20 | 76 | 0.3061 |
| P45 | 0.5880 | 5/5 | 3/5 | 0.8020 | 19.17 | 140 | 0.2457 |

Stronger caps reduce scaled exclusion among the highest-ranked candidates, but
they also create more rank movement and clipping. P15 changes four Top-10 sets
without changing any of the 5 Top-1 results. P30 changes four Top-10 sets and
three Top-1 results. P45 changes all 5 Top-10 sets and produces the largest
disruption.

Using scaled exclusion of at least 0.8 as a structural definition, the mean
high-exclusion share in the Top-10 falls from 0.14 under P00 to 0.06 under P30
and 0.02 under P45. All 5 queries have complete observed coverage, so the
denominator is the full Top-10. The 0.8 cutoff is query-relative and is not a
biological high-expression threshold.

## P30 by query

| Query | Top-10 Jaccard vs P00 | P00 Top-1 | P30 Top-1 | Mean absolute rank change | Maximum rank change | P30 zero-score rows |
|---|---:|---|---|---:|---:|---:|
| FGFR2--MET, gastric | 0.6667 | SNU16 | SNU16 | 3.40 | 16.0 | 10 |
| ERBB2--ESR1, breast | 0.8182 | EFM192A | SUM190PT | 6.38 | 23.0 | 1 |
| CD3E--CD86, blood | 0.6667 | MOLT3 | MOLT3 | 44.08 | 128.0 | 39 |
| ASGR1--AFP, liver | 1.0000 | JHH5 | HUH1 | 1.48 | 4.0 | 0 |
| MSLN--MUC1, lung | 0.4286 | NCIH322M | NCIH2052 | 20.66 | 114.5 | 26 |

The effect varies substantially by query. CD3E--CD86 and MSLN--MUC1 have the
largest rank movements and the most clipped rows. ASGR1--AFP keeps the same
Top-10 membership at P30, although its internal order and Top-1 change. Since
exclusion expression is scaled within each candidate set, a single cap can be
mild for one query and disruptive for another.

Clipping also creates large tie groups. The largest final-score group contains
39 candidates under P30 and 71 under P45. The production sorter uses
`(finalScore, confidenceScore)` as its full comparator; under that comparator,
the largest groups contain 17 and 33 candidates. Confidence can still order
candidates with equal final scores. None of the 5 real queries has
a tie at the Top-10 cutoff, but the runner keeps tie-aware selection for future
query sets.

## Interpretation

The audit supports the following restricted conclusions. The Penalty is
subtracted correctly and reproducibly in the tested real and synthetic paths.
Increasing the cap lowers candidate scores according to scaled exclusion and
reduces high-exclusion representation near the top of the ranking. Individual
ranks need not fall monotonically because every candidate moves relative to the
others and clipping creates ties. P30 has a material, query-dependent effect;
P45 is useful as a stress condition but is not supported as a production
candidate.

The analysis does not show that any reordered candidate is better for a real
experiment. Five purposively chosen queries cannot represent all targets or
diseases, and the current data cannot determine whether 0.15, 0.30, or another
cap gives the best trade-off between retaining target-positive lines and
suppressing an exclusion marker. That comparison requires pair-specific
evidence for both target suitability and exclusion expression.

## Reproduction

Run the synthetic checks:

```bash
python scripts/exclusion_penalty_ablation_runner.py --self-test
```

Run the frozen structural experiment:

```bash
python scripts/exclusion_penalty_ablation_runner.py
```

Existing results are protected from replacement. Use `--force` only to
reproduce the same frozen design intentionally. The 5 CSV files and JSON
manifest are written to `results/exclusion_penalty_structural_v1/`. The same
directory contains `exclusion_penalty_structural_v1.xlsx` as a formatted review
copy; the CSV files and manifest remain the machine-readable record.

## Possible labelled extension

A later pair-specific pool could be drawn from the frozen P00/P15/P30/P45
results. Reviewers would need to assess both target suitability and whether the
exclusion marker is high or low. Existing V5 target-only labels cannot answer
that question. If the new labels are used to revise the cap, the pair-level
benchmark becomes development evidence rather than an untouched holdout.
