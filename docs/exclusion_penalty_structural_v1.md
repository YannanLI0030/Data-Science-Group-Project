# Exclusion penalty structural ablation V1

## 1. Experimental question and claim boundary

This supplementary experiment tests whether the exclusion penalty is
implemented consistently and how strongly it changes rankings under fixed
candidate sets. It does not use a Gold standard and does not estimate NDCG,
Recall, MRR, or biological accuracy. It therefore cannot establish that the
production cap of 0.30 is optimal.

The experiment is separate from the frozen V5 target-gene holdout and was
designed after V5 was completed. None of the V5 labels, ranks, or internal audit
fields are inputs to this run. The runner records only the V5 design-file hash
as an integrity guard; V5 contents do not affect query selection or scoring.

## 2. Frozen design

The machine-readable design is
`config/exclusion_penalty_structural_v1_frozen_design.json`. The query table is
`benchmarks/penalty/penalty_structural_v1_queries.csv`.

All four profiles retain the same production target-scoring terms:

- RNA weight: 0.55;
- Protein weight: 0.30;
- Confidence weight: 0.15;
- Confidence composition: completeness 0.40, source support 0.35, and
  RNA-Protein consistency 0.25.

Only the maximum exclusion penalty changes:

| Profile | Maximum penalty | Role |
|---|---:|---|
| `P00_no_penalty` | 0.00 | Negative control |
| `P15_weak` | 0.15 | Weaker sensitivity |
| `P30_team_baseline` | 0.30 | Current production setting |
| `P45_stress` | 0.45 | Stress sensitivity only |

Five disease-context scenarios were purposively chosen from genes already
present in the local cache and frozen before the full P00/P15/P30/P45 output.
Selection used cache availability, candidate count, exclusion RNA coverage, and
distinct disease contexts. It did not use Gold labels or profile-comparison
outcomes. This is a convenience scenario set, not a random, deterministic, or
representative sample of genes; the exact pair choice contains human scenario
judgement.

| Query | Target | Exclusion | Disease | Candidates | Exclusion RNA observed |
|---|---|---|---|---:|---:|
| Q01 | FGFR2 | MET | Gastric | 43 | 43 |
| Q02 | ERBB2 | ESR1 | Breast | 73 | 73 |
| Q03 | CD3E | CD86 | Blood | 239 | 239 |
| Q04 | ASGR1 | AFP | Liver | 25 | 25 |
| Q05 | MSLN | MUC1 | Lung | 236 | 236 |

The pairs are mechanics-only scenarios rather than externally verified
experimental requirements. They do not assert that MET, ESR1, CD86, AFP, or
MUC1 is undesirable in every corresponding target-gene experiment. Reusing
familiar genes is acceptable here because the analysis measures score
mechanics, not held-out accuracy.

## 3. Production semantics retained by the runner

For the exclusion gene, the current production path calculates the available
mean of the DepMap, HPA, and GEO standardised RNA values. It then applies
Min-Max scaling within the disease-filtered target candidate set. The final
score is:

\[
\operatorname{Final}_\alpha =
\operatorname{clip}\left(
0.85\,\operatorname{Biology}
+0.15\,\operatorname{Confidence}
-\alpha\,\operatorname{ExclusionScaled},
0,1
\right).
\]

The scaled exclusion value is query-relative. A value of 0.8 in one disease
query is not an absolute expression threshold and is not directly comparable
with 0.8 in another query.

Missing exclusion RNA receives zero penalty in the production implementation.
The runner preserves that behaviour but records the observation as missing; it
does not reinterpret missing evidence as confirmed low expression.

Each query is fetched once. P00, P15, P30, and P45 then score copied versions of
the same candidate rows. The runner changes the module-level cap only in memory,
sequentially, and restores it after every call. The production source file and
all cache files remain unchanged.

## 4. Required implementation checks

The run completed 155 substantive check records, all with `PASS` status. They
include:

- identical candidate IDs across the four profiles;
- identical RNA, Protein, Biology, and Confidence components across profiles;
- zero Penalty for every P00 candidate;
- Penalty equal to the rounded profile cap multiplied by scaled exclusion;
- P30 equal to twice P15 and P45 equal to three times P15 within the required
  four-decimal rounding tolerance;
- candidate-level Penalty non-decreasing as the cap increases;
- candidate-level final score non-increasing as the cap increases;
- score loss equal to the reported Penalty unless clipping at zero limits the
  realised loss;
- P30 exact parity with an unmodified production scorer call for all five
  queries;
- repeated-run determinism;
- unchanged hashes for the production script, team data modules, sample info,
  frozen designs, query table, and ten gene caches.

All five real queries have complete observed exclusion-RNA coverage. They do not
empirically test how missing exclusion evidence affects real rankings. Synthetic
checks cover all-missing and partially missing exclusion values, lower-bound
clipping, constant positive/zero/negative exclusion values, and the
`GENE_MULTIOMICS`, `PROTEIN_ONLY`, and `COMBINED` subtraction paths. The
three query-mode checks verify both the reported Penalty field and the actual
subtraction from the final score. A separate synthetic cutoff case verifies
that exact final-score ties spanning ordinal positions 10 and 11 are both
retained. Checks that would be not applicable on a real query, such as a
missing-value assertion when no real candidate is missing exclusion RNA, are
not counted as automatic passes. The
constant-positive case documents an existing Min-Max edge behaviour:
when every observed exclusion value is the same positive value, the production
scaler assigns the full cap to every observed row. This experiment exposes but
does not change that behaviour.

## 5. Overall structural results

Values below are macro means across the five frozen queries. Top-10 sets include
all exact four-decimal final-score ties at the cutoff.

| Profile | Mean Top-10 Jaccard vs P00 | Queries with Top-10 change | Queries with Top-1 change | Mean Spearman vs P00 | Mean absolute rank change | Zero-score rows | Mean Top-10 scaled exclusion |
|---|---:|---:|---:|---:|---:|---:|---:|
| P00 | 1.0000 | 0/5 | 0/5 | 1.0000 | 0.00 | 0 | 0.4021 |
| P15 | 0.7463 | 4/5 | 0/5 | 0.9475 | 9.46 | 7 | 0.3160 |
| P30 | 0.7160 | 4/5 | 3/5 | 0.8778 | 15.20 | 76 | 0.3061 |
| P45 | 0.5880 | 5/5 | 3/5 | 0.8020 | 19.17 | 140 | 0.2457 |

Increasing the cap progressively reduces exclusion expression among the Top-10,
so the implementation acts in the intended direction. At the same time,
ranking disruption and clipping increase. P15 changes four Top-10 sets but
preserves all five Top-1 selections. P30 changes four Top-10 sets and three
Top-1 selections. P45 changes every Top-10 set and produces the largest rank
disruption.

The high-exclusion share in Top-10, defined structurally as scaled exclusion at
least 0.8 among candidates with observed exclusion evidence, falls from a
five-query mean of 0.14 under P00 to 0.06 under P30 and 0.02 under P45. All
five current queries have complete observed coverage, so the denominator is the
full Top-10 here. The 0.8 threshold is a query-relative structural diagnostic,
not a biological high-expression cutoff. This confirms stronger suppression,
not better biological recommendations.

## 6. Query-level P30 effects

| Query | Top-10 Jaccard vs P00 | P00 Top-1 | P30 Top-1 | Mean absolute rank change | Maximum rank change | P30 zero-score rows |
|---|---:|---|---|---:|---:|---:|
| FGFR2-MET, gastric | 0.6667 | SNU16 | SNU16 | 3.40 | 16.0 | 10 |
| ERBB2-ESR1, breast | 0.8182 | EFM192A | SUM190PT | 6.38 | 23.0 | 1 |
| CD3E-CD86, blood | 0.6667 | MOLT3 | MOLT3 | 44.08 | 128.0 | 39 |
| ASGR1-AFP, liver | 1.0000 | JHH5 | HUH1 | 1.48 | 4.0 | 0 |
| MSLN-MUC1, lung | 0.4286 | NCIH322M | NCIH2052 | 20.66 | 114.5 | 26 |

The strength of the response is heterogeneous. The CD3E-CD86 and MSLN-MUC1
queries show the largest rank movements and most zero-score clipping. The
ASGR1-AFP Top-10 membership remains unchanged at P30 even though its internal
order and Top-1 change. A single global cap can therefore be mild in one query
and disruptive in another because exclusion expression is scaled within each
candidate set.

Zero-score rows also create large tie groups. The largest P30 final-score tie is
39 candidates, and the largest P45 final-score tie is 71 candidates. The
production sorter uses `(finalScore, confidenceScore)` as its full comparator;
under that comparator the corresponding largest groups are 17 and 33. The
distinction matters because equal final scores may still be ordered by
confidence. The Top-10 cutoff itself has no ties in these five real runs, but
the runner keeps a tie-aware policy so future query sets do not silently
truncate equal final-score candidates.

## 7. Interpretation

The structural experiment supports four conclusions:

1. The current Penalty subtraction is executed correctly and reproducibly for
   the tested real queries and synthetic boundary paths.
2. Stronger caps lower every candidate's score monotonically according to its
   scaled exclusion value, and the aggregate Top-10 representation of
   high-exclusion candidates falls. An individual candidate's rank is not
   guaranteed to move downward because all candidates move relative to one
   another and clipping can create ties.
3. P30 is an active intervention, not a negligible term. Its impact varies
   markedly by query and can create substantial zero-score saturation.
4. P45 is useful as a stress condition but is not supported as a production
   candidate.

It does not show whether any changed recommendation is more suitable for a real
experiment, and the five convenience scenarios are not representative of all
genes or diseases. It also cannot determine whether 0.15, 0.30, or another cap
provides the best trade-off between retaining target-positive cell lines and
avoiding the exclusion marker. That question requires a separate pair-specific
evidence set that verifies both target suitability and exclusion expression.

## 8. Reproduction

Run the synthetic checks only:

```bash
cd "/Users/liyannan/Desktop/Data-Science-Group-Project"

/Users/liyannan/miniconda3/envs/cellline-v3/bin/python \
  scripts/exclusion_penalty_ablation_runner.py \
  --self-test
```

Run the frozen structural experiment:

```bash
cd "/Users/liyannan/Desktop/Data-Science-Group-Project"

/Users/liyannan/miniconda3/envs/cellline-v3/bin/python \
  scripts/exclusion_penalty_ablation_runner.py
```

The runner refuses to overwrite an existing result set. Use `--force` only when
intentionally reproducing the same frozen design. The five CSV outputs and JSON
manifest are written to `results/exclusion_penalty_structural_v1/`. The same
directory also contains `exclusion_penalty_structural_v1.xlsx`, a formatted
review copy of the CSV results. The CSV files and JSON manifest remain the
machine-readable source of truth.

## 9. Optional next stage

If time remains, create a new pair-specific candidate pool from the already
frozen P00/P15/P30/P45 results. Reviewers would need to verify two dimensions:
target suitability and whether the exclusion marker is high or low. Existing V5
target-only labels are insufficient for this purpose. If the pair-level labels
are then used to change the cap, that new benchmark must be described as
development evidence rather than an untouched holdout.
