#!/usr/bin/env python3
"""
ablation_runner.py
==================
CellLineSelector — controlled ablation over scoring configurations.

Answers the supervisor's requirement: "多组权重对照实验确定最终权重".

Design
------
Reuses the team's merged_cellline_selector.py for data access
(MergedDataRepository + fetch_candidate_evidence), then RE-SCORES the same
candidate rows under multiple configurations. This guarantees every config
sees byte-identical evidence — the only variable is the scoring rule.

Configurations
--------------
B0   Naive baseline: rank purely by mean standardized RNA. The "why bother"
     control — every ablation table needs it.
A0   Team formula as-is: RNA .55 / Protein .30 / Confidence .15
     (confidence = .40 completeness + .35 source support + .25 consistency).
A1   No protein: biological score = RNA only (redistributed).
A2   No confidence block (weight -> 0, renormalised biological terms).
A3   Equal biological weights RNA .425 / Protein .425 / Confidence .15.
A4   A0 + gene-adaptive RNA trust (OnCorr thresholds 0.50 / 0.25) — YL's
     contribution, ported from scoring_adaptive v3.
A5   v3 4-term structure: rna .444 / protein .278 / consistency .167 /
     completeness .111 (no composite confidence block).

Metrics (against benchmarks/gold_standard.csv, verified rows only)
------------------------------------------------------------------
Precision@K, Recall@K, NDCG@K, MRR for positives; mean rank percentile of
negatives (higher = better, negatives should sink); Jaccard overlap of
top-K with A0 (ranking stability across configs).

Usage
-----
    # 1. Put this next to merged_cellline_selector.py (renamed, importable)
    # 2. Build/obtain master_table.csv + cellline_annotations.csv
    python ablation_runner.py \
        --data-dir ./data_merged/merged \
        --gold benchmarks/gold_standard.csv \
        --corr gene_rna_protein_correlations.csv \
        --out results/ablation_results.csv

    # Smoke test with synthetic rows (no data needed):
    python ablation_runner.py --self-test
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------
# Gene-adaptive RNA trust (ported from scoring_adaptive.py v3)
# ----------------------------------------------------------------------
R_HIGH, R_LOW = 0.50, 0.25          # OnCorr (Nawaz et al. 2026) thresholds
TRUST_FLOOR, TRUST_DEFAULT = 0.30, 0.60


def rna_reliability(r: float) -> float:
    if r is None or (isinstance(r, float) and math.isnan(r)):
        return TRUST_DEFAULT
    if r >= R_HIGH:
        return 1.0
    if r < R_LOW:
        return TRUST_FLOOR
    return TRUST_FLOOR + (r - R_LOW) / (R_HIGH - R_LOW) * (1.0 - TRUST_FLOOR)


# ----------------------------------------------------------------------
# Configurable re-implementation of the team's score_candidates()
# ----------------------------------------------------------------------
DEFAULT_CFG = {
    "rna_w": 0.55,
    "protein_w": 0.30,
    "confidence_w": 0.15,
    "conf_completeness": 0.40,
    "conf_source_support": 0.35,
    "conf_consistency": 0.25,
    "exclusion_max_penalty": 0.30,
    "adaptive_trust": False,       # A4 switches this on
    "v3_structure": False,         # A5: 4-term flat structure
    "rna_only_rank": False,        # B0
}

# Two families:
#   Leave-one-out (A1..A4, A3b): each differs from A0 by exactly ONE design
#     decision -> attributes contribution to individual components.
#   Architecture comparison (B0, A5, A6): whole-structure model selection.
#     A0/A4/A5/A6 form a 2x2 factorial: {team, v3 structure} x {adaptive off, on}.
CONFIGS: dict[str, dict] = {
    "B0_naive_rna":       {**DEFAULT_CFG, "rna_only_rank": True},
    "A0_team_baseline":   {**DEFAULT_CFG},
    "A1_no_protein":      {**DEFAULT_CFG, "protein_w": 0.0, "rna_w": 0.85},
    "A2_no_confidence":   {**DEFAULT_CFG, "confidence_w": 0.0,
                           "rna_w": 0.647, "protein_w": 0.353},
    "A3_equal_bio":       {**DEFAULT_CFG, "rna_w": 0.425, "protein_w": 0.425},
    "A3b_reversed_bio":   {**DEFAULT_CFG, "rna_w": 0.30, "protein_w": 0.55},
    "A4_adaptive_trust":  {**DEFAULT_CFG, "adaptive_trust": True},
    "A5_v3_structure":    {**DEFAULT_CFG, "v3_structure": True},
    "A6_v3_full":         {**DEFAULT_CFG, "v3_structure": True,
                           "adaptive_trust": True},
}


def _minmax_masked(vals: list[float | None]) -> list[float | None]:
    xs = [v for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))]
    if not xs:
        return [None] * len(vals)
    lo, hi = min(xs), max(xs)
    if hi == lo:
        return [0.5 if v is not None else None for v in vals]
    return [((v - lo) / (hi - lo)) if v is not None and not (isinstance(v, float) and math.isnan(v))
            else None for v in vals]


def score_rows(rows: list[dict], cfg: dict, gene: str,
               corr_table: dict | None = None) -> list[dict]:
    """Re-scores fetch_candidate_evidence() rows under one configuration.

    Mirrors the team's score_candidates() semantics (row-wise redistribution
    when protein is missing; expression-based exclusion penalty; clip to
    [0,1]) with weights/structure taken from cfg.
    """
    if not rows:
        return []

    rna_s = _minmax_masked([r.get("rnaExpr") for r in rows])
    prot_s = _minmax_masked([r.get("protExpr") for r in rows])
    excl_s = _minmax_masked([r.get("exclusionExpr") for r in rows])

    trust = 1.0
    if cfg["adaptive_trust"] and corr_table is not None:
        trust = rna_reliability(corr_table.get(gene, float("nan")))

    out = []
    for i, row in enumerate(rows):
        has_rna, has_prot = rna_s[i] is not None, prot_s[i] is not None
        rna = float(rna_s[i]) if has_rna else 0.0
        prot = float(prot_s[i]) if has_prot else 0.0
        penalty = cfg["exclusion_max_penalty"] * float(excl_s[i]) if excl_s[i] is not None else 0.0

        completeness = (float(has_rna) + float(has_prot)) / 2.0
        support = (float(bool(row.get("hasDepMapRNA"))) + float(bool(row.get("hasHpaRNA")))
                   + float(bool(row.get("hasGeoRNA"))) + float(bool(row.get("hasProteomics")))) / 4.0
        if has_rna and has_prot:
            consistency = max(0.0, min(1.0, 1.0 - abs(rna - prot)))
        elif has_rna or has_prot:
            consistency = 0.5
        else:
            consistency = 0.0

        if cfg["rna_only_rank"]:                       # B0
            final = rna if has_rna else 0.0

        elif cfg["v3_structure"]:                      # A5: 4-term flat
            w = {"rna": 0.444, "protein": 0.278,
                 "consistency": 0.167, "completeness": 0.111}
            rna_w = w["rna"] * (trust if cfg["adaptive_trust"] else 1.0)
            freed = w["rna"] - rna_w
            prot_w = w["protein"] + (freed if has_prot else 0.0)
            if not has_prot:
                rna_w += freed + w["protein"] + w["consistency"]
                prot_w, cons_w = 0.0, 0.0
            else:
                cons_w = w["consistency"]
            final = (rna_w * rna + prot_w * prot
                     + cons_w * consistency + w["completeness"] * completeness) - penalty

        else:                                          # team structure
            rna_w = cfg["rna_w"] * trust
            freed = cfg["rna_w"] - rna_w
            prot_w = cfg["protein_w"] + (freed if has_prot else 0.0)
            if not has_prot:
                rna_w += freed
            avail, weighted = 0.0, 0.0
            if has_rna:
                avail += rna_w
                weighted += rna_w * rna
            if has_prot:
                avail += prot_w
                weighted += prot_w * prot
            if avail == 0.0:
                continue
            bio = weighted / avail
            conf = (cfg["conf_completeness"] * completeness
                    + cfg["conf_source_support"] * support
                    + cfg["conf_consistency"] * consistency)
            final = (cfg["rna_w"] + cfg["protein_w"]) * bio + cfg["confidence_w"] * conf - penalty

        out.append({"cellLine": row.get("cellLine"),
                    "finalScore": max(0.0, min(1.0, final))})

    out.sort(key=lambda d: d["finalScore"], reverse=True)
    for rank, d in enumerate(out, 1):
        d["rank"] = rank
    return out


# ----------------------------------------------------------------------
# Metrics
# ----------------------------------------------------------------------
def _norm_name(s: str) -> str:
    import re
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def evaluate_ranking(ranked: list[dict], positives: set[str],
                     negatives: set[str], k: int = 10) -> dict:
    names = [_norm_name(d["cellLine"]) for d in ranked]
    pos = {_norm_name(p) for p in positives}
    neg = {_norm_name(n) for n in negatives}
    n = len(names)
    if n == 0 or not pos:
        return {}

    topk = names[:k]
    hits = [nm for nm in topk if nm in pos]
    precision = len(hits) / k
    recall = len(hits) / len(pos & set(names)) if pos & set(names) else float("nan")

    dcg = sum(1.0 / math.log2(i + 2) for i, nm in enumerate(topk) if nm in pos)
    ideal = sum(1.0 / math.log2(i + 2) for i in range(min(k, len(pos & set(names)))))
    ndcg = dcg / ideal if ideal > 0 else float("nan")

    rr = 0.0
    for i, nm in enumerate(names):
        if nm in pos:
            rr = 1.0 / (i + 1)
            break

    neg_pcts = [1.0 - (names.index(nm) / (n - 1)) if n > 1 else 0.5
                for nm in neg if nm in names]
    neg_depth = float(np.mean(neg_pcts)) if neg_pcts else float("nan")

    return {"precision_at_k": precision, "recall_at_k": recall,
            "ndcg_at_k": ndcg, "mrr": rr,
            "neg_sink": 1.0 - neg_depth if not math.isnan(neg_depth) else float("nan")}


def jaccard_topk(a: list[dict], b: list[dict], k: int = 10) -> float:
    sa = {_norm_name(d["cellLine"]) for d in a[:k]}
    sb = {_norm_name(d["cellLine"]) for d in b[:k]}
    return len(sa & sb) / len(sa | sb) if sa | sb else float("nan")


# ----------------------------------------------------------------------
# Runner
# ----------------------------------------------------------------------
def run_ablation(repo, gold: pd.DataFrame, corr_table: dict,
                 k: int = 10, disease: str | None = None) -> pd.DataFrame:
    """repo: merged_cellline_selector.MergedDataRecommender instance."""
    import merged_cellline_selector as mcs

    # 组员的 _disease_filter 是强制的; [""] 使 `alias in text` 恒真 = 不过滤
    aliases = mcs.build_disease_aliases(disease) if disease else [""]

    genes = sorted(gold["gene"].unique())
    per_config_metrics: dict[str, list[dict]] = {c: [] for c in CONFIGS}
    per_config_rank_a0: dict[str, list[float]] = {c: [] for c in CONFIGS}

    for gene in genes:
        g = gold[gold["gene"] == gene]
        positives = set(g.loc[g["relation"] == "positive", "expected_cell_line"])
        negatives = set(g.loc[g["relation"] == "negative", "expected_cell_line"])
        if not positives:
            continue

        try:
            rows = repo.fetch_candidate_evidence(gene, aliases, None)
        except Exception as exc:
            print(f"  [skip] {gene}: {exc}")
            continue
        if not rows:
            print(f"  [skip] {gene}: no candidate rows (gene not in master table?)")
            continue

        ranked_a0 = score_rows(rows, CONFIGS["A0_team_baseline"], gene, corr_table)
        for cname, cfg in CONFIGS.items():
            ranked = score_rows(rows, cfg, gene, corr_table)
            m = evaluate_ranking(ranked, positives, negatives, k=k)
            if m:
                m["gene"] = gene
                per_config_metrics[cname].append(m)
                per_config_rank_a0[cname].append(jaccard_topk(ranked, ranked_a0, k))

    records = []
    for cname in CONFIGS:
        ms = per_config_metrics[cname]
        if not ms:
            continue
        df = pd.DataFrame(ms)
        records.append({
            "config": cname,
            "n_genes": len(ms),
            f"P@{k}": df["precision_at_k"].mean(),
            f"NDCG@{k}": df["ndcg_at_k"].mean(),
            "MRR": df["mrr"].mean(),
            "neg_sink": df["neg_sink"].mean(),
            f"Jaccard_vs_A0@{k}": float(np.mean(per_config_rank_a0[cname])),
        })
    return pd.DataFrame(records)




# ----------------------------------------------------------------------
# Self-test with synthetic rows (no data needed)
# ----------------------------------------------------------------------
def self_test() -> int:
    rng = np.random.default_rng(0)
    rows = []
    for i in range(60):
        expr = float(rng.normal(5, 2))
        rows.append({
            "cellLine": f"LINE-{i:03d}", "rnaExpr": expr, "nRna": 3,
            "protExpr": expr * 0.8 + float(rng.normal(0, 1)) if i % 3 else None,
            "nProt": 1 if i % 3 else 0, "exclusionExpr": None, "nExclusion": 0,
            "hasDepMapRNA": True, "hasHpaRNA": i % 2 == 0,
            "hasGeoRNA": i % 4 == 0, "hasProteomics": i % 3 != 0,
        })
    rows[5]["rnaExpr"], rows[5]["protExpr"] = 12.0, 10.0   # planted positive
    rows[7]["rnaExpr"] = 0.1                                # planted negative

    gold_pos, gold_neg = {"LINE-005"}, {"LINE-007"}
    corr = {"TESTGENE": 0.86}

    print(f"{'config':<20} {'P@10':>6} {'NDCG@10':>8} {'MRR':>6} {'neg_sink':>9} {'Jac_vs_A0':>10}")
    ranked_a0 = score_rows(rows, CONFIGS["A0_team_baseline"], "TESTGENE", corr)
    ok = True
    for cname, cfg in CONFIGS.items():
        ranked = score_rows(rows, cfg, "TESTGENE", corr)
        m = evaluate_ranking(ranked, gold_pos, gold_neg, k=10)
        j = jaccard_topk(ranked, ranked_a0, 10)
        print(f"{cname:<20} {m['precision_at_k']:>6.2f} {m['ndcg_at_k']:>8.2f} "
              f"{m['mrr']:>6.2f} {m['neg_sink']:>9.2f} {j:>10.2f}")
        if m["mrr"] < 0.2:
            ok = False
    print("\n" + ("self-test OK: planted positive ranks near top in all configs"
                  if ok else "self-test FAILED"))
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, help="folder with master_table.csv")
    ap.add_argument("--gold", type=Path, default=Path("benchmarks/gold_standard.csv"))
    ap.add_argument("--corr", type=Path, default=Path("gene_rna_protein_correlations.csv"))
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--disease", default=None)
    ap.add_argument("--out", type=Path, default=Path("results/ablation_results.csv"))
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    if not args.data_dir:
        ap.error("--data-dir is required (or use --self-test)")

    import merged_cellline_selector as mcs
    repo = mcs.MergedDataRecommender(args.data_dir)

    gold = pd.read_csv(args.gold)
    if "verified" in gold.columns:
        unverified = (gold["verified"].str.lower() != "yes").sum()
        gold = gold[gold["verified"].str.lower() == "yes"]
        print(f"gold standard: using {len(gold)} verified rows ({unverified} unverified excluded)")
    if gold.empty:
        print("No verified gold rows. Verify entries (set verified=yes) first.")
        return 1

    corr_df = pd.read_csv(args.corr)
    corr_table = dict(zip(corr_df["gene"], corr_df["correlation"]))

    results = run_ablation(repo, gold, corr_table, k=args.k, disease=args.disease)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.out, index=False)
    print("\n" + results.to_string(index=False))
    print(f"\nsaved -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
