#!/usr/bin/env python3
"""CellLineSelector controlled ablation runner.

Evidence retrieval is shared across configurations. Evaluation uses DepMap IDs,
penalises missing gold positives, handles exact-score ties, and writes both a
summary and a per-gene audit table.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import merged_cellline_selector as mcs


R_HIGH, R_LOW = 0.50, 0.25
TRUST_FLOOR, TRUST_DEFAULT = 0.30, 0.60


def rna_reliability(r: float | None) -> float:
    if r is None or pd.isna(r):
        return TRUST_DEFAULT
    if r >= R_HIGH:
        return 1.0
    if r < R_LOW:
        return TRUST_FLOOR
    return TRUST_FLOOR + (r - R_LOW) / (R_HIGH - R_LOW) * (1.0 - TRUST_FLOOR)


_PRODUCTION_BASELINE = mcs.DEFAULT_SCORING_CONFIG

DEFAULT_CFG = {
    "rna_w": _PRODUCTION_BASELINE.rna_weight,
    "protein_w": _PRODUCTION_BASELINE.protein_weight,
    "confidence_w": _PRODUCTION_BASELINE.confidence_weight,
    "conf_completeness": _PRODUCTION_BASELINE.confidence_completeness,
    "conf_source_support": _PRODUCTION_BASELINE.confidence_source_support,
    "conf_consistency": _PRODUCTION_BASELINE.confidence_consistency,
    "exclusion_max_penalty": _PRODUCTION_BASELINE.max_exclusion_penalty,
    "adaptive_trust": False,
    "v3_structure": False,
    "rna_only_rank": False,
    # full: protein affects biology and confidence
    # direct_off: direct term off, protein-derived confidence retained
    # none: protein and all protein-derived confidence features removed
    "protein_mode": "full",
}


CONFIGS: dict[str, dict] = {
    "B0_rna_mean": {
        **DEFAULT_CFG,
        "rna_only_rank": True,
        "protein_mode": "none",
    },
    "A0_team_baseline": {**DEFAULT_CFG},
    "A1a_no_direct_protein": {
        **DEFAULT_CFG,
        "protein_w": 0.0,
        "rna_w": 0.85,
        "protein_mode": "direct_off",
    },
    "A1b_no_protein_evidence": {
        **DEFAULT_CFG,
        "protein_w": 0.0,
        "rna_w": 0.85,
        "protein_mode": "none",
    },
    "A2_no_confidence": {
        **DEFAULT_CFG,
        "confidence_w": 0.0,
        "rna_w": 0.647,
        "protein_w": 0.353,
    },
    "A3_equal_bio": {
        **DEFAULT_CFG,
        "rna_w": 0.425,
        "protein_w": 0.425,
    },
    "A3b_reversed_bio": {
        **DEFAULT_CFG,
        "rna_w": 0.30,
        "protein_w": 0.55,
    },
    "A4_adaptive_trust": {
        **DEFAULT_CFG,
        "adaptive_trust": True,
    },
    "A5_v3_structure": {
        **DEFAULT_CFG,
        "v3_structure": True,
    },
    "A6_v3_full": {
        **DEFAULT_CFG,
        "v3_structure": True,
        "adaptive_trust": True,
    },
}


CONFIG_DESCRIPTIONS = {
    "B0_rna_mean": "Mean standardized DepMap/HPA/GEO RNA only",
    "A0_team_baseline": "Team formula: RNA .55, protein .30, confidence .15",
    "A1a_no_direct_protein": "Remove direct protein term; retain protein-derived confidence",
    "A1b_no_protein_evidence": "Remove protein and all protein-derived confidence features",
    "A2_no_confidence": "Remove composite confidence block",
    "A3_equal_bio": "Equal RNA/protein biological weights",
    "A3b_reversed_bio": "Protein biological weight greater than RNA",
    "A4_adaptive_trust": "Team formula plus gene-adaptive RNA trust",
    "A5_v3_structure": "Flat v3 four-term structure without adaptive trust",
    "A6_v3_full": "Flat v3 four-term structure plus adaptive trust",
}


def _validate_configs() -> None:
    for name, cfg in CONFIGS.items():
        if cfg["protein_mode"] not in {"full", "direct_off", "none"}:
            raise ValueError(f"{name}: invalid protein_mode")
        if cfg["rna_only_rank"] or cfg["v3_structure"]:
            continue
        total = cfg["rna_w"] + cfg["protein_w"] + cfg["confidence_w"]
        if not math.isclose(total, 1.0, abs_tol=1e-6):
            raise ValueError(f"{name}: outer weights sum to {total}, not 1")


def _is_present(value: object) -> bool:
    return value is not None and not pd.isna(value)


def _flag(value: object) -> bool:
    return bool(value) if _is_present(value) else False


def _minmax_masked(vals: list[float | None]) -> list[float | None]:
    xs = [float(v) for v in vals if _is_present(v)]
    if not xs:
        return [None] * len(vals)
    lo, hi = min(xs), max(xs)
    if math.isclose(hi, lo):
        constant = 1.0 if hi > 0 else 0.0
        return [constant if _is_present(v) else None for v in vals]
    return [((float(v) - lo) / (hi - lo)) if _is_present(v) else None for v in vals]


def _norm_name(value: object) -> str:
    import re

    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def _candidate_id(row: dict) -> str:
    depmap_id = row.get("DepMap_ID")
    if _is_present(depmap_id) and str(depmap_id).strip():
        return str(depmap_id).strip().upper()
    return "NAME:" + _norm_name(row.get("cellLine", ""))


def _confidence_components(
    row: dict,
    has_rna: bool,
    has_protein: bool,
    rna_score: float,
    protein_score: float,
    protein_mode: str,
) -> tuple[float, float, float | None]:
    if protein_mode == "none":
        completeness = float(has_rna)
        support = (
            float(_flag(row.get("hasDepMapRNA")))
            + float(_flag(row.get("hasHpaRNA")))
            + float(_flag(row.get("hasGeoRNA")))
        ) / 3.0
        return completeness, support, None

    completeness = (float(has_rna) + float(has_protein)) / 2.0
    support = (
        float(_flag(row.get("hasDepMapRNA")))
        + float(_flag(row.get("hasHpaRNA")))
        + float(_flag(row.get("hasGeoRNA")))
        + float(_flag(row.get("hasProteomics")))
    ) / 4.0
    if has_rna and has_protein:
        consistency = max(0.0, min(1.0, 1.0 - abs(rna_score - protein_score)))
    elif has_rna or has_protein:
        consistency = 0.5
    else:
        consistency = 0.0
    return completeness, support, consistency


def _composite_confidence(
    cfg: dict,
    completeness: float,
    support: float,
    consistency: float | None,
) -> float:
    parts = [
        (cfg["conf_completeness"], completeness),
        (cfg["conf_source_support"], support),
    ]
    if consistency is not None:
        parts.append((cfg["conf_consistency"], consistency))
    active_weight = sum(weight for weight, _ in parts)
    return sum(weight * value for weight, value in parts) / active_weight


def _assign_tie_ranks(ranked: list[dict]) -> None:
    start = 0
    while start < len(ranked):
        end = start + 1
        score = ranked[start]["finalScore"]
        while end < len(ranked) and ranked[end]["finalScore"] == score:
            end += 1
        midrank = ((start + 1) + end) / 2.0
        for i in range(start, end):
            ranked[i]["rank"] = i + 1
            ranked[i]["rankMid"] = midrank
        start = end


def score_rows(
    rows: list[dict],
    cfg: dict,
    gene: str,
    corr_table: dict[str, float] | None = None,
) -> list[dict]:
    """Re-score one gene's candidate evidence under one configuration."""
    if not rows:
        return []

    rna_s = _minmax_masked([row.get("rnaExpr") for row in rows])
    prot_s = _minmax_masked([row.get("protExpr") for row in rows])
    excl_s = _minmax_masked([row.get("exclusionExpr") for row in rows])
    trust = 1.0
    if cfg["adaptive_trust"]:
        correlation = corr_table.get(gene, np.nan) if corr_table is not None else np.nan
        trust = rna_reliability(correlation)

    out: list[dict] = []
    for i, row in enumerate(rows):
        has_rna = rna_s[i] is not None
        has_protein = prot_s[i] is not None and cfg["protein_mode"] != "none"
        rna = float(rna_s[i]) if has_rna else 0.0
        protein = float(prot_s[i]) if has_protein else 0.0
        penalty = (
            cfg["exclusion_max_penalty"] * float(excl_s[i])
            if excl_s[i] is not None
            else 0.0
        )
        completeness, support, consistency = _confidence_components(
            row, has_rna, has_protein, rna, protein, cfg["protein_mode"]
        )

        if cfg["rna_only_rank"]:
            final = rna if has_rna else 0.0

        elif cfg["v3_structure"]:
            weights = {
                "rna": 0.444,
                "protein": 0.278,
                "consistency": 0.167,
                "completeness": 0.111,
            }
            rna_w = weights["rna"] * (trust if cfg["adaptive_trust"] else 1.0)
            freed = weights["rna"] - rna_w
            protein_w = weights["protein"] + (freed if has_protein else 0.0)
            if not has_protein:
                rna_w += freed + weights["protein"] + weights["consistency"]
                protein_w, consistency_w = 0.0, 0.0
            else:
                consistency_w = weights["consistency"]
            final = (
                rna_w * rna
                + protein_w * protein
                + consistency_w * float(consistency)
                + weights["completeness"] * completeness
                - penalty
            )

        else:
            rna_w = cfg["rna_w"] * trust
            freed = cfg["rna_w"] - rna_w
            protein_w = cfg["protein_w"] + (freed if has_protein else 0.0)
            if not has_protein:
                rna_w += freed
            available_weight, weighted_biology = 0.0, 0.0
            if has_rna:
                available_weight += rna_w
                weighted_biology += rna_w * rna
            if has_protein and protein_w > 0:
                available_weight += protein_w
                weighted_biology += protein_w * protein
            if available_weight == 0:
                continue
            biological_score = weighted_biology / available_weight
            confidence_score = _composite_confidence(
                cfg, completeness, support, consistency
            )
            final = (
                (cfg["rna_w"] + cfg["protein_w"]) * biological_score
                + cfg["confidence_w"] * confidence_score
                - penalty
            )

        out.append(
            {
                "DepMap_ID": _candidate_id(row),
                "cellLine": row.get("cellLine"),
                "finalScore": max(0.0, min(1.0, float(final))),
            }
        )

    out.sort(key=lambda item: (-item["finalScore"], item["DepMap_ID"]))
    candidate_ids = [item["DepMap_ID"] for item in out]
    if len(candidate_ids) != len(set(candidate_ids)):
        raise ValueError(f"{gene}: duplicate candidate DepMap IDs")
    _assign_tie_ranks(out)
    return out


def _norm_ids(values: Iterable[object]) -> set[str]:
    return {str(value).strip().upper() for value in values if _is_present(value)}


def _score_groups(ranked: list[dict]) -> list[list[dict]]:
    groups: list[list[dict]] = []
    for item in ranked:
        if not groups or groups[-1][0]["finalScore"] != item["finalScore"]:
            groups.append([item])
        else:
            groups[-1].append(item)
    return groups


def _tie_aware_at_k(
    groups: list[list[dict]], positives: set[str], k: int
) -> tuple[float, float, float]:
    expected_hits, dcg, position = 0.0, 0.0, 1
    for group in groups:
        if position > k:
            break
        group_size = len(group)
        group_positive = sum(item["DepMap_ID"] in positives for item in group)
        slots = min(group_size, k - position + 1)
        positive_fraction = group_positive / group_size
        expected_hits += slots * positive_fraction
        dcg += positive_fraction * sum(
            1.0 / math.log2(rank + 1)
            for rank in range(position, position + slots)
        )
        position += group_size
    precision = expected_hits / k
    recall = expected_hits / len(positives) if positives else np.nan
    ideal = sum(
        1.0 / math.log2(rank + 1)
        for rank in range(1, min(k, len(positives)) + 1)
    )
    return precision, recall, dcg / ideal if ideal else np.nan


def evaluate_ranking(
    ranked: list[dict],
    positives: set[str],
    negatives: set[str],
    ks: Iterable[int] = (3, 5, 10),
) -> dict:
    positives, negatives = _norm_ids(positives), _norm_ids(negatives)
    if not ranked or not positives:
        return {}
    candidate_ids = {item["DepMap_ID"] for item in ranked}
    groups = _score_groups(ranked)
    result: dict[str, float | int] = {
        "n_candidates": len(ranked),
        "n_pos_gold": len(positives),
        "n_pos_found": len(positives & candidate_ids),
        "n_neg_gold": len(negatives),
        "n_neg_found": len(negatives & candidate_ids),
    }
    for k in sorted(set(ks)):
        precision, recall, ndcg = _tie_aware_at_k(groups, positives, k)
        result[f"P@{k}"] = precision
        result[f"Recall@{k}"] = recall
        result[f"NDCG@{k}"] = ndcg
    positive_midrank = [
        item["rankMid"] for item in ranked if item["DepMap_ID"] in positives
    ]
    result["MRR"] = 1.0 / min(positive_midrank) if positive_midrank else 0.0
    negative_midrank = [
        item["rankMid"] for item in ranked if item["DepMap_ID"] in negatives
    ]
    n = len(ranked)
    result["neg_sink"] = (
        float(np.mean([(rank - 1.0) / (n - 1.0) for rank in negative_midrank]))
        if negative_midrank and n > 1
        else np.nan
    )
    return result


def _topk_with_ties(ranked: list[dict], k: int) -> set[str]:
    if not ranked:
        return set()
    cutoff_score = ranked[min(k, len(ranked)) - 1]["finalScore"]
    return {
        item["DepMap_ID"] for item in ranked if item["finalScore"] >= cutoff_score
    }


def jaccard_topk(a: list[dict], b: list[dict], k: int = 10) -> float:
    set_a, set_b = _topk_with_ties(a, k), _topk_with_ties(b, k)
    return len(set_a & set_b) / len(set_a | set_b) if set_a | set_b else np.nan


def _bootstrap_mean_delta(
    values: np.ndarray, iterations: int, seed: int
) -> tuple[float, float, float]:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan, np.nan, np.nan
    delta = float(values.mean())
    if iterations <= 0:
        return delta, np.nan, np.nan
    rng = np.random.default_rng(seed)
    samples = rng.choice(values, size=(iterations, values.size), replace=True).mean(axis=1)
    low, high = np.quantile(samples, [0.025, 0.975])
    return delta, float(low), float(high)


def _summarise(
    detail: pd.DataFrame,
    ks: tuple[int, ...],
    bootstrap_iterations: int,
    seed: int,
) -> pd.DataFrame:
    metric_columns = [
        column
        for column in detail.columns
        if column.startswith(("P@", "Recall@", "NDCG@"))
        or column in {"MRR", "neg_sink", "Jaccard_vs_A0@10_ties"}
    ]
    records: list[dict] = []
    for config_name in CONFIGS:
        sub = detail[detail["config"] == config_name]
        if sub.empty:
            continue
        records.append(
            {
                "config": config_name,
                "description": CONFIG_DESCRIPTIONS[config_name],
                "n_genes": sub["gene"].nunique(),
                "n_gold_positive": int(sub["n_pos_gold"].sum()),
                "n_positive_found": int(sub["n_pos_found"].sum()),
                "positive_coverage": sub["n_pos_found"].sum()
                / sub["n_pos_gold"].sum(),
                "n_gold_negative": int(sub["n_neg_gold"].sum()),
                "n_negative_found": int(sub["n_neg_found"].sum()),
                **{column: sub[column].mean() for column in metric_columns},
            }
        )
    summary = pd.DataFrame(records)
    if summary.empty:
        return summary

    primary_k = 5 if 5 in ks else ks[0]
    for metric in (f"NDCG@{primary_k}", "MRR"):
        baseline = detail.loc[
            detail["config"] == "A0_team_baseline", ["gene", metric]
        ].rename(columns={metric: "baseline"})
        delta_rows = []
        for offset, config_name in enumerate(summary["config"]):
            current = detail.loc[
                detail["config"] == config_name, ["gene", metric]
            ].rename(columns={metric: "current"})
            paired = current.merge(baseline, on="gene", how="inner")
            values = (paired["current"] - paired["baseline"]).to_numpy(dtype=float)
            delta_rows.append(
                _bootstrap_mean_delta(values, bootstrap_iterations, seed + offset)
            )
        summary[f"delta_{metric}_vs_A0"] = [row[0] for row in delta_rows]
        summary[f"delta_{metric}_ci_low"] = [row[1] for row in delta_rows]
        summary[f"delta_{metric}_ci_high"] = [row[2] for row in delta_rows]
    return summary


def run_ablation(
    repo,
    gold: pd.DataFrame,
    corr_table: dict[str, float],
    ks: tuple[int, ...] = (3, 5, 10),
    disease: str | None = None,
    bootstrap_iterations: int = 2000,
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    import merged_cellline_selector as mcs

    aliases = mcs.build_disease_aliases(disease) if disease else [""]
    detail_records: list[dict] = []
    for gene in sorted(gold["gene"].unique()):
        gene_gold = gold[gold["gene"] == gene]
        positives = _norm_ids(
            gene_gold.loc[
                gene_gold["relation"] == "positive", "expected_depmap_id"
            ]
        )
        negatives = _norm_ids(
            gene_gold.loc[
                gene_gold["relation"] == "negative", "expected_depmap_id"
            ]
        )
        if not positives:
            print(f"  [skip] {gene}: no positive gold entries")
            continue
        try:
            rows = repo.fetch_candidate_evidence(gene, aliases, None)
        except Exception as exc:
            print(f"  [skip] {gene}: {exc}")
            continue
        if not rows:
            print(f"  [skip] {gene}: no candidate rows")
            continue

        ranked_a0 = score_rows(rows, CONFIGS["A0_team_baseline"], gene, corr_table)
        correlation = corr_table.get(gene, np.nan)
        trust = rna_reliability(correlation)
        for config_name, cfg in CONFIGS.items():
            ranked = score_rows(rows, cfg, gene, corr_table)
            metrics = evaluate_ranking(ranked, positives, negatives, ks=ks)
            if metrics:
                detail_records.append(
                    {
                        "config": config_name,
                        "gene": gene,
                        "rna_protein_correlation": correlation,
                        "rna_trust": trust,
                        **metrics,
                        "Jaccard_vs_A0@10_ties": jaccard_topk(
                            ranked, ranked_a0, 10
                        ),
                    }
                )
    detail = pd.DataFrame(detail_records)
    if detail.empty:
        return pd.DataFrame(), detail
    return (
        _summarise(detail, ks, bootstrap_iterations, seed),
        detail,
    )


def _prepare_gold(gold: pd.DataFrame) -> pd.DataFrame:
    required = {"gene", "expected_depmap_id", "relation"}
    missing = sorted(required - set(gold.columns))
    if missing:
        raise KeyError("gold standard missing columns: " + ", ".join(missing))
    out = gold.copy()
    if "verified" in out.columns:
        verified = out["verified"].fillna("").astype(str).str.lower().eq("yes")
        print(
            f"gold standard: using {int(verified.sum())} verified rows "
            f"({int((~verified).sum())} excluded)"
        )
        out = out.loc[verified].copy()
    if out["expected_depmap_id"].isna().any():
        raise ValueError("verified gold rows contain missing expected_depmap_id values")
    out["gene"] = out["gene"].astype(str).str.upper().str.strip()
    out["relation"] = out["relation"].astype(str).str.lower().str.strip()
    out["expected_depmap_id"] = (
        out["expected_depmap_id"].astype(str).str.upper().str.strip()
    )
    invalid = sorted(set(out["relation"]) - {"positive", "negative"})
    if invalid:
        raise ValueError(f"invalid gold relations: {invalid}")
    duplicate_mask = out.duplicated(
        ["gene", "expected_depmap_id", "relation"], keep=False
    )
    if duplicate_mask.any():
        duplicates = out.loc[
            duplicate_mask, ["gene", "expected_depmap_id", "relation"]
        ]
        raise ValueError("duplicate gold entries:\n" + duplicates.to_string(index=False))
    conflicting = out.groupby(["gene", "expected_depmap_id"])["relation"].nunique()
    if (conflicting > 1).any():
        raise ValueError("the same gene/DepMap ID is labelled both positive and negative")
    if out.empty:
        raise ValueError("no verified gold rows")
    return out


def self_test() -> int:
    _validate_configs()
    assert rna_reliability(np.nan) == TRUST_DEFAULT
    assert rna_reliability(0.10) == TRUST_FLOOR
    assert rna_reliability(0.50) == 1.0
    rows = [
        {
            "DepMap_ID": "ACH-000001",
            "cellLine": "RNA-high",
            "rnaExpr": 10.0,
            "protExpr": 0.0,
            "exclusionExpr": None,
            "hasDepMapRNA": True,
            "hasHpaRNA": True,
            "hasGeoRNA": True,
            "hasProteomics": True,
            "nRna": 3,
            "nProt": 1,
            "nExclusion": 0,
        },
        {
            "DepMap_ID": "ACH-000002",
            "cellLine": "protein-high",
            "rnaExpr": 3.0,
            "protExpr": 10.0,
            "exclusionExpr": None,
            "hasDepMapRNA": True,
            "hasHpaRNA": True,
            "hasGeoRNA": True,
            "hasProteomics": True,
            "nRna": 3,
            "nProt": 1,
            "nExclusion": 0,
        },
        {
            "DepMap_ID": "ACH-000003",
            "cellLine": "negative",
            "rnaExpr": 0.0,
            "protExpr": 0.0,
            "exclusionExpr": None,
            "hasDepMapRNA": True,
            "hasHpaRNA": False,
            "hasGeoRNA": False,
            "hasProteomics": True,
            "nRna": 1,
            "nProt": 1,
            "nExclusion": 0,
        },
    ]
    high_corr, low_corr = {"TEST": 0.80}, {"TEST": 0.10}
    a0 = score_rows(rows, CONFIGS["A0_team_baseline"], "TEST", high_corr)
    production_a0 = mcs.score_candidates(rows, config=mcs.DEFAULT_SCORING_CONFIG)
    assert [x["DepMap_ID"] for x in a0] == [x["DepMap_ID"] for x in production_a0]
    assert all(
        math.isclose(prod["finalScore"], round(ablation["finalScore"], 4))
        for prod, ablation in zip(production_a0, a0)
    )
    a4_high = score_rows(rows, CONFIGS["A4_adaptive_trust"], "TEST", high_corr)
    a4_low = score_rows(rows, CONFIGS["A4_adaptive_trust"], "TEST", low_corr)
    assert [(x["DepMap_ID"], x["finalScore"]) for x in a0] == [
        (x["DepMap_ID"], x["finalScore"]) for x in a4_high
    ]
    assert [x["DepMap_ID"] for x in a0] != [x["DepMap_ID"] for x in a4_low]

    no_protein_cfg = CONFIGS["A1b_no_protein_evidence"]
    before = score_rows(rows, no_protein_cfg, "TEST", low_corr)
    changed = [dict(row) for row in rows]
    changed[0]["protExpr"], changed[0]["hasProteomics"] = 1000.0, False
    changed[1]["protExpr"] = -1000.0
    after = score_rows(changed, no_protein_cfg, "TEST", low_corr)
    assert [(x["DepMap_ID"], x["finalScore"]) for x in before] == [
        (x["DepMap_ID"], x["finalScore"]) for x in after
    ]

    tie_ranked = [
        {"DepMap_ID": "ACH-POS", "finalScore": 1.0, "rankMid": 1.5},
        {"DepMap_ID": "ACH-TIE", "finalScore": 1.0, "rankMid": 1.5},
        {"DepMap_ID": "ACH-LOW", "finalScore": 0.0, "rankMid": 3.0},
    ]
    metrics = evaluate_ranking(
        tie_ranked, {"ACH-POS", "ACH-MISSING"}, {"ACH-LOW"}, ks=(1, 2)
    )
    assert math.isclose(metrics["P@1"], 0.5)
    assert math.isclose(metrics["Recall@1"], 0.25)
    assert metrics["n_pos_found"] == 1
    assert math.isclose(metrics["MRR"], 2.0 / 3.0)
    assert math.isclose(metrics["neg_sink"], 1.0)
    print("self-test OK")
    print("  production scorer A0 == ablation A0")
    print("  high-correlation adaptive == baseline")
    print("  low-correlation adaptive changes ranking")
    print("  full no-protein ablation ignores all protein fields")
    print("  missing positives and exact-score ties are handled")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path)
    parser.add_argument(
        "--gold", type=Path, default=Path("benchmarks/gold_standard_v2.csv")
    )
    parser.add_argument(
        "--corr", type=Path, default=Path("gene_rna_protein_correlations.csv")
    )
    parser.add_argument("--ks", type=int, nargs="+", default=[3, 5, 10])
    parser.add_argument("--k", type=int, default=None, help="legacy single-K mode")
    parser.add_argument("--disease", default=None)
    parser.add_argument(
        "--out", type=Path, default=Path("results/ablation_results.csv")
    )
    parser.add_argument("--detail-out", type=Path, default=None)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        return self_test()
    if not args.data_dir:
        parser.error("--data-dir is required (or use --self-test)")
    _validate_configs()
    ks = tuple(sorted(set([args.k] if args.k is not None else args.ks)))
    if not ks or any(k <= 0 for k in ks):
        parser.error("all K values must be positive")

    repo = mcs.MergedDataRecommender(args.data_dir)
    gold = _prepare_gold(pd.read_csv(args.gold))
    corr_df = pd.read_csv(args.corr)
    if not {"gene", "correlation"}.issubset(corr_df.columns):
        raise KeyError("correlation table must contain gene and correlation columns")
    corr_table = dict(
        zip(
            corr_df["gene"].astype(str).str.upper().str.strip(),
            pd.to_numeric(corr_df["correlation"], errors="coerce"),
        )
    )
    summary, detail = run_ablation(
        repo,
        gold,
        corr_table,
        ks=ks,
        disease=args.disease,
        bootstrap_iterations=args.bootstrap,
        seed=args.seed,
    )
    if summary.empty:
        print("No evaluable genes were found.")
        return 1

    detail_out = args.detail_out or args.out.with_name(
        f"{args.out.stem}_by_gene{args.out.suffix}"
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    detail_out.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.out, index=False)
    detail.to_csv(detail_out, index=False)
    display_columns = [
        "config",
        "n_genes",
        "positive_coverage",
        *[column for k in ks for column in (f"P@{k}", f"Recall@{k}", f"NDCG@{k}")],
        "MRR",
        "neg_sink",
    ]
    print("\n" + summary[display_columns].to_string(index=False))
    print(f"\nsummary -> {args.out}")
    print(f"per-gene audit -> {detail_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
