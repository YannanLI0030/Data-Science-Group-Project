#!/usr/bin/env python3
"""Search the V4 outer-weight grid with leave-one-gene-out evaluation.

Only RNA, direct-Protein, and Confidence weights vary. Confidence internals,
the penalty cap, and adaptive trust remain fixed. Selection uses NDCG@5 and a
one-standard-error rule anchored to A1a.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

import ablation_runner as ar


PRIMARY_METRIC = "NDCG@5"
SECONDARY_METRICS = ["MRR", "Recall@5", "NDCG@10"]
A0_WEIGHTS = (0.55, 0.30, 0.15)
A1A_WEIGHTS = (0.85, 0.00, 0.15)
A0_ID = "W_rna0.55_protein0.30_confidence0.15"
A1A_ID = "W_rna0.85_protein0.00_confidence0.15"
CHALLENGER_ID = "W_rna0.70_protein0.05_confidence0.25"


def generate_grid() -> list[dict[str, Any]]:
    """Return the predeclared 0.05 weight grid."""
    configs: list[dict[str, Any]] = []
    for protein_pct in range(0, 31, 5):
        for confidence_pct in range(10, 31, 5):
            rna_pct = 100 - protein_pct - confidence_pct
            if not 50 <= rna_pct <= 90:
                continue
            rna = rna_pct / 100.0
            protein = protein_pct / 100.0
            confidence = confidence_pct / 100.0
            config_id = (
                f"W_rna{rna:.2f}_protein{protein:.2f}_confidence{confidence:.2f}"
            )
            configs.append(
                {
                    **ar.DEFAULT_CFG,
                    "config": config_id,
                    "rna_w": rna,
                    "protein_w": protein,
                    "confidence_w": confidence,
                    "adaptive_trust": False,
                    "protein_mode": "direct_off" if protein_pct == 0 else "full",
                }
            )
    ids = {config["config"] for config in configs}
    required = {A0_ID, A1A_ID, CHALLENGER_ID}
    if not required.issubset(ids):
        raise AssertionError("weight grid is missing required anchor configuration(s)")
    if len(ids) != len(configs):
        raise AssertionError("weight grid contains duplicate configuration IDs")
    return configs


def _distance_to_anchor(
    row: pd.Series | dict[str, Any], anchor: tuple[float, float, float]
) -> float:
    return float(
        abs(float(row["rna_w"]) - anchor[0])
        + abs(float(row["protein_w"]) - anchor[1])
        + abs(float(row["confidence_w"]) - anchor[2])
    )


def _bootstrap_delta(
    detail: pd.DataFrame,
    config: str,
    baseline: str,
    metric: str,
    iterations: int,
    seed: int,
) -> tuple[float, float, float]:
    current = detail.loc[
        detail["config"].eq(config), ["gene", metric]
    ].rename(columns={metric: "current"})
    reference = detail.loc[
        detail["config"].eq(baseline), ["gene", metric]
    ].rename(columns={metric: "baseline"})
    paired = current.merge(reference, on="gene", validate="one_to_one")
    values = (paired["current"] - paired["baseline"]).to_numpy(float)
    return ar._bootstrap_mean_delta(values, iterations, seed)


def _summarise(
    detail: pd.DataFrame, bootstrap: int, seed: int
) -> pd.DataFrame:
    metric_columns = [
        "P@3",
        "Recall@3",
        "NDCG@3",
        "P@5",
        "Recall@5",
        "NDCG@5",
        "P@10",
        "Recall@10",
        "NDCG@10",
        "MRR",
        "neg_sink",
    ]
    group_columns = [
        "config",
        "rna_w",
        "protein_w",
        "confidence_w",
        "protein_mode",
    ]
    summary = (
        detail.groupby(group_columns, sort=False, dropna=False)[metric_columns]
        .mean()
        .reset_index()
    )
    grouped = detail.groupby("config", sort=False)
    summary["n_genes"] = summary["config"].map(grouped["gene"].nunique())
    summary["n_gold_positive"] = summary["config"].map(
        grouped["n_pos_gold"].sum()
    )
    summary["n_positive_found"] = summary["config"].map(
        grouped["n_pos_found"].sum()
    )
    summary["positive_coverage"] = (
        summary["n_positive_found"] / summary["n_gold_positive"]
    )
    summary["n_gold_negative"] = summary["config"].map(
        grouped["n_neg_gold"].sum()
    )
    summary["n_negative_found"] = summary["config"].map(
        grouped["n_neg_found"].sum()
    )
    summary["distance_to_A0"] = summary.apply(
        lambda row: _distance_to_anchor(row, A0_WEIGHTS), axis=1
    )
    summary["distance_to_A1a"] = summary.apply(
        lambda row: _distance_to_anchor(row, A1A_WEIGHTS), axis=1
    )

    for baseline_name, baseline_id in (("A0", A0_ID), ("A1a", A1A_ID)):
        for metric_index, metric in enumerate((PRIMARY_METRIC, "MRR", "NDCG@10")):
            deltas = [
                _bootstrap_delta(
                    detail,
                    str(config),
                    baseline_id,
                    metric,
                    bootstrap,
                    seed + 1000 * metric_index + offset,
                )
                for offset, config in enumerate(summary["config"])
            ]
            summary[f"delta_{metric}_vs_{baseline_name}"] = [x[0] for x in deltas]
            summary[f"delta_{metric}_vs_{baseline_name}_ci_low"] = [
                x[1] for x in deltas
            ]
            summary[f"delta_{metric}_vs_{baseline_name}_ci_high"] = [
                x[2] for x in deltas
            ]
    return summary


def _one_se_choice(detail: pd.DataFrame) -> tuple[pd.Series, float, float, str]:
    """Choose the A1a-nearest model within one standard error of the best."""
    grouped = (
        detail.groupby(
            ["config", "rna_w", "protein_w", "confidence_w", "protein_mode"],
            sort=False,
            dropna=False,
        )[[PRIMARY_METRIC, *SECONDARY_METRICS]]
        .agg(["mean", "std", "count"])
    )
    grouped.columns = ["_".join(column) for column in grouped.columns]
    grouped = grouped.reset_index()
    ranked = grouped.sort_values(
        [
            f"{PRIMARY_METRIC}_mean",
            "MRR_mean",
            "Recall@5_mean",
            "NDCG@10_mean",
            "config",
        ],
        ascending=[False, False, False, False, True],
        kind="stable",
    )
    best = ranked.iloc[0]
    best_se = float(best[f"{PRIMARY_METRIC}_std"]) / math.sqrt(
        float(best[f"{PRIMARY_METRIC}_count"])
    )
    if not math.isfinite(best_se):
        best_se = 0.0
    threshold = float(best[f"{PRIMARY_METRIC}_mean"]) - best_se
    eligible = grouped[grouped[f"{PRIMARY_METRIC}_mean"] >= threshold].copy()
    eligible["distance_to_A1a"] = eligible.apply(
        lambda row: _distance_to_anchor(row, A1A_WEIGHTS), axis=1
    )
    eligible = eligible.sort_values(
        [
            "distance_to_A1a",
            "protein_w",
            f"{PRIMARY_METRIC}_mean",
            "MRR_mean",
            "Recall@5_mean",
            "NDCG@10_mean",
            "config",
        ],
        ascending=[True, True, False, False, False, False, True],
        kind="stable",
    )
    return eligible.iloc[0], threshold, best_se, str(best["config"])


def run_search(
    repository: Any,
    gold: pd.DataFrame,
    correlation_table: dict[str, float],
    bootstrap: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    configs = generate_grid()
    records: list[dict[str, Any]] = []

    for gene in sorted(gold["gene"].unique()):
        gene_gold = gold[gold["gene"].eq(gene)]
        positives = ar._norm_ids(
            gene_gold.loc[
                gene_gold["relation"].eq("positive"), "expected_depmap_id"
            ]
        )
        negatives = ar._norm_ids(
            gene_gold.loc[
                gene_gold["relation"].eq("negative"), "expected_depmap_id"
            ]
        )
        if not positives:
            continue
        rows = repository.fetch_candidate_evidence(gene, [""], None)
        if not rows:
            continue
        for config in configs:
            ranked = ar.score_rows(rows, config, gene, correlation_table)
            metrics = ar.evaluate_ranking(
                ranked, positives, negatives, ks=(3, 5, 10)
            )
            records.append(
                {
                    "config": config["config"],
                    "gene": gene,
                    "rna_w": config["rna_w"],
                    "protein_w": config["protein_w"],
                    "confidence_w": config["confidence_w"],
                    "protein_mode": config["protein_mode"],
                    **metrics,
                }
            )

    detail = pd.DataFrame(records)
    if detail.empty:
        raise ValueError("no genes could be evaluated")
    expected_genes = int(
        gold.loc[gold["relation"].eq("positive"), "gene"].nunique()
    )
    observed_genes = int(detail["gene"].nunique())
    if observed_genes != expected_genes:
        raise ValueError(
            f"evaluated {observed_genes} genes; expected {expected_genes}"
        )
    expected_rows = len(configs) * expected_genes
    if len(detail) != expected_rows:
        raise ValueError(f"detail has {len(detail)} rows; expected {expected_rows}")
    if detail.duplicated(["config", "gene"]).any():
        raise ValueError("detail contains duplicate configuration-gene rows")

    summary = _summarise(detail, bootstrap=bootstrap, seed=seed)
    full_choice, threshold, best_se, unconstrained_best = _one_se_choice(detail)
    selected_config = str(full_choice["config"])
    summary["unconstrained_best"] = summary["config"].eq(unconstrained_best)
    summary["best_NDCG@5_standard_error"] = best_se
    summary["one_se_threshold"] = threshold
    summary["eligible_one_se"] = summary[PRIMARY_METRIC] >= threshold
    summary["selected_one_se"] = summary["config"].eq(selected_config)
    summary = summary.sort_values(
        [PRIMARY_METRIC, "MRR", "Recall@5", "NDCG@10", "distance_to_A1a"],
        ascending=[False, False, False, False, True],
        kind="stable",
    ).reset_index(drop=True)

    logo_records: list[dict[str, Any]] = []
    for held_out_gene in sorted(detail["gene"].unique()):
        training = detail[~detail["gene"].eq(held_out_gene)]
        choice, fold_threshold, fold_best_se, fold_best = _one_se_choice(training)
        chosen_config = str(choice["config"])
        test = detail[
            detail["gene"].eq(held_out_gene)
            & detail["config"].eq(chosen_config)
        ].iloc[0]
        unconstrained_test = detail[
            detail["gene"].eq(held_out_gene)
            & detail["config"].eq(fold_best)
        ].iloc[0]
        a0_test = detail[
            detail["gene"].eq(held_out_gene) & detail["config"].eq(A0_ID)
        ].iloc[0]
        a1a_test = detail[
            detail["gene"].eq(held_out_gene) & detail["config"].eq(A1A_ID)
        ].iloc[0]
        logo_records.append(
            {
                "held_out_gene": held_out_gene,
                "selected_config": chosen_config,
                "rna_w": float(choice["rna_w"]),
                "protein_w": float(choice["protein_w"]),
                "confidence_w": float(choice["confidence_w"]),
                "protein_mode": str(choice["protein_mode"]),
                "training_n_genes": int(training["gene"].nunique()),
                "training_NDCG@5": float(choice["NDCG@5_mean"]),
                "training_unconstrained_best": fold_best,
                "training_one_se_threshold": fold_threshold,
                "training_best_se": fold_best_se,
                "held_out_P@5": test["P@5"],
                "held_out_Recall@5": test["Recall@5"],
                "held_out_NDCG@5": test["NDCG@5"],
                "held_out_NDCG@10": test["NDCG@10"],
                "held_out_MRR": test["MRR"],
                "held_out_neg_sink": test["neg_sink"],
                "A0_NDCG@5": a0_test["NDCG@5"],
                "A0_NDCG@10": a0_test["NDCG@10"],
                "A0_MRR": a0_test["MRR"],
                "A1a_NDCG@5": a1a_test["NDCG@5"],
                "A1a_NDCG@10": a1a_test["NDCG@10"],
                "A1a_MRR": a1a_test["MRR"],
                "delta_NDCG@5_vs_A0": test["NDCG@5"] - a0_test["NDCG@5"],
                "delta_NDCG@5_vs_A1a": test["NDCG@5"] - a1a_test["NDCG@5"],
                "delta_MRR_vs_A0": test["MRR"] - a0_test["MRR"],
                "delta_MRR_vs_A1a": test["MRR"] - a1a_test["MRR"],
                "unconstrained_rna_w": unconstrained_test["rna_w"],
                "unconstrained_protein_w": unconstrained_test["protein_w"],
                "unconstrained_confidence_w": unconstrained_test["confidence_w"],
                "unconstrained_protein_mode": unconstrained_test["protein_mode"],
                "unconstrained_held_out_P@5": unconstrained_test["P@5"],
                "unconstrained_held_out_Recall@5": unconstrained_test["Recall@5"],
                "unconstrained_held_out_NDCG@5": unconstrained_test["NDCG@5"],
                "unconstrained_held_out_NDCG@10": unconstrained_test["NDCG@10"],
                "unconstrained_held_out_MRR": unconstrained_test["MRR"],
                "unconstrained_delta_NDCG@5_vs_A0": (
                    unconstrained_test["NDCG@5"] - a0_test["NDCG@5"]
                ),
                "unconstrained_delta_NDCG@5_vs_A1a": (
                    unconstrained_test["NDCG@5"] - a1a_test["NDCG@5"]
                ),
                "unconstrained_delta_MRR_vs_A0": (
                    unconstrained_test["MRR"] - a0_test["MRR"]
                ),
                "unconstrained_delta_MRR_vs_A1a": (
                    unconstrained_test["MRR"] - a1a_test["MRR"]
                ),
            }
        )
    logo = pd.DataFrame(logo_records)

    selected_row = summary.loc[summary["config"].eq(selected_config)].iloc[0]
    eligible_count = int(summary["eligible_one_se"].sum())
    logo_selected_counts = {
        str(key): int(value)
        for key, value in logo["selected_config"].value_counts().items()
    }
    logo_unconstrained_counts = {
        str(key): int(value)
        for key, value in logo["training_unconstrained_best"].value_counts().items()
    }
    caution = (
        "Selected on 10 development genes with only 2 negative labels; "
        "do not describe the selected weights as a statistically unique optimum."
    )
    if eligible_count == len(summary):
        caution += (
            " The one-standard-error band contains every grid configuration, "
            "so the conservative selection is determined by proximity to the "
            "predeclared A1a anchor rather than clear separation in NDCG@5."
        )
    production_config = {
        "name": "selected_v4_development_logo_one_se_candidate",
        "status": "candidate_not_yet_frozen",
        "selection_dataset": "gold_standard_v4_development.csv",
        "selection_method": "predeclared_grid_plus_LOGO_one_standard_error",
        "primary_metric": PRIMARY_METRIC,
        "rna_weight": float(selected_row["rna_w"]),
        "protein_weight": float(selected_row["protein_w"]),
        "confidence_weight": float(selected_row["confidence_w"]),
        "protein_mode": str(selected_row["protein_mode"]),
        "confidence_completeness": float(ar.DEFAULT_CFG["conf_completeness"]),
        "confidence_source_support": float(
            ar.DEFAULT_CFG["conf_source_support"]
        ),
        "confidence_consistency": float(ar.DEFAULT_CFG["conf_consistency"]),
        "max_exclusion_penalty": float(ar.DEFAULT_CFG["exclusion_max_penalty"]),
        "adaptive_trust": False,
        "n_development_genes": expected_genes,
        "n_gold_positive": int((gold["relation"] == "positive").sum()),
        "n_gold_negative": int((gold["relation"] == "negative").sum()),
        "classic_case_policy": (
            "independent_post_freeze_evaluation_only_no_weight_retuning"
        ),
        "unconstrained_best_config": unconstrained_best,
        "one_se_threshold": float(threshold),
        "best_NDCG@5_standard_error": float(best_se),
        "one_se_eligible_config_count": eligible_count,
        "grid_config_count": int(len(summary)),
        "logo_one_se_selected_counts": logo_selected_counts,
        "logo_unconstrained_best_counts": logo_unconstrained_counts,
        "caution": caution,
    }
    return summary, detail, logo, production_config


def _atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False, encoding="utf-8-sig")
    temporary.replace(path)


def _atomic_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as file:
        json.dump(payload, file, ensure_ascii=False, indent=2)
        file.write("\n")
    temporary.replace(path)


def self_test() -> int:
    grid = generate_grid()
    assert len(grid) == 32
    assert all(
        math.isclose(
            config["rna_w"] + config["protein_w"] + config["confidence_w"],
            1.0,
            abs_tol=1e-12,
        )
        for config in grid
    )
    assert all(0.50 <= config["rna_w"] <= 0.90 for config in grid)
    assert all(0.00 <= config["protein_w"] <= 0.30 for config in grid)
    assert all(0.10 <= config["confidence_w"] <= 0.30 for config in grid)
    assert all(
        config["protein_mode"]
        == ("direct_off" if math.isclose(config["protein_w"], 0.0) else "full")
        for config in grid
    )
    ids = {config["config"] for config in grid}
    assert {A0_ID, A1A_ID, CHALLENGER_ID}.issubset(ids)
    print("self-test OK")
    print(f"  valid constrained configurations: {len(grid)}")
    print("  A0, A1a and the reviewed challenger are included")
    print("  every outer-weight sum equals one")
    print("  zero direct-protein weight retains protein-derived confidence")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path)
    parser.add_argument(
        "--gold",
        type=Path,
        default=Path("benchmarks/gold_standard_v4_development.csv"),
    )
    parser.add_argument("--corr", type=Path)
    parser.add_argument(
        "--summary-out",
        type=Path,
        default=Path("results/weight_search_v4_development_summary.csv"),
    )
    parser.add_argument(
        "--detail-out",
        type=Path,
        default=Path("results/weight_search_v4_development_by_gene.csv"),
    )
    parser.add_argument(
        "--logo-out",
        type=Path,
        default=Path("results/weight_search_v4_development_logo.csv"),
    )
    parser.add_argument(
        "--config-out",
        type=Path,
        default=Path("config/scoring_selected_v4_development_candidate.json"),
    )
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        return self_test()
    if not args.data_dir:
        parser.error("--data-dir is required")
    if not args.corr:
        parser.error("--corr is required")
    outputs = [args.summary_out, args.detail_out, args.logo_out, args.config_out]
    existing = [str(path) for path in outputs if path.exists()]
    if existing and not args.force:
        parser.error("output(s) already exist; use --force: " + ", ".join(existing))

    ar._validate_configs()
    gold = ar._prepare_gold(pd.read_csv(args.gold))
    corr = pd.read_csv(args.corr)
    if not {"gene", "correlation"}.issubset(corr.columns):
        raise KeyError("correlation table must contain gene and correlation columns")
    correlation_table = dict(
        zip(
            corr["gene"].astype(str).str.upper().str.strip(),
            pd.to_numeric(corr["correlation"], errors="coerce"),
        )
    )
    repository = ar.mcs.MergedDataRecommender(args.data_dir)
    summary, detail, logo, config = run_search(
        repository,
        gold,
        correlation_table,
        bootstrap=args.bootstrap,
        seed=args.seed,
    )

    _atomic_csv(summary, args.summary_out)
    _atomic_csv(detail, args.detail_out)
    _atomic_csv(logo, args.logo_out)
    _atomic_json(config, args.config_out)

    selected = summary.loc[summary["selected_one_se"]].iloc[0]
    best = summary.loc[summary["unconstrained_best"]].iloc[0]
    print(f"grid configurations: {summary.shape[0]}")
    print("unconstrained best:")
    print(
        best[
            ["config", "NDCG@5", "Recall@5", "NDCG@10", "MRR"]
        ].to_string()
    )
    print("selected by one-SE rule:")
    print(
        selected[
            [
                "config",
                "rna_w",
                "protein_w",
                "confidence_w",
                "protein_mode",
                "NDCG@5",
                "Recall@5",
                "NDCG@10",
                "MRR",
                "delta_NDCG@5_vs_A0",
                "delta_NDCG@5_vs_A1a",
            ]
        ].to_string()
    )
    print(
        "one-SE eligible configurations: "
        f"{int(summary['eligible_one_se'].sum())}/{len(summary)}"
    )
    print("LOGO held-out macro means:")
    print(
        logo[
            [
                "held_out_Recall@5",
                "held_out_NDCG@5",
                "held_out_NDCG@10",
                "held_out_MRR",
                "delta_NDCG@5_vs_A0",
                "delta_NDCG@5_vs_A1a",
            ]
        ]
        .mean()
        .to_string()
    )
    print("LOGO selected configuration counts:")
    print(logo["selected_config"].value_counts().to_string())
    print("LOGO unconstrained-best held-out macro means:")
    print(
        logo[
            [
                "unconstrained_held_out_Recall@5",
                "unconstrained_held_out_NDCG@5",
                "unconstrained_held_out_NDCG@10",
                "unconstrained_held_out_MRR",
                "unconstrained_delta_NDCG@5_vs_A0",
                "unconstrained_delta_NDCG@5_vs_A1a",
            ]
        ]
        .mean()
        .to_string()
    )
    print("LOGO unconstrained-best configuration counts:")
    print(logo["training_unconstrained_best"].value_counts().to_string())
    print(f"summary -> {args.summary_out}")
    print(f"per-gene -> {args.detail_out}")
    print(f"LOGO -> {args.logo_out}")
    print(f"candidate config -> {args.config_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
