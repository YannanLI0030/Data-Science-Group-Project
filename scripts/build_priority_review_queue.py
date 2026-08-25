#!/usr/bin/env python3
"""Build the high-information manual-review queue for weight selection.

The queue contains unverified candidates that occur in exactly one of the
Top-K sets from the conservative A0 baseline and the exploratory challenger.
These disagreements are the smallest useful judgement set for deciding which
weighting scheme is better. Ranking output never determines the biological
label: all manual-review fields remain unknown/unverified.
"""

from __future__ import annotations

import argparse
import math
import os
import re
from pathlib import Path
from typing import Any

import pandas as pd

import ablation_runner as ar


DEFAULT_CHALLENGER = {
    **ar.DEFAULT_CFG,
    "config": "W_rna0.70_protein0.05_confidence0.25",
    "rna_w": 0.70,
    "protein_w": 0.05,
    "confidence_w": 0.25,
    "adaptive_trust": False,
    "protein_mode": "full",
}

MANUAL_COLUMNS = [
    "judgement",
    "benchmark_task",
    "evidence_type",
    "source_url",
    "evidence_summary",
    "verified",
    "review_notes",
]


def _safe_text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    if isinstance(value, (list, tuple, set)):
        return ";".join(str(item) for item in value)
    return str(value)


def _read_correlations(path: Path) -> dict[str, float]:
    frame = pd.read_csv(path)
    required = {"gene", "correlation"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise KeyError("correlation table missing columns: " + ", ".join(missing))
    frame = frame.copy()
    frame["gene"] = frame["gene"].astype(str).str.upper().str.strip()
    if frame["gene"].duplicated().any():
        duplicates = sorted(frame.loc[frame["gene"].duplicated(False), "gene"].unique())
        raise ValueError("duplicate correlation genes: " + ", ".join(duplicates))
    return dict(
        zip(frame["gene"], pd.to_numeric(frame["correlation"], errors="coerce"))
    )


def _read_candidate_pool(path: Path) -> pd.DataFrame:
    pool = pd.read_csv(path, keep_default_na=False)
    required = {"gene", "DepMap_ID", "verified", *MANUAL_COLUMNS}
    missing = sorted(required - set(pool.columns))
    if missing:
        raise KeyError("candidate pool missing columns: " + ", ".join(missing))
    pool = pool.copy()
    pool["gene"] = pool["gene"].astype(str).str.upper().str.strip()
    pool["DepMap_ID"] = pool["DepMap_ID"].astype(str).str.upper().str.strip()
    if pool.duplicated(["gene", "DepMap_ID"]).any():
        raise ValueError("candidate pool contains duplicate gene/DepMap-ID keys")
    return pool


def _rank_maps(ranked: list[dict[str, Any]]) -> tuple[dict[str, int], dict[str, float]]:
    ranks = {
        str(row["DepMap_ID"]).upper().strip(): int(row["rank"])
        for row in ranked
    }
    scores = {
        str(row["DepMap_ID"]).upper().strip(): float(row["finalScore"])
        for row in ranked
    }
    return ranks, scores


def build_priority_queue(
    repository: Any,
    pool: pd.DataFrame,
    correlation_table: dict[str, float],
    top_k: int,
    challenger: dict[str, Any],
) -> pd.DataFrame:
    if top_k <= 0:
        raise ValueError("top_k must be positive")

    pool_by_key = {
        (row["gene"], row["DepMap_ID"]): row for _, row in pool.iterrows()
    }
    verified_keys = {
        key
        for key, row in pool_by_key.items()
        if str(row["verified"]).lower().strip() == "yes"
    }
    records: list[dict[str, Any]] = []

    for gene in pool["gene"].drop_duplicates():
        if not repository.check_gene_exists(gene):
            raise KeyError(f"candidate-pool gene {gene!r} is not in merged data")

        evidence_rows = repository.fetch_candidate_evidence(gene, [], None)
        evidence_by_id = {
            str(row["DepMap_ID"]).upper().strip(): row for row in evidence_rows
        }
        if len(evidence_by_id) != len(evidence_rows):
            raise ValueError(f"{gene}: duplicate candidate DepMap IDs in evidence")

        ranked_a0 = ar.score_rows(
            evidence_rows,
            ar.CONFIGS["A0_team_baseline"],
            gene,
            correlation_table,
        )
        ranked_challenger = ar.score_rows(
            evidence_rows,
            challenger,
            gene,
            correlation_table,
        )
        top_a0 = ar._topk_with_ties(ranked_a0, top_k)
        top_challenger = ar._topk_with_ties(ranked_challenger, top_k)
        disagreement_ids = top_a0 ^ top_challenger

        a0_ranks, a0_scores = _rank_maps(ranked_a0)
        challenger_ranks, challenger_scores = _rank_maps(ranked_challenger)
        correlation = correlation_table.get(gene, float("nan"))

        for depmap_id in disagreement_ids:
            key = (gene, depmap_id)
            if key in verified_keys:
                continue

            evidence = evidence_by_id[depmap_id]
            pool_row = pool_by_key.get(key)
            in_pool = pool_row is not None
            a0_only = depmap_id in top_a0
            side = "A0_only" if a0_only else "challenger_only"
            rank_a0 = a0_ranks[depmap_id]
            rank_challenger = challenger_ranks[depmap_id]

            record: dict[str, Any] = {
                "review_priority": 1 if not in_pool else 2,
                "priority_reason": (
                    "Top-10 disagreement and absent from candidate_pool_v3"
                    if not in_pool
                    else "Top-10 disagreement; unverified in candidate_pool_v3"
                ),
                "gene": gene,
                "DepMap_ID": depmap_id,
                "cell_line": _safe_text(evidence.get("cellLine")),
                "lineage": _safe_text(evidence.get("lineage")),
                "disease": _safe_text(evidence.get("disease")),
                "review_side": side,
                "in_candidate_pool_v3": in_pool,
                "rank_A0": rank_a0,
                "score_A0": a0_scores[depmap_id],
                "rank_challenger": rank_challenger,
                "score_challenger": challenger_scores[depmap_id],
                "rank_change_challenger_minus_A0": rank_challenger - rank_a0,
                "rna_protein_correlation": (
                    float(correlation) if not pd.isna(correlation) else pd.NA
                ),
                "rna_expr": evidence.get("rnaExpr"),
                "protein_expr": evidence.get("protExpr"),
                "n_rna_sources": evidence.get("nRna"),
                "has_depmap_rna": evidence.get("hasDepMapRNA"),
                "has_hpa_rna": evidence.get("hasHpaRNA"),
                "has_geo_rna": evidence.get("hasGeoRNA"),
                "has_protein": evidence.get("hasProteomics"),
                "alteration_status": _safe_text(evidence.get("alterationStatus")),
                "assay_ready_score": evidence.get("assayReadyScore"),
                "risk_flags": _safe_text(evidence.get("riskFlags")),
            }

            for column in MANUAL_COLUMNS:
                if pool_row is not None and _safe_text(pool_row.get(column)):
                    record[column] = _safe_text(pool_row.get(column))
                else:
                    record[column] = ""
            record["judgement"] = "unknown"
            record["benchmark_task"] = (
                record["benchmark_task"] if record["benchmark_task"] else "unassigned"
            )
            record["verified"] = "no"
            records.append(record)

    output = pd.DataFrame(records)
    if output.empty:
        raise ValueError("priority review queue is empty")
    if output.duplicated(["gene", "DepMap_ID"]).any():
        raise ValueError("priority queue contains duplicate gene/DepMap-ID keys")
    if not output["DepMap_ID"].str.fullmatch(r"ACH-\d{6}").all():
        raise ValueError("priority queue contains invalid DepMap IDs")
    if not output["judgement"].eq("unknown").all():
        raise ValueError("priority queue must not infer biological judgements")
    if not output["verified"].eq("no").all():
        raise ValueError("priority queue must contain only unverified rows")

    output["_best_disagreement_rank"] = output[["rank_A0", "rank_challenger"]].min(
        axis=1
    )
    output = output.sort_values(
        ["review_priority", "gene", "_best_disagreement_rank", "DepMap_ID"],
        kind="stable",
    ).drop(columns="_best_disagreement_rank")
    return output.reset_index(drop=True)


def self_test() -> int:
    total = sum(DEFAULT_CHALLENGER[key] for key in ("rna_w", "protein_w", "confidence_w"))
    assert math.isclose(total, 1.0, abs_tol=1e-9)
    assert DEFAULT_CHALLENGER["protein_mode"] == "full"
    assert re.fullmatch(r"ACH-\d{6}", "ACH-000832")
    print("self-test OK")
    print("  challenger weights sum to one")
    print("  challenger retains full protein/confidence structure")
    print("  DepMap ID validation is active")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=os.getenv("CELLLINESELECTOR_DATA_DIR"),
        help="merged-data directory; or set CELLLINESELECTOR_DATA_DIR",
    )
    parser.add_argument(
        "--candidate-pool",
        type=Path,
        default=Path("benchmarks/candidate_pool_v3.csv"),
    )
    parser.add_argument("--corr", type=Path, required=False)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--expected-rows", type=int, default=34)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("benchmarks/priority_review_queue_v4.csv"),
    )
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        return self_test()
    if not args.data_dir:
        parser.error("--data-dir is required (or set CELLLINESELECTOR_DATA_DIR)")
    if not args.corr:
        parser.error("--corr is required")

    import merged_cellline_selector as mcs

    pool = _read_candidate_pool(args.candidate_pool)
    correlations = _read_correlations(args.corr)
    repository = mcs.MergedDataRecommender(args.data_dir)
    queue = build_priority_queue(
        repository,
        pool,
        correlations,
        args.top_k,
        DEFAULT_CHALLENGER,
    )
    if args.expected_rows > 0 and len(queue) != args.expected_rows:
        raise ValueError(
            f"expected {args.expected_rows} queue rows, generated {len(queue)}; "
            "review candidate-pool or ranking changes before exporting"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    queue.to_csv(args.out, index=False, lineterminator="\n")
    print(f"priority review queue -> {args.out}")
    print(f"rows = {len(queue)}; genes = {queue['gene'].nunique()}")
    print(
        "outside candidate_pool_v3 = "
        f"{int((~queue['in_candidate_pool_v3']).sum())}"
    )
    print("\nby gene and review side:")
    print(
        queue.groupby(["gene", "review_side"])
        .size()
        .unstack(fill_value=0)
        .to_string()
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
