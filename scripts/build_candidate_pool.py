#!/usr/bin/env python3
"""Build the pooled-judgement table used for gold-standard expansion.

For each review gene, this script rescales the same candidate evidence under
the selected ablation configurations, takes the union of their Top-K results
(including exact-score ties at the boundary), and appends any existing gold
entries. New review labels are deliberately initialised as unknown: model rank
must never be used as biological ground truth.

Example
-------
python scripts/build_candidate_pool.py \
  --data-dir "/absolute/path/to/merged" \
  --corr "/absolute/path/to/gene_rna_protein_correlations.csv" \
  --gold benchmarks/gold_standard_v2.csv \
  --out benchmarks/candidate_pool_v3.csv
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

import ablation_runner as ar


DEFAULT_CONFIGS = [
    "B0_rna_mean",
    "A0_team_baseline",
    "A1a_no_direct_protein",
    "A1b_no_protein_evidence",
    "A4_adaptive_trust",
]

# Six current benchmark genes plus six expression-oriented expansion genes.
DEFAULT_REVIEW_GENES = [
    "EGFR",
    "ERBB2",
    "ESR1",
    "AR",
    "MET",
    "ALB",
    "DES",
    "AFP",
    "GFAP",
    "PTPRC",
    "MSLN",
    "MITF",
]

MANUAL_REVIEW_COLUMNS = [
    "judgement",
    "benchmark_task",
    "evidence_type",
    "source_url",
    "evidence_summary",
    "verified",
    "review_notes",
]


def _normalise_genes(values: Iterable[str]) -> list[str]:
    genes: list[str] = []
    seen: set[str] = set()
    for value in values:
        gene = str(value).upper().strip()
        if gene and gene not in seen:
            genes.append(gene)
            seen.add(gene)
    return genes


def _top_with_ties(ranked: list[dict[str, Any]], top_k: int) -> list[dict[str, Any]]:
    """Return Top-K plus every row exactly tied with the Kth score."""
    if top_k <= 0:
        raise ValueError("top_k must be positive")
    if len(ranked) <= top_k:
        return list(ranked)

    cutoff = float(ranked[top_k - 1]["finalScore"])
    return [
        row
        for row in ranked
        if float(row["finalScore"]) > cutoff
        or math.isclose(float(row["finalScore"]), cutoff, abs_tol=1e-12)
    ]


def _read_correlations(path: Path) -> tuple[pd.DataFrame, dict[str, float]]:
    frame = pd.read_csv(path)
    required = {"gene", "correlation"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise KeyError("correlation table missing columns: " + ", ".join(missing))

    frame = frame.copy()
    frame["gene"] = frame["gene"].astype(str).str.upper().str.strip()
    if frame["gene"].duplicated().any():
        duplicates = sorted(frame.loc[frame["gene"].duplicated(False), "gene"].unique())
        raise ValueError("duplicate genes in correlation table: " + ", ".join(duplicates))

    table = dict(zip(frame["gene"], pd.to_numeric(frame["correlation"], errors="coerce")))
    return frame, table


def _prepare_reference_gold(path: Path) -> pd.DataFrame:
    gold = pd.read_csv(path)
    required = {
        "gene",
        "expected_cell_line",
        "expected_depmap_id",
        "relation",
    }
    missing = sorted(required - set(gold.columns))
    if missing:
        raise KeyError("gold standard missing columns: " + ", ".join(missing))

    gold = gold.copy()
    gold["gene"] = gold["gene"].astype(str).str.upper().str.strip()
    gold["expected_depmap_id"] = (
        gold["expected_depmap_id"].astype(str).str.upper().str.strip()
    )
    duplicated = gold.duplicated(["gene", "expected_depmap_id"], keep=False)
    if duplicated.any():
        rows = gold.loc[duplicated, ["gene", "expected_depmap_id"]]
        raise ValueError("duplicate gold entries:\n" + rows.to_string(index=False))
    return gold


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


def build_candidate_pool(
    repository: Any,
    gold: pd.DataFrame,
    correlation_table: dict[str, float],
    genes: list[str],
    config_names: list[str],
    top_k: int,
) -> pd.DataFrame:
    """Create one flat, auditable review table keyed by gene and DepMap ID."""
    records: list[dict[str, Any]] = []
    config_order = {name: index for index, name in enumerate(config_names)}

    for gene in genes:
        if not repository.check_gene_exists(gene):
            raise KeyError(f"review gene {gene!r} is not in the merged gene panel")

        evidence_rows = repository.fetch_candidate_evidence(
            target_gene=gene,
            disease_aliases=[],
            exclusion_gene=None,
        )
        evidence_by_id = {
            str(row["DepMap_ID"]).upper().strip(): row for row in evidence_rows
        }
        if len(evidence_by_id) != len(evidence_rows):
            raise ValueError(f"{gene}: duplicate candidate DepMap IDs in evidence")

        ranked_by_config: dict[str, list[dict[str, Any]]] = {}
        top_ids_by_config: dict[str, set[str]] = {}
        for config_name in config_names:
            ranked = ar.score_rows(
                evidence_rows,
                ar.CONFIGS[config_name],
                gene,
                correlation_table,
            )
            ranked_by_config[config_name] = ranked
            top_ids_by_config[config_name] = {
                str(row["DepMap_ID"]).upper().strip()
                for row in _top_with_ties(ranked, top_k)
            }

        gene_gold = gold.loc[gold["gene"] == gene].copy()
        gold_by_id = {
            str(row["expected_depmap_id"]).upper().strip(): row
            for _, row in gene_gold.iterrows()
        }

        selected_ids: set[str] = set(gold_by_id)
        for ids in top_ids_by_config.values():
            selected_ids.update(ids)

        rank_maps = {
            config_name: {
                str(row["DepMap_ID"]).upper().strip(): int(row["rank"])
                for row in ranked
            }
            for config_name, ranked in ranked_by_config.items()
        }
        score_maps = {
            config_name: {
                str(row["DepMap_ID"]).upper().strip(): float(row["finalScore"])
                for row in ranked
            }
            for config_name, ranked in ranked_by_config.items()
        }

        correlation = correlation_table.get(gene, float("nan"))
        trust = ar.rna_reliability(correlation)

        for depmap_id in selected_ids:
            evidence = evidence_by_id.get(depmap_id, {})
            reference = gold_by_id.get(depmap_id)
            retrieved_by = [
                name for name in config_names if depmap_id in top_ids_by_config[name]
            ]
            is_current_gold = reference is not None

            origin_parts: list[str] = []
            if retrieved_by:
                origin_parts.append("top_k_union")
            if is_current_gold:
                origin_parts.append("existing_gold")

            available_ranks = [
                rank_maps[name][depmap_id]
                for name in config_names
                if depmap_id in rank_maps[name]
            ]

            record: dict[str, Any] = {
                "gene": gene,
                "DepMap_ID": depmap_id,
                "cell_line": _safe_text(
                    evidence.get("cellLine")
                    or (reference.get("expected_cell_line") if reference is not None else "")
                ),
                "lineage": _safe_text(evidence.get("lineage")),
                "disease": _safe_text(evidence.get("disease")),
                "candidate_origin": "+".join(origin_parts),
                "configs_retrieved": ";".join(
                    sorted(retrieved_by, key=lambda name: config_order[name])
                ),
                "best_rank": min(available_ranks) if available_ranks else pd.NA,
                "rna_protein_correlation": (
                    float(correlation) if not pd.isna(correlation) else pd.NA
                ),
                "rna_trust": float(trust),
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
                "current_gold_relation": (
                    _safe_text(reference.get("relation")) if reference is not None else ""
                ),
                "current_gold_evidence_type": (
                    _safe_text(reference.get("evidence_type"))
                    if reference is not None
                    else ""
                ),
                "current_gold_verified": (
                    _safe_text(reference.get("verified")) if reference is not None else ""
                ),
                "current_gold_source_hint": (
                    _safe_text(reference.get("source_hint")) if reference is not None else ""
                ),
                "current_gold_notes": (
                    _safe_text(reference.get("notes")) if reference is not None else ""
                ),
                # Manual-review fields. Never infer these from ranking.
                "judgement": "unknown",
                "benchmark_task": "unassigned",
                "evidence_type": "",
                "source_url": "",
                "evidence_summary": "",
                "verified": "no",
                "review_notes": "",
            }

            for config_name in config_names:
                record[f"rank_{config_name}"] = rank_maps[config_name].get(
                    depmap_id, pd.NA
                )
                record[f"score_{config_name}"] = score_maps[config_name].get(
                    depmap_id, pd.NA
                )

            records.append(record)

    output = pd.DataFrame(records)
    if output.empty:
        raise ValueError("candidate pool is empty")

    if output.duplicated(["gene", "DepMap_ID"]).any():
        duplicates = output.loc[
            output.duplicated(["gene", "DepMap_ID"], keep=False),
            ["gene", "DepMap_ID"],
        ]
        raise ValueError("duplicate candidate-pool keys:\n" + duplicates.to_string(index=False))

    gene_order = {gene: index for index, gene in enumerate(genes)}
    output["_gene_order"] = output["gene"].map(gene_order)
    output = output.sort_values(
        ["_gene_order", "best_rank", "DepMap_ID"],
        na_position="last",
        kind="stable",
    ).drop(columns="_gene_order")

    fixed_columns = [
        "gene",
        "DepMap_ID",
        "cell_line",
        "lineage",
        "disease",
        "candidate_origin",
        "configs_retrieved",
        "best_rank",
        "rna_protein_correlation",
        "rna_trust",
        "rna_expr",
        "protein_expr",
        "n_rna_sources",
        "has_depmap_rna",
        "has_hpa_rna",
        "has_geo_rna",
        "has_protein",
        "alteration_status",
        "assay_ready_score",
        "risk_flags",
    ]
    rank_score_columns = [
        column
        for config_name in config_names
        for column in (f"rank_{config_name}", f"score_{config_name}")
    ]
    reference_columns = [
        "current_gold_relation",
        "current_gold_evidence_type",
        "current_gold_verified",
        "current_gold_source_hint",
        "current_gold_notes",
    ]
    return output[
        fixed_columns + rank_score_columns + reference_columns + MANUAL_REVIEW_COLUMNS
    ].reset_index(drop=True)


def self_test() -> int:
    ranked = [
        {"DepMap_ID": "A", "finalScore": 1.0},
        {"DepMap_ID": "B", "finalScore": 0.8},
        {"DepMap_ID": "C", "finalScore": 0.8},
        {"DepMap_ID": "D", "finalScore": 0.5},
    ]
    assert [row["DepMap_ID"] for row in _top_with_ties(ranked, 2)] == ["A", "B", "C"]
    assert _normalise_genes(["egfr", " EGFR ", "met"]) == ["EGFR", "MET"]
    print("self-test OK")
    print("  Top-K boundary ties are retained")
    print("  gene inputs are normalised and de-duplicated")
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
        "--gold",
        type=Path,
        default=Path("benchmarks/gold_standard_v2.csv"),
    )
    parser.add_argument("--corr", type=Path)
    parser.add_argument("--genes", nargs="+", default=DEFAULT_REVIEW_GENES)
    parser.add_argument("--configs", nargs="+", default=DEFAULT_CONFIGS)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("benchmarks/candidate_pool_v3.csv"),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite an existing output; this discards manual review fields",
    )
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        return self_test()
    if not args.data_dir:
        parser.error("--data-dir is required unless CELLLINESELECTOR_DATA_DIR is set")
    if not args.corr:
        parser.error("--corr is required because A4 needs gene correlations")
    if args.top_k <= 0:
        parser.error("--top-k must be positive")
    if args.out.exists() and not args.force:
        parser.error(
            f"output already exists: {args.out}; use a new --out path or --force"
        )

    genes = _normalise_genes(args.genes)
    if not genes:
        parser.error("at least one review gene is required")

    config_names = list(dict.fromkeys(args.configs))
    unknown_configs = sorted(set(config_names) - set(ar.CONFIGS))
    if unknown_configs:
        parser.error("unknown configuration(s): " + ", ".join(unknown_configs))

    ar._validate_configs()
    _, correlation_table = _read_correlations(args.corr)
    gold = _prepare_reference_gold(args.gold)
    repository = ar.mcs.MergedDataRecommender(args.data_dir)
    pool = build_candidate_pool(
        repository=repository,
        gold=gold,
        correlation_table=correlation_table,
        genes=genes,
        config_names=config_names,
        top_k=args.top_k,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = args.out.with_suffix(args.out.suffix + ".tmp")
    pool.to_csv(temporary_path, index=False, encoding="utf-8-sig")
    temporary_path.replace(args.out)

    summary = (
        pool.groupby("gene", sort=False)
        .agg(
            candidates=("DepMap_ID", "size"),
            existing_gold=("current_gold_relation", lambda values: (values != "").sum()),
            protein_available=("has_protein", lambda values: values.fillna(False).sum()),
        )
        .reset_index()
    )
    print("candidate pool created")
    print(summary.to_string(index=False))
    print(f"rows: {len(pool)}")
    print(f"output: {args.out.resolve()}")
    print("manual labels initialised as judgement=unknown, benchmark_task=unassigned")
    return 0


if __name__ == "__main__":
    sys.exit(main())
