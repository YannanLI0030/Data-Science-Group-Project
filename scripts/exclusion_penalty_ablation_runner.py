#!/usr/bin/env python3
"""Run a label-free structural ablation of the exclusion penalty.

This supplementary experiment is deliberately separate from the frozen V5
holdout.  It imports the current dynamic production scorer, fetches each
pre-specified target/exclusion query once, and re-scores the identical candidate
rows under P00/P15/P30/P45.  The outputs quantify implementation correctness,
clipping, ties, and ranking sensitivity.  No Gold file is read and no accuracy
metric is calculated.

The production scorer currently exposes the penalty cap as a module constant.
This runner changes that constant only in memory, sequentially, and restores it
with ``try/finally`` after every call.  It never edits the production script.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import subprocess
import sys
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DESIGN = PROJECT_ROOT / "config" / "exclusion_penalty_structural_v1_frozen_design.json"
DEFAULT_QUERIES = PROJECT_ROOT / "benchmarks" / "penalty" / "penalty_structural_v1_queries.csv"
DEFAULT_OUT_DIR = PROJECT_ROOT / "results" / "exclusion_penalty_structural_v1"
DEFAULT_PRODUCTION_SCRIPT = Path(
    "/Users/liyannan/Desktop/cellline_selector_v3/dynamic_cellline_selector_gene_protein.py"
)
DEFAULT_RAW_DATA_DIR = Path("/Users/liyannan/Desktop/cellline_selector_v3/data_s3")
DEFAULT_CACHE_DIR = Path("/Users/liyannan/Desktop/cellline_selector_v3/gene_cache")
DEFAULT_V5_FREEZE = PROJECT_ROOT / "config" / "ablation_v5_frozen_design.json"

PROFILE_ORDER = (
    "P00_no_penalty",
    "P15_weak",
    "P30_team_baseline",
    "P45_stress",
)
SCORE_COMPONENT_FIELDS = (
    "rnaScore",
    "proteinScore",
    "biologicalScore",
    "completenessScore",
    "sourceSupportScore",
    "rnaProteinConsistencyScore",
    "confidenceScore",
)
P30_PARITY_FIELDS = (*SCORE_COMPONENT_FIELDS, "exclusionPenalty", "finalScore")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", type=Path, default=DEFAULT_DESIGN)
    parser.add_argument("--queries", type=Path, default=DEFAULT_QUERIES)
    parser.add_argument("--production-script", type=Path, default=DEFAULT_PRODUCTION_SCRIPT)
    parser.add_argument("--raw-data-dir", type=Path, default=DEFAULT_RAW_DATA_DIR)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--v5-freeze", type=Path, default=DEFAULT_V5_FREEZE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def git_state(repo: Path) -> dict[str, Any]:
    def run(*args: str) -> str | None:
        try:
            result = subprocess.run(
                ["git", *args],
                cwd=repo,
                check=True,
                capture_output=True,
                text=True,
            )
            return result.stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            return None

    return {
        "commit": run("rev-parse", "HEAD"),
        "branch": run("branch", "--show-current"),
        "status_short": run("status", "--short"),
    }


def load_python_module(path: Path):
    spec = importlib.util.spec_from_file_location("penalty_ablation_production", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import production script: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_design(path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    design = json.loads(path.read_text(encoding="utf-8"))
    profiles = design.get("profiles")
    if not isinstance(profiles, list) or not profiles:
        raise ValueError("Frozen design must contain a non-empty profiles list")
    observed_names = tuple(str(item.get("name")) for item in profiles)
    if observed_names != PROFILE_ORDER:
        raise ValueError(
            f"Profile order must be {PROFILE_ORDER}; observed {observed_names}"
        )
    caps = []
    for profile in profiles:
        cap = float(profile["max_exclusion_penalty"])
        if not math.isfinite(cap) or cap < 0:
            raise ValueError(f"Invalid penalty cap in {profile}")
        profile["max_exclusion_penalty"] = cap
        caps.append(cap)
    if caps != sorted(caps) or caps != [0.0, 0.15, 0.30, 0.45]:
        raise ValueError(f"Frozen caps must be [0, .15, .30, .45]; observed {caps}")
    if design.get("uses_gold_standard") is not False:
        raise ValueError("Structural design must explicitly disable Gold-standard use")
    return design, profiles


def load_queries(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    required = {
        "query_id",
        "target_gene",
        "exclusion_gene",
        "disease",
        "query_mode",
        "scenario_rationale",
        "selection_basis",
        "expected_candidates",
        "expected_exclusion_observed",
    }
    missing = sorted(required - set(table.columns))
    if missing:
        raise KeyError(f"Query table missing columns: {missing}")
    table = table.copy()
    for column in ("query_id", "target_gene", "exclusion_gene", "disease", "query_mode"):
        table[column] = table[column].astype(str).str.strip()
    table["target_gene"] = table["target_gene"].str.upper()
    table["exclusion_gene"] = table["exclusion_gene"].str.upper()
    if table["query_id"].duplicated().any():
        raise ValueError("Query IDs must be unique")
    if (table["target_gene"] == table["exclusion_gene"]).any():
        bad = table.loc[table["target_gene"] == table["exclusion_gene"], "query_id"].tolist()
        raise ValueError(f"Target and exclusion genes must differ: {bad}")
    if set(table["query_mode"]) != {"GENE_MULTIOMICS"}:
        raise ValueError("The frozen structural experiment supports GENE_MULTIOMICS only")
    for column in ("expected_candidates", "expected_exclusion_observed"):
        table[column] = pd.to_numeric(table[column], errors="raise").astype(int)
    if (table["expected_candidates"] < 10).any():
        raise ValueError("Every frozen query must have at least ten candidates")
    return table


def required_cache_files(queries: pd.DataFrame, cache_dir: Path) -> list[Path]:
    genes = sorted(set(queries["target_gene"]) | set(queries["exclusion_gene"]))
    paths = [cache_dir / f"{gene}.parquet" for gene in genes]
    missing = [path for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Frozen structural run requires existing caches and will not build them:\n  "
            + "\n  ".join(str(path) for path in missing)
        )
    return paths


def ensure_unique_finite_rows(rows: list[dict[str, Any]], query_id: str) -> None:
    ids = [str(row.get("DepMap_ID")) for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError(f"{query_id}: duplicate DepMap_ID rows in candidate set")
    for row in rows:
        if row.get("nExclusion", 0):
            value = row.get("exclusionExpr")
            if value is None or not math.isfinite(float(value)):
                raise ValueError(
                    f"{query_id}/{row.get('DepMap_ID')}: non-finite observed exclusionExpr"
                )


@contextmanager
def temporary_penalty_cap(production_module, cap: float):
    cap = float(cap)
    if not math.isfinite(cap) or cap < 0:
        raise ValueError("Penalty cap must be finite and non-negative")
    original = float(production_module.MAX_EXCLUSION_PENALTY)
    production_module.MAX_EXCLUSION_PENALTY = cap
    try:
        yield
    finally:
        production_module.MAX_EXCLUSION_PENALTY = original


def score_with_cap(
    production_module,
    rows: list[dict[str, Any]],
    query_mode: str,
    cap: float,
) -> list[dict[str, Any]]:
    with temporary_penalty_cap(production_module, cap):
        return production_module.score_candidates(rows, query_mode=query_mode)


def add_tie_ranks(ranked: list[dict[str, Any]]) -> None:
    """Add score-tie ranks without replacing the production ordinal rank."""
    position = 0
    while position < len(ranked):
        score = ranked[position]["finalScore"]
        end = position + 1
        while end < len(ranked) and ranked[end]["finalScore"] == score:
            end += 1
        rank_min = position + 1
        rank_max = end
        rank_mid = (rank_min + rank_max) / 2.0
        group_size = end - position
        for row in ranked[position:end]:
            row["scoreTieRankMin"] = rank_min
            row["scoreTieRankMax"] = rank_max
            row["scoreTieRankMid"] = rank_mid
            row["finalScoreTieSize"] = group_size
        position = end

    comparator_counts: dict[tuple[Any, Any], int] = {}
    for row in ranked:
        key = (row["finalScore"], row["confidenceScore"])
        comparator_counts[key] = comparator_counts.get(key, 0) + 1
    for row in ranked:
        row["fullComparatorTieSize"] = comparator_counts[
            (row["finalScore"], row["confidenceScore"])
        ]


def topk_with_ties(ranked: list[dict[str, Any]], top_k: int) -> list[dict[str, Any]]:
    if not ranked:
        return []
    cutoff = ranked[min(top_k, len(ranked)) - 1]["finalScore"]
    return [row for row in ranked if row["finalScore"] >= cutoff]


def rows_by_id(rows: Sequence[Mapping[str, Any]]) -> dict[str, Mapping[str, Any]]:
    return {str(row["DepMap_ID"]): row for row in rows}


def values_equal(left: Any, right: Any, tolerance: float = 0.0) -> bool:
    if left is None or right is None:
        return left is None and right is None
    if isinstance(left, (float, int, np.floating, np.integer)) and isinstance(
        right, (float, int, np.floating, np.integer)
    ):
        if math.isnan(float(left)) or math.isnan(float(right)):
            return math.isnan(float(left)) and math.isnan(float(right))
        return math.isclose(float(left), float(right), abs_tol=tolerance, rel_tol=0.0)
    return left == right


def exact_scored_equal(
    left: Sequence[Mapping[str, Any]],
    right: Sequence[Mapping[str, Any]],
    fields: Iterable[str],
) -> tuple[bool, str]:
    left_map = rows_by_id(left)
    right_map = rows_by_id(right)
    if left_map.keys() != right_map.keys():
        return False, "candidate ID sets differ"
    for depmap_id in sorted(left_map):
        for field in fields:
            if not values_equal(left_map[depmap_id].get(field), right_map[depmap_id].get(field)):
                return (
                    False,
                    f"{depmap_id}/{field}: {left_map[depmap_id].get(field)!r} != "
                    f"{right_map[depmap_id].get(field)!r}",
                )
    return True, "all compared fields match"


def rank_spearman(current: Sequence[Mapping[str, Any]], baseline: Sequence[Mapping[str, Any]]) -> float:
    current_map = {str(row["DepMap_ID"]): float(row["scoreTieRankMid"]) for row in current}
    baseline_map = {str(row["DepMap_ID"]): float(row["scoreTieRankMid"]) for row in baseline}
    keys = sorted(set(current_map) & set(baseline_map))
    if len(keys) < 2:
        return float("nan")
    left = np.asarray([current_map[key] for key in keys], dtype=float)
    right = np.asarray([baseline_map[key] for key in keys], dtype=float)
    if np.array_equal(left, right):
        return 1.0
    if left.std(ddof=0) == 0 or right.std(ddof=0) == 0:
        return float("nan")
    return float(np.corrcoef(left, right)[0, 1])


def jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else float("nan")


def tie_metrics(ranked: list[dict[str, Any]], top_k: int) -> dict[str, Any]:
    counts = pd.Series([row["finalScore"] for row in ranked], dtype=float).value_counts()
    top = topk_with_ties(ranked, top_k)
    cutoff = ranked[min(top_k, len(ranked)) - 1]["finalScore"]
    cutoff_tie = sum(row["finalScore"] == cutoff for row in ranked)
    return {
        "n_unique_final_scores": int(counts.size),
        "largest_final_score_tie": int(counts.max()),
        "tied_candidate_fraction": float(counts[counts > 1].sum() / len(ranked)),
        f"top{top_k}_with_ties_n": len(top),
        f"top{top_k}_cutoff_score": float(cutoff),
        f"top{top_k}_cutoff_tie_size": int(cutoff_tie),
        "top1_score_tie_size": int(ranked[0]["finalScoreTieSize"]),
        "largest_full_comparator_tie": int(
            max(int(row["fullComparatorTieSize"]) for row in ranked)
        ),
    }


def comparison_metrics(
    current: list[dict[str, Any]],
    baseline: list[dict[str, Any]],
    top_k: int,
) -> dict[str, Any]:
    current_map = rows_by_id(current)
    baseline_map = rows_by_id(baseline)
    current_top = {str(row["DepMap_ID"]) for row in topk_with_ties(current, top_k)}
    baseline_top = {str(row["DepMap_ID"]) for row in topk_with_ties(baseline, top_k)}
    rank_deltas = np.asarray(
        [
            float(current_map[key]["scoreTieRankMid"])
            - float(baseline_map[key]["scoreTieRankMid"])
            for key in sorted(current_map)
        ],
        dtype=float,
    )
    current_top1 = {
        str(row["DepMap_ID"])
        for row in current
        if row["finalScore"] == current[0]["finalScore"]
    }
    baseline_top1 = {
        str(row["DepMap_ID"])
        for row in baseline
        if row["finalScore"] == baseline[0]["finalScore"]
    }
    return {
        "candidate_set_same_as_P00": current_map.keys() == baseline_map.keys(),
        "ranking_identical_to_P00": [row["DepMap_ID"] for row in current]
        == [row["DepMap_ID"] for row in baseline],
        f"top{top_k}_jaccard_vs_P00": jaccard(current_top, baseline_top),
        f"top{top_k}_entered_vs_P00": len(current_top - baseline_top),
        f"top{top_k}_left_vs_P00": len(baseline_top - current_top),
        "top1_identical_to_P00": current_top1 == baseline_top1,
        "rank_spearman_vs_P00": rank_spearman(current, baseline),
        "mean_abs_rank_change_vs_P00": float(np.abs(rank_deltas).mean()),
        "max_abs_rank_change_vs_P00": float(np.abs(rank_deltas).max()),
        "n_ranked_down_vs_P00": int((rank_deltas > 0).sum()),
        "n_ranked_up_vs_P00": int((rank_deltas < 0).sum()),
        "n_rank_unchanged_vs_P00": int((rank_deltas == 0).sum()),
    }


def add_check(
    checks: list[dict[str, Any]],
    *,
    scope: str,
    check_name: str,
    passed: bool,
    observed: Any,
    expected: Any,
    query_id: str = "",
    profile: str = "",
    details: str = "",
) -> None:
    checks.append(
        {
            "scope": scope,
            "query_id": query_id,
            "profile": profile,
            "check_name": check_name,
            "status": "PASS" if passed else "FAIL",
            "observed": observed,
            "expected": expected,
            "details": details,
        }
    )


def synthetic_rows(exclusion_values: Sequence[float | None]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for index, exclusion in enumerate(exclusion_values, start=1):
        rows.append(
            {
                "DepMap_ID": f"ACH-9{index:05d}",
                "cellLine": f"synthetic-{index}",
                "rnaExpr": 2.0,
                "protExpr": 2.0,
                "exclusionExpr": exclusion,
                "nRna": 3,
                "nProt": 1,
                "nExclusion": int(exclusion is not None),
                "hasDepMapRNA": True,
                "hasHpaRNA": True,
                "hasGeoRNA": True,
                "hasProteomics": True,
            }
        )
    return rows


def run_self_test(production_module, verbose: bool = True) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    original_cap = float(production_module.MAX_EXCLUSION_PENALTY)
    if not math.isclose(original_cap, 0.30, abs_tol=1e-12):
        raise AssertionError(f"Production default penalty must be 0.30; observed {original_cap}")

    known = synthetic_rows([0.0, 5.0, 10.0, None])
    expected_scaled = [0.0, 0.5, 1.0, None]
    observed_scaled = production_module.scale_available_values(
        known, "exclusionExpr", "nExclusion"
    )
    scaled_ok = all(
        (expected is None and observed is None)
        or (expected is not None and observed is not None and math.isclose(expected, observed))
        for expected, observed in zip(expected_scaled, observed_scaled)
    )
    add_check(
        checks,
        scope="synthetic",
        check_name="known_vector_minmax",
        passed=scaled_ok,
        observed=observed_scaled,
        expected=expected_scaled,
    )

    for mode in ("GENE_MULTIOMICS", "PROTEIN_ONLY", "COMBINED"):
        p00_scored = score_with_cap(production_module, known, mode, 0.0)
        scored = score_with_cap(production_module, known, mode, 0.30)
        p00_map = rows_by_id(p00_scored)
        penalties = rows_by_id(scored)
        expected_penalties = [0.0, 0.15, 0.30, 0.0]
        observed = [penalties[row["DepMap_ID"]]["exclusionPenalty"] for row in known]
        ok = observed == expected_penalties
        add_check(
            checks,
            scope="synthetic",
            profile="P30_team_baseline",
            check_name=f"known_penalties_{mode}",
            passed=ok,
            observed=observed,
            expected=expected_penalties,
        )
        subtraction_ok = True
        subtraction_detail = "final = max(0, P00 final - reported penalty)"
        for row in known:
            depmap_id = row["DepMap_ID"]
            p00_final = float(p00_map[depmap_id]["finalScore"])
            penalty = float(penalties[depmap_id]["exclusionPenalty"])
            observed_final = float(penalties[depmap_id]["finalScore"])
            expected_final = round(max(0.0, p00_final - penalty), 4)
            if not math.isclose(observed_final, expected_final, abs_tol=0.00011):
                subtraction_ok = False
                subtraction_detail = (
                    f"{depmap_id}: observed_final={observed_final}, "
                    f"P00_final={p00_final}, penalty={penalty}, "
                    f"expected_final={expected_final}"
                )
                break
        add_check(
            checks,
            scope="synthetic",
            profile="P30_team_baseline",
            check_name=f"final_score_subtraction_{mode}",
            passed=subtraction_ok,
            observed=subtraction_detail,
            expected="final = max(0, P00 final - reported penalty)",
        )

    all_missing = synthetic_rows([None, None, None, None])
    all_missing_results = {
        cap: score_with_cap(production_module, all_missing, "GENE_MULTIOMICS", cap)
        for cap in (0.0, 0.15, 0.30, 0.45)
    }
    missing_reference = all_missing_results[0.0]
    missing_ok = all(
        exact_scored_equal(
            missing_reference,
            result,
            (*P30_PARITY_FIELDS, "rank"),
        )[0]
        for result in all_missing_results.values()
    )
    add_check(
        checks,
        scope="synthetic",
        check_name="all_missing_profiles_identical",
        passed=missing_ok,
        observed="identical" if missing_ok else "different",
        expected="identical",
    )

    clipping_rows = synthetic_rows([10.0, 0.0])
    clipping_rows[0]["rnaExpr"] = 0.0
    clipping_rows[0]["protExpr"] = None
    clipping_rows[0]["nProt"] = 0
    clipping_rows[0]["hasProteomics"] = False
    clipping_rows[1]["rnaExpr"] = 10.0
    clipping_rows[1]["protExpr"] = None
    clipping_rows[1]["nProt"] = 0
    clipping_rows[1]["hasProteomics"] = False
    clip_p0 = rows_by_id(score_with_cap(production_module, clipping_rows, "GENE_MULTIOMICS", 0.0))
    clip_p45 = rows_by_id(score_with_cap(production_module, clipping_rows, "GENE_MULTIOMICS", 0.45))
    low_id = clipping_rows[0]["DepMap_ID"]
    realised_drop = clip_p0[low_id]["finalScore"] - clip_p45[low_id]["finalScore"]
    clipping_ok = (
        clip_p45[low_id]["finalScore"] == 0.0
        and clip_p45[low_id]["exclusionPenalty"] > realised_drop
    )
    add_check(
        checks,
        scope="synthetic",
        profile="P45_stress",
        check_name="lower_clipping_reported_penalty_exceeds_realised_drop",
        passed=clipping_ok,
        observed=f"final={clip_p45[low_id]['finalScore']}, penalty={clip_p45[low_id]['exclusionPenalty']}, realised_drop={realised_drop:.4f}",
        expected="final=0 and reported penalty > realised score drop",
    )

    same_base = synthetic_rows([0.0, 10.0])
    same_base_scored = score_with_cap(production_module, same_base, "GENE_MULTIOMICS", 0.30)
    same_map = rows_by_id(same_base_scored)
    low_ex_id = same_base[0]["DepMap_ID"]
    high_ex_id = same_base[1]["DepMap_ID"]
    order_ok = same_map[low_ex_id]["finalScore"] >= same_map[high_ex_id]["finalScore"]
    add_check(
        checks,
        scope="synthetic",
        profile="P30_team_baseline",
        check_name="same_base_lower_exclusion_not_ranked_below",
        passed=order_ok,
        observed=[same_map[low_ex_id]["finalScore"], same_map[high_ex_id]["finalScore"]],
        expected="first score >= second score",
    )

    for label, values, expected_penalty in (
        ("constant_positive", [2.0, 2.0], 0.30),
        ("constant_zero", [0.0, 0.0], 0.0),
        ("constant_negative", [-2.0, -2.0], 0.0),
    ):
        result = score_with_cap(
            production_module, synthetic_rows(values), "GENE_MULTIOMICS", 0.30
        )
        observed = [row["exclusionPenalty"] for row in result]
        ok = observed == [expected_penalty, expected_penalty]
        add_check(
            checks,
            scope="synthetic",
            profile="P30_team_baseline",
            check_name=f"legacy_{label}_behaviour",
            passed=ok,
            observed=observed,
            expected=[expected_penalty, expected_penalty],
            details="Legacy Min-Max discontinuity is exposed, not corrected by this runner.",
        )

    tie_ranked: list[dict[str, Any]] = []
    tie_scores = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.1, 0.0]
    for index, score in enumerate(tie_scores, start=1):
        tie_ranked.append(
            {
                "DepMap_ID": f"ACH-8{index:05d}",
                "finalScore": score,
                "confidenceScore": round(1.0 - index / 100.0, 4),
            }
        )
    add_tie_ranks(tie_ranked)
    top_with_ties = topk_with_ties(tie_ranked, 10)
    cutoff_rows = tie_ranked[9:11]
    tie_ok = (
        len(top_with_ties) == 11
        and all(row["scoreTieRankMin"] == 10 for row in cutoff_rows)
        and all(row["scoreTieRankMax"] == 11 for row in cutoff_rows)
        and all(row["scoreTieRankMid"] == 10.5 for row in cutoff_rows)
        and all(row["finalScoreTieSize"] == 2 for row in cutoff_rows)
        and all(row["fullComparatorTieSize"] == 1 for row in cutoff_rows)
    )
    add_check(
        checks,
        scope="synthetic",
        check_name="top10_cutoff_score_tie_included",
        passed=tie_ok,
        observed=(
            f"top_n={len(top_with_ties)}, score_tie_rank="
            f"{cutoff_rows[0]['scoreTieRankMin']}-"
            f"{cutoff_rows[0]['scoreTieRankMax']}, "
            f"final_tie={cutoff_rows[0]['finalScoreTieSize']}, "
            f"full_comparator_tie={cutoff_rows[0]['fullComparatorTieSize']}"
        ),
        expected=(
            "11 returned; tied cutoff rows have ranks 10-11 (mid 10.5), "
            "final-score tie size 2, full-comparator tie size 1"
        ),
    )

    deterministic_left = score_with_cap(
        production_module, known, "GENE_MULTIOMICS", 0.30
    )
    deterministic_right = score_with_cap(
        production_module, known, "GENE_MULTIOMICS", 0.30
    )
    deterministic_ok = deterministic_left == deterministic_right
    add_check(
        checks,
        scope="synthetic",
        profile="P30_team_baseline",
        check_name="repeat_run_determinism",
        passed=deterministic_ok,
        observed="identical" if deterministic_ok else "different",
        expected="identical",
    )

    default_result = production_module.score_candidates(known, query_mode="GENE_MULTIOMICS")
    wrapped_result = score_with_cap(production_module, known, "GENE_MULTIOMICS", 0.30)
    parity_ok, parity_detail = exact_scored_equal(
        default_result, wrapped_result, (*P30_PARITY_FIELDS, "rank")
    )
    add_check(
        checks,
        scope="synthetic",
        profile="P30_team_baseline",
        check_name="p30_default_wrapper_parity",
        passed=parity_ok,
        observed=parity_detail,
        expected="exact match",
    )

    negative_rejected = False
    try:
        score_with_cap(production_module, known, "GENE_MULTIOMICS", -0.1)
    except ValueError:
        negative_rejected = True
    add_check(
        checks,
        scope="synthetic",
        check_name="negative_cap_rejected",
        passed=negative_rejected,
        observed=negative_rejected,
        expected=True,
    )

    restored = math.isclose(
        float(production_module.MAX_EXCLUSION_PENALTY), original_cap, abs_tol=0.0
    )
    add_check(
        checks,
        scope="synthetic",
        check_name="production_cap_restored",
        passed=restored,
        observed=production_module.MAX_EXCLUSION_PENALTY,
        expected=original_cap,
    )

    failed = [check for check in checks if check["status"] != "PASS"]
    if failed:
        raise AssertionError(f"Synthetic self-test failed: {failed[0]}")
    if verbose:
        print("self-test OK")
        print(f"  {len(checks)} synthetic checks passed")
        print("  P30 wrapper matches the unmodified production scorer")
        print("  missing, clipping, constant-value, and tie-sensitive cases are covered")
        print("  GENE_MULTIOMICS, PROTEIN_ONLY, and COMBINED subtraction paths passed")
    return checks


def input_hashes(
    production_path: Path,
    design_path: Path,
    query_path: Path,
    v5_freeze_path: Path,
    cache_paths: Sequence[Path],
    raw_data_dir: Path,
) -> dict[str, dict[str, Any]]:
    paths = [
        Path(__file__).resolve(),
        production_path,
        design_path,
        query_path,
        v5_freeze_path,
        *cache_paths,
    ]
    module_dir = production_path.parent / "src"
    for name in ("data_loader.py", "data_merger.py"):
        candidate = module_dir / name
        if candidate.exists():
            paths.append(candidate)
    sample_info = raw_data_dir / "nomenclature" / "9_DepMap_sample_info.csv"
    if sample_info.exists():
        paths.append(sample_info)
    records: dict[str, dict[str, Any]] = {}
    for path in paths:
        resolved = path.resolve()
        records[str(resolved)] = {
            "sha256": sha256_file(resolved),
            "size_bytes": resolved.stat().st_size,
        }
    return records


def check_row_monotonicity(
    ranked_by_profile: Mapping[str, list[dict[str, Any]]],
    profiles: Sequence[Mapping[str, Any]],
) -> tuple[bool, bool, str]:
    maps = {name: rows_by_id(rows) for name, rows in ranked_by_profile.items()}
    caps_and_names = [
        (float(profile["max_exclusion_penalty"]), str(profile["name"]))
        for profile in profiles
    ]
    ids = sorted(next(iter(maps.values())))
    penalty_ok = True
    score_ok = True
    detail = "all candidates monotonic"
    for depmap_id in ids:
        penalties = [float(maps[name][depmap_id]["exclusionPenalty"]) for _, name in caps_and_names]
        scores = [float(maps[name][depmap_id]["finalScore"]) for _, name in caps_and_names]
        if any(right + 1e-12 < left for left, right in zip(penalties, penalties[1:])):
            penalty_ok = False
            detail = f"{depmap_id}: penalties={penalties}"
            break
        if any(right - 1e-12 > left for left, right in zip(scores, scores[1:])):
            score_ok = False
            detail = f"{depmap_id}: scores={scores}"
            break
    return penalty_ok, score_ok, detail


def build_outputs(
    production_module,
    recommender,
    queries: pd.DataFrame,
    profiles: list[dict[str, Any]],
    top_k: int,
    synthetic_checks: list[dict[str, Any]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    checks = list(synthetic_checks)
    query_profile_rows: list[dict[str, Any]] = []
    candidate_rows: list[dict[str, Any]] = []
    original_cap = float(production_module.MAX_EXCLUSION_PENALTY)

    for query in queries.to_dict(orient="records"):
        query_id = str(query["query_id"])
        target = str(query["target_gene"])
        exclusion = str(query["exclusion_gene"])
        disease = str(query["disease"])
        query_mode = str(query["query_mode"])
        aliases = production_module.build_disease_aliases(disease)
        rows = recommender.fetch_candidate_evidence(
            target_gene=target,
            target_protein=None,
            disease_aliases=aliases,
            exclusion_gene=exclusion,
        )
        ensure_unique_finite_rows(rows, query_id)
        n_exclusion = sum(int(row.get("nExclusion", 0) > 0) for row in rows)
        expected_candidates = int(query["expected_candidates"])
        expected_exclusion = int(query["expected_exclusion_observed"])
        add_check(
            checks,
            scope="query",
            query_id=query_id,
            check_name="frozen_candidate_count",
            passed=len(rows) == expected_candidates,
            observed=len(rows),
            expected=expected_candidates,
        )
        add_check(
            checks,
            scope="query",
            query_id=query_id,
            check_name="frozen_exclusion_observed_count",
            passed=n_exclusion == expected_exclusion,
            observed=n_exclusion,
            expected=expected_exclusion,
        )
        if len(rows) < top_k:
            raise ValueError(f"{query_id}: only {len(rows)} candidates for Top-{top_k}")

        exclusion_scaled_values = production_module.scale_available_values(
            rows, "exclusionExpr", "nExclusion"
        )
        exclusion_scaled = {
            str(row["DepMap_ID"]): value
            for row, value in zip(rows, exclusion_scaled_values)
        }
        observed_expr = np.asarray(
            [float(row["exclusionExpr"]) for row in rows if row.get("nExclusion", 0)],
            dtype=float,
        )
        observed_scaled = np.asarray(
            [float(value) for value in exclusion_scaled_values if value is not None],
            dtype=float,
        )
        nonconstant = len(observed_expr) >= 2 and float(observed_expr.max() - observed_expr.min()) > 0
        add_check(
            checks,
            scope="query",
            query_id=query_id,
            check_name="exclusion_expression_nonconstant",
            passed=nonconstant,
            observed=float(observed_expr.max() - observed_expr.min()) if len(observed_expr) else None,
            expected="> 0",
        )

        if not math.isclose(original_cap, 0.30, abs_tol=1e-12):
            raise AssertionError(
                f"{query_id}: unmodified production cap is {original_cap}, expected 0.30"
            )
        production_default = production_module.score_candidates(rows, query_mode=query_mode)
        add_tie_ranks(production_default)

        ranked_by_profile: dict[str, list[dict[str, Any]]] = {}
        for profile in profiles:
            name = str(profile["name"])
            cap = float(profile["max_exclusion_penalty"])
            ranked = score_with_cap(production_module, rows, query_mode, cap)
            add_tie_ranks(ranked)
            ranked_by_profile[name] = ranked

        p30 = ranked_by_profile["P30_team_baseline"]
        parity_ok, parity_detail = exact_scored_equal(
            production_default, p30, (*P30_PARITY_FIELDS, "rank")
        )
        add_check(
            checks,
            scope="query",
            query_id=query_id,
            profile="P30_team_baseline",
            check_name="p30_production_parity",
            passed=parity_ok,
            observed=parity_detail,
            expected="exact match",
        )
        repeated_p30 = score_with_cap(production_module, rows, query_mode, 0.30)
        deterministic_ok = repeated_p30 == [
            {key: value for key, value in row.items() if key not in {
                "scoreTieRankMin", "scoreTieRankMax", "scoreTieRankMid",
                "finalScoreTieSize", "fullComparatorTieSize"
            }}
            for row in p30
        ]
        add_check(
            checks,
            scope="query",
            query_id=query_id,
            profile="P30_team_baseline",
            check_name="repeat_run_determinism",
            passed=deterministic_ok,
            observed="identical" if deterministic_ok else "different",
            expected="identical",
        )

        baseline = ranked_by_profile["P00_no_penalty"]
        baseline_map = rows_by_id(baseline)
        candidate_sets_same = all(
            rows_by_id(ranked).keys() == baseline_map.keys()
            for ranked in ranked_by_profile.values()
        )
        add_check(
            checks,
            scope="query",
            query_id=query_id,
            check_name="candidate_sets_same_across_profiles",
            passed=candidate_sets_same,
            observed=candidate_sets_same,
            expected=True,
        )

        non_penalty_ok = True
        non_penalty_detail = "all non-penalty components match"
        for name, ranked in ranked_by_profile.items():
            ok, detail = exact_scored_equal(baseline, ranked, SCORE_COMPONENT_FIELDS)
            if not ok:
                non_penalty_ok = False
                non_penalty_detail = f"{name}: {detail}"
                break
        add_check(
            checks,
            scope="query",
            query_id=query_id,
            check_name="non_penalty_components_fixed",
            passed=non_penalty_ok,
            observed=non_penalty_detail,
            expected="exact match",
        )

        penalty_monotonic, score_monotonic, monotonic_detail = check_row_monotonicity(
            ranked_by_profile, profiles
        )
        add_check(
            checks,
            scope="query",
            query_id=query_id,
            check_name="penalty_non_decreasing_by_cap",
            passed=penalty_monotonic,
            observed=monotonic_detail,
            expected="non-decreasing for every candidate",
        )
        add_check(
            checks,
            scope="query",
            query_id=query_id,
            check_name="final_score_non_increasing_by_cap",
            passed=score_monotonic,
            observed=monotonic_detail,
            expected="non-increasing for every candidate",
        )

        profile_maps = {
            name: rows_by_id(ranked) for name, ranked in ranked_by_profile.items()
        }
        ratio_ok = True
        ratio_detail = "P30=2*P15 and P45=3*P15 within rounding tolerance"
        for depmap_id in sorted(baseline_map):
            p15_penalty = float(profile_maps["P15_weak"][depmap_id]["exclusionPenalty"])
            p30_penalty = float(
                profile_maps["P30_team_baseline"][depmap_id]["exclusionPenalty"]
            )
            p45_penalty = float(profile_maps["P45_stress"][depmap_id]["exclusionPenalty"])
            if not (
                math.isclose(p30_penalty, 2.0 * p15_penalty, abs_tol=0.00011)
                and math.isclose(p45_penalty, 3.0 * p15_penalty, abs_tol=0.00021)
            ):
                ratio_ok = False
                ratio_detail = (
                    f"{depmap_id}: P15={p15_penalty}, P30={p30_penalty}, "
                    f"P45={p45_penalty}"
                )
                break
        add_check(
            checks,
            scope="query",
            query_id=query_id,
            check_name="penalty_profile_ratios",
            passed=ratio_ok,
            observed=ratio_detail,
            expected="P30=2*P15 and P45=3*P15 within four-decimal rounding",
        )

        for profile in profiles:
            name = str(profile["name"])
            cap = float(profile["max_exclusion_penalty"])
            ranked = ranked_by_profile[name]
            current_map = rows_by_id(ranked)
            penalties = np.asarray([float(row["exclusionPenalty"]) for row in ranked])
            finals = np.asarray([float(row["finalScore"]) for row in ranked])
            cap_ok = bool((penalties >= -1e-12).all() and (penalties <= cap + 0.00011).all())
            add_check(
                checks,
                scope="query_profile",
                query_id=query_id,
                profile=name,
                check_name="penalty_within_profile_cap",
                passed=cap_ok,
                observed=float(penalties.max()) if len(penalties) else None,
                expected=f"0 <= penalty <= {cap}",
            )
            expected_penalty_ok = True
            expected_penalty_detail = "all rows equal round(cap * scaled exclusion, 4)"
            realised_drop_ok = True
            realised_drop_detail = (
                "unclipped rows lose the reported penalty; clipped rows lose no more"
            )
            for depmap_id, row in current_map.items():
                scaled = exclusion_scaled[depmap_id]
                expected_penalty = round(cap * float(scaled), 4) if scaled is not None else 0.0
                actual_penalty = float(row["exclusionPenalty"])
                if not math.isclose(actual_penalty, expected_penalty, abs_tol=0.0):
                    expected_penalty_ok = False
                    expected_penalty_detail = (
                        f"{depmap_id}: actual={actual_penalty}, expected={expected_penalty}"
                    )
                    break
                baseline_score = float(baseline_map[depmap_id]["finalScore"])
                current_score = float(row["finalScore"])
                realised_drop = baseline_score - current_score
                if realised_drop > actual_penalty + 0.00021:
                    realised_drop_ok = False
                    realised_drop_detail = (
                        f"{depmap_id}: realised_drop={realised_drop}, "
                        f"reported_penalty={actual_penalty}"
                    )
                    break
                if current_score > 0 and not math.isclose(
                    realised_drop, actual_penalty, abs_tol=0.00021
                ):
                    realised_drop_ok = False
                    realised_drop_detail = (
                        f"{depmap_id}: unclipped realised_drop={realised_drop}, "
                        f"reported_penalty={actual_penalty}"
                    )
                    break
            add_check(
                checks,
                scope="query_profile",
                query_id=query_id,
                profile=name,
                check_name="penalty_equals_cap_times_scaled_exclusion",
                passed=expected_penalty_ok,
                observed=expected_penalty_detail,
                expected="exact equality after four-decimal rounding",
            )
            add_check(
                checks,
                scope="query_profile",
                query_id=query_id,
                profile=name,
                check_name="realised_score_drop_matches_penalty_or_clipping",
                passed=realised_drop_ok,
                observed=realised_drop_detail,
                expected="drop equals penalty unless the lower bound clips at zero",
            )
            if name == "P00_no_penalty":
                add_check(
                    checks,
                    scope="query_profile",
                    query_id=query_id,
                    profile=name,
                    check_name="p00_penalties_zero",
                    passed=bool((penalties == 0).all()),
                    observed=float(penalties.max()) if len(penalties) else None,
                    expected=0.0,
                )
            missing_ids = [
                str(row["DepMap_ID"]) for row in rows if not row.get("nExclusion", 0)
            ]
            if missing_ids:
                missing_zero = all(
                    current_map[key]["exclusionPenalty"] == 0 for key in missing_ids
                )
                add_check(
                    checks,
                    scope="query_profile",
                    query_id=query_id,
                    profile=name,
                    check_name="missing_exclusion_penalty_zero",
                    passed=missing_zero,
                    observed=len(missing_ids),
                    expected="all missing rows have penalty 0",
                )
            score_bounds = bool((finals >= 0).all() and (finals <= 1).all())
            add_check(
                checks,
                scope="query_profile",
                query_id=query_id,
                profile=name,
                check_name="final_scores_bounded",
                passed=score_bounds,
                observed=f"[{finals.min():.4f}, {finals.max():.4f}]",
                expected="[0, 1]",
            )

            top = topk_with_ties(ranked, top_k)
            top_ids = {str(row["DepMap_ID"]) for row in top}
            comparison = comparison_metrics(ranked, baseline, top_k)
            tie_info = tie_metrics(ranked, top_k)
            top_exclusion = [
                exclusion_scaled[str(row["DepMap_ID"])]
                for row in top
                if exclusion_scaled[str(row["DepMap_ID"])] is not None
            ]
            query_profile_rows.append(
                {
                    "query_id": query_id,
                    "target_gene": target,
                    "exclusion_gene": exclusion,
                    "disease": disease,
                    "query_mode": query_mode,
                    "profile": name,
                    "max_exclusion_penalty": cap,
                    "n_candidates": len(ranked),
                    "n_exclusion_observed": n_exclusion,
                    "n_exclusion_missing": len(rows) - n_exclusion,
                    "exclusion_observed_fraction": n_exclusion / len(rows),
                    "exclusion_expr_min": float(observed_expr.min()) if len(observed_expr) else np.nan,
                    "exclusion_expr_max": float(observed_expr.max()) if len(observed_expr) else np.nan,
                    "exclusion_expr_range": float(observed_expr.max() - observed_expr.min()) if len(observed_expr) else np.nan,
                    "n_unique_exclusion_expr": int(pd.Series(observed_expr).nunique()),
                    "scaled_exclusion_mean": float(observed_scaled.mean()) if len(observed_scaled) else np.nan,
                    "scaled_exclusion_median": float(np.median(observed_scaled)) if len(observed_scaled) else np.nan,
                    "scaled_exclusion_iqr": float(np.quantile(observed_scaled, 0.75) - np.quantile(observed_scaled, 0.25)) if len(observed_scaled) else np.nan,
                    "mean_final_score": float(finals.mean()),
                    "final_score_std": float(finals.std(ddof=0)),
                    "min_final_score": float(finals.min()),
                    "max_final_score": float(finals.max()),
                    "zero_score_count": int((finals == 0).sum()),
                    "zero_score_fraction": float((finals == 0).mean()),
                    f"top{top_k}_mean_scaled_exclusion": float(np.mean(top_exclusion)) if top_exclusion else np.nan,
                    f"top{top_k}_high_exclusion_fraction_ge_0_8": float(np.mean(np.asarray(top_exclusion) >= 0.8)) if top_exclusion else np.nan,
                    **tie_info,
                    **comparison,
                }
            )

            raw_map = {str(row["DepMap_ID"]): row for row in rows}
            baseline_top_ids = {
                str(row["DepMap_ID"]) for row in topk_with_ties(baseline, top_k)
            }
            for row in ranked:
                depmap_id = str(row["DepMap_ID"])
                base = baseline_map[depmap_id]
                raw = raw_map[depmap_id]
                scaled = exclusion_scaled[depmap_id]
                candidate_rows.append(
                    {
                        "query_id": query_id,
                        "target_gene": target,
                        "exclusion_gene": exclusion,
                        "disease_context": disease,
                        "profile": name,
                        "max_exclusion_penalty": cap,
                        "DepMap_ID": depmap_id,
                        "cell_line": row.get("cellLine"),
                        "lineage": row.get("lineage"),
                        "primary_disease": row.get("disease"),
                        "exclusion_evidence_observed": bool(raw.get("nExclusion", 0)),
                        "exclusion_expression_std_mean": raw.get("exclusionExpr"),
                        "exclusion_scaled_within_query": scaled,
                        "rna_score": row.get("rnaScore"),
                        "protein_score": row.get("proteinScore"),
                        "biological_score": row.get("biologicalScore"),
                        "confidence_score": row.get("confidenceScore"),
                        "exclusion_penalty": row.get("exclusionPenalty"),
                        "final_score": row.get("finalScore"),
                        "production_ordinal_rank": row.get("rank"),
                        "score_tie_rank_min": row.get("scoreTieRankMin"),
                        "score_tie_rank_max": row.get("scoreTieRankMax"),
                        "score_tie_rank_mid": row.get("scoreTieRankMid"),
                        "final_score_tie_size": row.get("finalScoreTieSize"),
                        "full_comparator_tie_size": row.get("fullComparatorTieSize"),
                        "p00_final_score": base.get("finalScore"),
                        "p00_score_tie_rank_mid": base.get("scoreTieRankMid"),
                        "realised_score_drop_vs_P00": round(
                            float(base["finalScore"]) - float(row["finalScore"]), 4
                        ),
                        "rank_change_vs_P00_positive_is_down": float(row["scoreTieRankMid"])
                        - float(base["scoreTieRankMid"]),
                        f"in_top{top_k}_with_ties": depmap_id in top_ids,
                        f"in_P00_top{top_k}_with_ties": depmap_id in baseline_top_ids,
                        "score_clipped_to_zero": row.get("finalScore") == 0,
                    }
                )

    return (
        pd.DataFrame(query_profile_rows),
        pd.DataFrame(candidate_rows),
        pd.DataFrame(checks),
        pd.DataFrame(
            [
                {
                    "production_default_cap_after_run": production_module.MAX_EXCLUSION_PENALTY,
                    "expected_production_default_cap": original_cap,
                }
            ]
        ),
    )


def summarise_profiles(
    by_query: pd.DataFrame,
    checks: pd.DataFrame,
    top_k: int,
) -> pd.DataFrame:
    rows = []
    for profile, group in by_query.groupby("profile", sort=False):
        profile_checks = checks[(checks["profile"].isin(["", profile]))]
        rows.append(
            {
                "profile": profile,
                "max_exclusion_penalty": float(group["max_exclusion_penalty"].iloc[0]),
                "n_queries": group["query_id"].nunique(),
                "total_candidate_profile_rows": int(group["n_candidates"].sum()),
                "mean_exclusion_observed_fraction": float(group["exclusion_observed_fraction"].mean()),
                f"mean_top{top_k}_jaccard_vs_P00": float(group[f"top{top_k}_jaccard_vs_P00"].mean()),
                f"queries_with_top{top_k}_change_vs_P00": int((group[f"top{top_k}_jaccard_vs_P00"] < 1).sum()),
                "queries_with_top1_change_vs_P00": int((~group["top1_identical_to_P00"]).sum()),
                "mean_rank_spearman_vs_P00": float(group["rank_spearman_vs_P00"].mean()),
                "mean_of_query_mean_abs_rank_change_vs_P00": float(group["mean_abs_rank_change_vs_P00"].mean()),
                "max_abs_rank_change_vs_P00": float(group["max_abs_rank_change_vs_P00"].max()),
                "mean_zero_score_fraction": float(group["zero_score_fraction"].mean()),
                "total_zero_score_count": int(group["zero_score_count"].sum()),
                f"mean_top{top_k}_scaled_exclusion": float(group[f"top{top_k}_mean_scaled_exclusion"].mean()),
                f"mean_top{top_k}_high_exclusion_fraction_ge_0_8": float(group[f"top{top_k}_high_exclusion_fraction_ge_0_8"].mean()),
                f"mean_top{top_k}_with_ties_n": float(group[f"top{top_k}_with_ties_n"].mean()),
                "largest_final_score_tie_across_queries": int(group["largest_final_score_tie"].max()),
                "largest_full_comparator_tie_across_queries": int(
                    group["largest_full_comparator_tie"].max()
                ),
                "all_applicable_invariant_checks_pass": bool((profile_checks["status"] == "PASS").all()),
            }
        )
    return pd.DataFrame(rows)


def ensure_outputs(paths: Sequence[Path], force: bool) -> None:
    existing = [path for path in paths if path.exists()]
    if existing and not force:
        raise FileExistsError(
            "Refusing to overwrite existing Penalty results; use --force:\n  "
            + "\n  ".join(str(path) for path in existing)
        )


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=str) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()
    production_path = args.production_script.expanduser().resolve()
    if not production_path.exists():
        raise FileNotFoundError(f"Production script not found: {production_path}")
    production_module = load_python_module(production_path)

    if args.self_test:
        run_self_test(production_module)
        return 0

    if args.top_k < 1:
        raise ValueError("--top-k must be positive")
    design_path = args.design.expanduser().resolve()
    query_path = args.queries.expanduser().resolve()
    raw_data_dir = args.raw_data_dir.expanduser().resolve()
    cache_dir = args.cache_dir.expanduser().resolve()
    v5_freeze_path = args.v5_freeze.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    for path, label in (
        (design_path, "frozen design"),
        (query_path, "frozen query table"),
        (raw_data_dir, "raw data directory"),
        (cache_dir, "gene cache directory"),
        (v5_freeze_path, "V5 frozen design"),
    ):
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")

    design, profiles = load_design(design_path)
    queries = load_queries(query_path)
    expected_n_queries = int(design["query_selection_policy"]["n_queries"])
    if len(queries) != expected_n_queries:
        raise ValueError(
            f"Frozen design expects {expected_n_queries} queries; observed {len(queries)}"
        )
    cache_paths = required_cache_files(queries, cache_dir)

    output_paths = {
        "summary": out_dir / "penalty_ablation_summary.csv",
        "by_query_profile": out_dir / "penalty_ablation_by_query_profile.csv",
        "candidate_movements": out_dir / "penalty_ablation_candidate_movements.csv",
        "top10": out_dir / f"penalty_ablation_top{args.top_k}.csv",
        "invariant_checks": out_dir / "penalty_ablation_invariant_checks.csv",
        "manifest": out_dir / "penalty_ablation_manifest.json",
    }
    ensure_outputs(list(output_paths.values()), args.force)
    out_dir.mkdir(parents=True, exist_ok=True)

    before_hashes = input_hashes(
        production_path,
        design_path,
        query_path,
        v5_freeze_path,
        cache_paths,
        raw_data_dir,
    )
    synthetic_checks = run_self_test(production_module, verbose=True)

    recommender = production_module.DynamicMultiOmicsRecommender(
        raw_data_dir=raw_data_dir,
        cache_dir=cache_dir,
        refresh_cache=False,
    )
    by_query, candidates, checks, cap_audit = build_outputs(
        production_module,
        recommender,
        queries,
        profiles,
        args.top_k,
        synthetic_checks,
    )
    production_cap_restored = math.isclose(
        float(cap_audit["production_default_cap_after_run"].iloc[0]),
        float(cap_audit["expected_production_default_cap"].iloc[0]),
        abs_tol=0.0,
    )
    add_check(
        checks_records := checks.to_dict(orient="records"),
        scope="global",
        check_name="production_cap_restored_after_real_queries",
        passed=production_cap_restored,
        observed=cap_audit["production_default_cap_after_run"].iloc[0],
        expected=cap_audit["expected_production_default_cap"].iloc[0],
    )
    checks = pd.DataFrame(checks_records)

    after_hashes = input_hashes(
        production_path,
        design_path,
        query_path,
        v5_freeze_path,
        cache_paths,
        raw_data_dir,
    )
    hashes_unchanged = before_hashes == after_hashes
    checks_records = checks.to_dict(orient="records")
    add_check(
        checks_records,
        scope="global",
        check_name="all_input_and_cache_hashes_unchanged",
        passed=hashes_unchanged,
        observed=hashes_unchanged,
        expected=True,
        details="Includes production script, team loader/merger, sample info, frozen designs, query table, and all target/exclusion caches.",
    )
    checks = pd.DataFrame(checks_records)
    summary = summarise_profiles(by_query, checks, args.top_k)
    top10 = candidates[candidates[f"in_top{args.top_k}_with_ties"]].copy()

    failed = checks[checks["status"] != "PASS"]
    if not failed.empty:
        first = failed.iloc[0].to_dict()
        raise AssertionError(f"Penalty structural invariant failed: {first}")

    summary.to_csv(output_paths["summary"], index=False)
    by_query.to_csv(output_paths["by_query_profile"], index=False)
    candidates.to_csv(output_paths["candidate_movements"], index=False)
    top10.to_csv(output_paths["top10"], index=False)
    checks.to_csv(output_paths["invariant_checks"], index=False)

    output_hashes = {
        name: {
            "path": str(path),
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        }
        for name, path in output_paths.items()
        if name != "manifest"
    }
    manifest = {
        "schema_version": 1,
        "experiment_id": design["experiment_id"],
        "status": "complete_all_invariants_passed",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "claim_boundary": design["claim_boundary"],
        "command": [sys.executable, *sys.argv],
        "python_version": sys.version,
        "pandas_version": pd.__version__,
        "numpy_version": np.__version__,
        "top_k_with_ties": args.top_k,
        "n_queries": int(queries["query_id"].nunique()),
        "n_profiles": len(profiles),
        "n_query_profile_rows": len(by_query),
        "n_candidate_profile_rows": len(candidates),
        "n_topk_with_ties_rows": len(top10),
        "n_invariant_checks": len(checks),
        "all_invariant_checks_passed": True,
        "profiles": profiles,
        "queries": queries.to_dict(orient="records"),
        "paths": {
            "project_root": str(PROJECT_ROOT),
            "production_script": str(production_path),
            "raw_data_dir": str(raw_data_dir),
            "cache_dir": str(cache_dir),
            "frozen_design": str(design_path),
            "frozen_queries": str(query_path),
            "v5_frozen_design_hash_only_not_input": str(v5_freeze_path),
            "output_dir": str(out_dir),
        },
        "input_hashes_before": before_hashes,
        "input_hashes_after": after_hashes,
        "input_hashes_unchanged": hashes_unchanged,
        "output_hashes": output_hashes,
        "git": git_state(PROJECT_ROOT),
    }
    write_json(output_paths["manifest"], manifest)

    print("\nPenalty structural ablation complete")
    print(f"  queries: {len(queries)}")
    print(f"  profiles: {len(profiles)}")
    print(f"  candidate-profile rows: {len(candidates)}")
    print(f"  invariant checks: {len(checks)} PASS")
    print(f"  output: {out_dir}")
    print("\nProfile summary:")
    print(summary.to_string(index=False))
    print("\nConclusion boundary: structural/implementation evidence only; no accuracy or optimal-cap claim.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
