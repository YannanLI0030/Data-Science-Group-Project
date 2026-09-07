#!/usr/bin/env python3
"""Unlabelled dynamic ablation and missing-data robustness runner.

The runner consumes the versioned long-table panels under
``benchmarks/panels`` and deliberately does not consume a gold standard.  It
answers whether a controlled component removal changes scores/rankings, and
whether configurations collapse when modalities are missing.  Accuracy metrics
such as NDCG, Recall, and MRR belong to the later independently judged V5
holdout and are intentionally absent here.

Duplicate (gene, source, DepMap_ID) measurements are aggregated by arithmetic
mean for both ``value_raw`` and ``value_std``.  This matches the current dynamic
``_build_gene_table_from_long`` implementation.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import re
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PANELS_DIR = PROJECT_ROOT / "benchmarks" / "panels"
DEFAULT_OUT_DIR = PROJECT_ROOT / "results" / "dynamic_ablation_unlabelled"

PANEL_FILES = {
    "unfiltered_random": (
        "ablation_100genes_unfiltered_random_list.csv",
        "ablation_100genes_unfiltered_random.parquet",
    ),
    "coverage_stratified": (
        "ablation_100genes_coverage_stratified_list.csv",
        "ablation_100genes_coverage_stratified.parquet",
    ),
}

RNA_SOURCES = ("depmap_rna", "hpa_rna", "geo_rna")
PROTEIN_SOURCE = "ccle_gygi_protein"
R_LOW = 0.25
R_HIGH = 0.50
TRUST_FLOOR = 0.30
TRUST_DEFAULT = 0.60


@dataclass(frozen=True)
class ControlledConfig:
    name: str
    description: str
    rna_weight: float = 0.55
    protein_weight: float = 0.30
    confidence_weight: float = 0.15
    confidence_completeness: float = 0.40
    confidence_source_support: float = 0.35
    confidence_consistency: float = 0.25
    protein_for_biology: bool = True
    protein_for_confidence: bool = True
    adaptive_trust: bool = False
    structure: str = "standard"  # standard, rna_only, v3

    def validate(self) -> None:
        if self.structure not in {"standard", "rna_only", "v3"}:
            raise ValueError(f"{self.name}: unsupported structure {self.structure!r}")
        values = (
            self.rna_weight,
            self.protein_weight,
            self.confidence_weight,
            self.confidence_completeness,
            self.confidence_source_support,
            self.confidence_consistency,
        )
        if any(not math.isfinite(value) or value < 0 for value in values):
            raise ValueError(f"{self.name}: weights must be finite and non-negative")
        if self.structure == "standard" and not math.isclose(
            self.rna_weight + self.protein_weight + self.confidence_weight,
            1.0,
            abs_tol=1e-9,
        ):
            raise ValueError(f"{self.name}: outer weights must sum to 1")
        if not math.isclose(
            self.confidence_completeness
            + self.confidence_source_support
            + self.confidence_consistency,
            1.0,
            abs_tol=1e-9,
        ):
            raise ValueError(f"{self.name}: confidence component weights must sum to 1")
        if not self.protein_for_biology and self.protein_weight != 0:
            raise ValueError(
                f"{self.name}: protein weight must be zero when direct protein is disabled"
            )


def _renormalise(*weights: float) -> tuple[float, ...]:
    total = sum(weights)
    if total <= 0:
        raise ValueError("At least one retained weight is required")
    return tuple(weight / total for weight in weights)


_RNA_NO_CONF, _PROTEIN_NO_CONF = _renormalise(0.55, 0.30)
_NO_COMPLETENESS_SUPPORT, _NO_COMPLETENESS_CONSISTENCY = _renormalise(0.35, 0.25)
_NO_SUPPORT_COMPLETENESS, _NO_SUPPORT_CONSISTENCY = _renormalise(0.40, 0.25)
_NO_CONSISTENCY_COMPLETENESS, _NO_CONSISTENCY_SUPPORT = _renormalise(0.40, 0.35)

CONFIGS: dict[str, ControlledConfig] = {
    "B0_rna_mean": ControlledConfig(
        name="B0_rna_mean",
        description="RNA-only mean ranking; no protein or confidence in final score",
        rna_weight=1.0,
        protein_weight=0.0,
        confidence_weight=0.0,
        protein_for_biology=False,
        protein_for_confidence=False,
        structure="rna_only",
    ),
    "A0_team_baseline": ControlledConfig(
        name="A0_team_baseline",
        description="Team baseline: RNA .55, protein .30, confidence .15",
    ),
    "A1a_no_direct_protein": ControlledConfig(
        name="A1a_no_direct_protein",
        description="Remove direct protein; retain protein-derived confidence",
        rna_weight=0.85,
        protein_weight=0.0,
        protein_for_biology=False,
        protein_for_confidence=True,
    ),
    "A1b_no_protein_evidence": ControlledConfig(
        name="A1b_no_protein_evidence",
        description="Remove direct protein and all protein-derived confidence",
        rna_weight=0.85,
        protein_weight=0.0,
        protein_for_biology=False,
        protein_for_confidence=False,
    ),
    "A1c_no_protein_confidence": ControlledConfig(
        name="A1c_no_protein_confidence",
        description="Retain direct protein; remove protein-derived confidence only",
        protein_for_biology=True,
        protein_for_confidence=False,
    ),
    "A2_no_confidence": ControlledConfig(
        name="A2_no_confidence",
        description="Remove the full confidence block; renormalise RNA/protein",
        rna_weight=_RNA_NO_CONF,
        protein_weight=_PROTEIN_NO_CONF,
        confidence_weight=0.0,
    ),
    "A2a_no_conf_completeness": ControlledConfig(
        name="A2a_no_conf_completeness",
        description="Remove completeness from confidence; renormalise retained components",
        confidence_completeness=0.0,
        confidence_source_support=_NO_COMPLETENESS_SUPPORT,
        confidence_consistency=_NO_COMPLETENESS_CONSISTENCY,
    ),
    "A2b_no_conf_source_support": ControlledConfig(
        name="A2b_no_conf_source_support",
        description="Remove source support from confidence; renormalise retained components",
        confidence_completeness=_NO_SUPPORT_COMPLETENESS,
        confidence_source_support=0.0,
        confidence_consistency=_NO_SUPPORT_CONSISTENCY,
    ),
    "A2c_no_conf_consistency": ControlledConfig(
        name="A2c_no_conf_consistency",
        description="Remove RNA-protein consistency; renormalise retained components",
        confidence_completeness=_NO_CONSISTENCY_COMPLETENESS,
        confidence_source_support=_NO_CONSISTENCY_SUPPORT,
        confidence_consistency=0.0,
    ),
    "A3_equal_bio": ControlledConfig(
        name="A3_equal_bio",
        description="Equal RNA/protein biological weights",
        rna_weight=0.425,
        protein_weight=0.425,
    ),
    "A3b_reversed_bio": ControlledConfig(
        name="A3b_reversed_bio",
        description="Protein biological weight greater than RNA",
        rna_weight=0.30,
        protein_weight=0.55,
    ),
    "A4_adaptive_trust": ControlledConfig(
        name="A4_adaptive_trust",
        description="Team baseline plus gene-level adaptive RNA trust",
        adaptive_trust=True,
    ),
    "A5_v3_structure": ControlledConfig(
        name="A5_v3_structure",
        description="Flat v3 RNA/protein/consistency/completeness structure",
        structure="v3",
    ),
    "A6_v3_full": ControlledConfig(
        name="A6_v3_full",
        description="Flat v3 structure plus adaptive trust",
        adaptive_trust=True,
        structure="v3",
    ),
}

for _config in CONFIGS.values():
    _config.validate()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--panels-dir", type=Path, default=DEFAULT_PANELS_DIR
    )
    parser.add_argument(
        "--panel",
        choices=("both", "unfiltered_random", "coverage_stratified"),
        default="both",
    )
    parser.add_argument(
        "--sample-info",
        type=Path,
        help="Matching local 9_DepMap_sample_info.csv.",
    )
    parser.add_argument(
        "--corr",
        type=Path,
        help="RNA-protein correlation CSV with gene and correlation columns.",
    )
    parser.add_argument(
        "--production-script",
        type=Path,
        default=None,
        help="Latest dynamic scorer used to assert real-data A0 numeric parity.",
    )
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument(
        "--configs",
        nargs="+",
        choices=tuple(CONFIGS),
        default=list(CONFIGS),
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if not args.self_test and (args.sample_info is None or args.corr is None):
        parser.error("--sample-info and --corr are required (or use --self-test)")
    if args.top_k <= 0:
        parser.error("--top-k must be positive")
    return args


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def git_state(repo: Path) -> dict[str, object]:
    try:
        commit = subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
        ).strip()
        dirty = bool(
            subprocess.check_output(
                ["git", "-C", str(repo), "status", "--porcelain"], text=True
            ).strip()
        )
        return {"commit": commit, "dirty": dirty}
    except (OSError, subprocess.CalledProcessError):
        return {"commit": None, "dirty": None}


def safe_float(value: object, default: float = 0.0) -> float:
    try:
        if value is None or pd.isna(value):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def is_present(value: object) -> bool:
    return value is not None and not pd.isna(value)


def minmax_masked(values: Sequence[float | None]) -> list[float | None]:
    present = [float(value) for value in values if is_present(value)]
    if not present:
        return [None] * len(values)
    low, high = min(present), max(present)
    if math.isclose(low, high):
        constant = 1.0 if high > 0 else 0.0
        return [constant if is_present(value) else None for value in values]
    return [
        (float(value) - low) / (high - low) if is_present(value) else None
        for value in values
    ]


def rna_reliability(correlation: float | None) -> float:
    if correlation is None or pd.isna(correlation):
        return TRUST_DEFAULT
    value = float(correlation)
    if value >= R_HIGH:
        return 1.0
    if value < R_LOW:
        return TRUST_FLOOR
    return TRUST_FLOOR + (value - R_LOW) / (R_HIGH - R_LOW) * (1.0 - TRUST_FLOOR)


def confidence_components(
    row: Mapping[str, object],
    has_rna: bool,
    has_protein_input: bool,
    rna_score: float,
    protein_score: float,
    config: ControlledConfig,
) -> tuple[float, float, float | None]:
    if not config.protein_for_confidence:
        completeness = float(has_rna)
        source_support = (
            float(bool(row.get("hasDepMapRNA")))
            + float(bool(row.get("hasHpaRNA")))
            + float(bool(row.get("hasGeoRNA")))
        ) / 3.0
        return completeness, source_support, None

    completeness = (float(has_rna) + float(has_protein_input)) / 2.0
    source_support = (
        float(bool(row.get("hasDepMapRNA")))
        + float(bool(row.get("hasHpaRNA")))
        + float(bool(row.get("hasGeoRNA")))
        + float(bool(row.get("hasProteomics")))
    ) / 4.0
    if has_rna and has_protein_input:
        consistency = 1.0 - abs(rna_score - protein_score)
    elif has_rna or has_protein_input:
        consistency = 0.5
    else:
        consistency = 0.0
    return completeness, source_support, max(0.0, min(1.0, consistency))


def composite_confidence(
    config: ControlledConfig,
    completeness: float,
    source_support: float,
    consistency: float | None,
) -> float:
    parts = [
        (config.confidence_completeness, completeness),
        (config.confidence_source_support, source_support),
    ]
    if consistency is not None:
        parts.append((config.confidence_consistency, consistency))
    active = [(weight, value) for weight, value in parts if weight > 0]
    total = sum(weight for weight, _ in active)
    if total <= 0:
        return 0.0
    return sum(weight * value for weight, value in active) / total


def assign_tie_ranks(ranked: list[dict[str, object]]) -> None:
    start = 0
    while start < len(ranked):
        end = start + 1
        score = ranked[start]["finalScore"]
        while end < len(ranked) and ranked[end]["finalScore"] == score:
            end += 1
        midrank = ((start + 1) + end) / 2.0
        for position in range(start, end):
            ranked[position]["rank"] = position + 1
            ranked[position]["rankMid"] = midrank
        start = end


def score_candidates_controlled(
    rows: list[dict[str, object]],
    config: ControlledConfig,
    correlation: float | None,
) -> list[dict[str, object]]:
    if not rows:
        return []
    rna_scaled = minmax_masked(
        [safe_float(row.get("rnaExpr")) if safe_float(row.get("nRna")) > 0 else None for row in rows]
    )
    protein_scaled = minmax_masked(
        [safe_float(row.get("protExpr")) if safe_float(row.get("nProt")) > 0 else None for row in rows]
    )
    trust = rna_reliability(correlation) if config.adaptive_trust else 1.0
    scored: list[dict[str, object]] = []

    for index, row in enumerate(rows):
        has_rna = rna_scaled[index] is not None
        has_protein_input = protein_scaled[index] is not None
        has_protein_biology = has_protein_input and config.protein_for_biology
        rna_score = float(rna_scaled[index]) if has_rna else 0.0
        protein_score = float(protein_scaled[index]) if has_protein_input else 0.0
        completeness, support, consistency = confidence_components(
            row,
            has_rna,
            has_protein_input,
            rna_score,
            protein_score,
            config,
        )
        confidence = composite_confidence(
            config, completeness, support, consistency
        )

        if config.structure == "rna_only":
            if not has_rna:
                continue
            biological = rna_score
            final = rna_score
            effective_rna_weight = 1.0
            effective_protein_weight = 0.0
        elif config.structure == "v3":
            if not has_rna and not has_protein_biology:
                continue
            base_rna, base_protein, base_consistency, completeness_weight = (
                0.444,
                0.278,
                0.167,
                0.111,
            )
            effective_rna_weight = base_rna * trust
            freed = base_rna - effective_rna_weight
            effective_protein_weight = base_protein + (
                freed if has_protein_biology else 0.0
            )
            consistency_weight = base_consistency
            if not has_protein_biology:
                effective_rna_weight += freed + base_protein + base_consistency
                effective_protein_weight = 0.0
                consistency_weight = 0.0
            biological = (
                effective_rna_weight * rna_score
                + effective_protein_weight * protein_score
            )
            final = (
                biological
                + consistency_weight * float(consistency or 0.0)
                + completeness_weight * completeness
            )
        else:
            effective_rna_weight = config.rna_weight * trust
            freed = config.rna_weight - effective_rna_weight
            effective_protein_weight = config.protein_weight + (
                freed if has_protein_biology else 0.0
            )
            if not has_protein_biology:
                effective_rna_weight += freed
            available_weight = 0.0
            weighted_biology = 0.0
            if has_rna:
                available_weight += effective_rna_weight
                weighted_biology += effective_rna_weight * rna_score
            if has_protein_biology and effective_protein_weight > 0:
                available_weight += effective_protein_weight
                weighted_biology += effective_protein_weight * protein_score
            if available_weight <= 0:
                continue
            biological = weighted_biology / available_weight
            final = (
                (config.rna_weight + config.protein_weight) * biological
                + config.confidence_weight * confidence
            )

        final = round(max(0.0, min(1.0, float(final))), 4)
        confidence = round(max(0.0, min(1.0, float(confidence))), 4)
        enriched = dict(row)
        enriched.update(
            {
                "config": config.name,
                "rnaScore": round(rna_score, 4) if has_rna else None,
                "proteinScore": round(protein_score, 4) if has_protein_input else None,
                "biologicalScore": round(float(biological), 4),
                "completenessScore": round(completeness, 4),
                "sourceSupportScore": round(support, 4),
                "rnaProteinConsistencyScore": (
                    round(consistency, 4) if consistency is not None else None
                ),
                "confidenceScore": confidence,
                "finalScore": final,
                "rnaTrust": round(trust, 4),
                "effectiveRnaWeight": round(effective_rna_weight, 6),
                "effectiveProteinWeight": round(effective_protein_weight, 6),
            }
        )
        scored.append(enriched)

    scored.sort(
        key=lambda item: (
            -float(item["finalScore"]),
            -float(item["confidenceScore"]),
            str(item["DepMap_ID"]),
        )
    )
    assign_tie_ranks(scored)
    return scored


def detect_column(columns: Iterable[object], names: Iterable[str], label: str) -> str:
    norm = {
        re.sub(r"[^a-z0-9]", "", str(column).lower()): str(column)
        for column in columns
    }
    for name in names:
        key = re.sub(r"[^a-z0-9]", "", name.lower())
        if key in norm:
            return norm[key]
    raise KeyError(f"Could not detect {label}; columns={list(columns)}")


def load_correlations(path: Path) -> dict[str, float]:
    table = pd.read_csv(path)
    gene_col = detect_column(table.columns, ("gene", "gene_symbol", "symbol"), "gene")
    corr_col = detect_column(
        table.columns,
        ("correlation", "rna_protein_correlation", "pearson_r", "corr"),
        "correlation",
    )
    values = pd.to_numeric(table[corr_col], errors="coerce")
    return {
        str(gene).strip().upper(): float(value)
        for gene, value in zip(table[gene_col], values)
        if pd.notna(gene) and pd.notna(value)
    }


def prepare_sample_info(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, low_memory=False)
    if "DepMap_ID" not in table.columns:
        raise KeyError("sample info must contain DepMap_ID")
    table = table.dropna(subset=["DepMap_ID"]).copy()
    table["DepMap_ID"] = table["DepMap_ID"].astype(str).str.strip().str.upper()
    if table["DepMap_ID"].duplicated().any():
        table = table.drop_duplicates("DepMap_ID", keep="first")
    return table.sort_values("DepMap_ID").reset_index(drop=True)


def panel_specs(panels_dir: Path, panel: str) -> list[tuple[str, Path, Path]]:
    names = list(PANEL_FILES) if panel == "both" else [panel]
    specs = []
    for name in names:
        list_name, parquet_name = PANEL_FILES[name]
        list_path = panels_dir / list_name
        parquet_path = panels_dir / parquet_name
        if not list_path.exists() or not parquet_path.exists():
            raise FileNotFoundError(
                f"Panel {name!r} requires {list_path} and {parquet_path}"
            )
        specs.append((name, list_path, parquet_path))
    return specs


def load_panel(
    panel_name: str,
    list_path: Path,
    parquet_path: Path,
    valid_ids: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    panel_list = pd.read_csv(list_path)
    gene_col = detect_column(panel_list.columns, ("gene", "gene_symbol"), "panel gene")
    panel_list = panel_list.rename(columns={gene_col: "gene"}).copy()
    panel_list["gene"] = panel_list["gene"].astype(str).str.strip()
    panel_list["gene_key"] = panel_list["gene"].str.upper()
    if panel_list["gene_key"].duplicated().any():
        raise ValueError(f"{panel_name}: panel list contains duplicate genes")
    if "selection_stratum" not in panel_list.columns:
        panel_list["selection_stratum"] = panel_name

    long_table = pd.read_parquet(parquet_path)
    required = {
        "DepMap_ID",
        "gene_symbol",
        "omics_layer",
        "value_raw",
        "source",
        "value_std",
    }
    missing = sorted(required - set(long_table.columns))
    if missing:
        raise KeyError(f"{panel_name}: long table missing columns {missing}")
    long_table = long_table.copy()
    long_table["DepMap_ID"] = long_table["DepMap_ID"].astype(str).str.strip().str.upper()
    long_table["gene_key"] = long_table["gene_symbol"].astype(str).str.strip().str.upper()
    requested = set(panel_list["gene_key"])
    long_table = long_table[long_table["gene_key"].isin(requested)]

    ids_before = set(long_table["DepMap_ID"])
    rows_before = len(long_table)
    long_table = long_table[long_table["DepMap_ID"].isin(valid_ids)].copy()
    audit = {
        "rows_before_sample_filter": rows_before,
        "rows_after_sample_filter": len(long_table),
        "depmap_ids_before_sample_filter": len(ids_before),
        "unmatched_depmap_ids": len(ids_before - valid_ids),
        "rows_removed_by_sample_filter": rows_before - len(long_table),
    }
    observed = set(long_table["gene_key"])
    missing_genes = sorted(requested - observed)
    if missing_genes:
        raise ValueError(f"{panel_name}: genes absent after ID filtering: {missing_genes}")
    return panel_list, long_table, audit


def aggregate_long_table(long_table: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    keys = ["gene_key", "source", "DepMap_ID"]
    sizes = long_table.groupby(keys, dropna=False).size()
    values = long_table.groupby(keys, dropna=False)["value_raw"].nunique(dropna=True)
    audit = {
        "input_rows": len(long_table),
        "duplicate_groups": int((sizes > 1).sum()),
        "conflicting_value_groups": int((values > 1).sum()),
    }
    aggregated = (
        long_table.groupby(keys, as_index=False, dropna=False)[["value_raw", "value_std"]]
        .mean()
    )
    audit["aggregated_rows"] = len(aggregated)
    audit["rows_collapsed"] = len(long_table) - len(aggregated)
    return aggregated, audit


def preferred_name(row: pd.Series) -> str:
    for column in (
        "stripped_cell_line_name",
        "cell_line_name",
        "CCLE_Name",
        "ModelID",
    ):
        if column in row.index and pd.notna(row[column]) and str(row[column]).strip():
            return str(row[column]).strip()
    return str(row["DepMap_ID"])


def build_candidate_rows(
    gene: str,
    aggregated: pd.DataFrame,
    sample_info: pd.DataFrame,
) -> list[dict[str, object]]:
    gene_key = gene.upper()
    gene_long = aggregated[aggregated["gene_key"] == gene_key]
    base = sample_info.copy()
    for source in sorted(gene_long["source"].dropna().unique()):
        source_rows = gene_long[gene_long["source"] == source][
            ["DepMap_ID", "value_raw", "value_std"]
        ]
        base = base.merge(
            source_rows.rename(
                columns={
                    "value_raw": f"{source}__raw",
                    "value_std": f"{source}__std",
                }
            ),
            on="DepMap_ID",
            how="left",
            validate="one_to_one",
        )
    for source in (*RNA_SOURCES, PROTEIN_SOURCE):
        for suffix in ("raw", "std"):
            column = f"{source}__{suffix}"
            if column not in base.columns:
                base[column] = np.nan

    rows: list[dict[str, object]] = []
    for _, row in base.iterrows():
        rna_values = [row[f"{source}__std"] for source in RNA_SOURCES]
        available_rna = [float(value) for value in rna_values if pd.notna(value)]
        protein = row[f"{PROTEIN_SOURCE}__std"]
        if not available_rna and pd.isna(protein):
            continue
        rows.append(
            {
                "DepMap_ID": str(row["DepMap_ID"]),
                "cellLine": preferred_name(row),
                "lineage": row.get("lineage"),
                "disease": row.get("primary_disease"),
                "rnaExpr": float(np.mean(available_rna)) if available_rna else None,
                "protExpr": float(protein) if pd.notna(protein) else None,
                "exclusionExpr": None,
                "nRna": len(available_rna),
                "nProt": int(pd.notna(protein)),
                "nExclusion": 0,
                "hasDepMapRNA": pd.notna(row["depmap_rna__std"]),
                "hasHpaRNA": pd.notna(row["hpa_rna__std"]),
                "hasGeoRNA": pd.notna(row["geo_rna__std"]),
                "hasProteomics": pd.notna(protein),
            }
        )
    return rows


def topk_with_ties(ranked: list[dict[str, object]], k: int) -> list[dict[str, object]]:
    if not ranked:
        return []
    cutoff = float(ranked[min(k, len(ranked)) - 1]["finalScore"])
    return [row for row in ranked if float(row["finalScore"]) >= cutoff]


def jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else np.nan


def tie_diagnostics(ranked: list[dict[str, object]]) -> tuple[int, float]:
    if not ranked:
        return 0, np.nan
    counts = pd.Series([row["finalScore"] for row in ranked]).value_counts()
    tied = int(counts[counts > 1].sum())
    return int(counts.max()), tied / len(ranked)


def compare_to_baseline(
    ranked: list[dict[str, object]],
    baseline: list[dict[str, object]],
    top_k: int,
) -> dict[str, object]:
    current_order = [str(row["DepMap_ID"]) for row in ranked]
    baseline_order = [str(row["DepMap_ID"]) for row in baseline]
    current_top = {str(row["DepMap_ID"]) for row in topk_with_ties(ranked, top_k)}
    baseline_top = {str(row["DepMap_ID"]) for row in topk_with_ties(baseline, top_k)}
    current_top1 = {
        str(row["DepMap_ID"]) for row in ranked if row["finalScore"] == ranked[0]["finalScore"]
    }
    baseline_top1 = {
        str(row["DepMap_ID"])
        for row in baseline
        if row["finalScore"] == baseline[0]["finalScore"]
    }
    current_ranks = {str(row["DepMap_ID"]): float(row["rankMid"]) for row in ranked}
    baseline_ranks = {
        str(row["DepMap_ID"]): float(row["rankMid"]) for row in baseline
    }
    current_scores = {
        str(row["DepMap_ID"]): float(row["finalScore"]) for row in ranked
    }
    baseline_scores = {
        str(row["DepMap_ID"]): float(row["finalScore"]) for row in baseline
    }
    shared = sorted(set(current_ranks) & set(baseline_ranks))
    rank_delta = np.array(
        [abs(current_ranks[key] - baseline_ranks[key]) for key in shared], dtype=float
    )
    score_delta = np.array(
        [abs(current_scores[key] - baseline_scores[key]) for key in shared], dtype=float
    )
    # rankMid already contains tie-aware ranks, so Spearman correlation is the
    # ordinary Pearson correlation of these rank vectors.  Computing it here
    # avoids Pandas' optional SciPy dependency.
    if len(shared) >= 2:
        current_vector = np.asarray([current_ranks[key] for key in shared], dtype=float)
        baseline_vector = np.asarray([baseline_ranks[key] for key in shared], dtype=float)
        current_std = float(current_vector.std(ddof=0))
        baseline_std = float(baseline_vector.std(ddof=0))
        if current_std > 0 and baseline_std > 0:
            spearman = float(np.corrcoef(current_vector, baseline_vector)[0, 1])
        elif np.array_equal(current_vector, baseline_vector):
            spearman = 1.0
        else:
            spearman = np.nan
    else:
        spearman = np.nan
    return {
        "candidate_set_same_as_A0": set(current_order) == set(baseline_order),
        "ranking_identical_to_A0": current_order == baseline_order,
        f"top{top_k}_identical_to_A0": current_top == baseline_top,
        f"top{top_k}_jaccard_vs_A0": jaccard(current_top, baseline_top),
        "top1_identical_to_A0": current_top1 == baseline_top1,
        "rank_spearman_vs_A0": spearman,
        "mean_abs_rank_change_vs_A0": float(rank_delta.mean()) if len(rank_delta) else np.nan,
        "max_abs_rank_change_vs_A0": float(rank_delta.max()) if len(rank_delta) else np.nan,
        "mean_abs_score_change_vs_A0": float(score_delta.mean()) if len(score_delta) else np.nan,
        "n_shared_candidates_with_A0": len(shared),
    }


def load_python_module(path: Path):
    spec = importlib.util.spec_from_file_location("cellline_dynamic_production", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import production script: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def assert_production_a0_parity(
    production_module,
    rows: list[dict[str, object]],
    controlled_a0: list[dict[str, object]],
    gene: str,
) -> None:
    production = production_module.score_candidates(
        rows, query_mode="GENE_MULTIOMICS"
    )
    expected = {
        str(row["DepMap_ID"]): (
            float(row["finalScore"]),
            float(row["confidenceScore"]),
        )
        for row in production
    }
    observed = {
        str(row["DepMap_ID"]): (
            float(row["finalScore"]),
            float(row["confidenceScore"]),
        )
        for row in controlled_a0
    }
    if expected.keys() != observed.keys():
        raise AssertionError(f"{gene}: production A0 candidate-set parity failed")
    mismatches = [key for key in expected if expected[key] != observed[key]]
    if mismatches:
        first = mismatches[0]
        raise AssertionError(
            f"{gene}: production A0 numeric parity failed for {first}: "
            f"expected={expected[first]}, observed={observed[first]}"
        )


def summarise_detail(detail: pd.DataFrame, top_k: int) -> pd.DataFrame:
    if detail.empty:
        return pd.DataFrame()
    overall = detail.assign(panel_stratum="ALL")
    data = pd.concat([detail, overall], ignore_index=True)
    group_cols = ["panel", "panel_stratum", "config"]
    rows = []
    for keys, group in data.groupby(group_cols, sort=True):
        panel, stratum, config = keys
        rows.append(
            {
                "panel": panel,
                "panel_stratum": stratum,
                "config": config,
                "description": CONFIGS[config].description,
                "n_genes": group["gene"].nunique(),
                "genes_with_protein": int(group["gene_has_protein"].sum()),
                "mean_candidates": group["n_candidates"].mean(),
                "ranking_identical_rate_vs_A0": group[
                    "ranking_identical_to_A0"
                ].mean(),
                f"top{top_k}_identical_rate_vs_A0": group[
                    f"top{top_k}_identical_to_A0"
                ].mean(),
                f"mean_top{top_k}_jaccard_vs_A0": group[
                    f"top{top_k}_jaccard_vs_A0"
                ].mean(),
                "top1_change_rate_vs_A0": 1.0 - group["top1_identical_to_A0"].mean(),
                "mean_rank_spearman_vs_A0": group["rank_spearman_vs_A0"].mean(),
                "mean_abs_rank_change_vs_A0": group[
                    "mean_abs_rank_change_vs_A0"
                ].mean(),
                "mean_abs_score_change_vs_A0": group[
                    "mean_abs_score_change_vs_A0"
                ].mean(),
                "confidence_constant_rate": group["confidence_is_constant"].mean(),
                "mean_final_score_std": group["final_score_std"].mean(),
                "mean_confidence_std": group["confidence_std"].mean(),
                "mean_largest_tie_size": group["largest_final_score_tie"].mean(),
                "mean_tied_candidate_fraction": group[
                    "tied_candidate_fraction"
                ].mean(),
            }
        )
    return pd.DataFrame(rows)


def run_panel(
    panel_name: str,
    panel_list: pd.DataFrame,
    long_table: pd.DataFrame,
    sample_info: pd.DataFrame,
    correlations: Mapping[str, float],
    config_names: Sequence[str],
    top_k: int,
    production_module,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    aggregated, aggregation_audit = aggregate_long_table(long_table)
    detail_rows: list[dict[str, object]] = []
    top_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    parity_genes = 0

    for panel_row in panel_list.itertuples(index=False):
        gene = str(panel_row.gene)
        gene_key = gene.upper()
        stratum = str(panel_row.selection_stratum)
        correlation = correlations.get(gene_key, np.nan)
        candidates = build_candidate_rows(gene, aggregated, sample_info)
        if not candidates:
            raise ValueError(f"{panel_name}/{gene}: no candidates after aggregation")
        ranked_by_config = {
            name: score_candidates_controlled(candidates, CONFIGS[name], correlation)
            for name in config_names
        }
        if "A0_team_baseline" not in ranked_by_config:
            baseline = score_candidates_controlled(
                candidates, CONFIGS["A0_team_baseline"], correlation
            )
        else:
            baseline = ranked_by_config["A0_team_baseline"]

        if production_module is not None:
            assert_production_a0_parity(production_module, candidates, baseline, gene)
            parity_genes += 1

        full_fingerprints = set()
        top_fingerprints = set()
        changed_top = 0
        baseline_top_ids = tuple(
            sorted(str(row["DepMap_ID"]) for row in topk_with_ties(baseline, top_k))
        )
        for config_name, ranked in ranked_by_config.items():
            if not ranked:
                raise ValueError(f"{panel_name}/{gene}/{config_name}: no ranked candidates")
            comparison = compare_to_baseline(ranked, baseline, top_k)
            scores = pd.Series([float(row["finalScore"]) for row in ranked], dtype=float)
            confidences = pd.Series(
                [float(row["confidenceScore"]) for row in ranked], dtype=float
            )
            largest_tie, tied_fraction = tie_diagnostics(ranked)
            gene_has_protein = any(bool(row["hasProteomics"]) for row in candidates)
            source_counts = {
                source: sum(bool(row[field]) for row in candidates)
                for source, field in (
                    ("depmap_rna", "hasDepMapRNA"),
                    ("hpa_rna", "hasHpaRNA"),
                    ("geo_rna", "hasGeoRNA"),
                    ("protein", "hasProteomics"),
                )
            }
            detail_rows.append(
                {
                    "panel": panel_name,
                    "panel_stratum": stratum,
                    "gene": gene,
                    "config": config_name,
                    "description": CONFIGS[config_name].description,
                    "rna_protein_correlation": correlation,
                    "adaptive_rna_trust": rna_reliability(correlation),
                    "n_candidates": len(ranked),
                    "gene_has_protein": gene_has_protein,
                    "n_candidates_depmap_rna": source_counts["depmap_rna"],
                    "n_candidates_hpa_rna": source_counts["hpa_rna"],
                    "n_candidates_geo_rna": source_counts["geo_rna"],
                    "n_candidates_protein": source_counts["protein"],
                    "final_score_std": float(scores.std(ddof=0)),
                    "final_score_range": float(scores.max() - scores.min()),
                    "n_unique_final_scores": int(scores.nunique()),
                    "largest_final_score_tie": largest_tie,
                    "tied_candidate_fraction": tied_fraction,
                    "confidence_std": float(confidences.std(ddof=0)),
                    "n_unique_confidence_scores": int(confidences.nunique()),
                    "confidence_is_constant": confidences.nunique() <= 1,
                    **comparison,
                }
            )
            full_ids = tuple(str(row["DepMap_ID"]) for row in ranked)
            top_ids = tuple(
                sorted(str(row["DepMap_ID"]) for row in topk_with_ties(ranked, top_k))
            )
            full_fingerprints.add(full_ids)
            top_fingerprints.add(top_ids)
            changed_top += int(top_ids != baseline_top_ids)
            for row in topk_with_ties(ranked, top_k):
                top_rows.append(
                    {
                        "panel": panel_name,
                        "panel_stratum": stratum,
                        "gene": gene,
                        "config": config_name,
                        "rna_protein_correlation": correlation,
                        "DepMap_ID": row["DepMap_ID"],
                        "cellLine": row.get("cellLine"),
                        "rank": row["rank"],
                        "rankMid": row["rankMid"],
                        "finalScore": row["finalScore"],
                        "rnaScore": row["rnaScore"],
                        "proteinScore": row["proteinScore"],
                        "confidenceScore": row["confidenceScore"],
                        "completenessScore": row["completenessScore"],
                        "sourceSupportScore": row["sourceSupportScore"],
                        "rnaProteinConsistencyScore": row[
                            "rnaProteinConsistencyScore"
                        ],
                        "hasDepMapRNA": row["hasDepMapRNA"],
                        "hasHpaRNA": row["hasHpaRNA"],
                        "hasGeoRNA": row["hasGeoRNA"],
                        "hasProteomics": row["hasProteomics"],
                    }
                )

        gene_rows.append(
            {
                "panel": panel_name,
                "panel_stratum": stratum,
                "gene": gene,
                "rna_protein_correlation": correlation,
                "adaptive_rna_trust": rna_reliability(correlation),
                "n_configs": len(config_names),
                "unique_full_rankings": len(full_fingerprints),
                f"unique_top{top_k}_sets": len(top_fingerprints),
                "all_config_rankings_identical": len(full_fingerprints) == 1,
                f"configs_with_changed_top{top_k}_vs_A0": changed_top,
            }
        )

    audit = {
        "aggregation": aggregation_audit,
        "production_a0_parity_checked": production_module is not None,
        "production_a0_parity_genes": parity_genes,
    }
    return (
        pd.DataFrame(detail_rows),
        pd.DataFrame(gene_rows),
        pd.DataFrame(top_rows),
        audit,
    )


def ensure_outputs(paths: Sequence[Path], force: bool) -> None:
    existing = [path for path in paths if path.exists()]
    if existing and not force:
        raise FileExistsError(
            "Refusing to overwrite existing results; use --force:\n  "
            + "\n  ".join(str(path) for path in existing)
        )


def run_self_test() -> None:
    rows = [
        {
            "DepMap_ID": "ACH-000001",
            "cellLine": "RNA-high",
            "rnaExpr": 10.0,
            "protExpr": 0.0,
            "nRna": 3,
            "nProt": 1,
            "hasDepMapRNA": True,
            "hasHpaRNA": True,
            "hasGeoRNA": True,
            "hasProteomics": True,
        },
        {
            "DepMap_ID": "ACH-000002",
            "cellLine": "protein-high",
            "rnaExpr": 3.0,
            "protExpr": 10.0,
            "nRna": 3,
            "nProt": 1,
            "hasDepMapRNA": True,
            "hasHpaRNA": True,
            "hasGeoRNA": True,
            "hasProteomics": True,
        },
        {
            "DepMap_ID": "ACH-000003",
            "cellLine": "low",
            "rnaExpr": 0.0,
            "protExpr": 0.0,
            "nRna": 1,
            "nProt": 1,
            "hasDepMapRNA": True,
            "hasHpaRNA": False,
            "hasGeoRNA": False,
            "hasProteomics": True,
        },
    ]
    a0 = score_candidates_controlled(rows, CONFIGS["A0_team_baseline"], 0.80)
    a1a = score_candidates_controlled(rows, CONFIGS["A1a_no_direct_protein"], 0.80)
    a0_scores = {row["DepMap_ID"]: row["finalScore"] for row in a0}
    a1a_scores = {row["DepMap_ID"]: row["finalScore"] for row in a1a}
    assert a0_scores != a1a_scores

    without_protein = []
    for row in rows:
        changed = dict(row)
        changed.update({"protExpr": None, "nProt": 0, "hasProteomics": False})
        without_protein.append(changed)
    no_p_a0 = score_candidates_controlled(
        without_protein, CONFIGS["A0_team_baseline"], np.nan
    )
    no_p_a1a = score_candidates_controlled(
        without_protein, CONFIGS["A1a_no_direct_protein"], np.nan
    )
    assert [(r["DepMap_ID"], r["finalScore"]) for r in no_p_a0] == [
        (r["DepMap_ID"], r["finalScore"]) for r in no_p_a1a
    ]

    a1b_before = score_candidates_controlled(
        rows, CONFIGS["A1b_no_protein_evidence"], 0.10
    )
    changed_rows = [dict(row) for row in rows]
    changed_rows[0]["protExpr"] = 1000.0
    changed_rows[1]["protExpr"] = -1000.0
    a1b_after = score_candidates_controlled(
        changed_rows, CONFIGS["A1b_no_protein_evidence"], 0.10
    )
    assert [(r["DepMap_ID"], r["finalScore"]) for r in a1b_before] == [
        (r["DepMap_ID"], r["finalScore"]) for r in a1b_after
    ]

    a1c_before = score_candidates_controlled(
        rows, CONFIGS["A1c_no_protein_confidence"], 0.10
    )
    a1c_after = score_candidates_controlled(
        changed_rows, CONFIGS["A1c_no_protein_confidence"], 0.10
    )
    before_confidence = {
        r["DepMap_ID"]: r["confidenceScore"] for r in a1c_before
    }
    after_confidence = {r["DepMap_ID"]: r["confidenceScore"] for r in a1c_after}
    before_final = {r["DepMap_ID"]: r["finalScore"] for r in a1c_before}
    after_final = {r["DepMap_ID"]: r["finalScore"] for r in a1c_after}
    assert before_confidence == after_confidence
    assert before_final != after_final

    a4_high = score_candidates_controlled(rows, CONFIGS["A4_adaptive_trust"], 0.80)
    a4_low = score_candidates_controlled(rows, CONFIGS["A4_adaptive_trust"], 0.10)
    assert [(r["DepMap_ID"], r["finalScore"]) for r in a0] == [
        (r["DepMap_ID"], r["finalScore"]) for r in a4_high
    ]
    assert [r["DepMap_ID"] for r in a4_high] != [r["DepMap_ID"] for r in a4_low]

    sys.path.insert(0, str(PROJECT_ROOT / "src"))
    import merged_cellline_selector as production

    prod = production.score_candidates(rows, config=production.DEFAULT_SCORING_CONFIG)
    expected = {r["DepMap_ID"]: (r["finalScore"], r["confidenceScore"]) for r in prod}
    observed = {r["DepMap_ID"]: (r["finalScore"], r["confidenceScore"]) for r in a0}
    assert expected == observed

    duplicate_long = pd.DataFrame(
        {
            "gene_key": ["TEST", "TEST", "TEST"],
            "source": ["geo_rna", "geo_rna", "depmap_rna"],
            "DepMap_ID": ["ACH-000001", "ACH-000001", "ACH-000001"],
            "value_raw": [1.0, 3.0, 4.0],
            "value_std": [-1.0, 1.0, 0.5],
        }
    )
    aggregated, audit = aggregate_long_table(duplicate_long)
    geo = aggregated[aggregated["source"] == "geo_rna"].iloc[0]
    assert math.isclose(geo["value_raw"], 2.0)
    assert math.isclose(geo["value_std"], 0.0)
    assert audit["duplicate_groups"] == 1 and audit["rows_collapsed"] == 1

    print("self-test OK")
    print("  A0 numeric output matches production scorer")
    print("  no-protein A0 and A1a collapse as expected")
    print("  A1b ignores all protein values")
    print("  A1c retains direct protein but removes protein-derived confidence")
    print("  low-correlation adaptive trust changes ranking")
    print("  duplicate source measurements use arithmetic mean")


def main() -> int:
    args = parse_args()
    if args.self_test:
        run_self_test()
        return 0

    panels_dir = args.panels_dir.expanduser().resolve()
    sample_path = args.sample_info.expanduser().resolve()
    corr_path = args.corr.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    production_path = (
        args.production_script.expanduser().resolve()
        if args.production_script is not None
        else None
    )
    for path, label in ((sample_path, "sample info"), (corr_path, "correlation")):
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")
    if production_path is not None and not production_path.exists():
        raise FileNotFoundError(f"production script not found: {production_path}")

    output_paths = {
        "summary": out_dir / "dynamic_ablation_unlabelled_summary.csv",
        "by_gene_config": out_dir / "dynamic_ablation_unlabelled_by_gene_config.csv",
        "gene_diagnostics": out_dir / "dynamic_ablation_unlabelled_gene_diagnostics.csv",
        "topk": out_dir / f"dynamic_ablation_unlabelled_top{args.top_k}.csv",
        "manifest": out_dir / "dynamic_ablation_unlabelled_manifest.json",
    }
    ensure_outputs(list(output_paths.values()), args.force)
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_info = prepare_sample_info(sample_path)
    valid_ids = set(sample_info["DepMap_ID"])
    correlations = load_correlations(corr_path)
    production_module = load_python_module(production_path) if production_path else None

    all_detail = []
    all_gene = []
    all_top = []
    panel_audits: dict[str, object] = {}
    input_records = []
    for panel_name, list_path, parquet_path in panel_specs(panels_dir, args.panel):
        print(f"\nPanel: {panel_name}")
        panel_list, long_table, load_audit = load_panel(
            panel_name, list_path, parquet_path, valid_ids
        )
        detail, genes, top, run_audit = run_panel(
            panel_name=panel_name,
            panel_list=panel_list,
            long_table=long_table,
            sample_info=sample_info,
            correlations=correlations,
            config_names=args.configs,
            top_k=args.top_k,
            production_module=production_module,
        )
        all_detail.append(detail)
        all_gene.append(genes)
        all_top.append(top)
        panel_audits[panel_name] = {"load": load_audit, **run_audit}
        input_records.append(
            {
                "panel": panel_name,
                "list": {"path": str(list_path), "sha256": sha256_file(list_path)},
                "long_table": {
                    "path": str(parquet_path),
                    "sha256": sha256_file(parquet_path),
                },
            }
        )

    detail = pd.concat(all_detail, ignore_index=True)
    gene_diagnostics = pd.concat(all_gene, ignore_index=True)
    topk = pd.concat(all_top, ignore_index=True)
    summary = summarise_detail(detail, args.top_k)

    summary.to_csv(output_paths["summary"], index=False)
    detail.to_csv(output_paths["by_gene_config"], index=False)
    gene_diagnostics.to_csv(output_paths["gene_diagnostics"], index=False)
    topk.to_csv(output_paths["topk"], index=False)

    script_path = Path(__file__).resolve()
    manifest = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "experiment_role": "unlabelled component identifiability and missing-data robustness",
        "gold_standard_used": False,
        "accuracy_metrics_computed": [],
        "prohibited_interpretation": [
            "Do not interpret rank changes as improved recommendation accuracy.",
            "Do not report NDCG, Recall, MRR, or a final winning configuration from this run.",
        ],
        "candidate_scope": "all cell lines in the matching sample-info snapshot with RNA or protein evidence",
        "duplicate_aggregation": "arithmetic mean of value_raw and value_std by gene, source, DepMap_ID",
        "adaptive_trust": {
            "low_boundary": R_LOW,
            "high_boundary": R_HIGH,
            "floor": TRUST_FLOOR,
            "missing_default": TRUST_DEFAULT,
        },
        "top_k": args.top_k,
        "configs": [asdict(CONFIGS[name]) for name in args.configs],
        "generator": {
            "script": str(script_path),
            "script_sha256": sha256_file(script_path),
            "git": git_state(PROJECT_ROOT),
            "python": sys.version,
            "pandas": pd.__version__,
            "numpy": np.__version__,
        },
        "inputs": {
            "sample_info": {"path": str(sample_path), "sha256": sha256_file(sample_path)},
            "correlation": {"path": str(corr_path), "sha256": sha256_file(corr_path)},
            "production_script": (
                {"path": str(production_path), "sha256": sha256_file(production_path)}
                if production_path
                else None
            ),
            "panels": input_records,
        },
        "panel_audits": panel_audits,
        "outputs": {},
    }
    for key, path in output_paths.items():
        if key == "manifest":
            continue
        manifest["outputs"][key] = {
            "path": str(path),
            "rows": int(len(pd.read_csv(path))),
            "sha256": sha256_file(path),
        }
    output_paths["manifest"].write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    print("\nDynamic unlabelled ablation complete")
    print(f"  panel genes: {gene_diagnostics.groupby('panel')['gene'].nunique().to_dict()}")
    print(f"  configurations: {len(args.configs)}")
    print(f"  summary: {output_paths['summary']}")
    print(f"  per-gene/config: {output_paths['by_gene_config']}")
    print(f"  gene diagnostics: {output_paths['gene_diagnostics']}")
    print(f"  Top-{args.top_k} audit: {output_paths['topk']}")
    print(f"  manifest: {output_paths['manifest']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
