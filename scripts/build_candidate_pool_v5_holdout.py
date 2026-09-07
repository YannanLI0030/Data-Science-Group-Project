#!/usr/bin/env python3
"""Freeze the V5 ablation design and build a blinded holdout review pool.

The input is the unlabelled dynamic-ablation Top-10 audit.  The script freezes
configuration definitions before any V5 labels exist, selects a pre-declared
12-gene subset from the coverage-stratified panel, and writes two deliberately
separate files:

* a reviewer-facing CSV without configuration names, ranks, or scores; and
* an internal audit CSV containing the retrieval provenance.

The resulting pool is not a gold standard.  It becomes a holdout benchmark
only after every retrieved candidate has been independently reviewed and the
labels have been frozen.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DYNAMIC_DIR = PROJECT_ROOT / "results" / "dynamic_ablation_unlabelled"
DEFAULT_PANEL = (
    PROJECT_ROOT
    / "benchmarks"
    / "panels"
    / "ablation_100genes_coverage_stratified_list.csv"
)
DEFAULT_DESIGN = PROJECT_ROOT / "config" / "ablation_v5_frozen_design.json"
DEFAULT_GENE_MANIFEST = (
    PROJECT_ROOT
    / "benchmarks"
    / "panels"
    / "ablation_v5_holdout_gene_manifest.csv"
)
DEFAULT_REVIEW = PROJECT_ROOT / "benchmarks" / "candidate_pool_v5_holdout_review.csv"
DEFAULT_AUDIT = (
    PROJECT_ROOT
    / "results"
    / "v5_holdout_internal"
    / "candidate_pool_v5_holdout_audit.csv"
)
DEFAULT_POOL_MANIFEST = (
    PROJECT_ROOT / "benchmarks" / "candidate_pool_v5_holdout_manifest.json"
)

PRIMARY_POOL_CONFIGS = (
    "B0_rna_mean",
    "A0_team_baseline",
    "A1a_no_direct_protein",
    "A1b_no_protein_evidence",
    "A1c_no_protein_confidence",
    "A2_no_confidence",
    "A2a_no_conf_completeness",
    "A2b_no_conf_source_support",
    "A2c_no_conf_consistency",
    "A4_adaptive_trust",
)

SECONDARY_UNLABELLED_CONFIGS = (
    "A3_equal_bio",
    "A3b_reversed_bio",
    "A5_v3_structure",
    "A6_v3_full",
)

# Chosen before V5 candidate labels are inspected.  Selection intentionally
# favours protein-coding, human-interpretable genes whose frozen configurations
# have a chance to differ.  DepMap-RNA-only genes remain in the missing-data
# stress test, but are not suitable for a literature-intensive performance set.
SELECTED_GENES = (
    (
        "MFN2",
        "protein_corr_low",
        "Low RNA-protein correlation; broadly observed protein and interpretable mitochondrial function.",
    ),
    (
        "CALML3",
        "protein_corr_low",
        "Low RNA-protein correlation; interpretable epithelial calcium-binding marker.",
    ),
    (
        "OGDH",
        "protein_corr_mid",
        "Mid RNA-protein correlation; interpretable metabolic enzyme with broad measurements.",
    ),
    (
        "PSME4",
        "protein_corr_mid",
        "Mid RNA-protein correlation; interpretable proteostasis-related protein with broad measurements.",
    ),
    (
        "ERBB3",
        "protein_corr_high",
        "High RNA-protein correlation; interpretable receptor and cancer-model marker.",
    ),
    (
        "CHKA",
        "protein_corr_high",
        "High RNA-protein correlation; interpretable phospholipid-metabolism enzyme.",
    ),
    (
        "ARID3A",
        "protein_corr_high",
        "High RNA-protein correlation; interpretable lineage-associated transcription factor.",
    ),
    (
        "GALNT7",
        "protein_no_reliable_corr",
        "Broad protein coverage but no reliable correlation; informative correlation-missing case.",
    ),
    (
        "SIGLEC9",
        "protein_no_reliable_corr",
        "Protein present without reliable correlation; interpretable immune-surface marker.",
    ),
    (
        "CALCB",
        "protein_no_reliable_corr",
        "Protein present without reliable correlation; interpretable neuroendocrine-associated marker.",
    ),
    (
        "POU2F3",
        "multisource_rna_no_protein",
        "Three RNA sources and no proteomics; interpretable lineage transcription factor.",
    ),
    (
        "KRT76",
        "multisource_rna_no_protein",
        "Three RNA sources and no proteomics; interpretable epithelial keratin.",
    ),
)

MANUAL_FIELDS = (
    "judgement",
    "benchmark_task",
    "evidence_type",
    "source_url",
    "evidence_summary",
    "verified",
    "review_notes",
)

REVIEW_COLUMNS = (
    "review_id",
    "benchmark_role",
    "review_priority",
    "priority_reason",
    "gene",
    "DepMap_ID",
    "cell_line",
    "lineage",
    "disease",
    "subtype",
    "query_mode",
    "query_scope",
    *MANUAL_FIELDS,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dynamic-dir", type=Path, default=DEFAULT_DYNAMIC_DIR)
    parser.add_argument("--panel-list", type=Path, default=DEFAULT_PANEL)
    parser.add_argument("--sample-info", type=Path, required=False)
    parser.add_argument("--design-out", type=Path, default=DEFAULT_DESIGN)
    parser.add_argument("--gene-manifest-out", type=Path, default=DEFAULT_GENE_MANIFEST)
    parser.add_argument("--review-out", type=Path, default=DEFAULT_REVIEW)
    parser.add_argument("--audit-out", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--manifest-out", type=Path, default=DEFAULT_POOL_MANIFEST)
    parser.add_argument("--seed", type=int, default=20260906)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if not args.self_test and args.sample_info is None:
        parser.error("--sample-info is required (or use --self-test)")
    return args


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def ensure_files(paths: Iterable[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required file(s):\n  " + "\n  ".join(missing))


def load_dynamic_sources(dynamic_dir: Path) -> tuple[dict, Path, pd.DataFrame]:
    manifest_path = dynamic_dir / "dynamic_ablation_unlabelled_manifest.json"
    top10_path = dynamic_dir / "dynamic_ablation_unlabelled_top10.csv"
    ensure_files((manifest_path, top10_path))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("gold_standard_used") is not False:
        raise ValueError("Expected an unlabelled dynamic-ablation manifest")
    if manifest.get("top_k") != 10:
        raise ValueError("V5 design requires a frozen Top-10 source")
    recorded = manifest.get("outputs", {}).get("topk", {}).get("sha256")
    observed = sha256_file(top10_path)
    if recorded != observed:
        raise ValueError("Top-10 CSV hash does not match its dynamic manifest")
    parity = manifest.get("panel_audits", {}).get("coverage_stratified", {})
    if not parity.get("production_a0_parity_checked") or parity.get(
        "production_a0_parity_genes"
    ) != 100:
        raise ValueError("Coverage-stratified A0 production parity is incomplete")
    top10 = pd.read_csv(top10_path)
    return manifest, manifest_path, top10


def build_expected_design(
    dynamic_manifest: dict,
    dynamic_manifest_path: Path,
    top10_path: Path,
    panel_path: Path,
    sample_info_path: Path,
) -> dict:
    configs = {item["name"]: item for item in dynamic_manifest["configs"]}
    required = set(PRIMARY_POOL_CONFIGS) | set(SECONDARY_UNLABELLED_CONFIGS)
    missing = sorted(required - set(configs))
    if missing:
        raise ValueError(f"Dynamic manifest is missing frozen configs: {missing}")
    return {
        "schema_version": 1,
        "status": "frozen_before_v5_candidate_labelling",
        "experiment_role": "targeted independent-label component holdout",
        "not_a_weight_search": True,
        "query_mode": "GENE_MULTIOMICS",
        "query_scope": "all_cell_lines_no_disease_filter",
        "candidate_pool_rule": (
            "Union of Top-10-with-exact-score-ties from every primary pool "
            "configuration for each selected gene"
        ),
        "primary_pool_configs": list(PRIMARY_POOL_CONFIGS),
        "secondary_unlabelled_sensitivity_configs": list(
            SECONDARY_UNLABELLED_CONFIGS
        ),
        "secondary_config_policy": (
            "Frozen for structural sensitivity reporting only; excluded from "
            "the V5 candidate union and from V5 accuracy claims"
        ),
        "config_definitions": [configs[name] for name in PRIMARY_POOL_CONFIGS],
        "secondary_config_definitions": [
            configs[name] for name in SECONDARY_UNLABELLED_CONFIGS
        ],
        "gene_selection_policy": {
            "n_genes": len(SELECTED_GENES),
            "basis": (
                "Coverage-stratum balance, human interpretability, and unlabelled "
                "structural discriminability; no V5 candidate labels were used"
            ),
            "excluded_from_labelled_holdout": (
                "DepMap-RNA-only stratum retained for missing-data stress testing "
                "because literature verification is impractical and configurations "
                "mostly collapse"
            ),
            "selected": [
                {"order": order, "gene": gene, "stratum": stratum, "reason": reason}
                for order, (gene, stratum, reason) in enumerate(SELECTED_GENES, 1)
            ],
        },
        "source_snapshot": {
            "dynamic_manifest": {
                "path": str(dynamic_manifest_path.resolve()),
                "sha256": sha256_file(dynamic_manifest_path),
            },
            "dynamic_top10": {
                "path": str(top10_path.resolve()),
                "sha256": sha256_file(top10_path),
            },
            "coverage_stratified_panel": {
                "path": str(panel_path.resolve()),
                "sha256": sha256_file(panel_path),
            },
            "sample_info": {
                "path": str(sample_info_path.resolve()),
                "sha256": sha256_file(sample_info_path),
            },
        },
        "interpretation_limits": [
            "The V5 gene subset is targeted and coverage-stratified, not a random population sample.",
            "The pool is not a gold standard until all rows are externally reviewed and labels are frozen.",
            "Unknown candidates must not be treated as negatives.",
            "No configuration may be modified after V5 review begins.",
            "A post-label configuration change converts V5 from holdout to development evidence.",
        ],
    }


def freeze_or_verify_design(path: Path, expected: dict) -> dict:
    if path.exists():
        observed = json.loads(path.read_text(encoding="utf-8"))
        comparable = dict(observed)
        comparable.pop("frozen_utc", None)
        if comparable != expected:
            raise ValueError(
                f"Frozen design differs from current inputs or policy: {path}. "
                "Create a new version instead of overwriting V5."
            )
        return observed
    path.parent.mkdir(parents=True, exist_ok=True)
    frozen = {**expected, "frozen_utc": datetime.now(timezone.utc).isoformat()}
    path.write_text(json.dumps(frozen, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return frozen


def prepare_panel_manifest(panel_path: Path, top10: pd.DataFrame) -> pd.DataFrame:
    panel = pd.read_csv(panel_path)
    required = {
        "panel_order",
        "gene",
        "selection_stratum",
        "correlation",
        "protein_n_cell_lines",
        "rna_source_count",
    }
    missing = sorted(required - set(panel.columns))
    if missing:
        raise KeyError(f"Panel list missing columns: {missing}")
    panel = panel.copy()
    panel["gene"] = panel["gene"].astype(str).str.strip().str.upper()
    panel = panel.set_index("gene", drop=False)
    rows = []
    filtered = top10[
        (top10["panel"] == "coverage_stratified")
        & top10["config"].isin(PRIMARY_POOL_CONFIGS)
    ]
    for order, (gene, expected_stratum, reason) in enumerate(SELECTED_GENES, 1):
        if gene not in panel.index:
            raise ValueError(f"Selected V5 gene absent from panel: {gene}")
        source = panel.loc[gene]
        if str(source["selection_stratum"]) != expected_stratum:
            raise ValueError(
                f"{gene}: expected {expected_stratum}, observed {source['selection_stratum']}"
            )
        gene_top = filtered[filtered["gene"].astype(str).str.upper() == gene]
        observed_configs = set(gene_top["config"])
        if observed_configs != set(PRIMARY_POOL_CONFIGS):
            raise ValueError(f"{gene}: incomplete primary configuration retrieval")
        rows.append(
            {
                "holdout_order": order,
                "panel_order": int(source["panel_order"]),
                "gene": gene,
                "selection_stratum": expected_stratum,
                "rna_protein_correlation": source["correlation"],
                "rna_source_count": int(source["rna_source_count"]),
                "protein_n_cell_lines": int(source["protein_n_cell_lines"]),
                "top10_union_candidates": int(gene_top["DepMap_ID"].nunique()),
                "v5_role": "targeted_holdout_candidate",
                "selection_reason": reason,
            }
        )
    return pd.DataFrame(rows)


def prepare_sample_info(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, low_memory=False)
    needed = {
        "DepMap_ID",
        "cell_line_name",
        "lineage",
        "primary_disease",
        "Subtype",
    }
    missing = sorted(needed - set(table.columns))
    if missing:
        raise KeyError(f"Sample info missing columns: {missing}")
    table = table.dropna(subset=["DepMap_ID"]).copy()
    table["DepMap_ID"] = table["DepMap_ID"].astype(str).str.strip().str.upper()
    if table["DepMap_ID"].duplicated().any():
        raise ValueError("Sample info contains duplicate DepMap_ID values")
    optional = [
        name
        for name in ("parent_depmap_id", "Cellosaurus_issues", "RRID")
        if name in table.columns
    ]
    return table[
        ["DepMap_ID", "cell_line_name", "lineage", "primary_disease", "Subtype", *optional]
    ]


def stable_review_id(gene: str, depmap_id: str) -> str:
    token = sha256_text(f"v5|{gene}|{depmap_id}")[:12].upper()
    return f"V5-{token}"


def blinded_order_key(seed: int, gene: str, depmap_id: str) -> str:
    return sha256_text(f"{seed}|{gene}|{depmap_id}")


def build_pools(
    top10: pd.DataFrame,
    sample_info: pd.DataFrame,
    gene_manifest: pd.DataFrame,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {
        "panel",
        "panel_stratum",
        "gene",
        "config",
        "DepMap_ID",
        "cellLine",
        "rank",
        "rankMid",
        "finalScore",
        "rnaScore",
        "proteinScore",
        "confidenceScore",
        "completenessScore",
        "sourceSupportScore",
        "rnaProteinConsistencyScore",
        "hasDepMapRNA",
        "hasHpaRNA",
        "hasGeoRNA",
        "hasProteomics",
    }
    missing = sorted(required - set(top10.columns))
    if missing:
        raise KeyError(f"Top-10 audit missing columns: {missing}")
    selected = set(gene_manifest["gene"])
    work = top10[
        (top10["panel"] == "coverage_stratified")
        & top10["gene"].astype(str).str.upper().isin(selected)
        & top10["config"].isin(PRIMARY_POOL_CONFIGS)
    ].copy()
    work["gene"] = work["gene"].astype(str).str.strip().str.upper()
    work["DepMap_ID"] = work["DepMap_ID"].astype(str).str.strip().str.upper()

    rows = []
    order_map = {name: index for index, name in enumerate(PRIMARY_POOL_CONFIGS)}
    for (gene, depmap_id), group in work.groupby(["gene", "DepMap_ID"], sort=False):
        configs = sorted(set(group["config"]), key=order_map.__getitem__)
        first = group.iloc[0]
        row = {
            "review_id": stable_review_id(gene, depmap_id),
            "gene": gene,
            "DepMap_ID": depmap_id,
            "scorer_cell_line": first["cellLine"],
            "selection_stratum": first["panel_stratum"],
            "configs_retrieved": ";".join(configs),
            "n_configs_retrieved": len(configs),
            "best_rank": int(group["rank"].min()),
            "best_rank_mid": float(group["rankMid"].min()),
            "score_min_across_retrieving_configs": float(group["finalScore"].min()),
            "score_max_across_retrieving_configs": float(group["finalScore"].max()),
            "rnaScore": first["rnaScore"],
            "proteinScore": first["proteinScore"],
            "hasDepMapRNA": bool(first["hasDepMapRNA"]),
            "hasHpaRNA": bool(first["hasHpaRNA"]),
            "hasGeoRNA": bool(first["hasGeoRNA"]),
            "hasProteomics": bool(first["hasProteomics"]),
        }
        for config in PRIMARY_POOL_CONFIGS:
            match = group[group["config"] == config]
            row[f"rank_{config}"] = (
                int(match.iloc[0]["rank"]) if not match.empty else pd.NA
            )
            row[f"score_{config}"] = (
                float(match.iloc[0]["finalScore"]) if not match.empty else pd.NA
            )
            row[f"confidence_{config}"] = (
                float(match.iloc[0]["confidenceScore"]) if not match.empty else pd.NA
            )
        rows.append(row)

    audit = pd.DataFrame(rows).merge(
        sample_info,
        on="DepMap_ID",
        how="left",
        validate="many_to_one",
        indicator="_sample_info_match",
    )
    if audit["_sample_info_match"].ne("both").any():
        missing_ids = sorted(
            audit.loc[audit["_sample_info_match"].ne("both"), "DepMap_ID"].unique()
        )
        raise ValueError(f"Candidate IDs absent from sample info: {missing_ids[:10]}")
    audit["cell_line_name_source"] = "sample_info"
    missing_name = audit["cell_line_name"].isna() | audit[
        "cell_line_name"
    ].astype(str).str.strip().eq("")
    audit.loc[missing_name, "cell_line_name"] = audit.loc[
        missing_name, "scorer_cell_line"
    ]
    audit.loc[missing_name, "cell_line_name_source"] = "dynamic_scorer_fallback"
    if audit["cell_line_name"].isna().any():
        raise ValueError("Candidate has no usable cell-line name in either source")
    audit = audit.drop(columns="_sample_info_match")
    reasons = gene_manifest.set_index("gene")["selection_reason"]
    gene_orders = gene_manifest.set_index("gene")["holdout_order"]
    audit["selection_reason"] = audit["gene"].map(reasons)
    audit["holdout_order"] = audit["gene"].map(gene_orders)
    audit["blinded_order_key"] = [
        blinded_order_key(seed, gene, depmap_id)
        for gene, depmap_id in zip(audit["gene"], audit["DepMap_ID"])
    ]
    audit = audit.sort_values(
        ["holdout_order", "blinded_order_key", "DepMap_ID"]
    ).reset_index(drop=True)

    review = pd.DataFrame(
        {
            "review_id": audit["review_id"],
            "benchmark_role": "independent_holdout",
            "review_priority": 1,
            "priority_reason": "Required member of frozen V5 Top-10 union",
            "gene": audit["gene"],
            "DepMap_ID": audit["DepMap_ID"],
            "cell_line": audit["cell_line_name"],
            "lineage": audit["lineage"],
            "disease": audit["primary_disease"],
            "subtype": audit["Subtype"],
            "query_mode": "GENE_MULTIOMICS",
            "query_scope": "all_cell_lines_no_disease_filter",
            "judgement": "unknown",
            "benchmark_task": "expression_suitability",
            "evidence_type": "",
            "source_url": "",
            "evidence_summary": "",
            "verified": "no",
            "review_notes": "",
        },
        columns=REVIEW_COLUMNS,
    )
    validate_review(review, gene_manifest)
    return review, audit


def validate_review(review: pd.DataFrame, gene_manifest: pd.DataFrame) -> None:
    if list(review.columns) != list(REVIEW_COLUMNS):
        raise ValueError("Reviewer-facing schema changed unexpectedly")
    leaked = [
        column
        for column in review.columns
        if "config" in column.lower()
        or "score" in column.lower()
        or column.lower().startswith("rank")
    ]
    if leaked:
        raise ValueError(f"Reviewer file leaks model provenance: {leaked}")
    if review.duplicated(["gene", "DepMap_ID"]).any():
        raise ValueError("Reviewer pool contains duplicate gene-DepMap pairs")
    if review["review_id"].duplicated().any():
        raise ValueError("Reviewer IDs are not unique")
    if not review["DepMap_ID"].str.fullmatch(r"ACH-\d{6}").all():
        raise ValueError("Invalid DepMap ID format in reviewer pool")
    if set(review["gene"]) != set(gene_manifest["gene"]):
        raise ValueError("Reviewer genes differ from frozen gene manifest")
    if not review["judgement"].eq("unknown").all():
        raise ValueError("New reviewer rows must start with judgement=unknown")
    if not review["verified"].eq("no").all():
        raise ValueError("New reviewer rows must start with verified=no")


def protect_reviewed_file(path: Path) -> None:
    if not path.exists():
        return
    existing = pd.read_csv(path, encoding="utf-8-sig", keep_default_na=False)
    reviewed = (
        existing.get("judgement", pd.Series(dtype=str)).astype(str).str.lower().ne("unknown")
        | existing.get("verified", pd.Series(dtype=str)).astype(str).str.lower().ne("no")
        | existing.get("source_url", pd.Series(dtype=str)).astype(str).str.strip().ne("")
        | existing.get("evidence_summary", pd.Series(dtype=str)).astype(str).str.strip().ne("")
        | existing.get("review_notes", pd.Series(dtype=str)).astype(str).str.strip().ne("")
    )
    if reviewed.any():
        raise ValueError(
            f"Refusing to overwrite reviewed V5 labels: {path}. Create a new version instead."
        )


def write_outputs(
    args: argparse.Namespace,
    design: dict,
    gene_manifest: pd.DataFrame,
    review: pd.DataFrame,
    audit: pd.DataFrame,
) -> dict:
    output_paths = (
        args.gene_manifest_out,
        args.review_out,
        args.audit_out,
        args.manifest_out,
    )
    existing = [path for path in output_paths if path.exists()]
    if existing and not args.force:
        raise FileExistsError(
            "Refusing to overwrite V5 outputs; use --force only before review:\n  "
            + "\n  ".join(str(path) for path in existing)
        )
    if args.force:
        protect_reviewed_file(args.review_out)
    for path in output_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    gene_manifest.to_csv(args.gene_manifest_out, index=False, encoding="utf-8-sig")
    review.to_csv(args.review_out, index=False, encoding="utf-8-sig")
    audit.to_csv(args.audit_out, index=False, encoding="utf-8-sig")
    per_gene = (
        review.groupby("gene", sort=False)
        .size()
        .rename("candidate_rows")
        .reset_index()
        .to_dict("records")
    )
    manifest = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "design_status": design["status"],
        "design_path": str(args.design_out.resolve()),
        "design_sha256": sha256_file(args.design_out),
        "review_blinding": (
            "Configuration names, retrieval counts, ranks, scores, and coverage strata "
            "are excluded from the reviewer-facing CSV"
        ),
        "review_instructions": {
            "allowed_judgements": ["positive", "negative", "unknown"],
            "unknown_is_not_negative": True,
            "required_manual_fields": list(MANUAL_FIELDS),
            "do_not_open_before_label_freeze": str(args.audit_out.resolve()),
        },
        "counts": {
            "genes": int(review["gene"].nunique()),
            "candidate_rows": int(len(review)),
            "unique_depmap_ids": int(review["DepMap_ID"].nunique()),
            "primary_configs": len(PRIMARY_POOL_CONFIGS),
            "rows_retrieved_by_all_primary_configs": int(
                audit["n_configs_retrieved"].eq(len(PRIMARY_POOL_CONFIGS)).sum()
            ),
            "rows_from_top10_disagreement": int(
                audit["n_configs_retrieved"].lt(len(PRIMARY_POOL_CONFIGS)).sum()
            ),
            "cell_line_name_sources": {
                str(name): int(count)
                for name, count in audit["cell_line_name_source"]
                .value_counts()
                .items()
            },
            "per_gene": per_gene,
        },
        "generator": {
            "script": str(Path(__file__).resolve()),
            "script_sha256": sha256_file(Path(__file__).resolve()),
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "seed": args.seed,
        },
        "outputs": {},
    }
    for name, path, rows in (
        ("gene_manifest", args.gene_manifest_out, len(gene_manifest)),
        ("review_pool", args.review_out, len(review)),
        ("internal_audit", args.audit_out, len(audit)),
    ):
        manifest["outputs"][name] = {
            "path": str(path.resolve()),
            "rows": int(rows),
            "sha256": sha256_file(path),
        }
    args.manifest_out.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def verify_existing(args: argparse.Namespace, design: dict) -> None:
    ensure_files(
        (
            args.design_out,
            args.gene_manifest_out,
            args.review_out,
            args.audit_out,
            args.manifest_out,
        )
    )
    manifest = json.loads(args.manifest_out.read_text(encoding="utf-8"))
    if manifest.get("design_sha256") != sha256_file(args.design_out):
        raise ValueError("Frozen design hash mismatch")
    for output in manifest.get("outputs", {}).values():
        path = Path(output["path"])
        if not path.exists() or sha256_file(path) != output["sha256"]:
            raise ValueError(f"V5 output hash mismatch: {path}")
    review = pd.read_csv(args.review_out, encoding="utf-8-sig", keep_default_na=False)
    gene_manifest = pd.read_csv(args.gene_manifest_out, encoding="utf-8-sig")
    validate_review(review, gene_manifest)
    if design.get("status") != "frozen_before_v5_candidate_labelling":
        raise ValueError("Unexpected V5 freeze status")
    print("V5 freeze and candidate-pool verification OK")
    print(f"  genes: {review['gene'].nunique()}")
    print(f"  candidate rows: {len(review)}")


def run_self_test() -> None:
    assert set(PRIMARY_POOL_CONFIGS).isdisjoint(SECONDARY_UNLABELLED_CONFIGS)
    assert len(PRIMARY_POOL_CONFIGS) == len(set(PRIMARY_POOL_CONFIGS)) == 10
    assert len(SELECTED_GENES) == len({item[0] for item in SELECTED_GENES}) == 12
    assert stable_review_id("ERBB3", "ACH-000001") == stable_review_id(
        "ERBB3", "ACH-000001"
    )
    assert blinded_order_key(42, "ERBB3", "ACH-000001") != blinded_order_key(
        43, "ERBB3", "ACH-000001"
    )
    assert not any(
        "config" in name.lower() or "score" in name.lower() or name.startswith("rank")
        for name in REVIEW_COLUMNS
    )
    print("self-test OK")
    print("  10 primary pool configs and 4 secondary unlabelled configs are disjoint")
    print("  12 selected genes and blinded reviewer schema are valid")


def main() -> int:
    args = parse_args()
    if args.self_test:
        run_self_test()
        return 0
    dynamic_manifest, dynamic_manifest_path, top10 = load_dynamic_sources(
        args.dynamic_dir
    )
    top10_path = args.dynamic_dir / "dynamic_ablation_unlabelled_top10.csv"
    ensure_files((args.panel_list, args.sample_info))
    expected_design = build_expected_design(
        dynamic_manifest,
        dynamic_manifest_path,
        top10_path,
        args.panel_list,
        args.sample_info,
    )
    design = freeze_or_verify_design(args.design_out, expected_design)
    if args.verify_only:
        verify_existing(args, design)
        return 0
    gene_manifest = prepare_panel_manifest(args.panel_list, top10)
    sample_info = prepare_sample_info(args.sample_info)
    review, audit = build_pools(top10, sample_info, gene_manifest, args.seed)
    manifest = write_outputs(args, design, gene_manifest, review, audit)
    print("V5 frozen design and blinded review pool created")
    print(f"  primary configs: {len(PRIMARY_POOL_CONFIGS)}")
    print(f"  secondary unlabelled configs: {len(SECONDARY_UNLABELLED_CONFIGS)}")
    print(f"  selected genes: {manifest['counts']['genes']}")
    print(f"  review rows: {manifest['counts']['candidate_rows']}")
    print(f"  frozen design: {args.design_out}")
    print(f"  reviewer CSV: {args.review_out}")
    print(f"  internal audit (keep blinded): {args.audit_out}")
    print(f"  manifest: {args.manifest_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
