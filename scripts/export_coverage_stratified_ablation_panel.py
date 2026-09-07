#!/usr/bin/env python3
"""Build a reproducible, coverage-stratified 100-gene ablation panel.

This script intentionally does not replace the teammate's original purely
random exporter.  The old panel is useful for missing-data stress tests; this
panel makes the protein, confidence, and adaptive-trust ablations identifiable
by sampling genes from pre-declared data-coverage strata.

Selection uses coverage metadata only.  It never inspects model rankings or
gold labels, which prevents performance-driven gene cherry-picking.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import re
import subprocess
import sys
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]

DEFAULT_EXCLUDED_GENES = (
    "AFP", "ALB", "AR", "DES", "ESR1", "GFAP", "MET", "MITF", "MSLN",
    "PTPRC", "EGFR", "KRAS", "ERBB2", "MYC", "TP53", "FGFR2", "CD86",
    "KLK4", "GAPDH", "ASGR1", "MUC1", "CD3E",
)

# Total = 100.  The first three strata identify adaptive trust; the fourth
# identifies missing-correlation behaviour; the final two measure robustness
# when protein evidence is unavailable.
DEFAULT_QUOTAS = OrderedDict(
    (
        ("protein_corr_low", 10),
        ("protein_corr_mid", 10),
        ("protein_corr_high", 10),
        ("protein_no_reliable_corr", 20),
        ("multisource_rna_no_protein", 25),
        ("depmap_rna_only", 25),
    )
)

REQUIRED_DATA_FILES = (
    Path("gene expression/1_4_hpa_rna_celline.tsv"),
    Path("gene expression/2_DepMap_OmicsExpressionAllGenesTPMLogp1Profile.csv"),
    Path("gene expression/3_GEOexpression.txt"),
    Path("gene expression/4_Harmonized_MS_CCLE_Gygi_subsetted.csv"),
    Path("gene properties/5_OmicsFusionFilteredSupplementary.csv"),
    Path("gene properties/6_OmicsSomaticMutationsProfile.csv"),
    Path("nomenclature/7_cellosaurus.csv"),
    Path("nomenclature/8_DepMap_OmicsProfiles.csv"),
    Path("nomenclature/9_DepMap_sample_info.csv"),
    Path("nomenclature/10_GEOInfo.txt"),
    Path("nomenclature/11_hpa_rna_celline_description.tsv"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export a coverage-stratified 100-gene ablation panel."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        help="Local raw-data root containing gene expression/, gene properties/, and nomenclature/.",
    )
    parser.add_argument(
        "--corr",
        type=Path,
        help="CSV with gene, correlation, and preferably n_cell_lines columns.",
    )
    parser.add_argument(
        "--src-dir",
        type=Path,
        default=PROJECT_ROOT / "src",
        help="Directory containing data_loader.py and data_merger.py.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=PROJECT_ROOT / "benchmarks" / "panels",
        help="Output directory for list, long-table Parquet, catalog, and manifest.",
    )
    parser.add_argument("--seed", type=int, default=42)
    # Match the frozen A4 adaptive-trust breakpoints used by the controlled
    # ablation runner: <=0.25 low, 0.25-0.50 mid, >=0.50 high.
    parser.add_argument("--low-corr-max", type=float, default=0.25)
    parser.add_argument("--high-corr-min", type=float, default=0.50)
    parser.add_argument("--min-correlation-pairs", type=int, default=10)
    parser.add_argument("--min-protein-cell-lines", type=int, default=10)
    parser.add_argument(
        "--min-rich-rna-sources",
        type=int,
        default=2,
        choices=(1, 2, 3),
        help="Minimum RNA sources for low/mid/high reliable-correlation strata.",
    )
    parser.add_argument(
        "--exclude-gene",
        action="append",
        default=[],
        help="Additional gene to exclude; may be supplied repeatedly.",
    )
    parser.add_argument(
        "--exclude-file",
        type=Path,
        action="append",
        default=[],
        help="CSV of additional genes to exclude; may be supplied repeatedly.",
    )
    parser.add_argument(
        "--skip-long-table",
        action="store_true",
        help="Generate the selection list/catalog/manifest without rebuilding Parquet.",
    )
    parser.add_argument(
        "--verify-input-hashes",
        action="store_true",
        help="Compute SHA-256 for every large raw input (slower but strongest provenance).",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite this panel's outputs.")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if not args.self_test and (args.data_dir is None or args.corr is None):
        parser.error("--data-dir and --corr are required (or use --self-test)")
    return args


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def normalise_key(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip().upper()


def parse_depmap_symbol(column: object) -> str:
    value = str(column).strip()
    match = re.match(r"^(.*?)\s+\(ENSG\d+(?:\.\d+)?\)$", value)
    return (match.group(1) if match else value).strip()


def parse_protein_symbol(column: object) -> str:
    value = str(column).strip()
    match = re.search(r"\(([^)]+)\)\s*$", value)
    return (match.group(1) if match else value).strip()


def find_column(columns: Iterable[object], candidates: Iterable[str], label: str) -> str:
    by_normalised = {
        re.sub(r"[^a-z0-9]", "", str(column).lower()): str(column)
        for column in columns
    }
    for candidate in candidates:
        key = re.sub(r"[^a-z0-9]", "", candidate.lower())
        if key in by_normalised:
            return by_normalised[key]
    raise ValueError(f"Could not detect {label}; columns={list(columns)}")


def load_hpa_index(path: Path, chunksize: int = 1_000_000) -> tuple[set[str], pd.DataFrame]:
    head = pd.read_csv(path, sep="\t", nrows=3)
    ensg_col = find_column(head.columns, ("Gene", "ensembl_id"), "HPA Ensembl column")
    symbol_col = find_column(head.columns, ("Gene name", "gene_name"), "HPA gene-symbol column")

    hpa_keys: set[str] = set()
    mapping_parts: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        path,
        sep="\t",
        usecols=[ensg_col, symbol_col],
        dtype="string",
        chunksize=chunksize,
    ):
        part = chunk.dropna().drop_duplicates()
        if part.empty:
            continue
        part = part.rename(columns={ensg_col: "ensembl_id", symbol_col: "gene_symbol"})
        part["ensembl_id"] = part["ensembl_id"].str.split(".").str[0].str.strip()
        part["gene_symbol"] = part["gene_symbol"].str.strip()
        part = part[(part["ensembl_id"] != "") & (part["gene_symbol"] != "")]
        hpa_keys.update(part["gene_symbol"].map(normalise_key))
        mapping_parts.append(part)

    mapping = (
        pd.concat(mapping_parts, ignore_index=True).drop_duplicates()
        if mapping_parts
        else pd.DataFrame(columns=["ensembl_id", "gene_symbol"])
    )
    return hpa_keys, mapping


def load_geo_keys(path: Path, ensg_to_symbol: Mapping[str, str]) -> set[str]:
    gene_column = pd.read_csv(path, sep="\t", nrows=0).columns[0]
    ids = pd.read_csv(path, sep="\t", usecols=[gene_column], dtype="string")[gene_column]
    ids = ids.str.split(".").str[0].str.strip()
    if ids.str.startswith("ENSG", na=False).mean() > 0.5:
        symbols = ids.map(ensg_to_symbol)
    else:
        symbols = ids
    return {normalise_key(value) for value in symbols.dropna() if normalise_key(value)}


def load_correlation(path: Path) -> pd.DataFrame:
    corr = pd.read_csv(path)
    gene_col = find_column(corr.columns, ("gene", "gene_symbol", "symbol"), "correlation gene")
    corr_col = find_column(
        corr.columns,
        ("correlation", "pearson_r", "rna_protein_correlation", "corr"),
        "correlation value",
    )
    try:
        n_col = find_column(
            corr.columns,
            ("n_cell_lines", "n", "sample_size", "n_pairs"),
            "correlation sample size",
        )
    except ValueError:
        n_col = None

    out = pd.DataFrame(
        {
            "gene_key": corr[gene_col].map(normalise_key),
            "correlation": pd.to_numeric(corr[corr_col], errors="coerce"),
            "correlation_n_cell_lines": (
                pd.to_numeric(corr[n_col], errors="coerce") if n_col else np.nan
            ),
        }
    )
    out = out[(out["gene_key"] != "") & out["correlation"].notna()]
    # If duplicates exist, retain the estimate with the largest stated sample size.
    out = out.sort_values(
        ["gene_key", "correlation_n_cell_lines"], na_position="last"
    ).drop_duplicates("gene_key", keep="last")
    return out.set_index("gene_key")


def build_coverage_catalog(
    data_dir: Path,
    corr_path: Path,
    excluded: set[str],
    min_corr_pairs: int,
    min_protein_cell_lines: int,
    min_rich_rna_sources: int,
    low_corr_max: float,
    high_corr_min: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rna_path = data_dir / "gene expression" / "2_DepMap_OmicsExpressionAllGenesTPMLogp1Profile.csv"
    protein_path = data_dir / "gene expression" / "4_Harmonized_MS_CCLE_Gygi_subsetted.csv"
    hpa_path = data_dir / "gene expression" / "1_4_hpa_rna_celline.tsv"
    geo_path = data_dir / "gene expression" / "3_GEOexpression.txt"

    rna_columns = pd.read_csv(rna_path, nrows=0).columns[1:]
    symbol_by_key: dict[str, str] = {}
    for column in rna_columns:
        symbol = parse_depmap_symbol(column)
        key = normalise_key(symbol)
        if key and key not in excluded:
            symbol_by_key.setdefault(key, symbol)

    protein = pd.read_csv(protein_path, index_col=0)
    protein_counts: dict[str, int] = {}
    for column, count in protein.notna().sum(axis=0).items():
        key = normalise_key(parse_protein_symbol(column))
        if key:
            protein_counts[key] = max(protein_counts.get(key, 0), int(count))

    hpa_keys, ensembl_map = load_hpa_index(hpa_path)
    ensg_to_symbol = dict(zip(ensembl_map["ensembl_id"], ensembl_map["gene_symbol"]))
    geo_keys = load_geo_keys(geo_path, ensg_to_symbol)
    correlation = load_correlation(corr_path)

    rows = []
    for key, symbol in symbol_by_key.items():
        protein_n = protein_counts.get(key, 0)
        has_hpa = key in hpa_keys
        has_geo = key in geo_keys
        corr_value = correlation.at[key, "correlation"] if key in correlation.index else np.nan
        corr_n = (
            correlation.at[key, "correlation_n_cell_lines"]
            if key in correlation.index
            else np.nan
        )
        has_reliable_corr = bool(
            protein_n >= min_protein_cell_lines
            and pd.notna(corr_value)
            and (pd.isna(corr_n) or float(corr_n) >= min_corr_pairs)
        )
        rows.append(
            {
                "gene": symbol,
                "gene_key": key,
                "has_depmap_rna": True,
                "has_hpa_rna": has_hpa,
                "has_geo_rna": has_geo,
                "rna_source_count": 1 + int(has_hpa) + int(has_geo),
                "protein_n_cell_lines": protein_n,
                "has_protein": protein_n >= min_protein_cell_lines,
                "correlation": corr_value,
                "correlation_n_cell_lines": corr_n,
                "has_reliable_correlation": has_reliable_corr,
            }
        )

    catalog = pd.DataFrame(rows).sort_values("gene").reset_index(drop=True)
    rich = catalog["rna_source_count"] >= min_rich_rna_sources
    reliable = catalog["has_reliable_correlation"] & rich

    catalog["eligible_stratum"] = "not_selected_pool"
    catalog.loc[
        reliable & (catalog["correlation"] <= low_corr_max), "eligible_stratum"
    ] = "protein_corr_low"
    catalog.loc[
        reliable
        & (catalog["correlation"] > low_corr_max)
        & (catalog["correlation"] < high_corr_min),
        "eligible_stratum",
    ] = "protein_corr_mid"
    catalog.loc[
        reliable & (catalog["correlation"] >= high_corr_min), "eligible_stratum"
    ] = "protein_corr_high"
    catalog.loc[
        catalog["has_protein"] & ~catalog["has_reliable_correlation"],
        "eligible_stratum",
    ] = "protein_no_reliable_corr"
    catalog.loc[
        ~catalog["has_protein"] & (catalog["rna_source_count"] >= 2),
        "eligible_stratum",
    ] = "multisource_rna_no_protein"
    catalog.loc[
        ~catalog["has_protein"] & (catalog["rna_source_count"] == 1),
        "eligible_stratum",
    ] = "depmap_rna_only"
    return catalog, ensembl_map


def select_stratified(
    catalog: pd.DataFrame,
    quotas: Mapping[str, int],
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    pieces = []
    for stratum, quota in quotas.items():
        eligible = catalog.loc[catalog["eligible_stratum"] == stratum].copy()
        eligible = eligible.sort_values("gene").reset_index(drop=True)
        if len(eligible) < quota:
            raise ValueError(
                f"Stratum {stratum!r} has only {len(eligible)} eligible genes; "
                f"requested {quota}. Adjust thresholds/quotas before running the model."
            )
        positions = rng.choice(len(eligible), size=quota, replace=False)
        chosen = eligible.iloc[positions].copy()
        chosen["selection_stratum"] = stratum
        chosen["selection_rank_in_stratum"] = np.arange(1, quota + 1)
        pieces.append(chosen)

    selected = pd.concat(pieces, ignore_index=True)
    if selected["gene_key"].duplicated().any():
        duplicates = selected.loc[selected["gene_key"].duplicated(False), "gene"].tolist()
        raise AssertionError(f"Selection contains duplicate gene symbols: {duplicates}")
    selected.insert(0, "panel_order", np.arange(1, len(selected) + 1))
    return selected


def load_project_modules(src_dir: Path):
    src_dir = src_dir.resolve()
    if not (src_dir / "data_loader.py").exists() or not (src_dir / "data_merger.py").exists():
        raise FileNotFoundError(
            f"Expected data_loader.py and data_merger.py under --src-dir: {src_dir}"
        )
    sys.path.insert(0, str(src_dir))
    for name in ("data_loader", "data_merger"):
        sys.modules.pop(name, None)
    return importlib.import_module("data_loader"), importlib.import_module("data_merger")


def build_long_table(
    selected: pd.DataFrame,
    data_dir: Path,
    src_dir: Path,
    ensembl_map: pd.DataFrame,
) -> pd.DataFrame:
    data_loader, data_merger = load_project_modules(src_dir)
    loader = data_loader.CellLineDataLoader(data_dir)
    cell_resolver = data_merger.CellLineIDResolver(
        pd.read_csv(data_dir / "nomenclature" / "7_cellosaurus.csv"),
        loader.sample_info,
    )
    gene_resolver = data_merger.GeneIDResolver(ensembl_map)
    merger = (
        data_merger.MultiOmicsMerger(loader.sample_info["DepMap_ID"], cell_resolver, gene_resolver)
        .register(data_merger.DepMapRNASource(loader))
        .register(
            data_merger.HPARNASource(
                data_dir / "gene expression" / "1_4_hpa_rna_celline.tsv",
                data_dir / "nomenclature" / "11_hpa_rna_celline_description.tsv",
            )
        )
        .register(
            data_merger.GEOWideSource(
                data_dir / "gene expression" / "3_GEOexpression.txt",
                data_dir / "nomenclature" / "10_GEOInfo.txt",
            )
        )
        .register(data_merger.MatrixProteinSource(loader))
        .register(data_merger.MutationSource(loader))
        .register(data_merger.FusionSource(loader))
    )
    return merger.build_long_table(genes=set(selected["gene"]))


def add_observed_long_counts(selected: pd.DataFrame, long_table: pd.DataFrame) -> pd.DataFrame:
    result = selected.copy()
    if long_table.empty:
        return result
    counts = (
        long_table.groupby(["gene_symbol", "source"], dropna=False)["DepMap_ID"]
        .nunique()
        .unstack(fill_value=0)
    )
    counts.columns = [f"observed_{source}_cell_lines" for source in counts.columns]
    return result.merge(counts, left_on="gene", right_index=True, how="left").fillna(
        {column: 0 for column in counts.columns}
    )


def align_to_sample_info(
    long_table: pd.DataFrame, data_dir: Path
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Keep only IDs in the matching local sample-info snapshot.

    Cellosaurus can resolve historical ACH identifiers that are valid registry
    records but absent from this project's DepMap sample universe.  Retaining
    them in the exported evidence would make candidate counts depend on a
    different metadata version, so they are audited and removed here.
    """
    sample_path = data_dir / "nomenclature" / "9_DepMap_sample_info.csv"
    sample_info = pd.read_csv(sample_path, usecols=["DepMap_ID"])
    valid_ids = set(sample_info["DepMap_ID"].dropna().astype(str))
    observed_ids = set(long_table["DepMap_ID"].dropna().astype(str))
    unmatched_ids = observed_ids - valid_ids
    keep = long_table["DepMap_ID"].astype(str).isin(valid_ids)
    aligned = long_table.loc[keep].copy()
    return aligned, {
        "sample_info_depmap_ids": len(valid_ids),
        "long_table_depmap_ids_before_filter": len(observed_ids),
        "unmatched_depmap_ids_before_filter": len(unmatched_ids),
        "rows_dropped_for_unmatched_depmap_id": int((~keep).sum()),
        "long_table_depmap_ids_after_filter": int(aligned["DepMap_ID"].nunique()),
    }


def duplicate_audit(long_table: pd.DataFrame) -> dict[str, int]:
    if long_table.empty:
        return {"duplicate_rows": 0, "duplicate_groups": 0, "conflicting_value_groups": 0}
    keys = ["gene_symbol", "source", "DepMap_ID"]
    group_sizes = long_table.groupby(keys, dropna=False).size()
    duplicate_groups = int((group_sizes > 1).sum())
    duplicate_rows = int(long_table.duplicated(keys, keep=False).sum())
    value_counts = long_table.groupby(keys, dropna=False)["value_raw"].nunique(dropna=True)
    conflicts = int((value_counts > 1).sum())
    return {
        "duplicate_rows": duplicate_rows,
        "duplicate_groups": duplicate_groups,
        "conflicting_value_groups": conflicts,
    }


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


def raw_file_records(
    data_dir: Path,
    raw_manifest_path: Path,
    verify_hashes: bool,
) -> list[dict[str, object]]:
    recorded_hashes: dict[str, str] = {}
    if raw_manifest_path.exists():
        try:
            raw_manifest = json.loads(raw_manifest_path.read_text(encoding="utf-8"))
            recorded_hashes = {
                item["file"]: item.get("raw_sha256")
                for item in raw_manifest.get("files", [])
                if item.get("file")
            }
        except (OSError, ValueError, TypeError):
            recorded_hashes = {}

    records = []
    for relative in REQUIRED_DATA_FILES:
        path = data_dir / relative
        if not path.exists():
            raise FileNotFoundError(f"Required raw input missing: {path}")
        actual_hash = sha256_file(path) if verify_hashes else None
        expected_hash = recorded_hashes.get(path.name)
        if actual_hash and expected_hash and actual_hash != expected_hash:
            raise ValueError(f"SHA-256 mismatch for {path}: {actual_hash} != {expected_hash}")
        records.append(
            {
                "relative_path": relative.as_posix(),
                "bytes": path.stat().st_size,
                "sha256": actual_hash or expected_hash,
                "hash_source": "computed" if actual_hash else ("data_manifest" if expected_hash else None),
            }
        )
    return records


def ensure_writable_outputs(paths: Iterable[Path], force: bool) -> None:
    existing = [path for path in paths if path.exists()]
    if existing and not force:
        joined = "\n  ".join(str(path) for path in existing)
        raise FileExistsError(f"Refusing to overwrite existing output(s):\n  {joined}\nUse --force.")


def genes_from_exclusion_file(path: Path) -> set[str]:
    table = pd.read_csv(path)
    if table.empty and not len(table.columns):
        return set()
    try:
        gene_col = find_column(table.columns, ("gene", "gene_symbol", "symbol"), "gene column")
    except ValueError:
        if len(table.columns) != 1:
            raise
        gene_col = str(table.columns[0])
    return {normalise_key(value) for value in table[gene_col] if normalise_key(value)}


def run_self_test() -> None:
    rows = []
    quotas = OrderedDict((name, 2) for name in DEFAULT_QUOTAS)
    for stratum in quotas:
        for index in range(5):
            rows.append(
                {
                    "gene": f"{stratum}_{index}",
                    "gene_key": f"{stratum}_{index}".upper(),
                    "eligible_stratum": stratum,
                }
            )
    catalog = pd.DataFrame(rows)
    first = select_stratified(catalog, quotas, seed=42)
    second = select_stratified(catalog, quotas, seed=42)
    assert first["gene"].tolist() == second["gene"].tolist()
    assert len(first) == 12 and first["gene_key"].nunique() == 12
    assert first.groupby("selection_stratum").size().to_dict() == dict(quotas)
    print("self-test OK")
    print("  deterministic stratified sampling")
    print("  quotas and unique-gene invariant enforced")


def main() -> int:
    args = parse_args()
    if args.self_test:
        run_self_test()
        return 0

    data_dir = args.data_dir.expanduser().resolve()
    corr_path = args.corr.expanduser().resolve()
    src_dir = args.src_dir.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    raw_manifest_path = PROJECT_ROOT / "data_manifest" / "manifest.json"

    if args.low_corr_max >= args.high_corr_min:
        raise ValueError("--low-corr-max must be smaller than --high-corr-min")
    if not corr_path.exists():
        raise FileNotFoundError(f"Correlation file not found: {corr_path}")

    excluded = {normalise_key(gene) for gene in DEFAULT_EXCLUDED_GENES}
    excluded.update(normalise_key(gene) for gene in args.exclude_gene)
    exclusion_files = [path.expanduser().resolve() for path in args.exclude_file]
    for exclusion_file in exclusion_files:
        if not exclusion_file.exists():
            raise FileNotFoundError(f"Exclusion file not found: {exclusion_file}")
        excluded.update(genes_from_exclusion_file(exclusion_file))

    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = "ablation_100genes_coverage_stratified"
    list_path = out_dir / f"{prefix}_list.csv"
    parquet_path = out_dir / f"{prefix}.parquet"
    catalog_path = out_dir / f"{prefix}_coverage_catalog.csv"
    manifest_path = out_dir / f"{prefix}_manifest.json"
    outputs = [list_path, catalog_path, manifest_path]
    if not args.skip_long_table:
        outputs.append(parquet_path)
    ensure_writable_outputs(outputs, args.force)

    print("Building gene-level coverage catalog...")
    catalog, ensembl_map = build_coverage_catalog(
        data_dir=data_dir,
        corr_path=corr_path,
        excluded=excluded,
        min_corr_pairs=args.min_correlation_pairs,
        min_protein_cell_lines=args.min_protein_cell_lines,
        min_rich_rna_sources=args.min_rich_rna_sources,
        low_corr_max=args.low_corr_max,
        high_corr_min=args.high_corr_min,
    )
    eligible_counts = catalog["eligible_stratum"].value_counts().to_dict()
    print("Eligible genes by stratum:")
    for stratum, quota in DEFAULT_QUOTAS.items():
        print(f"  {stratum}: {eligible_counts.get(stratum, 0)} eligible; selecting {quota}")

    selected = select_stratified(catalog, DEFAULT_QUOTAS, args.seed)
    selected_keys = set(selected["gene_key"])
    if selected_keys & excluded:
        raise AssertionError(f"Excluded genes selected: {sorted(selected_keys & excluded)}")

    long_table = pd.DataFrame()
    duplicate_summary = {
        "duplicate_rows": None,
        "duplicate_groups": None,
        "conflicting_value_groups": None,
    }
    sample_alignment = {
        "sample_info_depmap_ids": None,
        "long_table_depmap_ids_before_filter": None,
        "unmatched_depmap_ids_before_filter": None,
        "rows_dropped_for_unmatched_depmap_id": None,
        "long_table_depmap_ids_after_filter": None,
    }
    if not args.skip_long_table:
        print("Building the selected 100-gene unified long table...")
        long_table = build_long_table(selected, data_dir, src_dir, ensembl_map)
        long_table, sample_alignment = align_to_sample_info(long_table, data_dir)
        observed_genes = set(long_table["gene_symbol"].dropna())
        missing_genes = sorted(set(selected["gene"]) - observed_genes)
        if missing_genes:
            raise AssertionError(f"Selected genes absent from long-table output: {missing_genes}")
        selected = add_observed_long_counts(selected, long_table)
        duplicate_summary = duplicate_audit(long_table)

    catalog_to_write = catalog.copy()
    catalog_to_write["selected"] = catalog_to_write["gene_key"].isin(selected_keys)
    catalog_to_write = catalog_to_write.merge(
        selected[["gene_key", "selection_stratum", "panel_order"]],
        on="gene_key",
        how="left",
    )

    selected.to_csv(list_path, index=False)
    catalog_to_write.to_csv(catalog_path, index=False)
    if not args.skip_long_table:
        long_table.to_parquet(parquet_path, index=False)

    script_path = Path(__file__).resolve()
    source_files = [src_dir / "data_loader.py", src_dir / "data_merger.py"]
    manifest = {
        "schema_version": 1,
        "panel_name": "ablation_100genes_coverage_stratified",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "intended_role": [
            "controlled_component_ablation",
            "coverage-stratified_missing-data_robustness",
        ],
        "not_yet_valid_for": [
            "NDCG/Recall/MRR without independently verified gold labels",
            "final holdout evaluation before configurations are frozen",
        ],
        "selection": {
            "method": "coverage-only stratified random sampling without replacement",
            "seed": args.seed,
            "quotas": dict(DEFAULT_QUOTAS),
            "thresholds": {
                "low_correlation_max": args.low_corr_max,
                "high_correlation_min": args.high_corr_min,
                "min_correlation_pairs": args.min_correlation_pairs,
                "min_protein_cell_lines": args.min_protein_cell_lines,
                "min_rich_rna_sources": args.min_rich_rna_sources,
            },
            "excluded_genes": sorted(excluded),
            "excluded_gene_files": [
                {
                    "path": str(path),
                    "sha256": sha256_file(path),
                }
                for path in exclusion_files
            ],
            "candidate_gene_count_after_exclusion": int(len(catalog)),
            "eligible_counts": {key: int(value) for key, value in eligible_counts.items()},
            "selected_gene_count": int(len(selected)),
        },
        "generator": {
            "script": str(script_path),
            "script_sha256": sha256_file(script_path),
            "project_git": git_state(PROJECT_ROOT),
            "src_dir": str(src_dir),
            "source_code": {
                path.name: sha256_file(path) for path in source_files if path.exists()
            },
            "python": sys.version,
            "pandas": pd.__version__,
            "numpy": np.__version__,
        },
        "inputs": {
            "data_dir": str(data_dir),
            "raw_files": raw_file_records(
                data_dir, raw_manifest_path, args.verify_input_hashes
            ),
            "correlation_file": {
                "path": str(corr_path),
                "bytes": corr_path.stat().st_size,
                "sha256": sha256_file(corr_path),
            },
        },
        "long_table_audit": {
            "generated": not args.skip_long_table,
            "rows": int(len(long_table)) if not args.skip_long_table else None,
            "genes": (
                int(long_table["gene_symbol"].nunique()) if not args.skip_long_table else None
            ),
            "depmap_ids": (
                int(long_table["DepMap_ID"].nunique()) if not args.skip_long_table else None
            ),
            "sample_info_alignment": sample_alignment,
            **duplicate_summary,
        },
        "outputs": {},
    }
    for name, path in (
        ("panel_list", list_path),
        ("coverage_catalog", catalog_path),
        ("long_table", parquet_path),
    ):
        if path.exists():
            manifest["outputs"][name] = {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }

    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    print("Coverage-stratified panel created")
    print(f"  genes: {len(selected)}")
    print(f"  list: {list_path}")
    print(f"  catalog: {catalog_path}")
    if not args.skip_long_table:
        print(f"  long table: {parquet_path} ({len(long_table):,} rows)")
        print(f"  sample-info alignment: {sample_alignment}")
        print(f"  duplicate audit: {duplicate_summary}")
    print(f"  manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
