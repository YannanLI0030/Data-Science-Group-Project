#!/usr/bin/env python3
"""Export the verified non-classic v4 development benchmark.

The editable source of truth remains ``candidate_pool_v4_review.csv``.  This
script selects only rows assigned to ``weight_tuning`` and then delegates the
formal gold-schema validation to ``export_gold_standard_v3``.  Classic case
study genes are excluded defensively even if a row is accidentally assigned
the wrong role.

The current reviewed v4 benchmark is frozen at 50 rows across 10 genes:
48 positive and 2 negative labels.  These expectations are checked by default
so a partially saved or accidentally edited review file cannot silently become
the benchmark used by ablation and weight-search experiments.
"""

from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path

import pandas as pd

import export_gold_standard_v3 as v3


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEVELOPMENT_ROLE = "weight_tuning"
CLASSIC_GENES = frozenset({"EGFR", "KRAS", "ERBB2", "MYC", "TP53"})
REQUIRED_V4_COLUMNS = v3.REQUIRED_POOL_COLUMNS | {"benchmark_role"}


def read_pool(path: Path, input_encoding: str = "auto") -> tuple[pd.DataFrame, str]:
    """Read the review CSV without changing it and report the chosen encoding."""
    raw = path.read_bytes()
    encodings = (
        [input_encoding]
        if input_encoding != "auto"
        else ["utf-8-sig", "cp1250", "cp1252"]
    )
    failures: list[str] = []
    for encoding in encodings:
        try:
            text = raw.decode(encoding)
            frame = pd.read_csv(io.StringIO(text), dtype=str, keep_default_na=False)
            return frame, encoding
        except (UnicodeDecodeError, pd.errors.ParserError) as exc:
            failures.append(f"{encoding}: {exc}")
    raise UnicodeError(
        f"could not decode/parse candidate pool {path}; " + " | ".join(failures)
    )


def build_development_gold(
    pool: pd.DataFrame,
    *,
    expected_count: int = 50,
    expected_genes: int = 10,
    expected_positive: int = 48,
    expected_negative: int = 2,
) -> pd.DataFrame:
    missing = sorted(REQUIRED_V4_COLUMNS - set(pool.columns))
    if missing:
        raise KeyError("candidate pool missing v4 columns: " + ", ".join(missing))

    work = pool.copy()
    work["benchmark_role"] = (
        work["benchmark_role"].fillna("").astype(str).str.strip().str.lower()
    )
    work["gene"] = work["gene"].fillna("").astype(str).str.strip().str.upper()

    development = work.loc[work["benchmark_role"].eq(DEVELOPMENT_ROLE)].copy()
    if development.empty:
        raise ValueError("candidate pool contains no weight_tuning rows")

    classic_overlap = sorted(set(development["gene"]) & CLASSIC_GENES)
    if classic_overlap:
        raise ValueError(
            "classic genes must not be assigned to weight_tuning: "
            + ", ".join(classic_overlap)
        )

    gold = v3.export_gold(development, expected_count=expected_count)
    v3.validate_export(gold, expected_count=expected_count)
    gold = gold.sort_values(
        ["gene", "relation", "expected_depmap_id"], kind="stable"
    ).reset_index(drop=True)

    observed_genes = int(gold["gene"].nunique())
    if observed_genes != expected_genes:
        raise ValueError(
            f"expected {expected_genes} development genes, found {observed_genes}"
        )
    relation_counts = gold["relation"].value_counts()
    observed_positive = int(relation_counts.get("positive", 0))
    observed_negative = int(relation_counts.get("negative", 0))
    if observed_positive != expected_positive or observed_negative != expected_negative:
        raise ValueError(
            "unexpected relation counts: "
            f"positive={observed_positive} (expected {expected_positive}), "
            f"negative={observed_negative} (expected {expected_negative})"
        )

    positive_genes = int(
        gold.loc[gold["relation"].eq("positive"), "gene"].nunique()
    )
    if positive_genes != expected_genes:
        raise ValueError(
            f"all development genes need a positive label; found {positive_genes} "
            f"positive-labelled genes out of {expected_genes}"
        )
    final_overlap = sorted(set(gold["gene"]) & CLASSIC_GENES)
    if final_overlap:
        raise AssertionError(
            "classic genes leaked into development gold: " + ", ".join(final_overlap)
        )
    return gold


def self_test() -> int:
    rows: list[dict[str, str]] = []
    base = {
        "benchmark_task": "expression_suitability",
        "evidence_type": "rna_expression",
        "source_url": "https://example.org/source",
        "evidence_summary": "External evidence summary.",
        "review_notes": "Reviewed for exporter self-test.",
    }
    for gene, depmap_id, judgement, role, verified in [
        ("AFP", "ACH-000001", "positive", "weight_tuning", "yes"),
        ("AFP", "ACH-000002", "unknown", "weight_tuning", "no"),
        ("EGFR", "ACH-000003", "positive", "classic_case_study", "yes"),
    ]:
        rows.append(
            {
                **base,
                "gene": gene,
                "DepMap_ID": depmap_id,
                "cell_line": depmap_id,
                "judgement": judgement,
                "benchmark_role": role,
                "verified": verified,
            }
        )
    gold = build_development_gold(
        pd.DataFrame(rows),
        expected_count=1,
        expected_genes=1,
        expected_positive=1,
        expected_negative=0,
    )
    assert len(gold) == 1
    assert gold.iloc[0]["gene"] == "AFP"
    assert gold.iloc[0]["expected_depmap_id"] == "ACH-000001"
    print("self-test OK")
    print("  unknown/unverified rows are excluded")
    print("  classic case-study rows are excluded")
    print("  formal gold schema remains v3-compatible")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pool",
        type=Path,
        default=REPOSITORY_ROOT / "benchmarks/candidate_pool_v4_review.csv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=REPOSITORY_ROOT / "benchmarks/gold_standard_v4_development.csv",
    )
    parser.add_argument(
        "--input-encoding",
        default="auto",
        help="source CSV encoding; auto tries UTF-8, cp1250, then cp1252",
    )
    parser.add_argument("--expected-count", type=int, default=50)
    parser.add_argument("--expected-genes", type=int, default=10)
    parser.add_argument("--expected-positive", type=int, default=48)
    parser.add_argument("--expected-negative", type=int, default=2)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        return self_test()
    if args.out.exists() and not args.force:
        parser.error(f"output already exists: {args.out}; use --force to replace it")

    pool, encoding = read_pool(args.pool, args.input_encoding)
    gold = build_development_gold(
        pool,
        expected_count=args.expected_count,
        expected_genes=args.expected_genes,
        expected_positive=args.expected_positive,
        expected_negative=args.expected_negative,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.out.with_suffix(args.out.suffix + ".tmp")
    gold.to_csv(temporary, index=False, encoding="utf-8-sig")
    temporary.replace(args.out)

    print("gold_standard_v4_development created")
    print(f"input encoding: {encoding}")
    print(f"rows: {len(gold)}")
    print(f"genes: {gold['gene'].nunique()}")
    print(gold.groupby(["gene", "relation"], sort=True).size().to_string())
    print(f"output: {args.out.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
