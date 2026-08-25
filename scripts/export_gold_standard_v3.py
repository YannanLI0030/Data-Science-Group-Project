#!/usr/bin/env python3
"""Export verified expression-review rows into the formal v3 gold standard.

The candidate pool remains the editable review source. Only rows satisfying all
of the following enter the formal benchmark:

* verified == yes
* judgement is positive or negative
* benchmark_task is expression_suitability or both
* identifiers, evidence type, source URL, summary, and notes are present

The output keeps the v2-compatible columns used by ablation_runner.py and adds
structured provenance columns for auditability.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


REQUIRED_POOL_COLUMNS = {
    "gene",
    "DepMap_ID",
    "cell_line",
    "judgement",
    "benchmark_task",
    "evidence_type",
    "source_url",
    "evidence_summary",
    "verified",
    "review_notes",
}

OUTPUT_COLUMNS = [
    "gene",
    "expected_cell_line",
    "expected_depmap_id",
    "relation",
    "benchmark_task",
    "evidence_type",
    "source_url",
    "evidence_summary",
    "source_hint",
    "verified",
    "notes",
]


def _normalise(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def export_gold(pool: pd.DataFrame, expected_count: int | None = None) -> pd.DataFrame:
    missing = sorted(REQUIRED_POOL_COLUMNS - set(pool.columns))
    if missing:
        raise KeyError("candidate pool missing columns: " + ", ".join(missing))

    work = pool.copy()
    for column in REQUIRED_POOL_COLUMNS:
        work[column] = _normalise(work[column])

    verified_mask = work["verified"].str.lower().eq("yes")
    verified = work.loc[verified_mask].copy()
    if expected_count is not None and len(verified) != expected_count:
        raise ValueError(
            f"expected {expected_count} verified rows, found {len(verified)}"
        )
    if verified.empty:
        raise ValueError("candidate pool contains no verified=yes rows")

    invalid_relation = ~verified["judgement"].str.lower().isin(
        {"positive", "negative"}
    )
    if invalid_relation.any():
        raise ValueError(
            "verified rows must have positive/negative judgement:\n"
            + verified.loc[
                invalid_relation, ["gene", "DepMap_ID", "judgement"]
            ].to_string(index=False)
        )

    invalid_task = ~verified["benchmark_task"].str.lower().isin(
        {"expression_suitability", "both"}
    )
    if invalid_task.any():
        raise ValueError(
            "verified rows outside the expression benchmark:\n"
            + verified.loc[
                invalid_task, ["gene", "DepMap_ID", "benchmark_task"]
            ].to_string(index=False)
        )

    required_nonblank = [
        "gene",
        "DepMap_ID",
        "cell_line",
        "evidence_type",
        "source_url",
        "evidence_summary",
        "review_notes",
    ]
    blank_messages: list[str] = []
    for column in required_nonblank:
        blank = verified[column].eq("")
        if blank.any():
            keys = verified.loc[blank, ["gene", "DepMap_ID"]]
            blank_messages.append(f"{column}:\n{keys.to_string(index=False)}")
    if blank_messages:
        raise ValueError("blank required review fields:\n" + "\n".join(blank_messages))

    verified["gene"] = verified["gene"].str.upper()
    verified["DepMap_ID"] = verified["DepMap_ID"].str.upper()
    valid_id = verified["DepMap_ID"].str.fullmatch(r"ACH-\d{6}")
    if not valid_id.all():
        raise ValueError(
            "invalid DepMap IDs:\n"
            + verified.loc[~valid_id, ["gene", "DepMap_ID", "cell_line"]]
            .to_string(index=False)
        )

    duplicates = verified.duplicated(["gene", "DepMap_ID"], keep=False)
    if duplicates.any():
        raise ValueError(
            "duplicate gene/DepMap IDs:\n"
            + verified.loc[duplicates, ["gene", "DepMap_ID", "cell_line"]]
            .to_string(index=False)
        )

    output = pd.DataFrame(
        {
            "gene": verified["gene"],
            "expected_cell_line": verified["cell_line"],
            "expected_depmap_id": verified["DepMap_ID"],
            "relation": verified["judgement"].str.lower(),
            "benchmark_task": verified["benchmark_task"].str.lower(),
            "evidence_type": verified["evidence_type"],
            "source_url": verified["source_url"],
            "evidence_summary": verified["evidence_summary"],
            "source_hint": (
                verified["evidence_summary"]
                + "; Sources: "
                + verified["source_url"]
            ),
            "verified": "yes",
            "notes": verified["review_notes"],
        }
    ).reset_index(drop=True)

    if list(output.columns) != OUTPUT_COLUMNS:
        raise AssertionError("unexpected output column order")
    return output


def validate_export(gold: pd.DataFrame, expected_count: int | None = None) -> None:
    if list(gold.columns) != OUTPUT_COLUMNS:
        raise ValueError("gold v3 field order does not match the formal schema")
    if expected_count is not None and len(gold) != expected_count:
        raise ValueError(f"expected {expected_count} rows, found {len(gold)}")
    if gold.isna().any().any():
        columns = gold.columns[gold.isna().any()].tolist()
        raise ValueError("gold v3 contains missing values in: " + ", ".join(columns))
    if not gold["expected_depmap_id"].astype(str).str.fullmatch(r"ACH-\d{6}").all():
        raise ValueError("gold v3 contains invalid DepMap IDs")
    if gold.duplicated(["gene", "expected_depmap_id"]).any():
        raise ValueError("gold v3 contains duplicate gene/DepMap IDs")
    if set(gold["relation"]) - {"positive", "negative"}:
        raise ValueError("gold v3 contains invalid relations")
    if set(gold["verified"]) != {"yes"}:
        raise ValueError("gold v3 must contain only verified=yes rows")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pool",
        type=Path,
        default=Path("benchmarks/candidate_pool_v3.csv"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("benchmarks/gold_standard_v3.csv"),
    )
    parser.add_argument("--expected-count", type=int, default=None)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    if args.out.exists() and not args.force:
        parser.error(f"output already exists: {args.out}; use --force to replace it")

    pool = pd.read_csv(args.pool, dtype=str, keep_default_na=False)
    gold = export_gold(pool, expected_count=args.expected_count)
    validate_export(gold, expected_count=args.expected_count)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.out.with_suffix(args.out.suffix + ".tmp")
    gold.to_csv(temporary, index=False, encoding="utf-8-sig")
    temporary.replace(args.out)

    print("gold_standard_v3 created")
    print(f"rows: {len(gold)}")
    print(f"genes: {gold['gene'].nunique()}")
    print(gold.groupby(["gene", "relation"], sort=False).size().to_string())
    print(f"output: {args.out.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
