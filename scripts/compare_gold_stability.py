#!/usr/bin/env python3
"""Compare selected ablation configurations under two gold standards."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


DEFAULT_CONFIGS = [
    "A0_team_baseline",
    "A1a_no_direct_protein",
    "A1b_no_protein_evidence",
    "A4_adaptive_trust",
]

DEFAULT_METRICS = [
    "P@5",
    "Recall@5",
    "NDCG@5",
    "P@10",
    "Recall@10",
    "NDCG@10",
    "MRR",
    "neg_sink",
]


def compare(
    old_summary: pd.DataFrame,
    new_summary: pd.DataFrame,
    old_detail: pd.DataFrame,
    new_detail: pd.DataFrame,
    configs: list[str],
    metrics: list[str],
) -> pd.DataFrame:
    for label, frame in {
        "old summary": old_summary,
        "new summary": new_summary,
        "old detail": old_detail,
        "new detail": new_detail,
    }.items():
        required = {"config"} | (set(metrics) if "detail" in label else set(metrics))
        if "detail" in label:
            required.add("gene")
        missing = sorted(required - set(frame.columns))
        if missing:
            raise KeyError(f"{label} missing columns: {', '.join(missing)}")

    old = old_summary.set_index("config")
    new = new_summary.set_index("config")
    absent = sorted(set(configs) - set(old.index) | (set(configs) - set(new.index)))
    if absent:
        raise KeyError("selected configs absent from summaries: " + ", ".join(absent))

    common_genes = sorted(set(old_detail["gene"]) & set(new_detail["gene"]))
    records: list[dict] = []
    for metric in metrics:
        old_values = old.loc[configs, metric].astype(float)
        new_values = new.loc[configs, metric].astype(float)
        old_ranks = old_values.rank(ascending=False, method="min")
        new_ranks = new_values.rank(ascending=False, method="min")
        rank_spearman = (
            old_ranks.corr(new_ranks, method="pearson")
            if old_ranks.nunique() > 1 and new_ranks.nunique() > 1
            else float("nan")
        )

        for config in configs:
            old_common = old_detail.loc[
                (old_detail["config"] == config)
                & (old_detail["gene"].isin(common_genes)),
                ["gene", metric],
            ].rename(columns={metric: "old_metric"})
            new_common = new_detail.loc[
                (new_detail["config"] == config)
                & (new_detail["gene"].isin(common_genes)),
                ["gene", metric],
            ].rename(columns={metric: "new_metric"})
            paired = old_common.merge(new_common, on="gene", validate="one_to_one")
            common_abs_delta = (
                (paired["new_metric"] - paired["old_metric"]).abs().mean()
                if not paired.empty
                else float("nan")
            )

            records.append(
                {
                    "config": config,
                    "metric": metric,
                    "old_value": old_values[config],
                    "new_value": new_values[config],
                    "expanded_benchmark_delta": new_values[config]
                    - old_values[config],
                    "old_rank_among_four": old_ranks[config],
                    "new_rank_among_four": new_ranks[config],
                    "rank_change": new_ranks[config] - old_ranks[config],
                    "rank_spearman_among_four": rank_spearman,
                    "n_common_genes": len(common_genes),
                    "common_gene_mean_abs_delta": common_abs_delta,
                }
            )
    return pd.DataFrame(records)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-summary", type=Path, required=True)
    parser.add_argument("--new-summary", type=Path, required=True)
    parser.add_argument("--old-detail", type=Path, required=True)
    parser.add_argument("--new-detail", type=Path, required=True)
    parser.add_argument("--configs", nargs="+", default=DEFAULT_CONFIGS)
    parser.add_argument("--metrics", nargs="+", default=DEFAULT_METRICS)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    result = compare(
        pd.read_csv(args.old_summary),
        pd.read_csv(args.new_summary),
        pd.read_csv(args.old_detail),
        pd.read_csv(args.new_detail),
        list(dict.fromkeys(args.configs)),
        list(dict.fromkeys(args.metrics)),
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.out.with_suffix(args.out.suffix + ".tmp")
    result.to_csv(temporary, index=False, encoding="utf-8-sig")
    temporary.replace(args.out)

    focus = result[result["metric"].isin(["NDCG@5", "NDCG@10", "MRR"])]
    print(focus.to_string(index=False))
    print(f"output: {args.out.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
