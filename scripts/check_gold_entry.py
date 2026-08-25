#!/usr/bin/env python3
"""
check_gold_entry.py
===================
Sanity-check a gold-standard entry against the project's own merged data.

WHAT THIS DOES AND DOES NOT DO
------------------------------
This does NOT establish ground truth. The biological claim ("HCC827 is
EGFR-dependent") must come from an external source: DepMap mutation records,
Cellosaurus, HPA, or a paper. Using your own ranking to define your own
benchmark would be circular reasoning.

What this DOES is check that the claim is *detectable in the data the system
can see*. A gold-standard entry only works as a benchmark if the evidence is
present. Three outcomes, all informative:

  - Cell line ranks high  -> entry is usable, keep it.
  - Cell line missing     -> not in the 1840-line panel, drop the entry.
  - Cell line ranks low   -> INVESTIGATE. Either the literature claim is about
                             something other than RNA abundance (e.g. protein
                             activity, or a mutation that doesn't raise
                             expression), or there is a data coverage gap.
                             Both are worth knowing before you rely on it.

Usage
-----
    python check_gold_entry.py --data-dir ./merged EGFR HCC827
    python check_gold_entry.py --data-dir ./merged EGFR            # top 20 only
    python check_gold_entry.py --data-dir ./merged --batch gold_standard.csv
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


def norm(s: str) -> str:
    """Cell-line names vary wildly: 'Hs 746T' / 'HS746T' / 'NCI-H1975'."""
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def load(data_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    master = pd.read_csv(data_dir / "master_table.csv")
    ann = pd.read_csv(data_dir / "cellline_annotations.csv")
    return master, ann


def gene_ranking(master: pd.DataFrame, ann: pd.DataFrame, gene: str,
                 value_col: str = "depmap_rna__raw") -> pd.DataFrame | None:
    gene = gene.upper().strip()
    sub = master[master["gene"].str.upper() == gene]
    if sub.empty:
        return None

    keep = ["DepMap_ID", "stripped_cell_line_name", "lineage", "primary_disease"]
    keep = [c for c in keep if c in ann.columns]
    df = sub.merge(ann[keep], on="DepMap_ID", how="left")

    df = df.dropna(subset=[value_col]).sort_values(value_col, ascending=False)
    df = df.reset_index(drop=True)
    df["rank"] = df.index + 1
    df["pct"] = 100.0 * (1.0 - df.index / max(len(df) - 1, 1))
    return df


def report_one(df: pd.DataFrame, gene: str, cell_line: str | None,
               value_col: str, top_n: int = 15) -> None:
    show = ["rank", "stripped_cell_line_name", value_col, "lineage"]
    show = [c for c in show if c in df.columns]

    print(f"\n{gene}: {len(df)} cell lines with {value_col}")
    print(f"--- top {top_n} ---")
    print(df.head(top_n)[show].to_string(index=False))

    if not cell_line:
        return

    target = norm(cell_line)
    hit = df[df["stripped_cell_line_name"].map(norm) == target]
    if hit.empty:
        loose = df[df["stripped_cell_line_name"].map(norm).str.contains(target, na=False)]
        if loose.empty:
            print(f"\n  >> '{cell_line}' NOT FOUND for {gene}.")
            print("     Either the name differs, or this line has no value here.")
            print("     Drop the entry, or fix the spelling and re-check.")
            return
        hit = loose
        print(f"\n  (matched loosely on '{cell_line}')")

    for _, r in hit.iterrows():
        verdict = ("USABLE" if r["pct"] >= 90 else
                   "WEAK — investigate" if r["pct"] >= 50 else
                   "CONTRADICTED — investigate before using")
        print(f"\n  >> {r['stripped_cell_line_name']}  ({r['DepMap_ID']})")
        print(f"     rank {int(r['rank'])} of {len(df)}   top {100 - r['pct']:.1f}%"
              f"   {value_col} = {r[value_col]:.3f}")
        print(f"     lineage: {r.get('lineage', 'n/a')}")
        print(f"     -> {verdict}")


def batch(master: pd.DataFrame, ann: pd.DataFrame, gold_path: Path,
          value_col: str) -> None:
    gold = pd.read_csv(gold_path)
    rows = []
    for _, e in gold.iterrows():
        gene, cl = str(e["gene"]).upper(), str(e["expected_cell_line"])
        expected_id = str(e.get("expected_depmap_id", "")).upper().strip()
        if expected_id in {"", "NAN", "NONE"}:
            expected_id = ""
        rel = e.get("relation", "positive")
        df = gene_ranking(master, ann, gene, value_col)
        if df is None:
            rows.append({"gene": gene, "cell_line": cl, "relation": rel,
                         "status": "GENE NOT IN PANEL", "rank": None,
                         "n": None, "top_pct": None, "depmap_id": None})
            continue
        match_by = "depmap_id" if expected_id else "cell_line"
        if expected_id:
            hit = df[df["DepMap_ID"].astype(str).str.upper() == expected_id]
        else:
            hit = df[df["stripped_cell_line_name"].map(norm) == norm(cl)]
            if hit.empty:
                hit = df[
                    df["stripped_cell_line_name"]
                    .map(norm)
                    .str.contains(norm(cl), na=False)
                ]
        if hit.empty:
            master_gene = master[master["gene"].astype(str).str.upper() == gene]
            if expected_id:
                master_hit = master_gene[
                    master_gene["DepMap_ID"].astype(str).str.upper() == expected_id
                ]
            else:
                master_hit = master_gene.iloc[0:0]
            status = "VALUE MISSING" if not master_hit.empty else "CELL LINE NOT FOUND"
            rows.append({"gene": gene, "cell_line": cl, "relation": rel,
                         "status": status, "rank": None,
                         "n": len(df), "top_pct": None,
                         "depmap_id": expected_id or None,
                         "match_by": match_by})
            continue
        r = hit.iloc[0]
        top_pct = 100 - r["pct"]
        if rel == "positive":
            status = ("OK" if top_pct <= 10 else
                      "WEAK" if top_pct <= 50 else "CONTRADICTED")
        else:  # negative controls should sit LOW
            status = ("OK" if top_pct >= 75 else
                      "WEAK" if top_pct >= 50 else "CONTRADICTED")
        rows.append({"gene": gene, "cell_line": cl, "relation": rel,
                     "status": status, "rank": int(r["rank"]), "n": len(df),
                     "top_pct": round(top_pct, 1),
                     "depmap_id": r["DepMap_ID"], "match_by": match_by})

    out = pd.DataFrame(rows)
    print(out.to_string(index=False))
    print("\nsummary:")
    print(out["status"].value_counts().to_string())
    print("\nOK = evidence is where the literature says it should be.")
    print("WEAK / CONTRADICTED = check the source claim before setting verified=yes.")
    print("VALUE MISSING = the cell line exists, but this modality has no value.")
    print("NOT FOUND = the DepMap ID is absent from this gene panel or is incorrect.")
    print("\nNOTE: 'OK' here means detectable in the data, NOT that the biology is")
    print("      verified. You still need an external citation for each entry.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gene", nargs="?", help="gene symbol, e.g. EGFR")
    ap.add_argument("cell_line", nargs="?", help="cell line name, e.g. HCC827")
    ap.add_argument("--data-dir", type=Path, required=True,
                    help="folder containing master_table.csv and cellline_annotations.csv")
    ap.add_argument("--batch", type=Path, help="check every row of a gold standard csv")
    ap.add_argument("--value-col", default="depmap_rna__raw",
                    help="evidence column (default: depmap_rna__raw). Try "
                         "ccle_gygi_protein__raw or hpa_rna__raw.")
    ap.add_argument("--top", type=int, default=15)
    args = ap.parse_args()

    master, ann = load(args.data_dir)

    if args.batch:
        batch(master, ann, args.batch, args.value_col)
        return 0

    if not args.gene:
        ap.error("give a gene, or use --batch")

    df = gene_ranking(master, ann, args.gene, args.value_col)
    if df is None:
        print(f"'{args.gene.upper()}' is not in the merged gene panel.")
        panel = sorted(master["gene"].unique())
        print(f"Panel has {len(panel)} genes: {', '.join(panel[:20])} ...")
        return 1

    report_one(df, args.gene.upper(), args.cell_line, args.value_col, args.top)
    return 0


if __name__ == "__main__":
    sys.exit(main())
