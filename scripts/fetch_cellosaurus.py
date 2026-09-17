#!/usr/bin/env python3
"""Fetch Cellosaurus records for cell lines in a benchmark.

The output supports manual evidence review; it does not assign benchmark
labels. Use the command-line help for single-line and batch options.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("pip install requests")

import pandas as pd

API = "https://api.cellosaurus.org"
PAUSE = 0.5  # Limit request frequency.

# Cellosaurus fields included in the review output.
INTERESTING = {
    "ID": "name",
    "AC": "accession (CVCL)",
    "SY": "synonyms",
    "CC": "comment",
    "DI": "disease",
    "DR": "cross-reference",
    "RX": "reference (PubMed)",
    "OX": "species",
    "CA": "category",
}


def query(cell_line: str, timeout: int = 30) -> str | None:
    """Search Cellosaurus by cell line name; return the flat-file text record."""
    params = {
        "q": f'id:"{cell_line}"',
        "format": "txt",
        "rows": "1",
        "fields": "id,ac,sy,cc,di,ox,ca,dr,rx,sequence-variation",
    }
    try:
        r = requests.get(f"{API}/search/cell-line", params=params, timeout=timeout)
        if r.status_code == 200 and r.text.strip():
            return r.text
        # Retry without the exact-name filter.
        params["q"] = cell_line
        r = requests.get(f"{API}/search/cell-line", params=params, timeout=timeout)
        if r.status_code == 200 and r.text.strip():
            return r.text
        print(f"    HTTP {r.status_code}", file=sys.stderr)
    except requests.RequestException as exc:
        print(f"    request failed: {exc}", file=sys.stderr)
    return None


def summarise(record: str, gene: str | None = None) -> str:
    """Return the review fields and optionally highlight a gene."""
    out, hits = [], []
    for line in record.splitlines():
        if len(line) < 2:
            continue
        code = line[:2]
        if code in INTERESTING or line.startswith((" ", "\t")):
            out.append(line)
        if gene and gene.upper() in line.upper():
            hits.append(line.strip())

    text = "\n".join(out) if out else record
    if gene:
        if hits:
            text += f"\n\n  >>> lines mentioning {gene.upper()}:\n"
            text += "\n".join(f"      {h}" for h in hits)
        else:
            text += (f"\n\n  >>> no mention of {gene.upper()} in this record.\n"
                     f"      Cellosaurus may not curate this feature. Try PubMed,\n"
                     f"      or change evidence_type to match what IS recorded.")
    return text


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gold", type=Path, help="gold standard csv")
    ap.add_argument("--cell-line", help="single cell line name")
    ap.add_argument("--gene", help="highlight lines mentioning this gene")
    ap.add_argument("--out", type=Path, help="also write output to this file")
    args = ap.parse_args()

    if args.cell_line:
        jobs = [(args.cell_line, args.gene)]
    elif args.gold:
        gold = pd.read_csv(args.gold)
        if args.gene:
            gold = gold[gold["gene"].str.upper() == args.gene.upper()]
        seen, jobs = set(), []
        for _, r in gold.iterrows():
            key = (str(r["expected_cell_line"]), str(r["gene"]).upper())
            if key not in seen:
                seen.add(key)
                jobs.append(key)
    else:
        ap.error("give --gold or --cell-line")

    chunks = []
    for i, (cl, gene) in enumerate(jobs, 1):
        header = f"\n{'=' * 72}\n[{i}/{len(jobs)}]  {cl}" + (f"   (gene: {gene})" if gene else "")
        header += f"\n{'=' * 72}"
        print(header)
        chunks.append(header)

        rec = query(cl)
        body = summarise(rec, gene) if rec else (
            f"  NOT FOUND in Cellosaurus under the name '{cl}'.\n"
            f"  Try a synonym (e.g. 'NCIN87' vs 'NCI-N87', 'HS746T' vs 'Hs 746T'),\n"
            f"  or look it up manually at https://www.cellosaurus.org")
        print(body)
        chunks.append(body)

        if i < len(jobs):
            time.sleep(PAUSE)

    if args.out:
        args.out.write_text("\n".join(chunks), encoding="utf-8")
        print(f"\nwritten to {args.out}")

    print("\n" + "-" * 72)
    print("Cite Cellosaurus as: Bairoch A. The Cellosaurus, a cell-line knowledge")
    print("resource. J Biomol Tech. 2018;29(2):25-38. Record accessed <today>.")
    print("Copy real PMIDs from the RX / sequence-variation lines above into")
    print("source_hint, then set verified=yes. Do not guess.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
