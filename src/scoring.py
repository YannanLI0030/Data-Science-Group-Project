"""
scoring.py
==========
Implements and tests the CellLineFinder confidence-score formula proposed
in `CellLineSelector_技术方案与需求分析.docx`:

    Final Score = 0.40 x RNA expression
                + 0.25 x Protein evidence
                + 0.15 x Data consistency
                + 0.10 x Tissue relevance
                + 0.10 x Data completeness
                - Exclusion-gene penalty

This module turns that formula into testable code on top of
`data_loader.CellLineDataLoader.build_gene_profile()`, and documents
where real data forced a design decision the original spec didn't
specify.

Key findings from testing against real data (see notebook for full
exploration):

  - RNA vs protein correlation across cell lines with BOTH measurements
    is ~0.86 -- strong enough that "data consistency" is a meaningful
    signal, not noise.
  - Only ~26% of candidate cell lines have BOTH RNA and protein data for
    a typical gene (EGFR test case: 370 / 1434). The other 74% have RNA
    only. This means protein evidence and data consistency CANNOT be
    scored as 0 for missing data without systematically penalising the
    74% of cell lines that simply weren't profiled by mass spec -- that
    would bias the ranking towards whichever cell lines happened to be
    chosen for the (smaller, separate) proteomics study, not towards
    biological relevance. See `_score_with_missing` below for the
    resolution used.
  - "Tissue relevance" is necessarily *relative to a target tissue*, not
    an absolute property of a cell line -- the original spec didn't
    define what it's relative to. This implementation takes the
    candidate's own lineage distribution as the relevance signal: cell
    lines from a lineage enriched in the high-expression tail score
    higher. A `target_lineage` argument lets the spec's idea of
    "tissue relevance" be made explicit by the user instead.

Usage
-----
    from data_loader import CellLineDataLoader
    from scoring import score_gene_profile

    loader = CellLineDataLoader("data")
    profile = loader.build_gene_profile("EGFR", exclude_genes=["KRAS"])
    scored = score_gene_profile(profile)
"""

from __future__ import annotations

import pandas as pd
import numpy as np


# Default weights, taken directly from the group's spec document.
DEFAULT_WEIGHTS = {
    "rna": 0.40,
    "protein": 0.25,
    "consistency": 0.15,
    "tissue": 0.10,
    "completeness": 0.10,
}


def _minmax(series: pd.Series) -> pd.Series:
    """Scales a numeric series to [0, 1]. NaNs are preserved (not filled)."""
    lo, hi = series.min(), series.max()
    if pd.isna(lo) or pd.isna(hi) or hi == lo:
        return series * 0.0  # degenerate case: everything ties at 0
    return (series - lo) / (hi - lo)


def _score_with_missing(raw_score: pd.Series, has_data: pd.Series) -> tuple[pd.Series, pd.Series]:
    """
    Splits a normalised [0,1] score into (score, weight_multiplier) so that
    missing data reduces a row's effective denominator instead of being
    treated as a real zero.

    Returns
    -------
    (score_filled, multiplier)
        score_filled : the input score with NaN replaced by 0 (safe to
            multiply by a weight without propagating NaN)
        multiplier : 1.0 where data exists, 0.0 where it's missing -- used
            to redistribute that row's weight across the remaining terms.
    """
    multiplier = has_data.astype(float)
    score_filled = raw_score.fillna(0.0)
    return score_filled, multiplier


def score_gene_profile(
    profile: pd.DataFrame,
    weights: dict[str, float] | None = None,
    target_lineage: str | None = None,
) -> pd.DataFrame:
    """
    Applies the CellLineFinder scoring formula to the output of
    `CellLineDataLoader.build_gene_profile()`.

    Parameters
    ----------
    profile : DataFrame
        Output of `build_gene_profile()`. Must contain columns:
        ['DepMap_ID', 'lineage', 'rna_expression', 'protein_expression',
         'excluded'].
    weights : dict, optional
        Override any of the five DEFAULT_WEIGHTS terms. Unspecified terms
        keep their default value. Weights are NOT required to sum to 1;
        they're applied as given (see "Design notes" below).
    target_lineage : str, optional
        If given, "tissue relevance" becomes 1.0 for cell lines matching
        this lineage and 0.0 otherwise (the literal reading of the spec).
        If omitted (default), relevance is computed as each lineage's
        share of the top RNA-expression quartile vs. its overall share in
        the candidate pool -- i.e. "is this tissue over-represented among
        the highest-expressing candidates for this gene".

    Returns
    -------
    DataFrame
        Input `profile` with added columns:
        ['rna_score', 'protein_score', 'consistency_score',
         'tissue_score', 'completeness_score', 'final_score']
        sorted by `final_score` descending. Excluded rows are kept but
        sorted to the bottom (final_score forced to -1) so the caller can
        see what was excluded and why.

    Design notes (decisions made where the spec was silent)
    ---------------------------------------------------------------
    1. PROTEIN COVERAGE: ~74% of candidates have no protein measurement at
       all (mass spec wasn't run on every cell line). Scoring missing
       protein as 0 would penalise cell lines for not being *measured*,
       not for low expression -- a confound, not a signal. Instead, this
       implementation re-normalises weights per-row: a cell line with only
       RNA data has its `rna` weight inflated to absorb the `protein` and
       `consistency` weight it cannot otherwise earn, so it is judged
       fairly on the evidence that DOES exist instead of being dragged
       down by evidence that was never collected.
    2. DATA CONSISTENCY: defined here as 1 - |rank(rna) - rank(protein)|
       (percentile-rank agreement), only computable when both values
       exist. This is undefined behaviour in the original spec; rank
       agreement was chosen over raw correlation because RNA and protein
       are on different scales (log2 TPM+1 vs log2 MS intensity) and rank
       agreement is scale-invariant.
    3. TISSUE RELEVANCE: the spec doesn't say what tissue relevance is
       *relative to*. Two modes are implemented -- see `target_lineage`
       above.
    4. EXCLUSION PENALTY: the spec says "exclusion criteria" without
       specifying whether exclusion should subtract from score or hard-
       filter the cell line. This implementation hard-excludes (consistent
       with `build_gene_profile`'s `excluded` flag) rather than applying a
       soft penalty, because a partial penalty risks an excluded cell line
       still ranking highly enough to be recommended -- which defeats the
       purpose of specifying an exclusion gene in the first place.
    """
    w = {**DEFAULT_WEIGHTS, **(weights or {})}
    df = profile.copy()

    has_rna = df["rna_expression"].notna()
    has_protein = df["protein_expression"].notna()
    has_both = has_rna & has_protein

    # --- 1. RNA score: straightforward min-max normalisation -----------
    df["rna_score"] = _minmax(df["rna_expression"]).fillna(0.0)

    # --- 2. Protein score: min-max normalisation, 0 where absent --------
    df["protein_score"] = _minmax(df["protein_expression"]).fillna(0.0)

    # --- 3. Consistency score: percentile-rank agreement -----------------
    rna_pct = df["rna_expression"].rank(pct=True)
    prot_pct = df["protein_expression"].rank(pct=True)
    consistency_raw = 1.0 - (rna_pct - prot_pct).abs()
    df["consistency_score"] = consistency_raw.where(has_both, other=np.nan).fillna(0.0)

    # --- 4. Tissue relevance ---------------------------------------------
    if target_lineage is not None:
        df["tissue_score"] = (df["lineage"].str.lower() == target_lineage.lower()).astype(float)
    else:
        # Relative enrichment: share of this lineage in the top RNA quartile
        # vs. its share in the whole candidate pool.
        top_quartile = df["rna_expression"] >= df["rna_expression"].quantile(0.75)
        overall_share = df["lineage"].value_counts(normalize=True)
        top_share = df.loc[top_quartile, "lineage"].value_counts(normalize=True)
        enrichment = (top_share / overall_share).reindex(df["lineage"].unique()).fillna(0.0)
        enrichment = (enrichment / enrichment.max()).clip(upper=1.0)  # normalise to [0,1]
        df["tissue_score"] = df["lineage"].map(enrichment).fillna(0.0)

    # --- 5. Completeness score: how many of the 2 core evidence types exist
    df["completeness_score"] = (has_rna.astype(int) + has_protein.astype(int)) / 2.0

    # --- Combine with per-row weight renormalisation for missing protein -
    # Where protein is missing, redistribute its weight + the consistency
    # weight (which is also necessarily 0 without protein) onto rna+tissue,
    # proportionally, so a protein-less row is judged on what it DOES have.
    base_weights = pd.DataFrame({
        "rna": w["rna"],
        "protein": w["protein"],
        "consistency": w["consistency"],
        "tissue": w["tissue"],
        "completeness": w["completeness"],
    }, index=df.index)

    redistribute = base_weights["protein"] + base_weights["consistency"]
    base_weights.loc[~has_protein, "protein"] = 0.0
    base_weights.loc[~has_protein, "consistency"] = 0.0
    # Push the freed-up weight onto rna and tissue, split proportionally
    # to their original share of (rna + tissue).
    rna_tissue_total = base_weights["rna"] + base_weights["tissue"]
    rna_share = base_weights["rna"] / rna_tissue_total
    tissue_share = base_weights["tissue"] / rna_tissue_total
    base_weights.loc[~has_protein, "rna"] += redistribute[~has_protein] * rna_share[~has_protein]
    base_weights.loc[~has_protein, "tissue"] += redistribute[~has_protein] * tissue_share[~has_protein]

    df["final_score"] = (
        base_weights["rna"] * df["rna_score"]
        + base_weights["protein"] * df["protein_score"]
        + base_weights["consistency"] * df["consistency_score"]
        + base_weights["tissue"] * df["tissue_score"]
        + base_weights["completeness"] * df["completeness_score"]
    )

    # Hard-exclude (see Design note 4): excluded rows sink to the bottom
    # but are NOT dropped, so the caller can audit what was excluded.
    if "excluded" in df.columns:
        df.loc[df["excluded"], "final_score"] = -1.0

    return df.sort_values("final_score", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------- #
# Self-test when run directly
# ---------------------------------------------------------------------- #
if __name__ == "__main__":
    import sys
    sys.path.insert(0, ".")
    from data_loader import CellLineDataLoader

    data_dir = sys.argv[1] if len(sys.argv) > 1 else "data"
    loader = CellLineDataLoader(data_dir)

    print("Building EGFR profile...")
    profile = loader.build_gene_profile("EGFR", exclude_genes=["KRAS"])

    print("Scoring...")
    scored = score_gene_profile(profile)

    cols = ["DepMap_ID", "cell_line_name", "lineage", "rna_expression",
            "protein_expression", "final_score", "excluded"]
    print("\nTop 15 by final_score:")
    print(scored[cols].head(15).to_string(index=False))

    print(f"\nTotal candidates: {len(scored)}")
    print(f"Excluded (KRAS mut/fusion): {(scored['final_score'] == -1.0).sum()}")
    print(f"Protein-less candidates re-weighted: {profile['protein_expression'].isna().sum()}")
