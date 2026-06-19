"""
data_loader.py
================
Unified data loader for the CellLineFinder project.

Resolves the 7-ID-system problem documented in the project Data Dictionary
by converging everything onto DepMap_ID (ACH-xxxxxx) as the universal key.

BATCH 1 sources (core query loop: expression + exclusion criteria):
    File 9  - DepMap_sample_info.csv      -> master cell line table
    File 8  - DepMap_OmicsProfiles.csv     -> ProfileID / SequencingID -> ACH-ID bridge
    File 2  - DepMap_OmicsExpression...csv -> mRNA expression (log2 TPM+1)
    File 4  - Harmonized_MS_CCLE_Gygi.csv  -> protein expression (log2 intensity)
    File 6  - OmicsSomaticMutationsProfile.csv -> somatic mutations (exclusion criteria)
    File 5  - OmicsFusionFilteredSupplementary.csv -> gene fusions (exclusion criteria)

BATCH 2 sources (additional omics layers):
    File 13 - CCLE_miRNA_20181103.gct      -> miRNA expression (RPM)
              Columns are CCLE_Name (e.g. "DMS53_LUNG"), resolved to ACH-ID
              via file 9's CCLE_Name column. ~4 of 954 columns have no exact
              match in file 9 (see get_mirna_expression docstring).
    File 12 - CCLE_metabolomics_20190502.csv -> metabolite levels (log2 intensity)
              Already ships with DepMap_ID directly -- no bridging needed.
    File 14 - OmicsGlobalSignatures.csv    -> genome-wide signatures (MSI, CIN,
              ploidy, etc). Has duplicate ModelID rows (repeat sequencing);
              filtered to IsDefaultEntryForModel == "Yes" for a 1-row-per-
              cell-line table.

Usage
-----
    from data_loader import CellLineDataLoader

    loader = CellLineDataLoader("data-2")
    rna    = loader.get_rna_expression("EGFR")
    prot   = loader.get_protein_expression("EGFR")
    muts   = loader.get_mutations("EGFR")
    fus    = loader.get_fusions("EGFR")
    profile = loader.build_gene_profile("EGFR")   # merges Batch 1 sources

    mirna  = loader.get_mirna_expression("hsa-miR-21")
    metab  = loader.get_metabolite_level("glutamine")
    sigs   = loader.get_global_signatures()        # whole-genome features, all cell lines
"""

from __future__ import annotations

import re
from pathlib import Path
from functools import lru_cache

import pandas as pd
import numpy as np


class CellLineDataLoader:
    """
    Lazily loads and harmonises the Batch 1 + Batch 2 CellLineFinder data
    sources.

    All public methods return DataFrames indexed by / containing
    `DepMap_ID` (ACH-xxxxxx) so results from different omics layers
    can be merged directly.

    Parameters
    ----------
    data_dir : str | Path
        Path to the root data folder (e.g. "data-2"), containing the
        subfolders: "gene expression", "gene properties", "nomenclature",
        "non gene expression".
    """

    # Map logical file IDs to (subfolder, filename, separator)
    _FILE_MAP = {
        9:  ("nomenclature",        "9_DepMap_sample_info.csv", ","),
        8:  ("nomenclature",        "8_DepMap_OmicsProfiles.csv", ","),
        2:  ("gene expression",     "2_DepMap_OmicsExpressionAllGenesTPMLogp1Profile.csv", ","),
        4:  ("gene expression",     "4_Harmonized_MS_CCLE_Gygi_subsetted.csv", ","),
        6:  ("gene properties",     "6_OmicsSomaticMutationsProfile.csv", ","),
        5:  ("gene properties",     "5_OmicsFusionFilteredSupplementary.csv", ","),
        13: ("non gene expression", "13_CCLE_miRNA_20181103.gct", "\t"),
        12: ("non gene expression", "12_CCLE_metabolomics_20190502.csv", ","),
        14: ("non gene expression", "14_OmicsGlobalSignatures.csv", ","),
    }

    def __init__(self, data_dir: str | Path):
        self.data_dir = Path(data_dir)
        if not self.data_dir.exists():
            raise FileNotFoundError(f"Data directory not found: {self.data_dir}")

        # Cached frames (populated lazily on first access)
        self._sample_info: pd.DataFrame | None = None
        self._profiles: pd.DataFrame | None = None
        self._rna: pd.DataFrame | None = None       # large -> loaded once, kept wide
        self._protein: pd.DataFrame | None = None
        self._mutations: pd.DataFrame | None = None
        self._fusions: pd.DataFrame | None = None
        self._mirna: pd.DataFrame | None = None
        self._metabolomics: pd.DataFrame | None = None
        self._global_signatures: pd.DataFrame | None = None
        self._mirna_unmatched_cols: list[str] | None = None  # populated on first mirna load

    # ------------------------------------------------------------------ #
    # Internal path / IO helpers
    # ------------------------------------------------------------------ #

    def _path(self, file_id: int) -> Path:
        subfolder, filename, _ = self._FILE_MAP[file_id]
        p = self.data_dir / subfolder / filename
        if not p.exists():
            raise FileNotFoundError(
                f"Expected file {file_id} at: {p}\n"
                f"Check that '{subfolder}/' exists under {self.data_dir} "
                f"and contains '{filename}'."
            )
        return p

    # ------------------------------------------------------------------ #
    # File 9: master sample info  (the hub of the ID system)
    # ------------------------------------------------------------------ #

    @property
    def sample_info(self) -> pd.DataFrame:
        """Master cell line table. DepMap_ID is a regular column (not the index) so it merges cleanly with other tables."""
        if self._sample_info is None:
            df = pd.read_csv(self._path(9))
            self._sample_info = df
        return self._sample_info

    # ------------------------------------------------------------------ #
    # File 8: ProfileID / SequencingID -> ModelID (ACH-ID) bridge
    # ------------------------------------------------------------------ #

    @property
    def profiles(self) -> pd.DataFrame:
        """Bridge table: ProfileID (PR-xxxxxx) <-> ModelID (ACH-xxxxxx)."""
        if self._profiles is None:
            df = pd.read_csv(self._path(8))
            self._profiles = df
        return self._profiles

    @lru_cache(maxsize=1)
    def _profile_to_ach(self) -> dict:
        """ProfileID -> ACH-ID lookup dict (for resolving file 2's row index)."""
        df = self.profiles
        return dict(zip(df["ProfileID"], df["ModelID"]))

    # ------------------------------------------------------------------ #
    # File 2: mRNA expression  (wide matrix, ProfileID rows)
    # ------------------------------------------------------------------ #

    def _load_rna_matrix(self) -> pd.DataFrame:
        """
        Loads the full RNA expression matrix once and re-indexes rows
        from ProfileID (PR-xxxxxx) to DepMap_ID (ACH-xxxxxx) via file 8.

        NOTE: this matrix is ~1,495 rows x 53,961 cols (~600MB in memory
        as float64). Loaded once and cached; avoid calling repeatedly
        in a loop over many genes -- query columns from the cached frame
        instead.
        """
        if self._rna is None:
            df = pd.read_csv(self._path(2), index_col=0)
            mapping = self._profile_to_ach()
            ach_ids = df.index.map(mapping)

            n_unmapped = ach_ids.isna().sum()
            if n_unmapped:
                print(
                    f"[data_loader] Warning: {n_unmapped} of {len(df)} "
                    f"RNA profile rows could not be mapped to an ACH-ID "
                    f"(missing from file 8). These rows are dropped."
                )
            df = df.loc[~ach_ids.isna()]
            df.index = ach_ids[~ach_ids.isna()]
            df.index.name = "DepMap_ID"
            self._rna = df
        return self._rna

    def get_rna_expression(self, gene: str) -> pd.DataFrame:
        """
        Returns mRNA expression (log2 TPM+1) for `gene` across all cell
        lines, indexed by DepMap_ID.

        `gene` may be a bare gene symbol ("EGFR") or match the column
        format used in file 2 ("EGFR (ENSG00000146648)") -- partial,
        case-insensitive matching on the symbol is used.

        Returns
        -------
        DataFrame with columns: ['DepMap_ID', 'rna_expression']
        """
        df = self._load_rna_matrix()
        col = self._match_column(df.columns, gene)
        out = df[[col]].reset_index()
        out.columns = ["DepMap_ID", "rna_expression"]
        return out

    # ------------------------------------------------------------------ #
    # File 4: protein expression (wide matrix, ACH-ID rows already)
    # ------------------------------------------------------------------ #

    def _load_protein_matrix(self) -> pd.DataFrame:
        if self._protein is None:
            df = pd.read_csv(self._path(4), index_col=0)
            df.index.name = "DepMap_ID"
            self._protein = df
        return self._protein

    def get_protein_expression(self, gene: str) -> pd.DataFrame:
        """
        Returns protein expression (log2 normalised MS intensity) for
        `gene`, indexed by DepMap_ID. NaN means the protein was not
        detected by mass spec in that cell line (not necessarily absent).

        Returns
        -------
        DataFrame with columns: ['DepMap_ID', 'protein_expression']
        """
        df = self._load_protein_matrix()
        col = self._match_column(df.columns, gene)
        out = df[[col]].reset_index()
        out.columns = ["DepMap_ID", "protein_expression"]
        return out

    # ------------------------------------------------------------------ #
    # File 6: somatic mutations (long/event table)
    # ------------------------------------------------------------------ #

    @property
    def mutations(self) -> pd.DataFrame:
        """
        Full somatic mutation table (~1M rows). Loaded once, cached.

        NOTE: file 6 identifies samples by ProfileID (PR-xxxxxx), NOT
        ModelID -- unlike files 5 and 14. We resolve ProfileID -> ACH-ID
        via file 8 here, same as we do for the file 2 RNA matrix.
        """
        if self._mutations is None:
            usecols = [
                "ProfileID", "HugoSymbol", "EnsemblGeneID", "ProteinChange",
                "VariantType", "VepImpact", "AF", "DP", "LofGeneName",
            ]
            df = pd.read_csv(self._path(6), usecols=lambda c: c in usecols,
                              low_memory=False)

            mapping = self._profile_to_ach()
            ach_ids = df["ProfileID"].map(mapping)
            n_unmapped = ach_ids.isna().sum()
            if n_unmapped:
                print(
                    f"[data_loader] Warning: {n_unmapped} of {len(df)} "
                    f"mutation records could not be mapped to an ACH-ID "
                    f"(ProfileID missing from file 8). These rows are dropped."
                )
            df = df.loc[~ach_ids.isna()].copy()
            df["DepMap_ID"] = ach_ids[~ach_ids.isna()]
            self._mutations = df
        return self._mutations

    def get_mutations(self, gene: str) -> pd.DataFrame:
        """
        Returns all somatic mutation records for `gene`, one row per
        mutation event per cell line (a cell line may have multiple
        mutations in the same gene).

        Returns
        -------
        DataFrame with columns:
            ['DepMap_ID', 'ProteinChange', 'VariantType', 'VepImpact', 'AF', 'DP', 'is_lof']
        """
        df = self.mutations
        mask = df["HugoSymbol"].str.upper() == gene.upper()
        out = df.loc[mask, [
            "DepMap_ID", "ProteinChange", "VariantType", "VepImpact", "AF", "DP", "LofGeneName"
        ]].copy()
        out["is_lof"] = out["LofGeneName"].notna()
        out = out.drop(columns="LofGeneName")
        return out.reset_index(drop=True)

    # ------------------------------------------------------------------ #
    # File 5: gene fusions (long/event table)
    # ------------------------------------------------------------------ #

    @property
    def fusions(self) -> pd.DataFrame:
        """Full gene fusion event table (~184k rows). Loaded once, cached."""
        if self._fusions is None:
            usecols = [
                "ModelID", "CanonicalFusionName", "FFPM", "confidence", "reading_frame",
            ]
            df = pd.read_csv(self._path(5), usecols=lambda c: c in usecols)
            self._fusions = df
        return self._fusions

    def get_fusions(self, gene: str) -> pd.DataFrame:
        """
        Returns all fusion events involving `gene` (as either the 5' or
        3' partner), one row per fusion event per cell line.

        Returns
        -------
        DataFrame with columns:
            ['DepMap_ID', 'CanonicalFusionName', 'FFPM', 'confidence', 'reading_frame']
        """
        df = self.fusions
        # CanonicalFusionName format: "GENE1--GENE2"
        pattern = rf"(?:^|--){re.escape(gene.upper())}(?:--|$)"
        mask = df["CanonicalFusionName"].str.upper().str.contains(pattern, regex=True, na=False)
        out = df.loc[mask].rename(columns={"ModelID": "DepMap_ID"})
        return out.reset_index(drop=True)

    # ------------------------------------------------------------------ #
    # File 13: miRNA expression (GCT format, CCLE_Name columns)
    # ------------------------------------------------------------------ #

    def _load_mirna_matrix(self) -> pd.DataFrame:
        """
        Loads the miRNA matrix (GCT v1.2 format) once and re-indexes
        columns from CCLE_Name (e.g. "DMS53_LUNG") to DepMap_ID
        (ACH-xxxxxx) via file 9's CCLE_Name column.

        A small number of file 13's 954 CCLE_Name columns have no exact
        match in file 9 (renamed/retired cell lines). These columns are
        dropped; the dropped names are stored in
        `self._mirna_unmatched_cols` and printed as a warning so you can
        inspect them (e.g. to resolve manually via file 7 Cellosaurus
        synonyms) before deciding whether they matter for your analysis.
        """
        if self._mirna is None:
            df = pd.read_csv(self._path(13), sep="\t", skiprows=2)
            df = df.set_index("Description")  # standard hsa-miR-xxx name
            df = df.drop(columns=["Name"])     # drop the internal nmiR code

            ccle_to_ach = dict(zip(self.sample_info["CCLE_Name"], self.sample_info["DepMap_ID"]))
            ach_cols = df.columns.map(ccle_to_ach)

            unmatched_mask = ach_cols.isna()
            n_unmatched = unmatched_mask.sum()
            if n_unmatched:
                self._mirna_unmatched_cols = df.columns[unmatched_mask].tolist()
                print(
                    f"[data_loader] Warning: {n_unmatched} of {len(df.columns)} "
                    f"miRNA columns (CCLE_Name) had no exact match in file 9's "
                    f"CCLE_Name column and were dropped:"
                )
                for name in self._mirna_unmatched_cols:
                    print(f"    - {name}")
                print(
                    "    These may be renamed/retired cell lines. Cross-check "
                    "manually via file 7 (cellosaurus.csv) synonyms if you need them."
                )
            else:
                self._mirna_unmatched_cols = []

            df = df.loc[:, ~unmatched_mask]
            df.columns = ach_cols[~unmatched_mask]
            df.columns.name = "DepMap_ID"
            self._mirna = df.T  # rows = cell lines (DepMap_ID), cols = miRNAs
            self._mirna.index.name = "DepMap_ID"
        return self._mirna

    def get_mirna_expression(self, mirna: str) -> pd.DataFrame:
        """
        Returns miRNA expression (RPM) for `mirna` across all cell lines
        with a resolvable DepMap_ID, indexed by DepMap_ID.

        Parameters
        ----------
        mirna : str
            Standard miRNA name, e.g. "hsa-miR-21" or "hsa-let-7a".
            Case-insensitive exact match against file 13's Description column.

        Returns
        -------
        DataFrame with columns: ['DepMap_ID', 'mirna_expression']
        """
        df = self._load_mirna_matrix()
        matches = [c for c in df.columns if c.upper() == mirna.upper()]
        if not matches:
            raise KeyError(
                f"miRNA '{mirna}' not found. Use the standard hsa-miR-xxx / "
                f"hsa-let-7x naming as it appears in file 13's Description column."
            )
        col = matches[0]
        out = df[[col]].reset_index()
        out.columns = ["DepMap_ID", "mirna_expression"]
        return out

    # ------------------------------------------------------------------ #
    # File 12: metabolomics (wide matrix, ships with DepMap_ID directly)
    # ------------------------------------------------------------------ #

    def _load_metabolomics_matrix(self) -> pd.DataFrame:
        if self._metabolomics is None:
            df = pd.read_csv(self._path(12))
            df = df.drop(columns=["CCLE_ID"])  # DepMap_ID is the cleaner key, already present
            df = df.dropna(subset=["DepMap_ID"]).set_index("DepMap_ID")
            self._metabolomics = df
        return self._metabolomics

    def get_metabolite_level(self, metabolite: str) -> pd.DataFrame:
        """
        Returns metabolite level (log2-normalised intensity) for
        `metabolite` across all 928 profiled cell lines, indexed by
        DepMap_ID.

        Parameters
        ----------
        metabolite : str
            Metabolite name as it appears in file 12's column headers,
            e.g. "glutamine", "2-aminoadipate". Case-insensitive exact match.

        Returns
        -------
        DataFrame with columns: ['DepMap_ID', 'metabolite_level']
        """
        df = self._load_metabolomics_matrix()
        matches = [c for c in df.columns if c.upper() == metabolite.upper()]
        if not matches:
            raise KeyError(
                f"Metabolite '{metabolite}' not found. Check spelling against "
                f"file 12's column headers (225 metabolites total)."
            )
        col = matches[0]
        out = df[[col]].reset_index()
        out.columns = ["DepMap_ID", "metabolite_level"]
        return out

    def list_metabolites(self) -> list[str]:
        """Returns all 225 available metabolite names (file 12 columns)."""
        return self._load_metabolomics_matrix().columns.tolist()

    # ------------------------------------------------------------------ #
    # File 14: global genomic signatures (MSI, CIN, ploidy, etc.)
    # ------------------------------------------------------------------ #

    def get_global_signatures(self) -> pd.DataFrame:
        """
        Returns one row per cell line of whole-genome-derived signatures:
        MSIScore, LoHFraction, WGD, CIN, Ploidy, Aneuploidy.

        File 14 has multiple rows per ModelID (repeat sequencing
        profiles); this method filters to IsDefaultEntryForModel == "Yes"
        so each cell line appears exactly once.

        Returns
        -------
        DataFrame with columns:
            ['DepMap_ID', 'MSIScore', 'LoHFraction', 'WGD', 'CIN', 'Ploidy', 'Aneuploidy']
        """
        if self._global_signatures is None:
            df = pd.read_csv(self._path(14))
            df = df.loc[df["IsDefaultEntryForModel"] == "Yes"].copy()

            n_dupes = df["ModelID"].duplicated().sum()
            if n_dupes:
                print(
                    f"[data_loader] Warning: {n_dupes} ModelIDs still have more "
                    f"than one 'default' entry after filtering -- keeping the first."
                )
                df = df.drop_duplicates(subset="ModelID", keep="first")

            df = df.rename(columns={"ModelID": "DepMap_ID"})
            cols = ["DepMap_ID", "MSIScore", "LoHFraction", "WGD", "CIN", "Ploidy", "Aneuploidy"]
            self._global_signatures = df[cols].reset_index(drop=True)
        return self._global_signatures

    def get_global_signature_for_gene_profile(self, depmap_ids: list[str]) -> pd.DataFrame:
        """
        Convenience filter: returns global signatures only for the given
        list of DepMap_IDs (e.g. to attach to a build_gene_profile() result
        manually -- merging is intentionally left to the caller for now).
        """
        sigs = self.get_global_signatures()
        return sigs.loc[sigs["DepMap_ID"].isin(depmap_ids)].reset_index(drop=True)

    # ------------------------------------------------------------------ #
    # Combined gene profile (the actual CellLineFinder query)
    # ------------------------------------------------------------------ #

    def build_gene_profile(self, gene: str, exclude_genes: list[str] | None = None) -> pd.DataFrame:
        """
        Builds the core CellLineFinder query result for a single target
        gene: merges RNA + protein expression for every cell line, and
        flags cell lines that should be excluded based on `exclude_genes`
        (mutations or fusions in those genes).

        Parameters
        ----------
        gene : str
            Target gene symbol, e.g. "EGFR".
        exclude_genes : list[str], optional
            Genes whose mutation/fusion presence should flag a cell line
            for exclusion (e.g. resistance mutations).

        Returns
        -------
        DataFrame, one row per cell line, columns:
            ['DepMap_ID', 'cell_line_name', 'lineage', 'primary_disease',
             'rna_expression', 'protein_expression',
             'has_target_mutation', 'has_target_fusion', 'excluded']
        """
        rna = self.get_rna_expression(gene)
        prot = self.get_protein_expression(gene)
        muts = self.get_mutations(gene)
        fus = self.get_fusions(gene)

        # Start from the master sample table so every known cell line
        # is represented even if it has no expression data point.
        base = self.sample_info[["DepMap_ID", "cell_line_name", "lineage", "primary_disease"]].copy()

        out = base.merge(rna, on="DepMap_ID", how="left")
        out = out.merge(prot, on="DepMap_ID", how="left")

        mut_ids = set(muts["DepMap_ID"])
        fus_ids = set(fus["DepMap_ID"])
        out["has_target_mutation"] = out["DepMap_ID"].isin(mut_ids)
        out["has_target_fusion"] = out["DepMap_ID"].isin(fus_ids)

        out["excluded"] = False
        if exclude_genes:
            excluded_ids = set()
            for ex_gene in exclude_genes:
                excluded_ids |= set(self.get_mutations(ex_gene)["DepMap_ID"])
                excluded_ids |= set(self.get_fusions(ex_gene)["DepMap_ID"])
            out["excluded"] = out["DepMap_ID"].isin(excluded_ids)

        # Drop cell lines with no expression evidence at all (no RNA, no protein)
        has_evidence = out["rna_expression"].notna() | out["protein_expression"].notna()
        out = out.loc[has_evidence].reset_index(drop=True)

        return out.sort_values("rna_expression", ascending=False, na_position="last").reset_index(drop=True)

    # ------------------------------------------------------------------ #
    # Shared helper: flexible column matching for wide matrices
    # ------------------------------------------------------------------ #

    @staticmethod
    def _match_column(columns: pd.Index, gene: str) -> str:
        """
        Finds the column in a wide expression matrix corresponding to
        `gene`. Column headers in files 2 and 4 embed the gene symbol
        inside a longer string, e.g.:
            file 2: "EGFR (ENSG00000146648)"
            file 4: "P00533 (EGFR)"

        Strategy: extract the bracketed/leading symbol token from each
        header and look for an EXACT (case-insensitive) match first,
        to avoid "EGFR" incorrectly matching "EGFR-AS1". Only if no
        exact match exists do we fall back to substring matching.
        """
        gene_upper = gene.upper()

        # file 2 style: "SYMBOL (ENSGxxxxx)" -> symbol is before " ("
        # file 4 style: "UNIPROTID (SYMBOL)" -> symbol is inside parens
        exact_matches = []
        for c in columns:
            c_upper = c.upper()
            if c_upper == gene_upper:
                exact_matches.append(c)
                continue
            # token before " (" 
            head = c_upper.split(" (")[0]
            # token inside parens, if any
            paren = re.search(r"\(([^)]+)\)", c_upper)
            paren_token = paren.group(1) if paren else None
            if head == gene_upper or paren_token == gene_upper:
                exact_matches.append(c)

        if exact_matches:
            if len(exact_matches) > 1:
                print(
                    f"[data_loader] Warning: multiple exact matches for '{gene}': "
                    f"{exact_matches}. Using the first."
                )
            return exact_matches[0]

        # Fallback: substring/word-boundary search
        pattern = re.compile(rf"(?<![A-Z0-9]){re.escape(gene_upper)}(?![A-Z0-9])")
        matches = [c for c in columns if pattern.search(c.upper())]

        if not matches:
            raise KeyError(
                f"Gene '{gene}' not found in this matrix. "
                f"Check spelling/symbol -- e.g. use the HGNC-approved symbol."
            )
        if len(matches) > 1:
            print(
                f"[data_loader] Warning: {len(matches)} columns partially matched '{gene}': "
                f"{matches[:5]}{'...' if len(matches) > 5 else ''}. Using the first match."
            )
        return matches[0]


# ---------------------------------------------------------------------- #
# Quick self-test when run directly (not on import)
# ---------------------------------------------------------------------- #
if __name__ == "__main__":
    import sys
    data_dir = sys.argv[1] if len(sys.argv) > 1 else "data-2"
    loader = CellLineDataLoader(data_dir)

    print("Loading sample_info...")
    print(loader.sample_info.shape)

    print("\nQuerying EGFR profile (Batch 1)...")
    profile = loader.build_gene_profile("EGFR")
    print(profile.head(10))
    print(f"\nTotal cell lines with EGFR evidence: {len(profile)}")

    print("\n--- Batch 2 sanity checks ---")

    print("\nmiRNA: hsa-miR-21...")
    mirna = loader.get_mirna_expression("hsa-miR-21")
    print(f"Rows: {len(mirna)}")
    print(mirna.sort_values("mirna_expression", ascending=False).head(5))

    print("\nMetabolomics: glutamine...")
    metab = loader.get_metabolite_level("glutamine")
    print(f"Rows: {len(metab)}")
    print(metab.sort_values("metabolite_level", ascending=False).head(5))

    print("\nGlobal signatures (all cell lines)...")
    sigs = loader.get_global_signatures()
    print(f"Rows: {len(sigs)} (should be ~1,955 unique cell lines)")
    print(sigs.head(5))
