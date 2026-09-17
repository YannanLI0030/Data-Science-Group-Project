"""
data_merger.py  (v3)
====================
CellLineSelector: unified multi-source integration layer (DepMap + HPA + GEO)

Problems exposed by v2 on the full datasets and addressed in v3
----------------------------------------------------------------
1. [BUG] Column matching was too rigid. The Cellosaurus column is named
   'Accession (CVCL_xxxx)', not 'Accession'. `_find_col` now matches normalized
   names by containment and accepts parenthetical suffixes.
2. [BUG] Iterating over 152,231 Cellosaurus rows with iterrows took several
   minutes. The implementation is now vectorized.
3. [CRITICAL] v2 melted the complete DepMap RNA matrix into 80 million rows.
   Adding the 24 million HPA rows caused an out-of-memory failure. v3 pushes
   the target-gene set down to each source, which loads only the required genes
   and reduces memory use from gigabytes to megabytes.
4. [BUG] GEO File 3 is a wide table (Gene plus one column per GSM), whereas v2
   assumed long format. v3 detects the format automatically.
5. [ADDED] Automatic detection of RRID/Cellosaurus columns in DepMap
   sample_info. When present, this provides a more direct CVCL-to-ACH bridge
   than Cellosaurus Cross-references.
6. [ADDED] Running this file directly performs a self-check
   (python data_merger.py <data_dir>).
7. [ADDED] v3.2 adds MutationSource and FusionSource adapters.

Mapping to the four implementation requirements
------------------------------------------------
1. Unified ID alignment across sources: CellLineIDResolver, with CVCL as hub
2. Unified structured multi-omics storage: UNIFIED_COLUMNS long table and
   build_gene_table wide table
3. Standardized missing-value handling: outer-join scaffold, has_<layer>
   masks, and redistribute_weights
4. Extensible data-source support: OmicsSource registry pattern
"""

from __future__ import annotations

__version__ = "3.2"

import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterable, Optional, Set

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# 0. Unified schema
# ---------------------------------------------------------------------------

UNIFIED_COLUMNS = [
    "DepMap_ID", "gene_symbol", "ensembl_id",
    "omics_layer", "value_raw", "value_unit",
    "source", "source_version",
]

OMICS_LAYERS = ("rna", "protein", "mutation", "fusion", "mirna", "metabolite")


# ---------------------------------------------------------------------------
# 1. General utility: flexible column-name matching
# ---------------------------------------------------------------------------

def _norm_colname(c: str) -> str:
    """Normalize a column name to lowercase alphanumeric characters."""
    return re.sub(r"[^a-z0-9]", "", str(c).lower())


def _find_col(df: pd.DataFrame, keywords: list[str], what: str = "column") -> str:
    """Find a normalized column name containing a keyword.

    Parenthetical suffixes, capitalization, and spaces are ignored. Exact
    matches are preferred, followed by prefix matches and then containment.
    """
    norm_map = {_norm_colname(c): c for c in df.columns}
    keys = [_norm_colname(k) for k in keywords]
    # 1) Exact match
    for k in keys:
        if k in norm_map:
            return norm_map[k]
    # 2) Prefix match
    for k in keys:
        for nc, orig in norm_map.items():
            if nc.startswith(k):
                return orig
    # 3) Containment match
    for k in keys:
        for nc, orig in norm_map.items():
            if k in nc:
                return orig
    raise KeyError(
        f"Could not find {what}: expected a column containing one of {keywords}.\n"
        f"Available columns: {list(df.columns)[:15]}"
    )


def _find_col_optional(df: pd.DataFrame, keywords: list[str]) -> Optional[str]:
    try:
        return _find_col(df, keywords)
    except KeyError:
        return None


# ---------------------------------------------------------------------------
# 2. CellLineIDResolver: vectorized ID alignment using CVCL as the hub
# ---------------------------------------------------------------------------

class CellLineIDResolver:
    """Resolve a cell-line identifier to a DepMap_ID.

    Accepted identifiers include CVCL, cell-line names, aliases, and CCLE_Name.
    Two lookup tables are constructed:
      * cvcl_to_ach  : {'CVCL_0023': 'ACH-000681', ...}
      * name_to_cvcl : {normalized name -> 'CVCL_0023', ...}, derived from
        Identifier and Synonyms

    CVCL-to-ACH mappings come from two sources, ordered by reliability:
      (A) the RRID or Cellosaurus column in DepMap sample_info, when present;
      (B) 'DepMap; ACH-xxxxxx' entries in Cellosaurus Cross-references.
    Both sources are used, with source A taking precedence.
    """

    def __init__(self, cellosaurus: pd.DataFrame,
                 sample_info: Optional[pd.DataFrame] = None,
                 verbose: bool = True):
        acc_col = _find_col(cellosaurus, ["accession", "cvcl"], "Cellosaurus CVCL column")
        id_col = _find_col(cellosaurus, ["identifier", "cellline", "name"],
                           "Cellosaurus cell-line name column")
        xref_col = _find_col_optional(cellosaurus, ["crossreferences", "crossreference", "xref"])
        syn_col = _find_col_optional(cellosaurus, ["synonyms", "synonym"])

        if verbose:
            print(f"[resolver] Cellosaurus columns: CVCL={acc_col!r}, "
                  f"name={id_col!r}, xref={xref_col!r}, syn={syn_col!r}")

        cvcl = cellosaurus[acc_col].astype("string").str.strip()
        valid = cvcl.str.startswith("CVCL_", na=False)

        # ---------- Path A: CVCL column in sample_info (most reliable) ----------
        self.cvcl_to_ach: dict[str, str] = {}
        n_from_sample_info = 0
        if sample_info is not None:
            si_cvcl_col = _find_col_optional(sample_info, ["rrid", "cellosaurus", "cvcl"])
            if si_cvcl_col is not None:
                si_cvcl = sample_info[si_cvcl_col].astype("string").str.strip()
                # RRIDs commonly appear as 'CVCL_0023' or 'RRID:CVCL_0023'.
                si_cvcl = si_cvcl.str.extract(r"(CVCL_[A-Za-z0-9]+)", expand=False)
                pairs = pd.DataFrame({"cvcl": si_cvcl,
                                      "ach": sample_info["DepMap_ID"].astype("string")}).dropna()
                self.cvcl_to_ach.update(dict(zip(pairs["cvcl"], pairs["ach"])))
                n_from_sample_info = len(self.cvcl_to_ach)
                if verbose:
                    print(f"[resolver] Path A: sample_info column {si_cvcl_col!r} provided "
                          f"{n_from_sample_info} CVCL-to-ACH mappings (highest authority)")
            elif verbose:
                print("[resolver] Path A: no RRID/Cellosaurus column found in sample_info; skipped")

        # ---------- Path B: Cellosaurus Cross-references ----------
        n_from_xref = 0
        if xref_col is not None:
            ach = (cellosaurus[xref_col].astype("string")
                   .str.extract(r"(?i)depmap\s*[;=:]?\s*(ACH-\d+)", expand=False))
            pairs = pd.DataFrame({"cvcl": cvcl.where(valid), "ach": ach}).dropna()
            before = len(self.cvcl_to_ach)
            for c, a in zip(pairs["cvcl"], pairs["ach"]):
                self.cvcl_to_ach.setdefault(c, a.upper())
            n_from_xref = len(pairs)
            if verbose:
                print(f"[resolver] Path B: Cellosaurus Cross-references provided "
                      f"{n_from_xref} CVCL-to-ACH mappings "
                      f"({len(self.cvcl_to_ach) - before} newly added)")

        if not self.cvcl_to_ach:
            print("[resolver] WARNING: no CVCL-to-ACH mappings were created.\n"
                  "    Check whether (a) sample_info has an RRID column and "
                  "(b) Cellosaurus Cross-references contains 'DepMap; ACH-xxx'.\n"
                  "    Run diagnose_cellosaurus() to inspect the xref format.")

        # ---------- Vectorized name-to-CVCL dictionary ----------
        names = cellosaurus.loc[valid, id_col].astype("string")
        cvcl_valid = cvcl[valid]
        self._name_to_cvcl: dict[str, str] = {}
        self._bulk_register(names, cvcl_valid)

        if syn_col is not None:
            syn = cellosaurus.loc[valid, syn_col].astype("string")
            # Split synonyms on ';' or '||', then explode to one alias per row.
            exploded = (syn.str.split(r"\s*(?:\|\||;)\s*", regex=True)
                           .explode())
            exploded_cvcl = cvcl_valid.reindex(exploded.index)
            self._bulk_register(exploded, exploded_cvcl)

        # ---------- Attach sample_info names to CVCL for GEO ----------
        if sample_info is not None:
            ach_to_cvcl = {a: c for c, a in self.cvcl_to_ach.items()}
            si_cvcl_back = sample_info["DepMap_ID"].astype("string").map(ach_to_cvcl)
            for col in ("cell_line_name", "stripped_cell_line_name", "CCLE_Name"):
                if col in sample_info.columns:
                    self._bulk_register(sample_info[col].astype("string"), si_cvcl_back)

        if verbose:
            print(f"[resolver] Complete: {len(self.cvcl_to_ach)} CVCL-to-ACH mappings, "
                  f"{len(self._name_to_cvcl)} name-to-CVCL mappings")

    # ---- Vectorized batch name registration; no iterrows ----
    def _bulk_register(self, names: pd.Series, cvcls: pd.Series) -> None:
        keys = names.map(self._norm)
        df = pd.DataFrame({"k": keys, "v": cvcls}).dropna()
        df = df[df["k"] != ""]
        for k, v in zip(df["k"], df["v"]):
            self._name_to_cvcl.setdefault(k, v)

    @staticmethod
    def _norm(name) -> Optional[str]:
        if not isinstance(name, str) or not name.strip():
            return None
        return re.sub(r"[^A-Z0-9]", "", name.upper())

    # ---- Public API ----
    def resolve_cvcl_to_ach(self, cvcl) -> Optional[str]:
        """HPA primary path: File 11 already provides CVCL identifiers."""
        if isinstance(cvcl, str) and cvcl.startswith("CVCL_"):
            return self.cvcl_to_ach.get(cvcl)
        return None

    def resolve_name(self, name) -> Optional[str]:
        """GEO fallback path: name to CVCL to ACH."""
        if isinstance(name, str) and name.startswith("ACH-"):
            return name
        cvcl = self._name_to_cvcl.get(self._norm(name))
        return self.cvcl_to_ach.get(cvcl) if cvcl else None

    def resolve_cvcl_series(self, s: pd.Series) -> pd.Series:
        return s.map(self.resolve_cvcl_to_ach)

    def resolve_name_series(self, s: pd.Series) -> pd.Series:
        return s.map(self.resolve_name)

    def stats(self) -> dict:
        return {"cvcl_to_ach": len(self.cvcl_to_ach),
                "name_to_cvcl": len(self._name_to_cvcl)}


def diagnose_cellosaurus(cellosaurus: pd.DataFrame, n: int = 10) -> None:
    """Print resource names in Cellosaurus Cross-references and check for DepMap."""
    xref_col = _find_col_optional(cellosaurus, ["crossreferences", "xref"])
    if xref_col is None:
        print("No Cross-references column")
        return
    s = cellosaurus[xref_col].dropna().astype(str)
    print(f"Non-empty Cross-references rows: {len(s):,} / {len(cellosaurus):,}")
    # Extract every 'Resource;' token and count its frequency.
    resources = (s.str.findall(r"([A-Za-z0-9_]+)\s*[;=]")
                  .explode().dropna())
    print(f"\nTop {n} cross-reference resources:")
    print(resources.value_counts().head(n).to_string())
    has_depmap = resources.str.lower().eq("depmap").any()
    print(f"\nContains DepMap references: {has_depmap}")
    if has_depmap:
        sample = s[s.str.contains("depmap", case=False, na=False)].head(3)
        print("\nExample DepMap references:")
        for v in sample:
            print("  ", v[:150])


# ---------------------------------------------------------------------------
# 3. GeneIDResolver
# ---------------------------------------------------------------------------

class GeneIDResolver:
    def __init__(self, ensembl_map: Optional[pd.DataFrame] = None):
        self._ensg_to_sym: dict[str, str] = {}
        self._sym_to_ensg: dict[str, str] = {}
        if ensembl_map is not None:
            self._bulk(ensembl_map["ensembl_id"], ensembl_map["gene_symbol"])

    def _bulk(self, ensg: pd.Series, sym: pd.Series) -> None:
        df = pd.DataFrame({"e": ensg.astype("string").str.split(".").str[0],
                           "s": sym.astype("string")}).dropna().drop_duplicates()
        for e, s in zip(df["e"], df["s"]):
            self._ensg_to_sym.setdefault(e, s)
            self._sym_to_ensg.setdefault(s, e)

    def build_from_hpa_file(self, path, ensg_kw=("gene",), sym_kw=("genename",),
                            verbose: bool = True) -> "GeneIDResolver":
        """Build the symbol-to-ENSG dictionary from HPA File 1.

        HPA File 1 is a long table with approximately 24 million rows. Only two
        columns are read and deduplicated. Using nrows=100000 would cover only
        about 83 genes.
        """
        head = pd.read_csv(path, sep="\t", nrows=5)
        ensg_col = _find_col(head, list(ensg_kw), "HPA Gene (ENSG) column")
        sym_col = _find_col(head, list(sym_kw), "HPA Gene name column")
        df = pd.read_csv(path, sep="\t", usecols=[ensg_col, sym_col],
                         dtype="string").drop_duplicates()
        self._bulk(df[ensg_col], df[sym_col])
        if verbose:
            print(f"[gene_res] Loaded {len(self._sym_to_ensg):,} gene symbols from HPA")
        return self

    def to_symbol(self, ensg):
        return self._ensg_to_sym.get(ensg.split(".")[0]) if isinstance(ensg, str) else None

    def to_ensembl(self, symbol):
        return self._sym_to_ensg.get(symbol) if isinstance(symbol, str) else None


# ---------------------------------------------------------------------------
# 4. Abstract data-source interface (v3: supports pushed-down gene filtering)
# ---------------------------------------------------------------------------

class OmicsSource(ABC):
    source: str
    layer: str
    version: str

    @abstractmethod
    def load_long(self, resolver: "CellLineIDResolver",
                  gene_resolver: "GeneIDResolver",
                  genes: Optional[Set[str]] = None) -> pd.DataFrame:
        """Return a long table with UNIFIED_COLUMNS.

        If genes is provided, the source only needs to return data for those
        genes. Implementations should filter while reading instead of loading
        the complete dataset first.
        """
        ...

    def _pack(self, df: pd.DataFrame) -> pd.DataFrame:
        out = df.reindex(columns=UNIFIED_COLUMNS)
        out["source"] = self.source
        out["source_version"] = self.version
        out["omics_layer"] = self.layer
        n_before = len(out)
        # Drop only rows with missing keys. Retain NaN measurements to represent
        # undetected proteomics values.
        out = out[out["DepMap_ID"].notna() & out["gene_symbol"].notna()]
        dropped = n_before - len(out)
        if dropped:
            print(f"  [{self.source}] Skipped {dropped:,} rows with unresolved IDs "
                  f"({dropped / max(n_before, 1):.1%})")
        return out.reset_index(drop=True)


# ---------------------------------------------------------------------------
# 5. Source adapters
# ---------------------------------------------------------------------------

def _split_symbol_ensg(columns: pd.Index) -> pd.DataFrame:
    """Split a header such as 'EGFR (ENSG00000146648)' into symbol and Ensembl ID."""
    meta = pd.DataFrame({"orig": list(columns)})
    ext = meta["orig"].str.extract(r"^\s*(?P<sym>[^(]+?)\s*\((?P<ensg>[^)]+)\)\s*$")
    meta["gene_symbol"] = ext["sym"].fillna(meta["orig"].str.strip())
    meta["ensembl_id"] = ext["ensg"]
    return meta


class DepMapRNASource(OmicsSource):
    """DepMap RNA using the data_loader wide matrix.

    Rows are ACH identifiers and columns use 'SYMBOL (ENSG)'. v3 melts only
    the requested gene columns, avoiding an 80-million-row full melt.
    """
    layer, source = "rna", "depmap_rna"

    def __init__(self, loader, version: str = "DepMap-24Q2"):
        self.loader, self.version = loader, version
        self._colmeta: Optional[pd.DataFrame] = None

    def _meta(self, wide) -> pd.DataFrame:
        if self._colmeta is None:
            self._colmeta = _split_symbol_ensg(wide.columns)
        return self._colmeta

    def load_long(self, resolver, gene_resolver, genes=None):
        wide = self.loader._load_rna_matrix()
        meta = self._meta(wide)
        if genes is not None:
            meta = meta[meta["gene_symbol"].isin(genes)]
            if meta.empty:
                return self._pack(pd.DataFrame(columns=UNIFIED_COLUMNS))
        sub = wide[meta["orig"].tolist()]
        long = (sub.reset_index()
                   .melt(id_vars="DepMap_ID", var_name="orig", value_name="value_raw"))
        long = long.merge(meta[["orig", "gene_symbol", "ensembl_id"]], on="orig", how="left")
        long = long.drop(columns="orig")
        long = long[long["value_raw"].notna()]
        long["value_unit"] = "log2(TPM+1)"
        return self._pack(long)


class MatrixProteinSource(OmicsSource):
    """CCLE-Gygi protein data with ACH rows and 'UNIPROT (SYMBOL)' columns.

    NaN values are retained.
    """
    layer, source = "protein", "ccle_gygi_protein"

    def __init__(self, loader, version: str = "CCLE-Gygi-2020"):
        self.loader, self.version = loader, version
        self._colmeta: Optional[pd.DataFrame] = None

    def _meta(self, wide) -> pd.DataFrame:
        if self._colmeta is None:
            meta = pd.DataFrame({"orig": list(wide.columns)})
            # 'P00533 (EGFR)': the symbol is inside parentheses.
            meta["gene_symbol"] = meta["orig"].str.extract(r"\(([^)]+)\)")
            meta["gene_symbol"] = meta["gene_symbol"].fillna(meta["orig"].str.strip())
            self._colmeta = meta
        return self._colmeta

    def load_long(self, resolver, gene_resolver, genes=None):
        wide = self.loader._load_protein_matrix()
        meta = self._meta(wide)
        if genes is not None:
            meta = meta[meta["gene_symbol"].isin(genes)]
            if meta.empty:
                return self._pack(pd.DataFrame(columns=UNIFIED_COLUMNS))
        sub = wide[meta["orig"].tolist()]
        long = (sub.reset_index()
                   .melt(id_vars="DepMap_ID", var_name="orig", value_name="value_raw"))
        long = long.merge(meta[["orig", "gene_symbol"]], on="orig", how="left")
        long = long.drop(columns="orig")
        long["ensembl_id"] = long["gene_symbol"].map(gene_resolver.to_ensembl)
        long["value_unit"] = "log2_MS_intensity"
        return self._pack(long)


class HPARNASource(OmicsSource):
    """HPA cell-line RNA: File 1 long table plus File 11 metadata with CVCL.

    The ID path uses two exact lookups without string normalization:
        File 1 'Cell line' --File 11--> CVCL_xxxx --Cellosaurus--> ACH-xxxxxx

    File 1 contains approximately 24 million rows. v3 reads it in chunks and
    filters by genes, so memory use scales with the number of target genes
    rather than file size.
    """
    layer, source = "rna", "hpa_rna"

    def __init__(self, expr_path, desc_path,
                 version: str = "HPA-v23",
                 value_col: str = "nTPM",
                 chunksize: int = 2_000_000):
        self.expr_path = Path(expr_path)
        self.desc_path = Path(desc_path)
        self.version = version
        self.value_col = value_col
        self.chunksize = chunksize
        self._name_to_cvcl: Optional[dict] = None

    def _load_desc(self, verbose=True) -> dict:
        if self._name_to_cvcl is None:
            desc = pd.read_csv(self.desc_path, sep="\t", dtype="string")
            name_col = _find_col(desc, ["cellline"], "File 11 Cell line column")
            cvcl_col = _find_col(desc, ["cellosaurusid", "cellosaurus", "cvcl"],
                                 "File 11 Cellosaurus ID column")
            d = pd.DataFrame({"n": desc[name_col].str.strip(),
                              "c": desc[cvcl_col].str.strip()}).dropna(subset=["n"])
            self._name_to_cvcl = dict(zip(d["n"], d["c"].fillna("")))
            if verbose:
                n_ok = sum(1 for v in self._name_to_cvcl.values()
                           if isinstance(v, str) and v.startswith("CVCL_"))
                print(f"  [hpa_rna] File 11: {len(self._name_to_cvcl)} cell lines, "
                      f"{n_ok} with CVCL_id")
        return self._name_to_cvcl

    def load_long(self, resolver, gene_resolver, genes=None):
        name_to_cvcl = self._load_desc()

        head = pd.read_csv(self.expr_path, sep="\t", nrows=5)
        name_col = _find_col(head, ["cellline"], "File 1 Cell line column")
        ensg_col = _find_col(head, ["gene"], "File 1 Gene (ENSG) column")
        sym_col = _find_col(head, ["genename"], "File 1 Gene name column")
        val_col = _find_col(head, [self.value_col], f"File 1 {self.value_col} column")

        usecols = [name_col, ensg_col, sym_col, val_col]
        pieces = []
        reader = pd.read_csv(self.expr_path, sep="\t", usecols=usecols,
                             chunksize=self.chunksize)
        for chunk in reader:
            if genes is not None:
                chunk = chunk[chunk[sym_col].isin(genes)]
            if not len(chunk):
                continue
            pieces.append(chunk)
        expr = (pd.concat(pieces, ignore_index=True) if pieces
                else pd.DataFrame(columns=usecols))

        if not len(expr):
            return self._pack(pd.DataFrame(columns=UNIFIED_COLUMNS))

        # Two-hop resolution: name -> CVCL -> ACH.
        cvcl = expr[name_col].astype(str).str.strip().map(name_to_cvcl)
        ach = resolver.resolve_cvcl_series(cvcl)
        # Fallback for File 11 entries without CVCL: resolve names directly
        # through Cellosaurus aliases.
        need_fb = ach.isna()
        if need_fb.any():
            fb = resolver.resolve_name_series(expr.loc[need_fb, name_col])
            ach = ach.copy()
            ach[need_fb] = fb
            if fb.notna().sum():
                print(f"  [hpa_rna] Fallback recovered {fb.notna().sum():,} rows by direct name lookup")

        # out = pd.DataFrame({
        #     "DepMap_ID": ach,
        #     "gene_symbol": expr[sym_col],
        #     "ensembl_id": expr[ensg_col].astype(str).str.split(".").str[0],
        #     "value_raw": pd.to_numeric(expr[val_col], errors="coerce"),
        #     "value_unit": self.value_col,
        # })

        out = pd.DataFrame({
            "DepMap_ID": ach,
            "gene_symbol": expr[sym_col],
            "ensembl_id": expr[ensg_col].astype(str).str.split(".").str[0],
            "value_raw": np.log2(pd.to_numeric(expr[val_col], errors="coerce").clip(lower=0) + 1),
            "value_unit": f"log2({self.value_col}+1)",
        })

        return self._pack(out)


class GEOWideSource(OmicsSource):
    """GEO expression: File 3 wide table plus File 10 sample metadata.

    File 3 uses genes/ENSG as rows and one GSM per column. v3.1 makes two
    corrections for the full datasets:

    1. The gene column in File 3 contains ENSG identifiers, not gene symbols
       (`ENSG00000000003`, not `TSPAN6`). Target symbols are converted to ENSG
       before filtering and converted back after loading.
    2. File 10 contains a `Cellosaurus_ID` column pre-aligned through
       Cellosaurus GEO cross-references (`Matching_Type` = 'Cello GEO GSM').
       Exact CVCL lookup is therefore as reliable as the HPA path; name matching
       is used only as a fallback.

    Exact ID path:
        GSM --File 10 Cellosaurus_ID--> CVCL_xxxx --RRID/xref--> ACH-xxxxxx

    Multiple GSM samples often map to one cell line because of technical
    replicates or distinct GSE records. build_gene_table averages these values
    by DepMap_ID.
    """
    layer, source = "rna", "geo_rna"

    def __init__(self, expr_path, info_path,
                 version: str = "GEO-import",
                 log_transform: bool = True,
                 chunksize: int = 5000):
        """
        log_transform : File 3 contains linear microarray intensities
            (33.6 / 553.2 / 2182.3), whereas DepMap and HPA use log scales.
            The default log2(x+1) transform makes the three RNA sources
            comparable. `value_unit` records 'log2(intensity+1)' for
            traceability. False retains the original linear values.
        """
        self.expr_path = Path(expr_path)
        self.info_path = Path(info_path)
        self.version = version
        self.log_transform = log_transform
        self.chunksize = chunksize
        self._gsm_to_ach: Optional[pd.Series] = None

    def _load_gsm_map(self, resolver, verbose=True) -> pd.Series:
        """Build GSM-to-ACH mappings from File 10.

        Prefer the Cellosaurus_ID column and use names only as a fallback.
        """
        if self._gsm_to_ach is not None:
            return self._gsm_to_ach

        info = pd.read_csv(self.info_path, sep="\t", dtype="string")
        gsm_col = _find_col(info, ["geoaccession", "gsm", "sampleid", "accession"],
                            "File 10 GSM column")
        cvcl_col = _find_col_optional(info, ["cellosaurusid", "cellosaurus", "cvcl"])
        cell_col = _find_col_optional(info, ["cellline", "sourcename", "title"])

        if verbose:
            print(f"  [geo_rna] File 10 columns: GSM={gsm_col!r}, "
                  f"CVCL={cvcl_col!r}, name={cell_col!r}")

        gsm = info[gsm_col].astype(str).str.strip()
        ach = pd.Series([None] * len(info), index=info.index, dtype="object")

        # --- Primary path: exact Cellosaurus_ID lookup ---
        n_cvcl = 0
        if cvcl_col is not None:
            cvcl = (info[cvcl_col].astype("string").str.strip()
                    .str.extract(r"(CVCL_[A-Za-z0-9]+)", expand=False))
            ach = resolver.resolve_cvcl_series(cvcl)
            n_cvcl = ach.notna().sum()
            if verbose:
                n_has_cvcl = cvcl.notna().sum()
                print(f"  [geo_rna] Primary CVCL path: {n_has_cvcl:,} samples have CVCL, "
                      f"of which {n_cvcl:,} map to DepMap")

        # --- Fallback: resolve cell-line names through Cellosaurus aliases ---
        if cell_col is not None:
            need = ach.isna()
            if need.any():
                names = info.loc[need, cell_col].astype("string")
                names = names.str.replace(r"(?i)^.*cell\s*line\s*[:=]\s*", "",
                                          regex=True)
                fb = resolver.resolve_name_series(names)
                ach = ach.copy()
                ach[need] = fb
                if verbose and fb.notna().sum():
                    print(f"  [geo_rna] Fallback path recovered "
                          f"{fb.notna().sum():,} samples by direct name lookup")

        m = pd.Series(ach.values, index=gsm.values)
        m = m[m.notna()]
        m = m[~m.index.duplicated(keep="first")]
        self._gsm_to_ach = m
        if verbose:
            print(f"  [geo_rna] File 10 total: {len(info):,} samples -> "
                  f"{len(m):,} GSM identifiers resolved to ACH ({len(m)/max(len(info),1):.1%}), "
                  f"covering {m.nunique():,} distinct cell lines")
        return self._gsm_to_ach

    def load_long(self, resolver, gene_resolver, genes=None):
        gsm_map = self._load_gsm_map(resolver)
        if not len(gsm_map):
            print("  [geo_rna] WARNING: no GSM identifiers resolved to ACH; returning an empty table")
            return self._pack(pd.DataFrame(columns=UNIFIED_COLUMNS))

        head = pd.read_csv(self.expr_path, sep="\t", nrows=5)
        gene_col = head.columns[0]        # The first column in File 3 is the gene column.

        # --- Determine whether the gene column contains ENSG IDs or symbols ---
        sample_vals = head[gene_col].astype(str)
        is_ensg = sample_vals.str.startswith("ENSG").mean() > 0.5
        if is_ensg:
            print(f"  [geo_rna] Gene column {gene_col!r} contains ENSG IDs; converting to symbols")
        else:
            print(f"  [geo_rna] Gene column {gene_col!r} contains gene symbols")

        # --- Convert target genes to the identifier space used by the file ---
        wanted = None
        if genes is not None:
            if is_ensg:
                wanted = {e for e in (gene_resolver.to_ensembl(g) for g in genes)
                          if e}
                missing = {g for g in genes if not gene_resolver.to_ensembl(g)}
                if missing:
                    print(f"  [geo_rna] WARNING: these genes have no ENSG mapping and cannot be queried in GEO: "
                          f"{sorted(missing)}")
                if not wanted:
                    return self._pack(pd.DataFrame(columns=UNIFIED_COLUMNS))
            else:
                wanted = set(genes)

        # --- Read only the gene column and resolvable GSM columns; reading all
        # 3,267 columns is slow. ---
        keep_gsm = [c for c in head.columns[1:] if c in gsm_map.index]
        if not keep_gsm:
            print("  [geo_rna] WARNING: no GSM column in File 3 resolves to ACH through File 10")
            return self._pack(pd.DataFrame(columns=UNIFIED_COLUMNS))
        print(f"  [geo_rna] File 3: {len(head.columns)-1:,} GSM columns; "
              f"{len(keep_gsm):,} are resolvable and will be read")

        pieces = []
        reader = pd.read_csv(self.expr_path, sep="\t",
                             usecols=[gene_col] + keep_gsm,
                             chunksize=self.chunksize)
        for chunk in reader:
            if wanted is not None:
                chunk = chunk[chunk[gene_col].astype(str).str.split(".").str[0]
                              .isin(wanted)]
            if len(chunk):
                pieces.append(chunk)
        expr = pd.concat(pieces, ignore_index=True) if pieces else None
        if expr is None or not len(expr):
            print("  [geo_rna] No rows in File 3 match the target genes")
            return self._pack(pd.DataFrame(columns=UNIFIED_COLUMNS))

        long = expr.melt(id_vars=gene_col, var_name="GSM", value_name="value_raw")
        long["DepMap_ID"] = long["GSM"].map(gsm_map)
        long = long[long["value_raw"].notna()]

        # --- Restore gene identifiers ---
        gene_key = long[gene_col].astype(str).str.split(".").str[0]
        if is_ensg:
            long["ensembl_id"] = gene_key
            long["gene_symbol"] = gene_key.map(gene_resolver.to_symbol)
            n_unmapped = long["gene_symbol"].isna().sum()
            if n_unmapped:
                print(f"  [geo_rna] Skipped {n_unmapped:,} rows whose ENSG IDs could not be converted to symbols")
        else:
            long["gene_symbol"] = gene_key
            long["ensembl_id"] = gene_key.map(gene_resolver.to_ensembl)

        # --- Align measurement scales ---
        if self.log_transform:
            v = pd.to_numeric(long["value_raw"], errors="coerce")
            long["value_raw"] = np.log2(v.clip(lower=0) + 1)
            long["value_unit"] = "log2(intensity+1)"
        else:
            long["value_unit"] = "linear_intensity"

        return self._pack(long)


# ---------------------------------------------------------------------------
# 5b. Event adapters: mutation and fusion
# ---------------------------------------------------------------------------
#
# Unlike RNA and protein measurements, mutations and fusions are events. A
# (gene, cell-line) pair may have zero, one, or multiple records. Records are
# therefore aggregated to one row per (DepMap_ID, gene), with value_raw storing
# the event count. Functional-impact flags outside UNIFIED_COLUMNS capture
# high-impact mutations and fusion confidence for use by the scoring model.
#
# _pack retains only UNIFIED_COLUMNS, so it would discard extra aggregate
# columns. These sources therefore represent the count in value_raw and expose
# functional flags as separate virtual sources. For mutation data, value_raw is
# the mutation count and high_impact (0/1) is stored in a companion source.


class MutationSource(OmicsSource):
    """DepMap somatic mutations (File 6).

    Cell lines in File 6 are identified by **ProfileID (PR-xxxxxx)** rather
    than ACH IDs. The ProfileID-to-ACH bridge in File 8 (OmicsProfiles) is
    required, so the existing `profiles` / `_profile_to_ach()` implementation
    in data_loader is reused.

    Produces two sources within layer='mutation':
      * mutation_count      : number of mutations per (cell line, gene)
      * mutation_highimpact : whether any high-impact mutation is present
                              (oncogene / TSG / LoF) -> 1.0 / 0.0

    This lets the scoring model use both mutation presence and the more
    informative presence of a potentially pathogenic mutation.
    """
    layer = "mutation"
    source = "depmap_mutation"  # Used by the has_mutation mask.

    def __init__(self, loader, path=None, version: str = "DepMap-24Q2",
                 chunksize: int = 500_000):
        self.loader = loader
        self.version = version
        self.chunksize = chunksize
        # File 6 path: use an explicit path when provided; otherwise ask the loader.
        self.path = Path(path) if path else Path(loader._path(6))

    def load_long(self, resolver, gene_resolver, genes=None):
        prof2ach = self.loader._profile_to_ach()  # Reuse the existing bridge.

        head = pd.read_csv(self.path, sep=",", nrows=5)
        sym_col = _find_col(head, ["hugosymbol", "genesymbol", "gene"], "File 6 gene-symbol column")
        ensg_col = _find_col_optional(head, ["ensemblgeneid", "ensembl"])
        prof_col = _find_col(head, ["profileid"], "File 6 ProfileID column")
        # High-impact indicators are used only when present.
        onco_col = _find_col_optional(head, ["oncogenehighimpact"])
        tsg_col = _find_col_optional(head, ["tumorsuppressorhighimpact"])
        lof_col = _find_col_optional(head, ["likelylof"])
        impact_col = _find_col_optional(head, ["vepimpact"])

        usecols = [c for c in [prof_col, sym_col, ensg_col, onco_col, tsg_col,
                               lof_col, impact_col] if c]

        pieces = []
        for chunk in pd.read_csv(self.path, sep=",", usecols=usecols,
                                 chunksize=self.chunksize, low_memory=False):
            if genes is not None:
                chunk = chunk[chunk[sym_col].isin(genes)]
            if len(chunk):
                pieces.append(chunk)
        df = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(columns=usecols)
        if not len(df):
            return self._pack(pd.DataFrame(columns=UNIFIED_COLUMNS))

        df["DepMap_ID"] = df[prof_col].map(prof2ach)

        # High-impact flag: true when any indicator is true.
        def _truthy(col):
            if col is None or col not in df.columns:
                return pd.Series(False, index=df.index)
            s = df[col].astype("string").str.strip().str.lower()
            return s.isin(["true", "1", "yes", "high"])

        high = _truthy(onco_col) | _truthy(tsg_col) | _truthy(lof_col)
        if impact_col in df.columns:
            high = high | df[impact_col].astype("string").str.upper().eq("HIGH")
        df["_high"] = high.astype(int)

        # Aggregate mutation count and high-impact status per (cell line, gene).
        agg = (df.dropna(subset=["DepMap_ID"])
               .groupby(["DepMap_ID", sym_col])
               .agg(count=("_high", "size"), high=("_high", "max"))
               .reset_index()
               .rename(columns={sym_col: "gene_symbol"}))
        if ensg_col and ensg_col in df.columns:
            ensg_map = (df[[sym_col, ensg_col]].dropna().drop_duplicates()
            .set_index(sym_col)[ensg_col].astype(str).str.split(".").str[0])
            agg["ensembl_id"] = agg["gene_symbol"].map(ensg_map)
        else:
            agg["ensembl_id"] = agg["gene_symbol"].map(gene_resolver.to_ensembl)

        # Primary source: mutation count.
        base = agg.rename(columns={"count": "value_raw"}).copy()
        base["value_unit"] = "mutation_count"
        out_count = self._pack(base[["DepMap_ID", "gene_symbol", "ensembl_id",
                                     "value_raw", "value_unit"]])

        # Companion source: expose the high-impact flag separately for scoring.
        hi = agg.copy()
        hi["value_raw"] = hi["high"].astype(float)
        hi["value_unit"] = "high_impact_flag"
        hi = hi[["DepMap_ID", "gene_symbol", "ensembl_id", "value_raw", "value_unit"]]
        hi_packed = hi.reindex(columns=UNIFIED_COLUMNS)
        hi_packed["source"] = "depmap_mutation_highimpact"
        hi_packed["source_version"] = self.version
        hi_packed["omics_layer"] = "mutation"
        hi_packed = hi_packed[hi_packed["DepMap_ID"].notna()
                              & hi_packed["gene_symbol"].notna()].reset_index(drop=True)

        return pd.concat([out_count, hi_packed], ignore_index=True)


class FusionSource(OmicsSource):
    """DepMap gene fusions (File 5).

    File 5 identifies cell lines with **ModelID = ACH-xxxxxx**, which already
    corresponds to DepMap_ID and requires no conversion.

    Each fusion involves two genes (gene1--gene2), so every record is split
    into two rows, one for each gene. A query for EGFR therefore matches EGFR
    whether it appears in gene1 or gene2.

    Produces two sources within layer='fusion':
      * fusion_count       : number of fusions per (cell line, gene)
      * fusion_confidence  : maximum confidence
                             (high=3 / medium=2 / low=1, normalized to 0-1)
    """
    layer = "fusion"
    source = "depmap_fusion"

    _CONF = {"high": 3, "medium": 2, "low": 1}

    def __init__(self, loader, path=None, version: str = "DepMap-24Q2",
                 chunksize: int = 500_000):
        self.loader = loader
        self.version = version
        self.chunksize = chunksize
        self.path = Path(path) if path else Path(loader._path(5))

    @staticmethod
    def _parse_gene(cell):
        """'DLG1 (ENSG00000075711.21)' -> ('DLG1', 'ENSG00000075711'); '(.)' -> (sym, None)"""
        if not isinstance(cell, str):
            return None, None
        m = re.match(r"\s*([^(]+?)\s*\(([^)]*)\)", cell)
        if not m:
            return cell.strip() or None, None
        sym = m.group(1).strip() or None
        ensg = m.group(2).split(".")[0]
        ensg = ensg if ensg.startswith("ENSG") else None
        return sym, ensg

    def load_long(self, resolver, gene_resolver, genes=None):
        head = pd.read_csv(self.path, sep=",", nrows=5)
        model_col = _find_col(head, ["modelid"], "File 5 ModelID column")
        g1_col = _find_col(head, ["gene1"], "File 5 gene1 column")
        g2_col = _find_col(head, ["gene2"], "File 5 gene2 column")
        conf_col = _find_col_optional(head, ["confidence"])

        usecols = [c for c in [model_col, g1_col, g2_col, conf_col] if c]
        df = pd.read_csv(self.path, sep=",", usecols=usecols, low_memory=False)

        # Split gene1 and gene2 into separate rows.
        g1 = df[g1_col].map(self._parse_gene)
        g2 = df[g2_col].map(self._parse_gene)
        rows = pd.DataFrame({
            "DepMap_ID": pd.concat([df[model_col], df[model_col]], ignore_index=True),
            "gene_symbol": pd.concat([g1.str[0], g2.str[0]], ignore_index=True),
            "ensembl_id": pd.concat([g1.str[1], g2.str[1]], ignore_index=True),
            "conf": pd.concat([df[conf_col], df[conf_col]] if conf_col
                              else [pd.Series("", index=df.index)] * 2,
                              ignore_index=True),
        })
        if genes is not None:
            rows = rows[rows["gene_symbol"].isin(genes)]
        rows = rows.dropna(subset=["DepMap_ID", "gene_symbol"])
        if not len(rows):
            return self._pack(pd.DataFrame(columns=UNIFIED_COLUMNS))

        rows["conf_score"] = (rows["conf"].astype("string").str.lower()
                              .map(self._CONF).fillna(0))

        agg = (rows.groupby(["DepMap_ID", "gene_symbol"])
               .agg(count=("conf_score", "size"),
                    max_conf=("conf_score", "max"),
                    ensembl_id=("ensembl_id", "first"))
               .reset_index())

        # Primary source: fusion count.
        c = agg.rename(columns={"count": "value_raw"}).copy()
        c["value_unit"] = "fusion_count"
        out_count = self._pack(c[["DepMap_ID", "gene_symbol", "ensembl_id",
                                  "value_raw", "value_unit"]])

        # Companion source: confidence normalized to 0-1.
        cf = agg.copy()
        cf["value_raw"] = cf["max_conf"] / 3.0
        cf["value_unit"] = "fusion_confidence"
        cf = cf[["DepMap_ID", "gene_symbol", "ensembl_id", "value_raw", "value_unit"]]
        cf_packed = cf.reindex(columns=UNIFIED_COLUMNS)
        cf_packed["source"] = "depmap_fusion_confidence"
        cf_packed["source_version"] = self.version
        cf_packed["omics_layer"] = "fusion"
        cf_packed = cf_packed[cf_packed["DepMap_ID"].notna()
                              & cf_packed["gene_symbol"].notna()].reset_index(drop=True)

        return pd.concat([out_count, cf_packed], ignore_index=True)


# ---------------------------------------------------------------------------
# 6. Merge engine
# ---------------------------------------------------------------------------

class MultiOmicsMerger:
    def __init__(self, all_cell_lines: Iterable[str],
                 cell_resolver: CellLineIDResolver,
                 gene_resolver: GeneIDResolver):
        self.all_cell_lines = pd.Index(sorted(set(all_cell_lines)), name="DepMap_ID")
        self.cell_resolver = cell_resolver
        self.gene_resolver = gene_resolver
        self._sources: list[OmicsSource] = []
        self._cache: dict[frozenset | None, pd.DataFrame] = {}

    def register(self, source: OmicsSource) -> "MultiOmicsMerger":
        self._sources.append(source)
        self._cache.clear()
        return self

    def build_long_table(self, genes: Optional[Iterable[str]] = None) -> pd.DataFrame:
        """Build the unified long-format table.

        Supplying genes (for example, {'EGFR', 'KRAS'}) is strongly
        recommended. Omitting this argument attempts to load every gene and
        may require tens of gigabytes of memory with the full dataset.
        """
        key = frozenset(genes) if genes is not None else None
        if key in self._cache:
            return self._cache[key]
        if key is None:
            print("⚠️  No genes were specified. Loading every gene may exhaust memory.\n"
                  "    Recommended: build_long_table(genes={'EGFR', 'KRAS', ...})")
        gene_set = set(genes) if genes is not None else None

        frames = []
        for src in self._sources:
            print(f"→ {src.source} ({src.layer})")
            frames.append(src.load_long(self.cell_resolver, self.gene_resolver, gene_set))
        long = (pd.concat(frames, ignore_index=True) if frames
                else pd.DataFrame(columns=UNIFIED_COLUMNS))

        if len(long):
            long["value_std"] = (long.groupby(["source", "gene_symbol"])["value_raw"]
                                     .transform(lambda s: (s - s.mean()) / s.std(ddof=0)
                                                if s.std(ddof=0) and s.notna().sum() > 1
                                                else np.nan))
        else:
            long["value_std"] = pd.Series(dtype=float)
        self._cache[key] = long
        return long

    def build_gene_table(self, gene: str) -> pd.DataFrame:
        """Build a single-gene wide table with one row per cell line for scoring."""
        long = self.build_long_table(genes={gene})
        g = long[long["gene_symbol"] == gene]
        base = pd.DataFrame(index=self.all_cell_lines).reset_index()

        for src in sorted(g["source"].unique()):
            sub = (g[g["source"] == src]
                   .groupby("DepMap_ID")[["value_raw", "value_std"]].mean())
            base = base.merge(
                sub.rename(columns={"value_raw": f"{src}__raw",
                                    "value_std": f"{src}__std"}),
                on="DepMap_ID", how="left")

        for layer in OMICS_LAYERS:
            srcs = [s.source for s in self._sources if s.layer == layer]
            raw_cols = [f"{s}__raw" for s in srcs if f"{s}__raw" in base.columns]
            base[f"has_{layer}"] = (base[raw_cols].notna().any(axis=1)
                                    if raw_cols else False)

        rna_std_cols = [f"{s.source}__std" for s in self._sources
                        if s.layer == "rna" and f"{s.source}__std" in base.columns]
        if len(rna_std_cols) >= 2:
            spread = base[rna_std_cols].std(axis=1, ddof=0)
            base["rna_consistency"] = (1 - spread.clip(0, 1)).where(
                base[rna_std_cols].notna().sum(axis=1) >= 2)
        else:
            base["rna_consistency"] = np.nan

        present = [f"has_{l}" for l in OMICS_LAYERS
                   if any(s.layer == l for s in self._sources)]
        base["data_completeness"] = base[present].mean(axis=1) if present else 0.0
        return base

# Build the master table.
    def build_master_table(self, genes: Iterable[str]) -> pd.DataFrame:
        """Stack the wide tables for multiple genes into one master table.

        Each gene is first processed with build_gene_table (one row per cell
        line), labelled in a gene column, and then concatenated. In the
        resulting long-format master table, (gene, DepMap_ID) uniquely
        identifies a row; the remaining columns contain source values, masks,
        and completeness.

        This is the single-file deliverable and contains only the supplied
        genes. A wide table for all approximately 20,000 genes would be
        impractically large and is unnecessary because the scoring model
        queries one gene at a time.
        """
        frames = []
        for g in genes:
            t = self.build_gene_table(g)
            t.insert(1, "gene", g)
            frames.append(t)
        if not frames:
            return pd.DataFrame()
        # Source columns may differ by gene; concat aligns them and fills with NaN.
        master = pd.concat(frames, ignore_index=True)
        # Keep source data first and move masks and completeness to the end.
        tail = [c for c in master.columns
                if c.startswith("has_") or c == "data_completeness"
                or c == "rna_consistency"]
        head = [c for c in master.columns if c not in tail]
        return master[head + tail]


# ---------------------------------------------------------------------------
# 7. Missing-data weight redistribution
# ---------------------------------------------------------------------------

def redistribute_weights(mask: dict[str, bool],
                         base_weights: dict[str, float]) -> dict[str, float]:
    """Redistribute missing-layer weight proportionally without removing samples."""
    available = {k: w for k, w in base_weights.items() if mask.get(k, False)}
    total = sum(available.values())
    if total == 0:
        return {k: 0.0 for k in base_weights}
    return {k: (available.get(k, 0.0) / total) for k in base_weights}


# ---------------------------------------------------------------------------
# 8. Self-check entry point: python data_merger.py <data_dir>
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    print("=" * 66)
    print(f"data_merger.py self-check   (__version__ = {__version__})")
    print("=" * 66)

    # Locate the data directory: command-line argument, current/data, then parent/data.
    if len(sys.argv) > 1:
        data_dir = Path(sys.argv[1])
    else:
        here = Path(__file__).resolve().parent
        candidates = [Path.cwd() / "data", here / "data", here.parent / "data"]
        data_dir = next((c for c in candidates if c.exists()), candidates[0])

    print(f"data_dir = {data_dir.resolve()}")
    if not data_dir.exists():
        print("✗ Directory not found. Usage: python src/data_merger.py <data directory>")
        sys.exit(1)

    cellosaurus_path = data_dir / "nomenclature" / "7_cellosaurus.csv"
    print(f"\n[1/3] Reading Cellosaurus: {cellosaurus_path.name}")
    cs = pd.read_csv(cellosaurus_path)
    print(f"      {cs.shape[0]:,} rows x {cs.shape[1]} columns")

    print("\n[2/3] Checking whether Cross-references contain DepMap entries")
    diagnose_cellosaurus(cs)

    print("\n[3/3] Building CellLineIDResolver")
    try:
        from data_loader import CellLineDataLoader
        loader = CellLineDataLoader(data_dir)
        si = loader.sample_info
        print(f"      sample_info: {si.shape}")
        print(f"      sample_info columns: {list(si.columns)}")
    except Exception as e:
        print(f"      (data_loader skipped: {e})")
        si = None

    res = CellLineIDResolver(cs, si)
    print(f"\nFinal mapping sizes: {res.stats()}")

    print("\nSpot checks:")
    for probe in ["CVCL_0023", "A549", "A-549", "HeLa", "MCF7"]:
        if probe.startswith("CVCL_"):
            print(f"  CVCL {probe:12} -> {res.resolve_cvcl_to_ach(probe)}")
        else:
            print(f"  name {probe:12} -> {res.resolve_name(probe)}")
    print("\n✓ Self-check complete")
