from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd

# ============================== CONFIG ==============================
PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_RAW_DATA_DIR = PROJECT_ROOT / "data_s3"
DEFAULT_CACHE_DIR = PROJECT_ROOT / "gene_cache"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "outputs"

RNA_WEIGHT = 0.55
PROTEIN_WEIGHT = 0.30
CONFIDENCE_WEIGHT = 0.15
MAX_EXCLUSION_PENALTY = 0.30

# Protein-only query mode keeps the same overall split:
# 85% biological evidence + 15% confidence.
PROTEIN_ONLY_PROTEIN_WEIGHT = 0.85
PROTEIN_ONLY_CONFIDENCE_WEIGHT = 0.15

RNA_STD_COLUMNS = ["depmap_rna__std", "hpa_rna__std", "geo_rna__std"]
PROTEIN_STD_COLUMN = "ccle_gygi_protein__std"
PROTEIN_RAW_COLUMN = "ccle_gygi_protein__raw"
MUTATION_RAW_COLUMN = "depmap_mutation__raw"
MUTATION_HIGH_IMPACT_RAW_COLUMN = "depmap_mutation_highimpact__raw"
FUSION_RAW_COLUMN = "depmap_fusion__raw"
FUSION_CONFIDENCE_RAW_COLUMN = "depmap_fusion_confidence__raw"

EXPECTED_EVIDENCE_COLUMNS = [
    "depmap_rna__raw", "depmap_rna__std",
    "hpa_rna__raw", "hpa_rna__std",
    "geo_rna__raw", "geo_rna__std",
    "ccle_gygi_protein__raw", "ccle_gygi_protein__std",
    "depmap_mutation__raw", "depmap_mutation__std",
    "depmap_mutation_highimpact__raw", "depmap_mutation_highimpact__std",
    "depmap_fusion__raw", "depmap_fusion__std",
    "depmap_fusion_confidence__raw", "depmap_fusion_confidence__std",
]

DISEASE_ALIASES: Dict[str, List[str]] = {
    "brain": ["brain", "central nervous system", "cns", "glioma", "glioblastoma", "astrocytoma", "medulloblastoma"],
    "glioblastoma": ["glioblastoma", "gbm"],
    "glioma": ["glioma", "diffuse glioma", "high grade glioma", "low grade glioma"],
    "medulloblastoma": ["medulloblastoma"],
    "astrocytoma": ["astrocytoma"],
    "lung": ["lung", "non-small cell lung", "nsclc", "small cell lung", "sclc"],
    "breast": ["breast", "mammary"],
    "cervical": ["cervical", "cervix"],
    "colorectal": ["colorectal", "colon", "rectal", "large intestine"],
    "prostate": ["prostate"],
    "pancreatic": ["pancreatic", "pancreas"],
    "ovarian": ["ovarian", "ovary"],
    "liver": ["liver", "hepatic", "hepatocellular"],
    "kidney": ["kidney", "renal"],
    "skin": ["skin", "melanoma"],
    "gastric": ["gastric", "stomach"],
    "esophageal": ["esophageal", "oesophageal", "esophagus", "oesophagus"],
    "blood": ["blood", "haematopoietic", "hematopoietic", "leukaemia", "leukemia", "lymphoma", "myeloma"],
}

# ============================== MODULE LOADING ==============================
def _import_module_from_file(module_name: str, file_path: Path):
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module from: {file_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def load_team_modules():
    here = Path(__file__).resolve().parent
    module_dir = None
    for folder in [here / "src", here]:
        if (folder / "data_loader.py").exists() and (folder / "data_merger.py").exists():
            module_dir = folder
            break
    if module_dir is None:
        raise FileNotFoundError(
            "Cannot find data_loader.py and data_merger.py.\n"
            f"Put them beside this script or in: {here / 'src'}"
        )

    print(f"Team module directory: {module_dir}")
    dl = _import_module_from_file("cellline_team_data_loader", module_dir / "data_loader.py")
    dm = _import_module_from_file("cellline_team_data_merger", module_dir / "data_merger.py")

    names = [
        "CellLineIDResolver", "GeneIDResolver", "MultiOmicsMerger",
        "DepMapRNASource", "HPARNASource", "GEOWideSource",
        "MatrixProteinSource", "MutationSource", "FusionSource", "OMICS_LAYERS"
    ]
    if not hasattr(dl, "CellLineDataLoader"):
        raise ImportError("data_loader.py does not define CellLineDataLoader")
    for name in names:
        if not hasattr(dm, name):
            raise ImportError(f"data_merger.py does not define {name}")

    return dl, dm

# ============================== UTILITIES ==============================
def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        if value is None or pd.isna(value):
            return default
        return float(value)
    except Exception:
        return default


def min_max_scale(values: Sequence[float]) -> List[float]:
    values = [float(v) for v in values]
    if not values:
        return []
    lo, hi = min(values), max(values)
    if hi == lo:
        return [1.0 if hi > 0 else 0.0 for _ in values]
    return [(v - lo) / (hi - lo) for v in values]


def normalise_disease_key(text: Optional[str]) -> str:
    value = "" if text is None else str(text).lower().strip()
    value = value.replace("-", " ").replace("_", " ").replace("/", " ")
    return re.sub(r"\s+", " ", value).strip()


def build_disease_aliases(disease: str) -> List[str]:
    query = normalise_disease_key(disease)
    if not query:
        raise ValueError("Disease/tissue is required.")
    stop = {"cancer", "tumour", "tumor", "carcinoma", "disease", "cell", "cells", "line", "lines"}
    compact = " ".join(x for x in query.split() if x not in stop).strip()
    aliases = {query}
    if compact:
        aliases.add(compact)
    key = compact if compact in DISEASE_ALIASES else query if query in DISEASE_ALIASES else None
    if key:
        aliases.update(DISEASE_ALIASES[key])
    else:
        for canonical, known in DISEASE_ALIASES.items():
            if query in known or compact in known:
                aliases.update([canonical, query, compact])
                break
    return sorted({normalise_disease_key(x) for x in aliases if x}, key=len, reverse=True)


def confidence_level(score: float) -> str:
    if score >= 0.80: return "High"
    if score >= 0.65: return "Medium-High"
    if score >= 0.50: return "Medium"
    if score >= 0.35: return "Medium-Low"
    return "Low"


def recommendation_level(score: float) -> str:
    if score >= 0.80: return "Strongly Recommended"
    if score >= 0.65: return "Recommended"
    if score >= 0.50: return "Conditionally Recommended"
    if score >= 0.35: return "Low Priority"
    return "Not Recommended"


def mean_available(row: pd.Series, columns: Sequence[str]) -> Optional[float]:
    vals = [float(row[c]) for c in columns if c in row.index and pd.notna(row[c])]
    return float(np.mean(vals)) if vals else None


def count_available(row: pd.Series, columns: Sequence[str]) -> int:
    return sum(1 for c in columns if c in row.index and pd.notna(row[c]))


def alteration_status(mut: bool, fusion: bool, high: bool = False) -> str:
    if mut and fusion: return "MIXED: MUTATION + FUSION"
    if mut: return "MUTATION (HIGH-IMPACT)" if high else "MUTATION"
    if fusion: return "FUSION"
    return "NO ALTERATION"


def gene_label(gene: str, mut: bool, fusion: bool, high: bool = False) -> str:
    status = alteration_status(mut, fusion, high)
    return gene if status == "NO ALTERATION" else f"{gene} [{status}]"


def determine_query_mode(
    target_gene: Optional[str],
    target_protein: Optional[str],
) -> str:
    """Return one of: GENE_MULTIOMICS, PROTEIN_ONLY, COMBINED."""
    has_gene = bool(target_gene and str(target_gene).strip())
    has_protein = bool(target_protein and str(target_protein).strip())

    if has_gene and has_protein:
        return "COMBINED"
    if has_gene:
        return "GENE_MULTIOMICS"
    if has_protein:
        return "PROTEIN_ONLY"

    raise ValueError("At least one of target gene or target protein must be provided.")


def find_required_file(root: Path, candidates: Sequence[str], label: str) -> Path:
    for rel in candidates:
        p = root / rel
        if p.exists():
            return p
    names = {Path(x).name.lower() for x in candidates}
    for p in root.rglob("*"):
        if p.is_file() and p.name.lower() in names:
            return p
    if "hpa expression" in label.lower():
        for p in root.rglob("*.tsv"):
            n = p.name.lower()
            if "hpa" in n and "rna" in n and ("cellline" in n or "celline" in n) and "description" not in n:
                return p
    if "hpa description" in label.lower():
        for p in root.rglob("*.tsv"):
            n = p.name.lower()
            if "hpa" in n and "description" in n:
                return p
    raise FileNotFoundError(f"Cannot find {label}")

# ============================== DYNAMIC DATA ACCESS ==============================
class DynamicMultiOmicsRecommender:
    def __init__(self, raw_data_dir: Path, cache_dir: Path, refresh_cache: bool = False):
        self.raw_data_dir = Path(raw_data_dir)
        self.cache_dir = Path(cache_dir)
        self.refresh_cache = refresh_cache
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        if not self.raw_data_dir.exists():
            raise FileNotFoundError(f"Raw data directory not found: {self.raw_data_dir}")

        dl, dm = load_team_modules()
        self.OMICS_LAYERS = dm.OMICS_LAYERS
        print(f"Raw data directory: {self.raw_data_dir}")
        print(f"Gene cache directory: {self.cache_dir}")

        self.loader = dl.CellLineDataLoader(self.raw_data_dir)
        self.sample_info = self.loader.sample_info.copy()

        cellosaurus_path = find_required_file(self.raw_data_dir, ["nomenclature/7_cellosaurus.csv"], "Cellosaurus")
        hpa_expr_path = find_required_file(
            self.raw_data_dir,
            ["gene expression/1_4_hpa_rna_cellline.tsv", "gene expression/1_4_hpa_rna_celline.tsv"],
            "HPA expression",
        )
        hpa_desc_path = find_required_file(
            self.raw_data_dir,
            ["nomenclature/11_hpa_rna_cellline_description.tsv", "nomenclature/11_hpa_rna_celline_description.tsv"],
            "HPA description",
        )
        geo_expr_path = find_required_file(self.raw_data_dir, ["gene expression/3_GEOexpression.txt"], "GEO expression")
        geo_info_path = find_required_file(self.raw_data_dir, ["nomenclature/10_GEOInfo.txt"], "GEO info")

        cellosaurus = pd.read_csv(cellosaurus_path, low_memory=False)
        self.cell_resolver = dm.CellLineIDResolver(cellosaurus, self.sample_info, verbose=True)
        self.gene_resolver = dm.GeneIDResolver().build_from_hpa_file(hpa_expr_path, verbose=True)

        self.merger = (
            dm.MultiOmicsMerger(self.sample_info["DepMap_ID"], self.cell_resolver, self.gene_resolver)
            .register(dm.DepMapRNASource(self.loader, version="DepMap"))
            .register(dm.HPARNASource(hpa_expr_path, hpa_desc_path, value_col="nTPM", version="HPA"))
            .register(dm.GEOWideSource(geo_expr_path, geo_info_path, version="GEO", log_transform=True))
            .register(dm.MatrixProteinSource(self.loader, version="CCLE-Gygi"))
            .register(dm.MutationSource(self.loader, version="DepMap"))
            .register(dm.FusionSource(self.loader, version="DepMap"))
        )
        self.annotations = self._prepare_annotations(self.sample_info)
        print(f"Dynamic data layer ready. Cell-line annotations = {self.annotations.shape}")

    @staticmethod
    def _prepare_annotations(sample_info: pd.DataFrame) -> pd.DataFrame:
        if "DepMap_ID" not in sample_info.columns:
            raise KeyError("sample_info must contain DepMap_ID")
        out = sample_info.drop_duplicates("DepMap_ID").copy()
        for c in ["stripped_cell_line_name", "lineage", "primary_disease"]:
            if c not in out.columns:
                out[c] = np.nan
        if out["stripped_cell_line_name"].isna().all():
            for c in ["cell_line_name", "CCLE_Name", "ModelID"]:
                if c in out.columns:
                    out["stripped_cell_line_name"] = out[c]
                    break
        return out

    def _cache_path(self, gene: str) -> Path:
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", gene.upper().strip())
        return self.cache_dir / f"{safe}.parquet"

    def _ensure_schema(self, table: pd.DataFrame, gene: str) -> pd.DataFrame:
        table = table.copy()
        if "DepMap_ID" not in table.columns:
            raise KeyError("Gene table does not contain DepMap_ID")
        table["gene"] = gene
        for c in EXPECTED_EVIDENCE_COLUMNS:
            if c not in table.columns:
                table[c] = np.nan
        for c in ["has_rna", "has_protein", "has_mutation", "has_fusion"]:
            if c not in table.columns:
                table[c] = False
            table[c] = table[c].fillna(False).astype(bool)
        if "rna_consistency" not in table.columns:
            table["rna_consistency"] = np.nan
        if "data_completeness" not in table.columns:
            table["data_completeness"] = np.nan
        return table

    def _build_gene_table_from_long(self, long: pd.DataFrame, gene: str) -> pd.DataFrame:
        gene = gene.upper().strip()
        g = long.loc[long["gene_symbol"].astype(str).str.upper() == gene].copy()
        base = pd.DataFrame({"DepMap_ID": list(self.merger.all_cell_lines)})

        for src in sorted(g["source"].dropna().unique()):
            sub = g[g["source"] == src].groupby("DepMap_ID")[["value_raw", "value_std"]].mean()
            base = base.merge(
                sub.rename(columns={"value_raw": f"{src}__raw", "value_std": f"{src}__std"}),
                on="DepMap_ID", how="left"
            )

        for layer in self.OMICS_LAYERS:
            srcs = [s.source for s in self.merger._sources if s.layer == layer]
            raw_cols = [f"{s}__raw" for s in srcs if f"{s}__raw" in base.columns]
            base[f"has_{layer}"] = base[raw_cols].notna().any(axis=1) if raw_cols else False

        rna_cols = [f"{s.source}__std" for s in self.merger._sources if s.layer == "rna" and f"{s.source}__std" in base.columns]
        if len(rna_cols) >= 2:
            spread = base[rna_cols].std(axis=1, ddof=0)
            base["rna_consistency"] = (1 - spread.clip(0, 1)).where(base[rna_cols].notna().sum(axis=1) >= 2)
        else:
            base["rna_consistency"] = np.nan

        present = [f"has_{l}" for l in self.OMICS_LAYERS if any(s.layer == l for s in self.merger._sources)]
        base["data_completeness"] = base[present].mean(axis=1) if present else 0.0
        return self._ensure_schema(base, gene)

    @staticmethod
    def _gene_has_evidence(table: pd.DataFrame) -> bool:
        cols = ["depmap_rna__raw", "hpa_rna__raw", "geo_rna__raw", "ccle_gygi_protein__raw", "depmap_mutation__raw", "depmap_fusion__raw"]
        return any(c in table.columns and table[c].notna().any() for c in cols)

    @staticmethod
    def _protein_has_evidence(table: pd.DataFrame) -> bool:
        return (
            PROTEIN_RAW_COLUMN in table.columns
            and table[PROTEIN_RAW_COLUMN].notna().any()
        )

    def load_gene_tables(self, genes: Iterable[str]) -> Dict[str, pd.DataFrame]:
        requested = []
        for gene in genes:
            if gene:
                g = str(gene).upper().strip()
                if g and g not in requested:
                    requested.append(g)

        tables, missing = {}, []
        for gene in requested:
            p = self._cache_path(gene)
            if p.exists() and not self.refresh_cache:
                print(f"[cache] Loading {gene}: {p}")
                tables[gene] = self._ensure_schema(pd.read_parquet(p), gene)
            else:
                missing.append(gene)

        if missing:
            print("\nOn-demand raw-data query for: " + ", ".join(missing))
            print("Only requested genes are merged; no all-gene master table is created.")
            long = self.merger.build_long_table(genes=set(missing))
            for gene in missing:
                table = self._build_gene_table_from_long(long, gene)
                tables[gene] = table
                if self._gene_has_evidence(table):
                    p = self._cache_path(gene)
                    table.to_parquet(p, index=False)
                    print(f"[cache] Saved {gene}: {p}")
                else:
                    print(f"[warning] No evidence found for {gene}")
        return tables

    def _disease_filter(self, df: pd.DataFrame, aliases: Sequence[str]) -> pd.DataFrame:
        cols = [c for c in ["primary_disease", "lineage"] if c in df.columns]
        if not cols:
            raise KeyError("No disease/lineage metadata columns found")
        text = df[cols].fillna("").astype(str).agg(" ".join, axis=1).map(normalise_disease_key)
        mask = text.map(lambda value: any(alias in value for alias in aliases))
        out = df.loc[mask].copy()
        out["matchedDiseaseText"] = text.loc[mask]
        return out

    @staticmethod
    def _expression_map(table: Optional[pd.DataFrame], candidate_ids: Sequence[str]) -> Dict[str, Optional[float]]:
        if table is None:
            return {}
        ids = set(map(str, candidate_ids))
        sub = table.loc[table["DepMap_ID"].astype(str).isin(ids)]
        return {str(r["DepMap_ID"]): mean_available(r, RNA_STD_COLUMNS) for _, r in sub.iterrows()}

    def fetch_candidate_evidence(
        self,
        target_gene: Optional[str],
        target_protein: Optional[str],
        disease_aliases: List[str],
        exclusion_gene: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """
        Build candidates for three supported input modes.

        GENE_MULTIOMICS:
            target_gene is provided, target_protein is blank.
            RNA + protein for the same gene are used, exactly like the previous model.

        PROTEIN_ONLY:
            target_protein is provided, target_gene is blank.
            Protein is the biological ranking evidence.
            RNA/mutation/fusion for the corresponding coding-gene symbol are supporting context only.

        COMBINED:
            both target_gene and target_protein are provided.
            RNA comes from target_gene; protein comes from target_protein.
            The two symbols may be the same or different.
        """
        target_gene = target_gene.upper().strip() if target_gene else None
        target_protein = target_protein.upper().strip() if target_protein else None
        exclusion_gene = exclusion_gene.upper().strip() if exclusion_gene else None
        query_mode = determine_query_mode(target_gene, target_protein)

        genes_to_load: List[str] = []
        for symbol in [target_gene, target_protein, exclusion_gene]:
            if symbol and symbol not in genes_to_load:
                genes_to_load.append(symbol)

        tables = self.load_gene_tables(genes_to_load)

        gene_table = tables.get(target_gene) if target_gene else None

        # If the user does not explicitly provide a protein target,
        # use the protein measurement for the target gene (legacy behaviour).
        protein_symbol = target_protein or target_gene
        protein_table = tables.get(protein_symbol) if protein_symbol else None

        if target_gene:
            if gene_table is None or not self._gene_has_evidence(gene_table):
                raise KeyError(
                    f"Target gene '{target_gene}' was not found in registered sources"
                )

        if query_mode == "PROTEIN_ONLY":
            if protein_table is None or not self._protein_has_evidence(protein_table):
                raise KeyError(
                    f"No CCLE-Gygi proteomics evidence was found for target protein "
                    f"'{target_protein}'. Protein-only recommendation cannot be produced."
                )

        # For disease filtering, use target-gene rows when available;
        # otherwise use the protein-symbol rows. Both are indexed by DepMap_ID.
        base_table = gene_table if gene_table is not None else protein_table
        if base_table is None:
            raise KeyError("No target table could be constructed.")

        base = base_table.merge(
            self.annotations,
            on="DepMap_ID",
            how="left",
            validate="many_to_one",
        )
        base = self._disease_filter(base, disease_aliases)

        if base.empty:
            return []

        # Fast row lookups for cases where target gene and target protein differ.
        gene_lookup = (
            gene_table.drop_duplicates("DepMap_ID").set_index("DepMap_ID")
            if gene_table is not None else None
        )
        protein_lookup = (
            protein_table.drop_duplicates("DepMap_ID").set_index("DepMap_ID")
            if protein_table is not None else None
        )

        ex_table = None
        if exclusion_gene:
            ex_table = tables.get(exclusion_gene)
            if ex_table is not None and not self._gene_has_evidence(ex_table):
                print(
                    f"Warning: exclusion gene {exclusion_gene} has no evidence; penalty = 0"
                )
                ex_table = None

        ex_map = self._expression_map(
            ex_table,
            base["DepMap_ID"].astype(str).tolist(),
        )

        rows: List[Dict[str, Any]] = []

        for _, base_row in base.iterrows():
            depmap_id = str(base_row["DepMap_ID"])

            # RNA / mutation / fusion row:
            # - target gene if supplied
            # - otherwise the coding-gene symbol corresponding to target protein
            if gene_lookup is not None and depmap_id in gene_lookup.index:
                gene_row = gene_lookup.loc[depmap_id]
            elif protein_lookup is not None and depmap_id in protein_lookup.index:
                gene_row = protein_lookup.loc[depmap_id]
            else:
                gene_row = base_row

            # Protein row comes from explicit target protein when supplied,
            # otherwise from the target-gene table.
            if protein_lookup is not None and depmap_id in protein_lookup.index:
                protein_row = protein_lookup.loc[depmap_id]
            else:
                protein_row = base_row

            rna_value = mean_available(gene_row, RNA_STD_COLUMNS)
            n_rna = count_available(gene_row, RNA_STD_COLUMNS)

            protein_value = (
                float(protein_row[PROTEIN_STD_COLUMN])
                if PROTEIN_STD_COLUMN in protein_row.index
                and pd.notna(protein_row.get(PROTEIN_STD_COLUMN))
                else None
            )

            # Gene/combined mode needs RNA or protein.
            # Protein-only mode requires actual protein evidence.
            if query_mode == "PROTEIN_ONLY":
                if protein_value is None:
                    continue
            elif n_rna == 0 and protein_value is None:
                continue

            mutation_raw = gene_row.get(MUTATION_RAW_COLUMN)
            fusion_raw = gene_row.get(FUSION_RAW_COLUMN)
            hi_raw = gene_row.get(MUTATION_HIGH_IMPACT_RAW_COLUMN)

            has_mut = pd.notna(mutation_raw) and safe_float(mutation_raw) > 0
            has_fusion = pd.notna(fusion_raw) and safe_float(fusion_raw) > 0
            has_hi = pd.notna(hi_raw) and safe_float(hi_raw) > 0

            name = None
            for c in [
                "cellosaurus_name",
                "stripped_cell_line_name",
                "cell_line_name",
                "CCLE_Name",
            ]:
                if c in base_row.index and pd.notna(base_row.get(c)) and str(base_row.get(c)).strip():
                    name = str(base_row.get(c)).strip()
                    break
            name = name or depmap_id

            associated_gene = target_gene or target_protein
            candidate = {
                "DepMap_ID": depmap_id,
                "cellLine": name,
                "lineage": str(base_row.get("lineage")) if pd.notna(base_row.get("lineage")) else None,
                "disease": str(base_row.get("primary_disease")) if pd.notna(base_row.get("primary_disease")) else None,

                "queryMode": query_mode,
                "targetGene": target_gene,
                "targetProtein": target_protein,
                "proteinSymbolUsed": protein_symbol,

                # RNA is target-gene RNA in gene/combined modes.
                # In protein-only mode it is supplementary RNA for the corresponding symbol.
                "rnaExpr": rna_value,
                "protExpr": protein_value,
                "exclusionExpr": ex_map.get(depmap_id),
                "nRna": n_rna,
                "nProt": int(protein_value is not None),
                "nExclusion": int(ex_map.get(depmap_id) is not None),

                "hasDepMapRNA": pd.notna(gene_row.get("depmap_rna__std")),
                "hasHpaRNA": pd.notna(gene_row.get("hpa_rna__std")),
                "hasGeoRNA": pd.notna(gene_row.get("geo_rna__std")),
                "hasProteomics": pd.notna(protein_row.get(PROTEIN_STD_COLUMN)),

                "depmapRnaRaw": None if pd.isna(gene_row.get("depmap_rna__raw")) else float(gene_row.get("depmap_rna__raw")),
                "hpaRnaRaw": None if pd.isna(gene_row.get("hpa_rna__raw")) else float(gene_row.get("hpa_rna__raw")),
                "geoRnaRaw": None if pd.isna(gene_row.get("geo_rna__raw")) else float(gene_row.get("geo_rna__raw")),
                "proteinRaw": None if pd.isna(protein_row.get(PROTEIN_RAW_COLUMN)) else float(protein_row.get(PROTEIN_RAW_COLUMN)),

                "rnaConsistency": None if pd.isna(gene_row.get("rna_consistency")) else float(gene_row.get("rna_consistency")),
                "mergedDataCompleteness": None if pd.isna(gene_row.get("data_completeness")) else float(gene_row.get("data_completeness")),

                "hasMutation": has_mut,
                "hasFusion": has_fusion,
                "hasHighImpactMutation": has_hi,
                "mutationRaw": None if pd.isna(mutation_raw) else float(mutation_raw),
                "mutationHighImpactRaw": None if pd.isna(hi_raw) else float(hi_raw),
                "fusionRaw": None if pd.isna(fusion_raw) else float(fusion_raw),
                "fusionConfidenceRaw": (
                    None
                    if pd.isna(gene_row.get(FUSION_CONFIDENCE_RAW_COLUMN))
                    else float(gene_row.get(FUSION_CONFIDENCE_RAW_COLUMN))
                ),
            }

            candidate["alterationStatus"] = alteration_status(
                has_mut,
                has_fusion,
                has_hi,
            )
            candidate["targetGeneLabel"] = gene_label(
                associated_gene or "",
                has_mut,
                has_fusion,
                has_hi,
            )

            rows.append(candidate)

        print("Detected dynamic evidence mapping:")
        if target_gene:
            print(f"- Target-gene RNA: {target_gene} from DepMap + HPA + GEO")
        if protein_symbol:
            print(f"- Target protein: {protein_symbol} from CCLE-Gygi")
        if target_gene or target_protein:
            print(f"- Mutation/Fusion context: {target_gene or target_protein} from DepMap")
        if exclusion_gene:
            print(f"- Exclusion penalty: {exclusion_gene} RNA")
        print(f"- Query mode: {query_mode}")

        return rows

    @staticmethod
    def fetch_evidence_trace(row: Dict[str, Any]) -> List[Dict[str, Any]]:
        trace: List[Dict[str, Any]] = []
        query_mode = row.get("queryMode")

        # Protein-only mode: protein is primary evidence.
        if query_mode == "PROTEIN_ONLY" and row.get("proteinRaw") is not None:
            trace.append({
                "dataset": "CCLE_Gygi_Proteomics",
                "type": "PRIMARY_target_protein_expression",
                "value": row["proteinRaw"],
                "unit": "log2 MS intensity",
            })

        # RNA is primary in gene/combined modes and supplementary in protein-only mode.
        rna_prefix = "SUPPORTING" if query_mode == "PROTEIN_ONLY" else "PRIMARY"
        if row.get("depmapRnaRaw") is not None:
            trace.append({
                "dataset": "DepMap_RNAseq",
                "type": f"{rna_prefix}_rna_expression",
                "value": row["depmapRnaRaw"],
                "unit": "log2(TPM+1)",
            })
        if row.get("hpaRnaRaw") is not None:
            trace.append({
                "dataset": "HPA_RNAseq",
                "type": f"{rna_prefix}_rna_expression",
                "value": row["hpaRnaRaw"],
                "unit": "log2(nTPM+1)",
            })
        if row.get("geoRnaRaw") is not None:
            trace.append({
                "dataset": "GEO_RNA",
                "type": f"{rna_prefix}_rna_expression",
                "value": row["geoRnaRaw"],
                "unit": "log2(intensity+1)",
            })

        # In gene/combined modes protein is also primary multi-omics evidence.
        if query_mode != "PROTEIN_ONLY" and row.get("proteinRaw") is not None:
            trace.append({
                "dataset": "CCLE_Gygi_Proteomics",
                "type": "PRIMARY_protein_expression",
                "value": row["proteinRaw"],
                "unit": "log2 MS intensity",
            })

        mutation_prefix = "SUPPORTING" if query_mode == "PROTEIN_ONLY" else "PRIMARY"
        if row.get("hasMutation"):
            trace.append({
                "dataset": "DepMap_Mutation",
                "type": (
                    f"{mutation_prefix}_high_impact_mutation"
                    if row.get("hasHighImpactMutation")
                    else f"{mutation_prefix}_mutation"
                ),
                "value": row.get("mutationRaw"),
                "unit": "event count",
            })

        if row.get("hasFusion"):
            trace.append({
                "dataset": "DepMap_Fusion",
                "type": f"{mutation_prefix}_fusion",
                "value": row.get("fusionRaw"),
                "unit": "event count",
            })
            if row.get("fusionConfidenceRaw") is not None:
                trace.append({
                    "dataset": "DepMap_Fusion",
                    "type": f"{mutation_prefix}_fusion_confidence",
                    "value": row.get("fusionConfidenceRaw"),
                    "unit": "source confidence",
                })

        return trace

    @staticmethod
    def _top_numeric_features(row: Any, top_n: int) -> List[Dict[str, Any]]:
        if isinstance(row, pd.DataFrame):
            s = row.apply(pd.to_numeric, errors="coerce").mean(axis=0)
        else:
            s = pd.to_numeric(pd.Series(row), errors="coerce")
        s = s.replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False).head(top_n)
        return [{"feature": str(k), "value": float(v)} for k, v in s.items()]

    def get_supplementary_context(self, depmap_id: str, top_n: int = 5) -> Dict[str, Any]:
        """Scheme A: Files 12/13/14 are output-only cell-line context, never scoring inputs."""
        context = {"global_signatures": {}, "top_mirna": [], "top_metabolites": [], "notes": []}

        try:
            sigs = self.loader.get_global_signatures()
            hit = sigs.loc[sigs["DepMap_ID"].astype(str) == str(depmap_id)]
            if not hit.empty:
                r = hit.iloc[0]
                for c in ["MSIScore", "LoHFraction", "WGD", "CIN", "Ploidy", "Aneuploidy"]:
                    if c in r.index and pd.notna(r[c]):
                        v = r[c].item() if isinstance(r[c], np.generic) else r[c]
                        context["global_signatures"][c] = v
            else:
                context["notes"].append("No File 14 global-signature record for this cell line.")
        except Exception as exc:
            context["notes"].append(f"File 14 unavailable: {type(exc).__name__}: {exc}")

        try:
            mirna = self.loader._load_mirna_matrix()
            if str(depmap_id) in mirna.index:
                context["top_mirna"] = self._top_numeric_features(mirna.loc[str(depmap_id)], top_n)
            else:
                context["notes"].append("No File 13 miRNA profile for this cell line.")
        except Exception as exc:
            context["notes"].append(f"File 13 unavailable: {type(exc).__name__}: {exc}")

        try:
            metab = self.loader._load_metabolomics_matrix()
            if str(depmap_id) in metab.index:
                context["top_metabolites"] = self._top_numeric_features(metab.loc[str(depmap_id)], top_n)
            else:
                context["notes"].append("No File 12 metabolomics profile for this cell line.")
        except Exception as exc:
            context["notes"].append(f"File 12 unavailable: {type(exc).__name__}: {exc}")

        return context

    @staticmethod
    def fetch_similar_cell_lines(scored_rows: List[Dict[str, Any]], top_k: int = 5) -> List[Dict[str, Any]]:
        if len(scored_rows) <= 1:
            return []
        top, out = scored_rows[0], []
        for row in scored_rows[1:]:
            diffs = []
            if top.get("rnaScore") is not None and row.get("rnaScore") is not None:
                diffs.append(abs(float(top["rnaScore"]) - float(row["rnaScore"])))
            if top.get("proteinScore") is not None and row.get("proteinScore") is not None:
                diffs.append(abs(float(top["proteinScore"]) - float(row["proteinScore"])))
            sim = max(0.0, min(1.0, 1.0 - sum(diffs) / len(diffs))) if diffs else 0.0
            out.append({
                "alternativeCellLine": row.get("cellLine"),
                "similarityScore": round(sim, 4),
                "finalScore": row.get("finalScore"),
                "targetGeneLabel": row.get("targetGeneLabel"),
            })
        out.sort(key=lambda x: (safe_float(x["similarityScore"]), safe_float(x["finalScore"])), reverse=True)
        return out[:top_k]

# ============================== SCORING ==============================
def scale_available_values(rows: List[Dict[str, Any]], value_key: str, count_key: str) -> List[Optional[float]]:
    idx = [i for i, r in enumerate(rows) if safe_float(r.get(count_key)) > 0 and r.get(value_key) is not None]
    out: List[Optional[float]] = [None] * len(rows)
    if not idx:
        return out
    scaled = min_max_scale([safe_float(rows[i][value_key]) for i in idx])
    for i, v in zip(idx, scaled):
        out[i] = v
    return out


def score_candidates(
    rows: List[Dict[str, Any]],
    query_mode: str,
) -> List[Dict[str, Any]]:
    if not rows:
        return []

    rna_scores = scale_available_values(rows, "rnaExpr", "nRna")
    protein_scores = scale_available_values(rows, "protExpr", "nProt")
    exclusion_scores = scale_available_values(rows, "exclusionExpr", "nExclusion")
    scored_rows: List[Dict[str, Any]] = []

    for i, row in enumerate(rows):
        has_rna = rna_scores[i] is not None
        has_protein = protein_scores[i] is not None
        has_exclusion = exclusion_scores[i] is not None

        rna = float(rna_scores[i]) if has_rna else 0.0
        protein = float(protein_scores[i]) if has_protein else None
        penalty = (
            MAX_EXCLUSION_PENALTY * float(exclusion_scores[i])
            if has_exclusion
            else 0.0
        )

        # Confidence still uses data availability / cross-layer agreement.
        # In protein-only mode RNA can support confidence, but it does not
        # directly contribute to biological ranking.
        completeness = (float(has_rna) + float(has_protein)) / 2.0
        support = (
            float(bool(row.get("hasDepMapRNA")))
            + float(bool(row.get("hasHpaRNA")))
            + float(bool(row.get("hasGeoRNA")))
            + float(bool(row.get("hasProteomics")))
        ) / 4.0

        if has_rna and has_protein:
            rp_consistency = 1.0 - abs(rna - float(protein))
        elif has_rna or has_protein:
            rp_consistency = 0.5
        else:
            rp_consistency = 0.0

        rp_consistency = max(0.0, min(1.0, rp_consistency))

        confidence = (
            0.40 * completeness
            + 0.35 * support
            + 0.25 * rp_consistency
        )
        confidence = max(0.0, min(1.0, confidence))

        if query_mode == "PROTEIN_ONLY":
            # Explicit protein query: do not silently substitute RNA if protein is absent.
            if not has_protein:
                continue

            biological = float(protein)
            final = (
                PROTEIN_ONLY_PROTEIN_WEIGHT * biological
                + PROTEIN_ONLY_CONFIDENCE_WEIGHT * confidence
                - penalty
            )
        else:
            # Legacy gene/multi-omics scoring.
            available_weight = 0.0
            weighted = 0.0

            if has_rna:
                available_weight += RNA_WEIGHT
                weighted += RNA_WEIGHT * rna

            if has_protein:
                available_weight += PROTEIN_WEIGHT
                weighted += PROTEIN_WEIGHT * float(protein)

            if available_weight == 0:
                continue

            biological = weighted / available_weight
            final = (
                (RNA_WEIGHT + PROTEIN_WEIGHT) * biological
                + CONFIDENCE_WEIGHT * confidence
                - penalty
            )

        final = max(0.0, min(1.0, final))

        s = dict(row)
        s.update({
            "rnaScore": round(rna, 4) if has_rna else None,
            "proteinScore": round(float(protein), 4) if has_protein else None,
            "biologicalScore": round(biological, 4),
            "exclusionPenalty": round(penalty, 4),
            "completenessScore": round(completeness, 4),
            "sourceSupportScore": round(support, 4),
            "rnaProteinConsistencyScore": round(rp_consistency, 4),
            "confidenceScore": round(confidence, 4),
            "confidenceLevel": confidence_level(confidence),
            "finalScore": round(final, 4),
            "recommendationLevel": recommendation_level(final),
        })
        scored_rows.append(s)

    scored_rows.sort(
        key=lambda x: (x["finalScore"], x["confidenceScore"]),
        reverse=True,
    )

    for rank, row in enumerate(scored_rows, start=1):
        row["rank"] = rank

    return scored_rows

# ============================== OUTPUT ==============================
def build_reason(
    row: Dict[str, Any],
    target_gene: Optional[str],
    target_protein: Optional[str],
    exclusion_gene: Optional[str],
    disease: str,
    query_mode: str,
) -> List[str]:
    reasons = [
        f"The cell line passed the disease-specific hard filter for {disease}."
    ]

    if query_mode == "PROTEIN_ONLY":
        protein = row.get("proteinScore")
        if protein is not None and protein >= 0.70:
            reasons.append(
                f"Strong protein expression evidence for {target_protein}."
            )
        elif protein is not None and protein >= 0.40:
            reasons.append(
                f"Moderate protein expression evidence for {target_protein}."
            )
        else:
            reasons.append(
                f"Limited protein expression evidence for {target_protein}."
            )

        if row.get("rnaScore") is not None:
            reasons.append(
                "RNA evidence is available as supporting context but is not "
                "used directly in the protein-only biological score."
            )
    else:
        rna = row.get("rnaScore")
        if rna is not None and rna >= 0.70:
            reasons.append(f"Strong RNA expression evidence for {target_gene}.")
        elif rna is not None and rna >= 0.40:
            reasons.append(f"Moderate RNA expression evidence for {target_gene}.")
        else:
            reasons.append(f"Limited RNA expression evidence for {target_gene}.")

        if row.get("proteinScore") is not None:
            reasons.append(
                f"Supporting protein-level evidence is available for "
                f"{target_protein or target_gene}."
            )
        else:
            reasons.append(
                "Protein evidence is missing; its weight was redistributed "
                "and confidence was reduced."
            )

    if row.get("alterationStatus") != "NO ALTERATION":
        symbol = target_gene or target_protein
        reasons.append(
            f"Genomic alteration detected for {symbol}: "
            f"{row.get('alterationStatus')}."
        )

    if exclusion_gene:
        reasons.append(
            f"The exclusion gene {exclusion_gene} does not cause a major penalty."
            if row.get("exclusionPenalty", 0.0) <= 0.05
            else (
                f"A penalty is applied because {exclusion_gene} shows "
                "relatively high RNA expression evidence."
            )
        )

    return reasons


def print_report(
    target_gene: Optional[str],
    target_protein: Optional[str],
    exclusion_gene: Optional[str],
    disease: str,
    query_mode: str,
    top: Dict[str, Any],
    alternatives: List[Dict[str, Any]],
    trace: List[Dict[str, Any]],
    supplementary: Dict[str, Any],
) -> None:
    print("\n" + "=" * 80)
    print("CellLineSelector Dynamic Multi-Omics Recommendation Demo")
    print("=" * 80)

    print("\nInput:")
    if target_gene:
        print(f"Target gene: {target_gene}")
    if target_protein:
        print(f"Target protein: {target_protein}")
    if exclusion_gene:
        print(f"Exclusion gene: {exclusion_gene}")
    print(f"Disease hard filter: {disease}")
    print(f"Query mode: {query_mode}")

    print("\nOutput:")
    print(f"Recommended cell line: {top.get('cellLine')}")
    print(f"Recommendation level: {top.get('recommendationLevel')}")
    print(f"Final score: {top.get('finalScore'):.2f} / 1.00")
    print(f"Confidence: {top.get('confidenceLevel')}")
    print(f"Confidence score: {top.get('confidenceScore'):.2f} / 1.00")

    if query_mode == "PROTEIN_ONLY":
        print(f"Target protein: {target_protein}")
        print(
            f"Associated gene status: "
            f"{top.get('targetGeneLabel', target_protein)}"
        )
    else:
        print(
            f"Target gene status: "
            f"{top.get('targetGeneLabel', target_gene)}"
        )

    print("\nScore breakdown:")
    if query_mode == "PROTEIN_ONLY":
        print(
            f"- Protein expression score: "
            f"{top.get('proteinScore') if top.get('proteinScore') is not None else 'Missing'}"
        )
        print(
            f"- Supporting RNA score (not directly weighted): "
            f"{top.get('rnaScore')}"
        )
        print(f"- Biological score (protein): {top.get('biologicalScore')}")
        print(
            f"- Protein-only weights: "
            f"protein={PROTEIN_ONLY_PROTEIN_WEIGHT:.2f}, "
            f"confidence={PROTEIN_ONLY_CONFIDENCE_WEIGHT:.2f}"
        )
    else:
        print(f"- RNA expression score: {top.get('rnaScore')}")
        print(
            f"- Protein score: "
            f"{top.get('proteinScore') if top.get('proteinScore') is not None else 'Missing / redistributed'}"
        )
        print(f"- Biological score: {top.get('biologicalScore')}")
        print(
            f"- Multi-omics weights: RNA={RNA_WEIGHT:.2f}, "
            f"protein={PROTEIN_WEIGHT:.2f}, "
            f"confidence={CONFIDENCE_WEIGHT:.2f}"
        )

    print(f"- Exclusion penalty: {top.get('exclusionPenalty')}")
    print(f"- Completeness score: {top.get('completenessScore')}")
    print(f"- Source support score: {top.get('sourceSupportScore')}")
    if top.get("rnaConsistency") is not None:
        print(
            f"- RNA cross-source consistency: "
            f"{top.get('rnaConsistency'):.4f}"
        )

    print("\nMain reason:")
    for x in build_reason(
        top,
        target_gene,
        target_protein,
        exclusion_gene,
        disease,
        query_mode,
    ):
        print(f"- {x}")

    print("\nEvidence trace:")
    if trace:
        for e in trace:
            print(
                f"- {e['dataset']}: {e['type']}, "
                f"value={e['value']}, unit={e['unit']}"
            )
    else:
        print("- No detailed target evidence found.")

    print("\nSupplementary cell-line omics context (NOT used in scoring):")
    sigs = supplementary.get("global_signatures", {})
    print(
        "- File 14 Global signatures: "
        + (
            ", ".join(f"{k}={v}" for k, v in sigs.items())
            if sigs else "not available"
        )
    )

    mirna = supplementary.get("top_mirna", [])
    if mirna:
        print("- File 13 miRNA context (top measured features, RPM):")
        for item in mirna:
            print(f"  - {item['feature']}: {item['value']:.4f}")
    else:
        print("- File 13 miRNA context: not available")

    metab = supplementary.get("top_metabolites", [])
    if metab:
        print(
            "- File 12 metabolomics context "
            "(top measured features, log2-normalised intensity):"
        )
        for item in metab:
            print(f"  - {item['feature']}: {item['value']:.4f}")
    else:
        print("- File 12 metabolomics context: not available")

    for note in supplementary.get("notes", []):
        print(f"- Supplementary-data note: {note}")

    print("\nData gaps:")
    gaps: List[str] = []

    if query_mode == "PROTEIN_ONLY":
        if not top.get("hasProteomics"):
            gaps.append("CCLE-Gygi target-protein evidence is missing.")
        if not top.get("hasDepMapRNA"):
            gaps.append("Supporting DepMap RNA-seq evidence is missing.")
        if not top.get("hasHpaRNA"):
            gaps.append("Supporting HPA RNA-seq evidence is missing.")
        if not top.get("hasGeoRNA"):
            gaps.append("Supporting GEO RNA evidence is missing.")
    else:
        if not top.get("hasDepMapRNA"):
            gaps.append("DepMap RNA-seq evidence is missing.")
        if not top.get("hasHpaRNA"):
            gaps.append("HPA RNA-seq evidence is missing.")
        if not top.get("hasGeoRNA"):
            gaps.append("GEO RNA expression evidence is missing.")
        if not top.get("hasProteomics"):
            gaps.append(
                f"CCLE-Gygi protein evidence is missing for "
                f"{target_protein or target_gene}."
            )

    if gaps:
        for g in gaps:
            print(f"- {g}")
    else:
        print("- No major expression-data gap detected for this query.")

    print("\nSimilar alternative cell lines:")
    if alternatives:
        for i, alt in enumerate(alternatives, 1):
            print(
                f"{i}. {alt['alternativeCellLine']} — "
                f"similarity score: {alt['similarityScore']:.2f} | "
                f"{alt.get('targetGeneLabel', target_gene or target_protein)}"
            )
    else:
        print("- No similar alternatives found.")

    print("\nNote:")
    print(f"- All recommendations are restricted to the disease filter: {disease}.")
    print(
        "- Target evidence is queried dynamically from original source files; "
        "no all-gene master table is created."
    )
    print(
        "- Files 12/13/14 are supplementary context only and do not affect ranking."
    )
    print(
        "- Scheme A does not claim miRNA/metabolite features are target-specific "
        "because no regulator/pathway mapping is applied."
    )
    print(
        "- Mutation/fusion are reported as evidence but are not automatically "
        "beneficial or harmful."
    )
    if query_mode == "PROTEIN_ONLY":
        print(
            "- In protein-only mode, protein is the biological ranking evidence; "
            "RNA is retained only as supporting context/confidence evidence."
        )


def save_ranked_results(rows: List[Dict[str, Any]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fieldnames = [
        "rank", "DepMap_ID", "cellLine", "lineage", "disease",
        "queryMode", "targetGene", "targetProtein", "proteinSymbolUsed",
        "targetGeneLabel", "alterationStatus",
        "hasMutation", "hasHighImpactMutation", "mutationRaw", "mutationHighImpactRaw",
        "hasFusion", "fusionRaw", "fusionConfidenceRaw",
        "rnaExpr", "protExpr", "exclusionExpr", "rnaScore", "proteinScore", "biologicalScore",
        "exclusionPenalty", "completenessScore", "sourceSupportScore", "rnaProteinConsistencyScore",
        "rnaConsistency", "mergedDataCompleteness", "confidenceScore", "confidenceLevel",
        "finalScore", "recommendationLevel", "hasDepMapRNA", "hasHpaRNA", "hasGeoRNA", "hasProteomics",
        "depmapRnaRaw", "hpaRnaRaw", "geoRnaRaw", "proteinRaw",
        "globalSignaturesContext", "topMiRNAContext", "topMetaboliteContext",
    ]
    with output_path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

# ============================== CLI ==============================
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--target_gene", default=None)
    p.add_argument("--target_protein", default=None)
    p.add_argument("--exclusion_gene", default=None)
    p.add_argument("--disease", "--context", dest="disease", default=None)
    p.add_argument("--raw_data_dir", default=str(DEFAULT_RAW_DATA_DIR))
    p.add_argument("--cache_dir", default=str(DEFAULT_CACHE_DIR))
    p.add_argument("--refresh_cache", action="store_true")
    p.add_argument("--top_n", type=int, default=10)
    p.add_argument("--top_k_alternatives", type=int, default=5)
    p.add_argument("--supplementary_top_n", type=int, default=5)
    p.add_argument("--output_dir", default=DEFAULT_OUTPUT_DIR)
    p.add_argument(
        "--no_agent", "--no-agent", action="store_true",
        help="Return the deterministic recommendation without the explanation Agent.",
    )
    p.add_argument(
        "--non_interactive", "--non-interactive", action="store_true",
        help="Use only command-line values and do not prompt for omitted optional inputs.",
    )
    p.add_argument(
        "--agent_backend", "--agent-backend", default="scripted",
        choices=["scripted", "openai", "anthropic"],
        help="Writer used only for post-ranking explanations (default: offline scripted).",
    )
    p.add_argument(
        "--agent_model", "--agent-model", default=None,
        help="Optional hosted model id; ignored by the scripted backend.",
    )
    p.add_argument(
        "--agent_grounding", "--agent-grounding", default="l1l2",
        choices=["none", "l1l2", "full"],
        help="Explanation checks: citation+numeric by default; full also checks entities.",
    )
    return p.parse_args()


def ask(value: Optional[str], prompt: str, required: bool = False) -> Optional[str]:
    if value is not None and str(value).strip():
        return str(value).strip()
    while True:
        x = input(prompt).strip()
        if x:
            return x
        if not required:
            return None
        print("This field is required.")


def main() -> None:
    args = parse_args()

    # Gene and protein are both optional individually, but at least one is required.
    if not args.non_interactive:
        args.target_gene = ask(
            args.target_gene,
            "Target gene (press Enter to skip): ",
            False,
        )
        args.target_protein = ask(
            args.target_protein,
            "Target protein (press Enter to skip): ",
            False,
        )

    if not args.target_gene and not args.target_protein:
        print("At least one of Target gene or Target protein must be provided.")
        return

    if args.target_gene:
        args.target_gene = args.target_gene.upper().strip()
    if args.target_protein:
        args.target_protein = args.target_protein.upper().strip()

    if not args.non_interactive:
        args.exclusion_gene = ask(
            args.exclusion_gene,
            "Exclusion gene (press Enter to skip): ",
            False,
        )
    if args.exclusion_gene:
        args.exclusion_gene = args.exclusion_gene.upper().strip()

    if not args.non_interactive:
        args.disease = ask(
            args.disease,
            "Disease or tissue (required; used as a hard filter): ",
            True,
        )
    elif not args.disease or not str(args.disease).strip():
        raise ValueError("--disease is required in --non-interactive mode")

    query_mode = determine_query_mode(
        args.target_gene,
        args.target_protein,
    )
    aliases = build_disease_aliases(args.disease)

    try:
        print("\n" + "=" * 80)
        print("Dynamic multi-omics data initialization")
        print("=" * 80)

        data = DynamicMultiOmicsRecommender(
            Path(args.raw_data_dir),
            Path(args.cache_dir),
            args.refresh_cache,
        )

        print("\nReading target evidence from original multi-omics data...")
        print(f"Target gene: {args.target_gene or 'None'}")
        print(f"Target protein: {args.target_protein or 'None'}")
        print(f"Exclusion gene: {args.exclusion_gene or 'None'}")
        print(f"Disease hard filter: {args.disease}")
        print(f"Query mode: {query_mode}")
        print("Disease aliases: " + ", ".join(aliases))

        candidates = data.fetch_candidate_evidence(
            target_gene=args.target_gene,
            target_protein=args.target_protein,
            disease_aliases=aliases,
            exclusion_gene=args.exclusion_gene,
        )

        if not candidates:
            print(
                "No disease-matched candidates with the required target evidence "
                "were found."
            )
            return

        scored = score_candidates(
            candidates,
            query_mode=query_mode,
        )

        if not scored:
            if query_mode == "PROTEIN_ONLY":
                print(
                    "No candidates could be scored because no disease-matched "
                    "target-protein measurements were available."
                )
            else:
                print("No candidates could be scored.")
            return

        top_rows = scored[:max(1, args.top_n)]
        top = top_rows[0]

        alternatives = data.fetch_similar_cell_lines(
            scored,
            max(0, args.top_k_alternatives),
        )
        trace = data.fetch_evidence_trace(top)

        print("\nLoading supplementary context from Files 12/13/14...")
        supplementary = data.get_supplementary_context(
            top["DepMap_ID"],
            max(0, args.supplementary_top_n),
        )

        top["globalSignaturesContext"] = json.dumps(
            supplementary["global_signatures"],
            ensure_ascii=False,
        )
        top["topMiRNAContext"] = json.dumps(
            supplementary["top_mirna"],
            ensure_ascii=False,
        )
        top["topMetaboliteContext"] = json.dumps(
            supplementary["top_metabolites"],
            ensure_ascii=False,
        )

        print_report(
            target_gene=args.target_gene,
            target_protein=args.target_protein,
            exclusion_gene=args.exclusion_gene,
            disease=args.disease,
            query_mode=query_mode,
            top=top,
            alternatives=alternatives,
            trace=trace,
            supplementary=supplementary,
        )

        print("\nTop ranked cell lines:")
        for row in top_rows:
            marker = row.get(
                "targetGeneLabel",
                args.target_gene or args.target_protein,
            )
            print(
                f"{row['rank']}. {row['cellLine']} | "
                f"finalScore={row['finalScore']:.2f} | "
                f"confidence={row['confidenceScore']:.2f} | "
                f"level={row['recommendationLevel']} | "
                f"{marker}"
            )

        out_dir = Path(args.output_dir)
        if not out_dir.is_absolute():
            out_dir = Path(__file__).resolve().parent / out_dir

        out_path = out_dir / "dynamic_ranked_recommendations.csv"
        save_ranked_results(scored, out_path)

        print(f"\nSaved ranked results to: {out_path.resolve()}")

        # The Agent is deliberately downstream of score_candidates(). It receives
        # a read-only representation of the authoritative ranking and can only
        # produce explanatory claims. OutputAgent fingerprints the ranking before
        # and after generation and raises if rank, score, or level changes.
        if not args.no_agent:
            try:
                from src.agentic.output_agent import (
                    OutputAgent,
                    print_agent_output,
                    save_agent_output,
                )

                agent_request = {
                    "target_gene": args.target_gene,
                    "target_protein": args.target_protein,
                    "exclusion_gene": args.exclusion_gene,
                    "disease": args.disease,
                    "disease_aliases": aliases,
                    "query_mode": query_mode,
                    "top_n": len(top_rows),
                }
                output_agent = OutputAgent(
                    backend=args.agent_backend,
                    model=args.agent_model,
                    grounding_mode=args.agent_grounding,
                )
                agent_output = output_agent.run(
                    request=agent_request,
                    ranked_rows=top_rows,
                    evidence_trace=trace,
                    alternatives=alternatives,
                    supplementary=supplementary,
                )
                print_agent_output(agent_output)
                agent_path = out_dir / "agent_assisted_output.json"
                save_agent_output(agent_output, agent_path)
                print(f"Saved Agent explanation and audit trail to: {agent_path.resolve()}")
            except Exception as agent_exc:
                # Explanation failure must never invalidate or hide the scientific
                # result. The CSV and deterministic console report remain complete.
                print("\nAgent-assisted explanation was unavailable.")
                print(f"{type(agent_exc).__name__}: {agent_exc}")
                print("The deterministic ranking above is unchanged and remains the final result.")

    except Exception as exc:
        print("\nProgram stopped because dynamic loading or recommendation failed.")
        print(f"{type(exc).__name__}: {exc}")
        print(
            "Check that data_loader.py and data_merger.py are in ./src and "
            "that data_s3 contains the original files 1-14."
        )
        raise


if __name__ == "__main__":
    main()
