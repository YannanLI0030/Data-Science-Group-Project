"""
CellLineSelector - merged-data version

This script reads an already merged project dataset supplied through
--data_dir or the CELLLINESELECTOR_DATA_DIR environment variable.

Required files:
    master_table.csv
    cellline_annotations.csv

No Neo4j, AWS, boto3, s3fs, or data_loader.py is required.

Workflow:
1. Accept a target gene, optional exclusion gene, and optional disease/tissue.
2. Load the merged master table + cell-line annotations.
3. Apply disease/tissue as a HARD FILTER when supplied; otherwise run pan-cancer.
4. Combine DepMap RNA + HPA RNA + GEO RNA as target RNA evidence.
5. Use CCLE-Gygi protein as protein evidence.
6. Use the exclusion gene's RNA expression to calculate exclusion penalty.
7. Use a validated, serialisable scoring configuration (provisional A0 default):
       RNA 55% + Protein 30% + Confidence 15%
   with row-wise redistribution when protein is missing.
8. Annotate mutation/fusion status for the target gene:
       [MUTATION]
       [FUSION]
       [MIXED: MUTATION + FUSION]
9. Print the same recommendation-report sections as the previous program.
10. Save the full ranked table to CSV and a UI/RAG-facing JSON payload.

Recommended environment:
    numpy 1.26.4
    pandas 2.2.x / 2.3.x
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# CONFIG
# =============================================================================

DEFAULT_DATA_DIR = os.getenv("CELLLINESELECTOR_DATA_DIR")
MASTER_FILENAME = "master_table.csv"
ANNOTATION_FILENAME = "cellline_annotations.csv"

DEFAULT_TOP_N = 10
DEFAULT_TOP_K_ALTERNATIVES = 5
DEFAULT_OUTPUT_DIR = "results"
OUTPUT_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class ScoringConfig:
    """Versioned, serialisable scoring configuration shared with ablations."""

    name: str = "A0_team_baseline_provisional"
    rna_weight: float = 0.55
    protein_weight: float = 0.30
    confidence_weight: float = 0.15
    confidence_completeness: float = 0.40
    confidence_source_support: float = 0.35
    confidence_consistency: float = 0.25
    max_exclusion_penalty: float = 0.30

    def __post_init__(self) -> None:
        numeric_fields = [
            "rna_weight",
            "protein_weight",
            "confidence_weight",
            "confidence_completeness",
            "confidence_source_support",
            "confidence_consistency",
            "max_exclusion_penalty",
        ]
        for field_name in numeric_fields:
            value = float(getattr(self, field_name))
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"{field_name} must be a finite non-negative number")

        outer_total = self.rna_weight + self.protein_weight + self.confidence_weight
        if not math.isclose(outer_total, 1.0, abs_tol=1e-9):
            raise ValueError(f"outer scoring weights must sum to 1.0, got {outer_total}")

        confidence_total = (
            self.confidence_completeness
            + self.confidence_source_support
            + self.confidence_consistency
        )
        if not math.isclose(confidence_total, 1.0, abs_tol=1e-9):
            raise ValueError(
                "confidence component weights must sum to 1.0, "
                f"got {confidence_total}"
            )

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


DEFAULT_SCORING_CONFIG = ScoringConfig()

# Backwards-compatible aliases. New code should pass ScoringConfig explicitly.
RNA_WEIGHT = DEFAULT_SCORING_CONFIG.rna_weight
PROTEIN_WEIGHT = DEFAULT_SCORING_CONFIG.protein_weight
CONFIDENCE_WEIGHT = DEFAULT_SCORING_CONFIG.confidence_weight
MAX_EXCLUSION_PENALTY = DEFAULT_SCORING_CONFIG.max_exclusion_penalty

# The merged table contains three RNA evidence sources.
RNA_STD_COLUMNS = [
    "depmap_rna__std",
    "hpa_rna__std",
    "geo_rna__std",
]

RNA_RAW_COLUMNS = [
    "depmap_rna__raw",
    "hpa_rna__raw",
    "geo_rna__raw",
]

PROTEIN_STD_COLUMN = "ccle_gygi_protein__std"
PROTEIN_RAW_COLUMN = "ccle_gygi_protein__raw"

MUTATION_RAW_COLUMN = "depmap_mutation__raw"
MUTATION_HIGH_IMPACT_RAW_COLUMN = "depmap_mutation_highimpact__raw"
FUSION_RAW_COLUMN = "depmap_fusion__raw"
FUSION_CONFIDENCE_RAW_COLUMN = "depmap_fusion_confidence__raw"

# Disease/tissue aliases used by the hard filter.
DISEASE_ALIASES: Dict[str, List[str]] = {
    "brain": [
        "brain",
        "central nervous system",
        "cns",
        "glioma",
        "glioblastoma",
        "astrocytoma",
        "medulloblastoma",
    ],
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
    "blood": [
        "blood",
        "haematopoietic",
        "hematopoietic",
        "leukaemia",
        "leukemia",
        "lymphoma",
        "myeloma",
    ],
}


# =============================================================================
# Utilities
# =============================================================================

def safe_float(value: Any, default: float = 0.0) -> float:
    if value is None:
        return default
    try:
        if pd.isna(value):
            return default
    except Exception:
        pass
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def min_max_scale(values: Sequence[float]) -> List[float]:
    values = [float(v) for v in values]
    if not values:
        return []

    min_v = min(values)
    max_v = max(values)

    if max_v == min_v:
        return [1.0 if max_v > 0 else 0.0 for _ in values]

    return [(v - min_v) / (max_v - min_v) for v in values]


def normalise_text(text: Optional[str]) -> str:
    if text is None:
        return ""
    return str(text).lower().strip()


def normalise_disease_key(text: Optional[str]) -> str:
    value = normalise_text(text)
    value = value.replace("-", " ").replace("_", " ").replace("/", " ")
    value = re.sub(r"\s+", " ", value).strip()
    return value


def build_disease_aliases(disease: Optional[str]) -> List[str]:
    query = normalise_disease_key(disease)
    if not query:
        return []

    stop_words = {
        "cancer", "tumour", "tumor", "carcinoma", "disease",
        "cell", "cells", "cell line", "cell lines",
    }

    compact = " ".join(
        token for token in query.split()
        if token not in stop_words
    ).strip()

    aliases = {query}
    if compact:
        aliases.add(compact)

    mapping_key = None
    if compact in DISEASE_ALIASES:
        mapping_key = compact
    elif query in DISEASE_ALIASES:
        mapping_key = query

    if mapping_key:
        aliases.add(mapping_key)
        aliases.update(DISEASE_ALIASES[mapping_key])
    else:
        for canonical, known_aliases in DISEASE_ALIASES.items():
            if query in known_aliases or compact in known_aliases:
                aliases.update({canonical, query, compact})
                break

    return sorted(
        {
            normalise_disease_key(alias)
            for alias in aliases
            if alias
        },
        key=len,
        reverse=True,
    )


def confidence_level(score: float) -> str:
    if score >= 0.80:
        return "High"
    if score >= 0.65:
        return "Medium-High"
    if score >= 0.50:
        return "Medium"
    if score >= 0.35:
        return "Medium-Low"
    return "Low"


def recommendation_level(score: float) -> str:
    if score >= 0.80:
        return "Strongly Recommended"
    if score >= 0.65:
        return "Recommended"
    if score >= 0.50:
        return "Conditionally Recommended"
    if score >= 0.35:
        return "Low Priority"
    return "Not Recommended"


def mean_available(row: pd.Series, columns: Sequence[str]) -> Optional[float]:
    values: List[float] = []
    for col in columns:
        if col in row.index and pd.notna(row[col]):
            values.append(float(row[col]))
    if not values:
        return None
    return float(np.mean(values))


def count_available(row: pd.Series, columns: Sequence[str]) -> int:
    return sum(
        1
        for col in columns
        if col in row.index and pd.notna(row[col])
    )


def bool_value(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if value is None or pd.isna(value):
        return False
    if isinstance(value, str):
        return value.strip().lower() in {
            "true", "1", "yes", "y", "t"
        }
    try:
        return bool(int(value))
    except Exception:
        return bool(value)


def alteration_status(
    has_mutation: bool,
    has_fusion: bool,
    high_impact_mutation: bool = False,
) -> str:
    """
    Classify the target-gene genomic alteration for one recommended cell line.
    """
    if has_mutation and has_fusion:
        return "MIXED: MUTATION + FUSION"

    if has_mutation:
        if high_impact_mutation:
            return "MUTATION (HIGH-IMPACT)"
        return "MUTATION"

    if has_fusion:
        return "FUSION"

    return "NO ALTERATION"


def gene_label(
    gene: str,
    has_mutation: bool,
    has_fusion: bool,
    high_impact_mutation: bool = False,
) -> str:
    """
    Add mutation/fusion marker after the gene only when an alteration exists.
    """
    status = alteration_status(
        has_mutation,
        has_fusion,
        high_impact_mutation,
    )

    if status == "NO ALTERATION":
        return gene

    return f"{gene} [{status}]"


# =============================================================================
# Merged dataset access
# =============================================================================

class MergedDataRecommender:
    """
    Reads the pre-merged master_table.csv and cellline_annotations.csv.

    No data_loader.py is needed.
    """

    REQUIRED_MASTER_COLUMNS = {
        "DepMap_ID",
        "gene",
        "depmap_rna__std",
        "hpa_rna__std",
        "geo_rna__std",
        "ccle_gygi_protein__std",
        "has_rna",
        "has_protein",
        "has_mutation",
        "has_fusion",
        "rna_consistency",
        "data_completeness",
    }

    REQUIRED_ANNOTATION_COLUMNS = {
        "DepMap_ID",
        "stripped_cell_line_name",
        "lineage",
        "primary_disease",
    }

    def __init__(self, data_dir: Path):
        self.data_dir = Path(data_dir)

        master_path = self.data_dir / MASTER_FILENAME
        annotation_path = self.data_dir / ANNOTATION_FILENAME

        if not master_path.exists():
            raise FileNotFoundError(
                f"master_table.csv not found: {master_path}"
            )

        if not annotation_path.exists():
            raise FileNotFoundError(
                f"cellline_annotations.csv not found: {annotation_path}"
            )

        print(f"Loading master table: {master_path}")
        self.master = pd.read_csv(master_path, low_memory=False)

        print(f"Loading cell-line annotations: {annotation_path}")
        self.annotations = pd.read_csv(
            annotation_path,
            low_memory=False,
        )

        self._validate_columns()

        self.master["gene"] = (
            self.master["gene"]
            .astype(str)
            .str.upper()
            .str.strip()
        )

        self.annotations = self.annotations.drop_duplicates(
            subset=["DepMap_ID"],
            keep="first",
        ).copy()

        self._available_genes = sorted(
            self.master["gene"].dropna().unique().tolist()
        )

        self._last_rows: List[Dict[str, Any]] = []
        self._last_target_gene: Optional[str] = None
        self._last_exclusion_gene: Optional[str] = None

        print(
            f"Merged dataset OK. master_table shape = "
            f"{self.master.shape}"
        )
        print(
            f"Cell-line annotations shape = "
            f"{self.annotations.shape}"
        )
        print(
            f"Available target-gene panel = "
            f"{len(self._available_genes)} genes"
        )

    def _validate_columns(self) -> None:
        missing_master = sorted(
            self.REQUIRED_MASTER_COLUMNS
            - set(self.master.columns)
        )
        if missing_master:
            raise KeyError(
                "master_table.csv is missing required columns: "
                + ", ".join(missing_master)
            )

        missing_annotation = sorted(
            self.REQUIRED_ANNOTATION_COLUMNS
            - set(self.annotations.columns)
        )
        if missing_annotation:
            raise KeyError(
                "cellline_annotations.csv is missing required columns: "
                + ", ".join(missing_annotation)
            )

    @property
    def available_genes(self) -> List[str]:
        return list(self._available_genes)

    def check_gene_exists(self, gene: str) -> bool:
        return gene.upper() in set(self._available_genes)

    def _gene_table(self, gene: str) -> pd.DataFrame:
        gene = gene.upper().strip()
        table = self.master.loc[
            self.master["gene"] == gene
        ].copy()

        if table.empty:
            raise KeyError(
                f"Gene '{gene}' is not available in the merged "
                f"{len(self._available_genes)}-gene panel."
            )

        return table

    def _merge_annotations(
        self,
        gene_table: pd.DataFrame,
    ) -> pd.DataFrame:
        return gene_table.merge(
            self.annotations,
            on="DepMap_ID",
            how="left",
            validate="many_to_one",
        )

    def _disease_filter(
        self,
        df: pd.DataFrame,
        disease_aliases: Sequence[str],
    ) -> pd.DataFrame:
        # HARD disease filter:
        # use disease/lineage only. Do NOT use derived_site here because values
        # such as "Cervical lymph node" describe a sampling/metastatic site and
        # would incorrectly classify non-cervical cancers as cervical cancer.
        metadata_columns = [
            col for col in [
                "primary_disease",
                "lineage",
            ]
            if col in df.columns
        ]

        if not metadata_columns:
            raise KeyError(
                "No disease/lineage metadata columns were found."
            )

        text = (
            df[metadata_columns]
            .fillna("")
            .astype(str)
            .agg(" ".join, axis=1)
            .map(normalise_disease_key)
        )

        if not disease_aliases:
            out = df.copy()
            out["matchedDiseaseText"] = text
            return out

        mask = text.map(
            lambda value: any(
                alias in value
                for alias in disease_aliases
            )
        )

        out = df.loc[mask].copy()
        out["matchedDiseaseText"] = text.loc[mask]
        return out

    def _prepare_exclusion_expression(
        self,
        exclusion_gene: Optional[str],
        candidate_ids: Sequence[str],
    ) -> Dict[str, Optional[float]]:
        """
        Return exclusion-gene RNA evidence by DepMap_ID.

        Unlike the old loader version, exclusion penalty here is based on
        exclusion-gene RNA expression, not mutation/fusion presence.
        """
        if not exclusion_gene:
            return {}

        exclusion_gene = exclusion_gene.upper().strip()

        if not self.check_gene_exists(exclusion_gene):
            print(
                f"Warning: exclusion gene '{exclusion_gene}' is not "
                "available in the merged gene panel. "
                "Exclusion penalty will be 0."
            )
            return {}

        ex = self._gene_table(exclusion_gene)
        ex = ex.loc[
            ex["DepMap_ID"].isin(candidate_ids)
        ].copy()

        result: Dict[str, Optional[float]] = {}
        for _, row in ex.iterrows():
            result[str(row["DepMap_ID"])] = mean_available(
                row,
                RNA_STD_COLUMNS,
            )

        return result

    def fetch_candidate_evidence(
        self,
        target_gene: str,
        disease_aliases: List[str],
        exclusion_gene: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        target_gene = target_gene.upper().strip()

        target = self._gene_table(target_gene)
        target = self._merge_annotations(target)
        target = self._disease_filter(
            target,
            disease_aliases,
        )

        if target.empty:
            self._last_rows = []
            return []

        exclusion_map = self._prepare_exclusion_expression(
            exclusion_gene,
            target["DepMap_ID"].astype(str).tolist(),
        )

        rows: List[Dict[str, Any]] = []

        for _, row in target.iterrows():
            rna_value = mean_available(
                row,
                RNA_STD_COLUMNS,
            )
            n_rna = count_available(
                row,
                RNA_STD_COLUMNS,
            )

            protein_value = (
                float(row[PROTEIN_STD_COLUMN])
                if pd.notna(row.get(PROTEIN_STD_COLUMN))
                else None
            )

            n_protein = int(protein_value is not None)

            # Require at least RNA or protein target evidence.
            if n_rna == 0 and n_protein == 0:
                continue

            depmap_id = str(row["DepMap_ID"])

            has_mutation = bool_value(
                row.get("has_mutation")
            )
            has_fusion = bool_value(
                row.get("has_fusion")
            )

            high_impact_mutation = (
                pd.notna(
                    row.get(
                        MUTATION_HIGH_IMPACT_RAW_COLUMN
                    )
                )
                and safe_float(
                    row.get(
                        MUTATION_HIGH_IMPACT_RAW_COLUMN
                    )
                ) > 0
            )

            display_name = None

            for name_col in [
                "cellosaurus_name",
                "stripped_cell_line_name",
            ]:
                if (
                    name_col in row.index
                    and pd.notna(row[name_col])
                    and str(row[name_col]).strip()
                ):
                    display_name = str(
                        row[name_col]
                    ).strip()
                    break

            if not display_name:
                display_name = depmap_id

            exclusion_value = exclusion_map.get(
                depmap_id
            )

            candidate = {
                "DepMap_ID": depmap_id,
                "cellLine": display_name,
                "lineage": (
                    str(row["lineage"])
                    if pd.notna(row.get("lineage"))
                    else None
                ),
                "disease": (
                    str(row["primary_disease"])
                    if pd.notna(
                        row.get("primary_disease")
                    )
                    else None
                ),
                "matchedDiseaseText": row.get(
                    "matchedDiseaseText"
                ),

                # Biological expression layers.
                "rnaExpr": rna_value,
                "protExpr": protein_value,
                "exclusionExpr": exclusion_value,
                "nRna": n_rna,
                "nProt": n_protein,
                "nExclusion": int(
                    exclusion_value is not None
                ),

                # RNA-source flags.
                "hasDepMapRNA": pd.notna(
                    row.get("depmap_rna__std")
                ),
                "hasHpaRNA": pd.notna(
                    row.get("hpa_rna__std")
                ),
                "hasGeoRNA": pd.notna(
                    row.get("geo_rna__std")
                ),
                "hasProteomics": pd.notna(
                    row.get(
                        PROTEIN_STD_COLUMN
                    )
                ),

                # Raw evidence values for traceability.
                "depmapRnaRaw": (
                    None
                    if pd.isna(
                        row.get("depmap_rna__raw")
                    )
                    else float(
                        row["depmap_rna__raw"]
                    )
                ),
                "hpaRnaRaw": (
                    None
                    if pd.isna(
                        row.get("hpa_rna__raw")
                    )
                    else float(
                        row["hpa_rna__raw"]
                    )
                ),
                "geoRnaRaw": (
                    None
                    if pd.isna(
                        row.get("geo_rna__raw")
                    )
                    else float(
                        row["geo_rna__raw"]
                    )
                ),
                "proteinRaw": (
                    None
                    if pd.isna(
                        row.get(
                            PROTEIN_RAW_COLUMN
                        )
                    )
                    else float(
                        row[
                            PROTEIN_RAW_COLUMN
                        ]
                    )
                ),

                # Existing merged-QC fields.
                "rnaConsistency": (
                    None
                    if pd.isna(
                        row.get(
                            "rna_consistency"
                        )
                    )
                    else float(
                        row[
                            "rna_consistency"
                        ]
                    )
                ),
                "mergedDataCompleteness": (
                    None
                    if pd.isna(
                        row.get(
                            "data_completeness"
                        )
                    )
                    else float(
                        row[
                            "data_completeness"
                        ]
                    )
                ),

                # Mutation/fusion evidence.
                "hasMutation": has_mutation,
                "hasFusion": has_fusion,
                "hasHighImpactMutation": (
                    high_impact_mutation
                ),
                "mutationRaw": (
                    None
                    if pd.isna(
                        row.get(
                            MUTATION_RAW_COLUMN
                        )
                    )
                    else float(
                        row[
                            MUTATION_RAW_COLUMN
                        ]
                    )
                ),
                "mutationHighImpactRaw": (
                    None
                    if pd.isna(
                        row.get(
                            MUTATION_HIGH_IMPACT_RAW_COLUMN
                        )
                    )
                    else float(
                        row[
                            MUTATION_HIGH_IMPACT_RAW_COLUMN
                        ]
                    )
                ),
                "fusionRaw": (
                    None
                    if pd.isna(
                        row.get(
                            FUSION_RAW_COLUMN
                        )
                    )
                    else float(
                        row[
                            FUSION_RAW_COLUMN
                        ]
                    )
                ),
                "fusionConfidenceRaw": (
                    None
                    if pd.isna(
                        row.get(
                            FUSION_CONFIDENCE_RAW_COLUMN
                        )
                    )
                    else float(
                        row[
                            FUSION_CONFIDENCE_RAW_COLUMN
                        ]
                    )
                ),

                # Practical cell-line metadata.
                "assayReadyScore": (
                    None
                    if pd.isna(
                        row.get(
                            "assay_ready_score"
                        )
                    )
                    else float(
                        row[
                            "assay_ready_score"
                        ]
                    )
                ),
                "riskFlags": (
                    None
                    if pd.isna(
                        row.get("risk_flags")
                    )
                    else str(
                        row["risk_flags"]
                    )
                ),
            }

            candidate["alterationStatus"] = (
                alteration_status(
                    candidate[
                        "hasMutation"
                    ],
                    candidate[
                        "hasFusion"
                    ],
                    candidate[
                        "hasHighImpactMutation"
                    ],
                )
            )

            candidate["targetGeneLabel"] = (
                gene_label(
                    target_gene,
                    candidate[
                        "hasMutation"
                    ],
                    candidate[
                        "hasFusion"
                    ],
                    candidate[
                        "hasHighImpactMutation"
                    ],
                )
            )

            rows.append(candidate)

        self._last_rows = rows
        self._last_target_gene = target_gene
        self._last_exclusion_gene = exclusion_gene

        print("Detected merged evidence mapping:")
        print(
            "- RNA source(s): "
            "DepMap RNA + HPA RNA + GEO RNA"
        )
        print(
            "- Protein source: "
            "CCLE-Gygi proteomics"
        )
        print(
            "- Mutation source: "
            "DepMap somatic mutation"
        )
        print(
            "- Fusion source: "
            "DepMap fusion"
        )
        if exclusion_gene:
            print(
                "- Exclusion penalty source: "
                f"{exclusion_gene.upper()} "
                "RNA expression"
            )

        return rows

    def fetch_evidence_trace(
        self,
        row: Dict[str, Any],
        target_gene: str,
    ) -> List[Dict[str, Any]]:
        trace: List[Dict[str, Any]] = []

        if row.get("depmapRnaRaw") is not None:
            trace.append({
                "dataset": "DepMap_RNAseq",
                "evidenceType": "rna_expression",
                "value": row.get(
                    "depmapRnaRaw"
                ),
                "unit": "log2(TPM+1)",
            })

        if row.get("hpaRnaRaw") is not None:
            trace.append({
                "dataset": "HPA_RNAseq",
                "evidenceType": "rna_expression",
                "value": row.get(
                    "hpaRnaRaw"
                ),
                "unit": "log2(nTPM+1)",
            })

        if row.get("geoRnaRaw") is not None:
            trace.append({
                "dataset": "GEO_RNA",
                "evidenceType": "rna_expression",
                "value": row.get(
                    "geoRnaRaw"
                ),
                "unit": "log2(intensity+1)",
            })

        if row.get("proteinRaw") is not None:
            trace.append({
                "dataset": "CCLE_Gygi_Proteomics",
                "evidenceType": "protein_expression",
                "value": row.get(
                    "proteinRaw"
                ),
                "unit": "log2 MS intensity",
            })

        if row.get("hasMutation"):
            trace.append({
                "dataset": "DepMap_Mutation",
                "evidenceType": (
                    "high_impact_mutation"
                    if row.get(
                        "hasHighImpactMutation"
                    )
                    else "mutation"
                ),
                "value": row.get(
                    "mutationRaw"
                ),
                "unit": "event count",
            })

        if row.get("hasFusion"):
            trace.append({
                "dataset": "DepMap_Fusion",
                "evidenceType": "fusion",
                "value": row.get(
                    "fusionRaw"
                ),
                "unit": "event evidence",
            })

            if (
                row.get(
                    "fusionConfidenceRaw"
                )
                is not None
            ):
                trace.append({
                    "dataset": "DepMap_Fusion",
                    "evidenceType": (
                        "fusion_confidence"
                    ),
                    "value": row.get(
                        "fusionConfidenceRaw"
                    ),
                    "unit": "source confidence",
                })

        return trace

    @staticmethod
    def fetch_similar_cell_lines(
        scored_rows: List[Dict[str, Any]],
        top_k: int = 5,
    ) -> List[Dict[str, Any]]:
        """
        Preserve the previous fallback similarity logic:
        compare normalized RNA/protein components among disease-matched
        candidates.
        """
        if len(scored_rows) <= 1:
            return []

        top = scored_rows[0]
        alternatives: List[Dict[str, Any]] = []

        for row in scored_rows[1:]:
            differences: List[float] = []

            if (
                top.get("rnaScore") is not None
                and row.get("rnaScore") is not None
            ):
                differences.append(
                    abs(
                        float(top["rnaScore"])
                        - float(row["rnaScore"])
                    )
                )

            if (
                top.get("proteinScore") is not None
                and row.get("proteinScore") is not None
            ):
                differences.append(
                    abs(
                        float(
                            top["proteinScore"]
                        )
                        - float(
                            row["proteinScore"]
                        )
                    )
                )

            similarity = (
                1.0
                - (
                    sum(differences)
                    / len(differences)
                )
                if differences
                else 0.0
            )

            similarity = max(
                0.0,
                min(1.0, similarity),
            )

            alternatives.append({
                "alternativeCellLine": (
                    row.get("cellLine")
                ),
                "lineage": row.get(
                    "lineage"
                ),
                "disease": row.get(
                    "disease"
                ),
                "similarityScore": round(
                    similarity,
                    4,
                ),
                "finalScore": row.get(
                    "finalScore"
                ),
                "targetGeneLabel": row.get(
                    "targetGeneLabel"
                ),
            })

        alternatives.sort(
            key=lambda x: (
                safe_float(
                    x.get("similarityScore")
                ),
                safe_float(
                    x.get("finalScore")
                ),
            ),
            reverse=True,
        )

        return alternatives[:top_k]


# =============================================================================
# Scoring
# =============================================================================

def scale_available_values(
    rows: List[Dict[str, Any]],
    value_key: str,
    count_key: str,
) -> List[Optional[float]]:
    present_indices = [
        i
        for i, row in enumerate(rows)
        if safe_float(
            row.get(count_key),
            0.0,
        ) > 0
        and row.get(value_key) is not None
    ]

    output: List[Optional[float]] = [
        None
    ] * len(rows)

    if not present_indices:
        return output

    scaled = min_max_scale([
        safe_float(
            rows[i].get(value_key)
        )
        for i in present_indices
    ])

    for i, value in zip(
        present_indices,
        scaled,
    ):
        output[i] = value

    return output


def score_candidates(
    rows: List[Dict[str, Any]],
    config: ScoringConfig = DEFAULT_SCORING_CONFIG,
) -> List[Dict[str, Any]]:
    """
    Preserve the previous biological scoring framework.

    RNA = average of the available standardized:
        DepMap RNA + HPA RNA + GEO RNA

    Protein = standardized CCLE-Gygi protein

    Mutation/fusion are annotations and evidence-trace features.
    They do NOT automatically increase/decrease final score because the
    biological meaning of a mutation/fusion is context dependent.

    Exclusion penalty uses the exclusion gene's RNA expression.
    """
    if not rows:
        return []

    rna_scores = scale_available_values(
        rows,
        "rnaExpr",
        "nRna",
    )

    protein_scores = scale_available_values(
        rows,
        "protExpr",
        "nProt",
    )

    exclusion_scores = scale_available_values(
        rows,
        "exclusionExpr",
        "nExclusion",
    )

    scored_rows: List[Dict[str, Any]] = []

    for i, row in enumerate(rows):
        has_rna = (
            rna_scores[i] is not None
        )
        has_protein = (
            protein_scores[i] is not None
        )
        has_exclusion = (
            exclusion_scores[i] is not None
        )

        rna_score = (
            float(rna_scores[i])
            if has_rna
            else 0.0
        )

        protein_score = (
            float(protein_scores[i])
            if has_protein
            else None
        )

        exclusion_penalty = (
            config.max_exclusion_penalty
            * float(
                exclusion_scores[i]
            )
            if has_exclusion
            else 0.0
        )

        # Preserve expression-layer completeness from the old code.
        completeness_score = (
            float(has_rna)
            + float(has_protein)
        ) / 2.0

        # Updated source support: the merged table exposes
        # DepMap RNA, HPA RNA, GEO RNA and Gygi protein separately.
        source_support_score = (
            float(
                bool(
                    row.get(
                        "hasDepMapRNA"
                    )
                )
            )
            + float(
                bool(
                    row.get(
                        "hasHpaRNA"
                    )
                )
            )
            + float(
                bool(
                    row.get(
                        "hasGeoRNA"
                    )
                )
            )
            + float(
                bool(
                    row.get(
                        "hasProteomics"
                    )
                )
            )
        ) / 4.0

        if has_rna and has_protein:
            rna_protein_consistency = (
                1.0
                - abs(
                    rna_score
                    - float(
                        protein_score
                    )
                )
            )
        elif has_rna or has_protein:
            rna_protein_consistency = 0.5
        else:
            rna_protein_consistency = 0.0

        rna_protein_consistency = max(
            0.0,
            min(
                1.0,
                rna_protein_consistency,
            ),
        )

        confidence_score = (
            config.confidence_completeness
            * completeness_score
            + config.confidence_source_support
            * source_support_score
            + config.confidence_consistency
            * rna_protein_consistency
        )

        available_weight = 0.0
        weighted_biology = 0.0

        if has_rna:
            available_weight += config.rna_weight
            weighted_biology += (
                config.rna_weight
                * rna_score
            )

        if has_protein:
            available_weight += (
                config.protein_weight
            )
            weighted_biology += (
                config.protein_weight
                * float(
                    protein_score
                )
            )

        if available_weight == 0:
            continue

        # Row-wise redistribution when one biological layer is missing.
        biological_score = (
            weighted_biology
            / available_weight
        )

        final_score = (
            (
                config.rna_weight
                + config.protein_weight
            )
            * biological_score
            + config.confidence_weight
            * confidence_score
            - exclusion_penalty
        )

        final_score = max(
            0.0,
            min(1.0, final_score),
        )

        confidence_score = max(
            0.0,
            min(
                1.0,
                confidence_score,
            ),
        )

        scored = dict(row)
        scored.update({
            "scoringConfigName": config.name,
            "rnaWeight": config.rna_weight,
            "proteinWeight": config.protein_weight,
            "confidenceWeight": config.confidence_weight,
            "confidenceCompletenessWeight": config.confidence_completeness,
            "confidenceSourceSupportWeight": config.confidence_source_support,
            "confidenceConsistencyWeight": config.confidence_consistency,
            "maxExclusionPenalty": config.max_exclusion_penalty,
            "rnaScore": (
                round(
                    rna_score,
                    4,
                )
                if has_rna
                else None
            ),
            "proteinScore": (
                round(
                    float(
                        protein_score
                    ),
                    4,
                )
                if has_protein
                else None
            ),
            "biologicalScore": round(
                biological_score,
                4,
            ),
            "exclusionPenalty": round(
                exclusion_penalty,
                4,
            ),
            "completenessScore": round(
                completeness_score,
                4,
            ),
            "sourceSupportScore": round(
                source_support_score,
                4,
            ),
            "rnaProteinConsistencyScore": round(
                rna_protein_consistency,
                4,
            ),
            "confidenceScore": round(
                confidence_score,
                4,
            ),
            "confidenceLevel": (
                confidence_level(
                    confidence_score
                )
            ),
            "finalScore": round(
                final_score,
                4,
            ),
            "recommendationLevel": (
                recommendation_level(
                    final_score
                )
            ),
            "_sortFinalScore": final_score,
        })

        scored_rows.append(scored)

    scored_rows.sort(
        key=lambda item: (
            -item["_sortFinalScore"],
            str(item.get("DepMap_ID", "")),
        )
    )

    for rank, row in enumerate(
        scored_rows,
        start=1,
    ):
        row["rank"] = rank
        row.pop("_sortFinalScore", None)

    return scored_rows


# =============================================================================
# Output formatting
# =============================================================================

def build_reason(
    row: Dict[str, Any],
    target_gene: str,
    exclusion_gene: Optional[str],
    disease: Optional[str],
) -> List[str]:
    reasons: List[str] = []
    if disease:
        reasons.append(
            f"The cell line passed the disease-specific hard filter for {disease}."
        )
    else:
        reasons.append("No disease or lineage hard filter was applied (pan-cancer query).")

    rna_score = row.get("rnaScore")

    if (
        rna_score is not None
        and rna_score >= 0.70
    ):
        reasons.append(
            f"Strong RNA expression evidence for {target_gene}."
        )
    elif (
        rna_score is not None
        and rna_score >= 0.40
    ):
        reasons.append(
            f"Moderate RNA expression evidence for {target_gene}."
        )
    else:
        reasons.append(
            f"Limited RNA expression evidence for {target_gene}."
        )

    if row.get("proteinScore") is not None:
        reasons.append(
            "Supporting protein-level evidence is available."
        )
    else:
        reasons.append(
            "Protein evidence is missing; its weight was redistributed and confidence was reduced."
        )

    alteration = row.get(
        "alterationStatus"
    )

    if alteration != "NO ALTERATION":
        reasons.append(
            f"Genomic alteration detected for {target_gene}: {alteration}."
        )

    if exclusion_gene:
        if (
            row.get(
                "exclusionPenalty",
                0.0,
            )
            <= 0.05
        ):
            reasons.append(
                f"The exclusion gene {exclusion_gene} does not cause a major penalty."
            )
        else:
            reasons.append(
                f"A penalty is applied because {exclusion_gene} shows relatively high RNA expression evidence."
            )

    return reasons


def print_recommendation_report(
    target_gene: str,
    exclusion_gene: Optional[str],
    disease: Optional[str],
    top_row: Dict[str, Any],
    alternatives: List[Dict[str, Any]],
    evidence_trace: List[Dict[str, Any]],
) -> None:
    """
    Keep the same report structure as the previous code.
    Mutation/fusion lines are ADDED without removing the previous sections.
    """
    cell_line = top_row.get(
        "cellLine",
        "Unknown",
    )

    print("\n" + "=" * 80)
    print(
        "CellLineSelector Merged Multi-Omics Recommendation Demo"
    )
    print("=" * 80)

    print("\nInput:")
    print(
        f"Target gene: {target_gene}"
    )

    if exclusion_gene:
        print(
            f"Exclusion gene: {exclusion_gene}"
        )

    print(f"Disease hard filter: {disease or 'None (pan-cancer)'}")

    print("\nOutput:")
    print(
        f"Recommended cell line: {cell_line}"
    )
    print(
        f"Recommendation level: {top_row.get('recommendationLevel')}"
    )
    print(
        f"Final score: {top_row.get('finalScore'):.2f} / 1.00"
    )
    print(
        f"Confidence: {top_row.get('confidenceLevel')}"
    )
    print(
        f"Confidence score: {top_row.get('confidenceScore'):.2f} / 1.00"
    )

    # New requested mutation/fusion marker.
    print(
        f"Target gene status: {top_row.get('targetGeneLabel', target_gene)}"
    )

    print("\nScore breakdown:")
    print(
        f"- RNA expression score: {top_row.get('rnaScore')}"
    )
    print(
        "- Protein score: "
        + (
            str(
                top_row.get(
                    "proteinScore"
                )
            )
            if top_row.get(
                "proteinScore"
            ) is not None
            else "Missing / redistributed"
        )
    )
    print(
        f"- Biological score: {top_row.get('biologicalScore')}"
    )
    print(
        f"- Exclusion penalty: {top_row.get('exclusionPenalty')}"
    )
    print(
        f"- Completeness score: {top_row.get('completenessScore')}"
    )
    print(
        f"- Source support score: {top_row.get('sourceSupportScore')}"
    )

    # Extra information from the merged dataset.
    if (
        top_row.get(
            "rnaConsistency"
        )
        is not None
    ):
        print(
            "- RNA cross-source consistency: "
            f"{top_row.get('rnaConsistency'):.4f}"
        )

    print("\nMain reason:")
    for reason in build_reason(
        top_row,
        target_gene,
        exclusion_gene,
        disease,
    ):
        print(f"- {reason}")

    print("\nEvidence trace:")

    if evidence_trace:
        for evidence in evidence_trace[:15]:
            print(
                f"- {evidence.get('dataset')}: "
                f"{evidence.get('evidenceType')}, "
                f"value={evidence.get('value')}, "
                f"unit={evidence.get('unit')}"
            )
    else:
        print(
            "- No detailed evidence trace found."
        )

    print("\nData gaps:")

    gaps: List[str] = []

    if not top_row.get("hasDepMapRNA"):
        gaps.append(
            "DepMap RNA-seq evidence is missing."
        )

    if not top_row.get("hasHpaRNA"):
        gaps.append(
            "HPA RNA-seq evidence is missing."
        )

    if not top_row.get("hasGeoRNA"):
        gaps.append(
            "GEO RNA expression evidence is missing."
        )

    if not top_row.get("hasProteomics"):
        gaps.append(
            "CCLE-Gygi protein evidence is missing."
        )

    if gaps:
        for gap in gaps:
            print(f"- {gap}")
    else:
        print(
            "- No major expression-data gap detected for this target."
        )

    print("\nSimilar alternative cell lines:")

    if alternatives:
        for i, alt in enumerate(
            alternatives,
            start=1,
        ):
            target_marker = alt.get(
                "targetGeneLabel",
                target_gene,
            )
            print(
                f"{i}. {alt.get('alternativeCellLine')} "
                f"— similarity score: "
                f"{safe_float(alt.get('similarityScore')):.2f} "
                f"| {target_marker}"
            )
    else:
        print(
            "- No similar alternative cell lines found after the disease filter."
        )

    print("\nNote:")
    if disease:
        print(
            f"- All recommendations and alternatives are restricted to the disease filter: {disease}."
        )
    else:
        print("- No disease/lineage filter was applied; candidates are pan-cancer.")
    print(
        "- Data are read directly from the pre-merged multi-omics master table."
    )
    print(
        "- RNA evidence combines standardized DepMap, HPA and GEO RNA values."
    )
    print(
        "- Ranking preserves the previous RNA/protein/confidence weighting scheme."
    )
    print(
        "- Mutation/fusion status is reported as genomic evidence but is not automatically treated as beneficial or harmful."
    )
    print(
        "- Alternative similarity uses the same normalized RNA/protein-component fallback as the previous non-Neo4j version."
    )
    print(
        "- This workflow supports project validation and does not constitute a final biological decision."
    )


def save_ranked_results(
    rows: List[Dict[str, Any]],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not rows:
        return

    fieldnames = [
        "rank",
        "DepMap_ID",
        "cellLine",
        "lineage",
        "disease",
        "scoringConfigName",
        "rnaWeight",
        "proteinWeight",
        "confidenceWeight",
        "confidenceCompletenessWeight",
        "confidenceSourceSupportWeight",
        "confidenceConsistencyWeight",
        "maxExclusionPenalty",
        "targetGeneLabel",
        "alterationStatus",
        "hasMutation",
        "hasHighImpactMutation",
        "mutationRaw",
        "mutationHighImpactRaw",
        "hasFusion",
        "fusionRaw",
        "fusionConfidenceRaw",
        "rnaExpr",
        "protExpr",
        "exclusionExpr",
        "rnaScore",
        "proteinScore",
        "biologicalScore",
        "exclusionPenalty",
        "completenessScore",
        "sourceSupportScore",
        "rnaProteinConsistencyScore",
        "rnaConsistency",
        "mergedDataCompleteness",
        "confidenceScore",
        "confidenceLevel",
        "finalScore",
        "recommendationLevel",
        "hasDepMapRNA",
        "hasHpaRNA",
        "hasGeoRNA",
        "hasProteomics",
        "depmapRnaRaw",
        "hpaRnaRaw",
        "geoRnaRaw",
        "proteinRaw",
        "assayReadyScore",
        "riskFlags",
    ]

    with output_path.open(
        "w",
        newline="",
        encoding="utf-8-sig",
    ) as file:
        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames,
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def load_scoring_config(config_path: Optional[str]) -> ScoringConfig:
    """Load an optional JSON configuration; otherwise use provisional A0."""
    if not config_path:
        return DEFAULT_SCORING_CONFIG

    path = Path(config_path).expanduser()
    with path.open("r", encoding="utf-8") as file:
        payload = json.load(file)
    if not isinstance(payload, dict):
        raise TypeError("scoring config JSON must contain an object")

    allowed = {field.name for field in fields(ScoringConfig)}
    unknown = sorted(set(payload) - allowed)
    if unknown:
        raise ValueError("unknown scoring config field(s): " + ", ".join(unknown))
    return ScoringConfig(**payload)


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    return value


def save_recommendation_payload(
    output_path: Path,
    *,
    target_gene: str,
    exclusion_gene: Optional[str],
    disease: Optional[str],
    scoring_config: ScoringConfig,
    recommendations: List[Dict[str, Any]],
    alternatives: List[Dict[str, Any]],
    evidence_trace: List[Dict[str, Any]],
    ranking_csv: Path,
) -> None:
    """Write a stable UI/RAG-facing JSON payload alongside the full CSV."""
    payload = {
        "schemaVersion": OUTPUT_SCHEMA_VERSION,
        "query": {
            "targetGene": target_gene,
            "exclusionGene": exclusion_gene,
            "diseaseFilter": disease,
            "queryScope": "disease_filtered" if disease else "pan_cancer",
        },
        "scoringConfig": scoring_config.to_dict(),
        "recommendations": recommendations,
        "alternatives": alternatives,
        "topRecommendationEvidenceTrace": evidence_trace,
        "artifacts": {"fullRankingCsv": ranking_csv.name},
        "knownLimitations": [
            "The default A0 weights are provisional until benchmark expansion and weight search are complete.",
            "Alternative similarity currently uses target RNA/protein component distance, not a full-omics embedding.",
            "Mutation and fusion evidence is annotated but is not assigned a universal positive or negative weight.",
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as file:
        json.dump(_json_safe(payload), file, ensure_ascii=False, indent=2)


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run CellLineSelector directly from the "
            "pre-merged multi-omics dataset."
        )
    )

    parser.add_argument(
        "--target_gene",
        type=str,
        default=None,
        help="Target gene, e.g. EGFR.",
    )

    parser.add_argument(
        "--exclusion_gene",
        type=str,
        default=None,
        help=(
            "Optional exclusion gene, "
            "e.g. ABCB1."
        ),
    )

    parser.add_argument(
        "--disease",
        "--context",
        dest="disease",
        type=str,
        default=None,
        help=(
            "Optional disease/tissue hard filter. Omit for pan-cancer."
        ),
    )

    parser.add_argument(
        "--data_dir",
        type=str,
        default=DEFAULT_DATA_DIR,
        help=(
            "Folder containing master_table.csv "
            "and cellline_annotations.csv. Can also be set with "
            "CELLLINESELECTOR_DATA_DIR."
        ),
    )

    parser.add_argument(
        "--top_n",
        type=int,
        default=DEFAULT_TOP_N,
    )

    parser.add_argument(
        "--top_k_alternatives",
        type=int,
        default=DEFAULT_TOP_K_ALTERNATIVES,
    )

    parser.add_argument(
        "--output_dir",
        type=str,
        default=DEFAULT_OUTPUT_DIR,
    )

    parser.add_argument(
        "--scoring_config",
        type=str,
        default=None,
        help="Optional JSON ScoringConfig file; defaults to provisional A0.",
    )

    parser.add_argument(
        "--output_prefix",
        type=str,
        default=None,
        help="Optional prefix for CSV/JSON outputs.",
    )

    parser.add_argument(
        "--no_prompt",
        action="store_true",
        help="Disable interactive prompts; --target_gene and --data_dir are then required.",
    )

    return parser.parse_args()


def ask_user_input(
    value: Optional[str],
    prompt_text: str,
    required: bool = False,
) -> Optional[str]:
    if (
        value is not None
        and str(value).strip()
    ):
        return str(value).strip()

    while True:
        user_value = input(
            prompt_text
        ).strip()

        if not user_value:
            if required:
                print(
                    "This field is required. "
                    "Please enter a value."
                )
                continue
            return None

        return user_value


def main() -> None:
    args = parse_args()

    if args.no_prompt:
        if not args.target_gene:
            raise ValueError("--target_gene is required with --no_prompt")
    else:
        args.target_gene = ask_user_input(
            args.target_gene,
            "Target gene: ",
            required=True,
        )
        args.exclusion_gene = ask_user_input(
            args.exclusion_gene,
            "Exclusion gene (press Enter to skip): ",
            required=False,
        )
        args.disease = ask_user_input(
            args.disease,
            "Disease or tissue (press Enter for pan-cancer): ",
            required=False,
        )

    if not args.data_dir:
        raise ValueError(
            "Data directory is required. Pass --data_dir or set "
            "CELLLINESELECTOR_DATA_DIR."
        )

    args.target_gene = (
        str(args.target_gene)
        .upper()
        .strip()
    )

    if args.exclusion_gene:
        args.exclusion_gene = (
            str(args.exclusion_gene)
            .upper()
            .strip()
        )

    disease_aliases = build_disease_aliases(
        args.disease
    )
    scoring_config = load_scoring_config(args.scoring_config)

    try:
        print("\n" + "=" * 80)
        print(
            "Merged multi-omics data initialization"
        )
        print("=" * 80)

        data_dir = Path(
            args.data_dir
        ).expanduser()

        print(
            f"Data directory: {data_dir}"
        )

        data = MergedDataRecommender(
            data_dir
        )

        if not data.check_gene_exists(
            args.target_gene
        ):
            print(
                "\nTarget gene check failed."
            )
            print(
                f"Target gene '{args.target_gene}' "
                "is not available in the merged gene panel."
            )
            print(
                f"Available genes ({len(data.available_genes)}):"
            )
            print(
                ", ".join(
                    data.available_genes
                )
            )
            return

        print("\n" + "=" * 80)
        print(
            "Reading evidence from merged multi-omics data..."
        )
        print(
            f"Target gene: {args.target_gene}"
        )
        print(
            "Exclusion gene: "
            f"{args.exclusion_gene or 'None'}"
        )
        print(
            f"Disease hard filter: {args.disease or 'None (pan-cancer)'}"
        )
        if disease_aliases:
            print("Disease aliases: " + ", ".join(disease_aliases))
        print(f"Scoring config: {scoring_config.name}")

        candidate_rows = (
            data.fetch_candidate_evidence(
                target_gene=args.target_gene,
                disease_aliases=disease_aliases,
                exclusion_gene=args.exclusion_gene,
            )
        )

        if not candidate_rows:
            print(
                "\nNo candidate cell lines "
                "with target-gene evidence were found."
            )
            print(
                "Check the disease wording and the merged "
                "cellline_annotations.csv fields."
            )
            return

        scored_rows = score_candidates(
            candidate_rows,
            config=scoring_config,
        )

        if not scored_rows:
            print(
                "No candidates could be scored after filtering."
            )
            return

        top_rows = scored_rows[
            : max(
                1,
                args.top_n,
            )
        ]

        top_row = top_rows[0]

        alternatives = (
            data.fetch_similar_cell_lines(
                scored_rows,
                top_k=max(
                    0,
                    args.top_k_alternatives,
                ),
            )
        )

        evidence_trace = (
            data.fetch_evidence_trace(
                top_row,
                args.target_gene,
            )
        )

        print_recommendation_report(
            target_gene=args.target_gene,
            exclusion_gene=args.exclusion_gene,
            disease=args.disease,
            top_row=top_row,
            alternatives=alternatives,
            evidence_trace=evidence_trace,
        )

        print("\nTop ranked cell lines:")

        for row in top_rows:
            # Preserve the previous line structure and add the
            # requested mutation/fusion marker at the end.
            marker = row.get(
                "targetGeneLabel",
                args.target_gene,
            )

            print(
                f"{row['rank']}. "
                f"{row['cellLine']} | "
                f"finalScore={row['finalScore']:.2f} | "
                f"confidence={row['confidenceScore']:.2f} | "
                f"level={row['recommendationLevel']} | "
                f"{marker}"
            )

        output_dir = Path(
            args.output_dir
        )

        if not output_dir.is_absolute():
            output_dir = Path.cwd() / output_dir

        output_prefix = args.output_prefix
        if not output_prefix:
            context = args.disease or "pan_cancer"
            output_prefix = re.sub(
                r"[^A-Za-z0-9]+", "_", f"{args.target_gene}_{context}"
            ).strip("_").lower()

        output_path = (
            output_dir
            / f"{output_prefix}_ranked_recommendations.csv"
        )

        save_ranked_results(
            scored_rows,
            output_path,
        )

        json_output_path = output_dir / f"{output_prefix}_recommendations.json"
        save_recommendation_payload(
            json_output_path,
            target_gene=args.target_gene,
            exclusion_gene=args.exclusion_gene,
            disease=args.disease,
            scoring_config=scoring_config,
            recommendations=top_rows,
            alternatives=alternatives,
            evidence_trace=evidence_trace,
            ranking_csv=output_path,
        )

        print(
            f"\nSaved ranked results to: "
            f"{output_path.resolve()}"
        )
        print(f"Saved structured payload to: {json_output_path.resolve()}")

    except KeyboardInterrupt:
        print(
            "\nProgram interrupted by user."
        )
        raise SystemExit(130)

    except Exception as exc:
        print(
            "\nProgram stopped because merged-data "
            "loading or recommendation failed."
        )
        print(
            f"{type(exc).__name__}: {exc}"
        )
        print("\nQuick checks:")
        print(
            "1. Confirm the data directory is:"
        )
        print(
            "   Pass --data_dir or set CELLLINESELECTOR_DATA_DIR."
        )
        print(
            "2. Confirm master_table.csv exists."
        )
        print(
            "3. Confirm cellline_annotations.csv exists."
        )
        print(
            "4. Confirm pandas and numpy import successfully."
        )
        print(
            "5. If the target gene is not found, "
            "remember that the current merged master table "
            "contains a selected gene panel rather than every human gene."
        )
        raise


if __name__ == "__main__":
    main()
