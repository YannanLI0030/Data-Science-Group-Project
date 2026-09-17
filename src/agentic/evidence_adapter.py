"""Convert the deterministic recommender output into cited evidence cards.

This adapter is the only bridge between the statistical backend and the agent.
It copies results into immutable cards; it performs no scoring or ranking.
"""

from __future__ import annotations

import hashlib
import json
from typing import Any, Iterable

from .schema import EvidenceCard


def _id(kind: str, *parts: Any) -> str:
    raw = json.dumps([kind, *parts], ensure_ascii=False, sort_keys=True, default=str)
    return hashlib.sha1(raw.encode("utf-8")).hexdigest()[:12]


def _float(value: Any) -> float | None:
    try:
        if value is None:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _values(*items: Any) -> tuple[float, ...]:
    return tuple(v for item in items if (v := _float(item)) is not None)


def _target(request: dict[str, Any]) -> str:
    gene = request.get("target_gene")
    protein = request.get("target_protein")
    if gene and protein:
        return f"gene {gene} and protein {protein}"
    if gene:
        return f"gene {gene}"
    return f"protein {protein}"


def _ranking_card(row: dict[str, Any], request: dict[str, Any]) -> EvidenceCard:
    rank = int(row["rank"])
    score = float(row["finalScore"])
    confidence = float(row["confidenceScore"])
    depmap_id = str(row.get("DepMap_ID", "")) or None
    cell_line = str(row.get("cellLine", depmap_id or "unknown"))
    target = _target(request)
    text = (
        f"The deterministic recommender ranks {cell_line}"
        f"{f' ({depmap_id})' if depmap_id else ''} at position {rank} for {target}, "
        f"with final score {score:.4f}, confidence score {confidence:.4f}, and "
        f"recommendation level {row.get('recommendationLevel', 'Unspecified')}."
    )
    return EvidenceCard(
        evidence_id=_id("RANKING", depmap_id, rank, score, confidence, target),
        card_type="RANKING",
        text=text,
        source="dynamic_multiomics_recommender",
        depmap_id=depmap_id,
        cell_line=cell_line,
        target=target,
        values=_values(rank, score, confidence),
        metadata={"rank": rank, "authoritative": True},
    )


def _score_card(row: dict[str, Any], request: dict[str, Any]) -> EvidenceCard:
    depmap_id = str(row.get("DepMap_ID", "")) or None
    cell_line = str(row.get("cellLine", depmap_id or "unknown"))
    rna = _float(row.get("rnaScore"))
    protein = _float(row.get("proteinScore"))
    biological = _float(row.get("biologicalScore"))
    penalty = _float(row.get("exclusionPenalty"))
    completeness = _float(row.get("completenessScore"))
    support = _float(row.get("sourceSupportScore"))
    parts = [
        f"RNA score {rna:.4f}" if rna is not None else "RNA score not measured",
        f"protein score {protein:.4f}" if protein is not None else "protein score not measured",
        f"biological score {biological:.4f}" if biological is not None else None,
        f"exclusion penalty {penalty:.4f}" if penalty is not None else None,
        f"data completeness {completeness:.4f}" if completeness is not None else None,
        f"source support {support:.4f}" if support is not None else None,
    ]
    details = ", ".join(x for x in parts if x)
    text = f"For {cell_line}, the deterministic score breakdown is: {details}."
    return EvidenceCard(
        evidence_id=_id("SCORE", depmap_id, details),
        card_type="SCORE",
        text=text,
        source="score_candidates",
        depmap_id=depmap_id,
        cell_line=cell_line,
        target=_target(request),
        values=_values(rna, protein, biological, penalty, completeness, support),
        metadata={"query_mode": request.get("query_mode")},
    )


def _gap_texts(top: dict[str, Any], request: dict[str, Any]) -> list[str]:
    gaps: list[str] = []
    mode = request.get("query_mode")
    target = request.get("target_protein") or request.get("target_gene")
    if not top.get("hasDepMapRNA"):
        gaps.append("DepMap RNA-seq evidence is not measured for the top recommendation.")
    if not top.get("hasHpaRNA"):
        gaps.append("HPA RNA-seq evidence is not measured for the top recommendation.")
    if not top.get("hasGeoRNA"):
        gaps.append("GEO RNA expression evidence is not measured for the top recommendation.")
    if not top.get("hasProteomics"):
        gaps.append(
            f"CCLE-Gygi protein evidence for {target} is not measured for the top recommendation."
        )
    if mode == "PROTEIN_ONLY" and not top.get("hasProteomics"):
        gaps.append("The protein-only query cannot substitute RNA for missing protein evidence.")
    return gaps


def build_evidence_cards(
    request: dict[str, Any],
    ranked_rows: list[dict[str, Any]],
    evidence_trace: list[dict[str, Any]],
    alternatives: list[dict[str, Any]],
    supplementary: dict[str, Any],
) -> list[EvidenceCard]:
    cards: list[EvidenceCard] = []
    for row in ranked_rows:
        cards.append(_ranking_card(row, request))
    if ranked_rows:
        cards.append(_score_card(ranked_rows[0], request))
        for text in _gap_texts(ranked_rows[0], request):
            cards.append(EvidenceCard(
                evidence_id=_id("GAP", text), card_type="GAP", text=text,
                source="data_availability", depmap_id=ranked_rows[0].get("DepMap_ID"),
                cell_line=ranked_rows[0].get("cellLine"), target=_target(request),
            ))

    top = ranked_rows[0] if ranked_rows else {}
    for item in evidence_trace:
        value = _float(item.get("value"))
        dataset = str(item.get("dataset", "unknown dataset"))
        text = (
            f"{dataset} reports {item.get('type', 'evidence')} for "
            f"{top.get('cellLine', 'the top recommendation')} with value "
            f"{item.get('value')} {item.get('unit', '')}."
        ).strip()
        cards.append(EvidenceCard(
            evidence_id=_id("MEASUREMENT", top.get("DepMap_ID"), dataset, item),
            card_type="MEASUREMENT", text=text, source=dataset,
            depmap_id=top.get("DepMap_ID"), cell_line=top.get("cellLine"),
            target=_target(request), values=_values(value), metadata=dict(item),
        ))

    for item in alternatives:
        score = _float(item.get("similarityScore"))
        text = (
            f"{item.get('alternativeCellLine')} is a similar alternative to the top "
            f"recommendation, with similarity score {score:.4f}."
            if score is not None else
            f"{item.get('alternativeCellLine')} is listed as a similar alternative."
        )
        cards.append(EvidenceCard(
            evidence_id=_id("ALTERNATIVE", item), card_type="ALTERNATIVE",
            text=text, source="deterministic_similarity",
            cell_line=str(item.get("alternativeCellLine")), target=_target(request),
            values=_values(score), metadata=dict(item),
        ))

    signatures = supplementary.get("global_signatures") or {}
    if signatures and ranked_rows:
        values = _values(*signatures.values())
        text = (
            f"Supplementary global signatures are available for {top.get('cellLine')}; "
            "they are context only and were not used in scoring."
        )
        cards.append(EvidenceCard(
            evidence_id=_id("SUPPLEMENTARY", top.get("DepMap_ID"), signatures),
            card_type="SUPPLEMENTARY", text=text, source="File 14",
            depmap_id=top.get("DepMap_ID"), cell_line=top.get("cellLine"),
            target=_target(request), values=values, metadata={"signatures": signatures},
        ))
    return cards


def serialise_context(
    request: dict[str, Any],
    ranked_rows: Iterable[dict[str, Any]],
    cards: Iterable[EvidenceCard],
) -> dict[str, Any]:
    return {
        "instruction": (
            "Explain the supplied deterministic ranking only. Never calculate a new score, "
            "change candidate order, or add an unsupported biological fact."
        ),
        "query": dict(request),
        "ranked_results": [dict(row) for row in ranked_rows],
        "cards": [card.to_dict() for card in cards],
    }
