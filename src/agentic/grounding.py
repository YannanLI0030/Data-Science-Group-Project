"""Citation, numeric and entity checks adapted from the team Agent project."""

from __future__ import annotations

import re
from typing import Iterable

from .schema import Claim, ClaimVerdict, EvidenceCard


NUMBER = re.compile(r"(?<![A-Za-z])[-+]?\d+(?:\.\d+)?")


def _mask(text: str, cards: Iterable[EvidenceCard]) -> str:
    values = []
    for card in cards:
        values.extend([card.depmap_id, card.cell_line, card.target, card.source])
    out = text
    for value in sorted({str(x) for x in values if x}, key=len, reverse=True):
        out = re.sub(re.escape(value), " ", out, flags=re.IGNORECASE)
    out = re.sub(r"log2\([^)]*\)", " ", out, flags=re.IGNORECASE)
    return out


def _numbers(text: str, cards: Iterable[EvidenceCard]) -> list[float]:
    return [float(x) for x in NUMBER.findall(_mask(text, cards))]


def _numeric_supported(number: float, cards: list[EvidenceCard]) -> bool:
    allowed = [value for card in cards for value in card.values]
    return any(abs(number - value) <= max(1e-4, abs(value) * 1e-3) for value in allowed)


def verify_claim(
    claim: Claim,
    index: dict[str, EvidenceCard],
    all_cards: list[EvidenceCard],
    mode: str = "l1l2",
) -> ClaimVerdict:
    if mode == "none":
        return ClaimVerdict(claim, True)
    if not claim.evidence_ids:
        return ClaimVerdict(claim, False, "L1_citation", "claim has no evidence id")
    unknown = [x for x in claim.evidence_ids if x not in index]
    if unknown:
        return ClaimVerdict(
            claim, False, "L1_citation", f"unknown evidence id(s): {', '.join(unknown)}"
        )
    cited = [index[x] for x in claim.evidence_ids]
    unsupported = [n for n in _numbers(claim.text, cited) if not _numeric_supported(n, cited)]
    if unsupported:
        return ClaimVerdict(
            claim, False, "L2_numeric",
            "number(s) not carried by cited evidence: " + ", ".join(map(str, unsupported)),
        )
    if mode == "full":
        cited_entities = {
            str(x).lower()
            for card in cited
            for x in (card.depmap_id, card.cell_line, card.target)
            if x
        }
        known_entities = {
            str(x).lower()
            for card in all_cards
            for x in (card.depmap_id, card.cell_line)
            if x
        }
        lowered = claim.text.lower()
        foreign = [x for x in known_entities if x in lowered and x not in cited_entities]
        if foreign:
            return ClaimVerdict(
                claim, False, "L3_entity",
                "claim names an entity not present in its cited evidence: " + ", ".join(foreign),
            )
    return ClaimVerdict(claim, True)


def verify(
    claims: list[Claim], cards: list[EvidenceCard], mode: str = "l1l2"
) -> list[ClaimVerdict]:
    if mode not in {"none", "l1l2", "full"}:
        raise ValueError("grounding mode must be none, l1l2 or full")
    index = {card.evidence_id: card for card in cards}
    return [verify_claim(claim, index, cards, mode) for claim in claims]
