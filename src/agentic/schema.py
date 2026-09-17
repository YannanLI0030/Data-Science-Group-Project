"""Small, serialisable contracts shared by the output agent and grounder."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass(frozen=True)
class EvidenceCard:
    evidence_id: str
    card_type: str
    text: str
    source: str
    depmap_id: str | None = None
    cell_line: str | None = None
    target: str | None = None
    values: tuple[float, ...] = ()
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class Claim:
    text: str
    claim_type: str
    evidence_ids: tuple[str, ...]

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "Claim":
        return cls(
            text=str(value.get("text", "")).strip(),
            claim_type=str(value.get("claim_type", "SUPPORT")).upper(),
            evidence_ids=tuple(str(x) for x in value.get("evidence_ids", [])),
        )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ClaimVerdict:
    claim: Claim
    accepted: bool
    failed_layer: str | None = None
    reason: str = ""

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["claim"] = self.claim.to_dict()
        return value


@dataclass
class AgentOutput:
    backend: str
    grounding_mode: str
    claims: list[Claim]
    dropped_claims: list[ClaimVerdict]
    evidence_cards: list[EvidenceCard]
    ranking_fingerprint: str
    run_manifest: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return {
            "backend": self.backend,
            "grounding_mode": self.grounding_mode,
            "claims": [x.to_dict() for x in self.claims],
            "dropped_claims": [x.to_dict() for x in self.dropped_claims],
            "evidence_cards": [x.to_dict() for x in self.evidence_cards],
            "ranking_fingerprint": self.ranking_fingerprint,
            "run_manifest": self.run_manifest,
        }
