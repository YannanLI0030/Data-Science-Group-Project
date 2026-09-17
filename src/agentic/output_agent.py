"""Read-only agent that explains results from the user's canonical ranker."""

from __future__ import annotations

import hashlib
import json
import platform
import time
from pathlib import Path
from typing import Any

from . import grounding, llm
from .evidence_adapter import build_evidence_cards, serialise_context
from .schema import AgentOutput, Claim


def ranking_fingerprint(rows: list[dict[str, Any]]) -> str:
    """Hash only fields the Agent is forbidden to change."""
    authority = [
        {
            "rank": row.get("rank"),
            "DepMap_ID": row.get("DepMap_ID"),
            "cellLine": row.get("cellLine"),
            "finalScore": row.get("finalScore"),
            "confidenceScore": row.get("confidenceScore"),
            "recommendationLevel": row.get("recommendationLevel"),
        }
        for row in rows
    ]
    raw = json.dumps(authority, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


class OutputAgent:
    """Generate grounded prose without owning a ranker or matrix store."""

    def __init__(
        self,
        backend: str = "scripted",
        model: str | None = None,
        grounding_mode: str = "l1l2",
        api_key: str | None = None,
        endpoint: str | None = None,
    ) -> None:
        self.writer = llm.get_writer(
            backend,
            model,
            api_key=api_key,
            endpoint=endpoint,
        )
        self.grounding_mode = grounding_mode

    def run(
        self,
        request: dict[str, Any],
        ranked_rows: list[dict[str, Any]],
        evidence_trace: list[dict[str, Any]],
        alternatives: list[dict[str, Any]],
        supplementary: dict[str, Any],
    ) -> AgentOutput:
        before = ranking_fingerprint(ranked_rows)
        cards = build_evidence_cards(
            request, ranked_rows, evidence_trace, alternatives, supplementary
        )
        context = serialise_context(request, ranked_rows, cards)
        draft = self.writer.write(context)
        claims = [
            Claim.from_dict(item)
            for item in draft.get("claims", [])
            if isinstance(item, dict) and str(item.get("text", "")).strip()
        ]
        verdicts = grounding.verify(claims, cards, self.grounding_mode)

        after = ranking_fingerprint(ranked_rows)
        if after != before:
            raise RuntimeError(
                "The Agent boundary was violated: authoritative rank, score, confidence, "
                "or recommendation level changed during explanation generation."
            )

        accepted = [verdict.claim for verdict in verdicts if verdict.accepted]
        dropped = [verdict for verdict in verdicts if not verdict.accepted]
        manifest = {
            "run_id": "explain_" + time.strftime("%Y%m%dT%H%M%S"),
            "created_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "python": platform.python_version(),
            "ranking_authority": "dynamic_cellline_selector_gene_protein.score_candidates",
            "model_backend": self.writer.name + (
                f":{self.writer.model}" if hasattr(self.writer, "model") else ""
            ),
            "grounding_mode": self.grounding_mode,
            "claims_drafted": len(claims),
            "claims_emitted": len(accepted),
            "claims_dropped": len(dropped),
            "ranking_fingerprint_before": before,
            "ranking_fingerprint_after": after,
            "ranking_unchanged": before == after,
            "query": dict(request),
        }
        return AgentOutput(
            backend=self.writer.name,
            grounding_mode=self.grounding_mode,
            claims=accepted,
            dropped_claims=dropped,
            evidence_cards=cards,
            ranking_fingerprint=before,
            run_manifest=manifest,
        )


def print_agent_output(output: AgentOutput) -> None:
    print("\n" + "=" * 80)
    print("Agent-assisted interpretation (post-ranking; scores and order unchanged)")
    print("=" * 80)
    if not output.claims:
        print("No Agent claim survived grounding. The deterministic ranking above remains valid.")
    else:
        labels = {
            "SUPPORT": "Evidence",
            "LIMITATION": "Caution",
            "ALTERNATIVE": "Alternative",
            "CONTEXT": "Context",
        }
        for claim in output.claims:
            citations = ", ".join(claim.evidence_ids)
            print(
                f"- {labels.get(claim.claim_type, claim.claim_type.title())}: "
                f"{claim.text} [{citations}]"
            )
    print(
        f"\nGrounding: {output.grounding_mode}; emitted={len(output.claims)}, "
        f"dropped={len(output.dropped_claims)}"
    )
    print(f"Ranking integrity: unchanged ({output.ranking_fingerprint[:12]}...).")
    print("Ranking authority: dynamic score_candidates(); Agent cannot change scores or order.")


def save_agent_output(output: AgentOutput, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(output.to_dict(), indent=2, ensure_ascii=False), encoding="utf-8"
    )
