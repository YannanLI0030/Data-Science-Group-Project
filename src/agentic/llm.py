"""Writers for post-ranking explanations.

Adapted from the teammate Agent module. The interface has deliberately only one
job: turn immutable deterministic results into cited prose.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Protocol
from urllib.parse import urlparse


PROJECT_ROOT = Path(__file__).resolve().parents[2]
ENV_FILE = PROJECT_ROOT / ".env"

SYSTEM_PROMPT = """You are the explanation layer of a cell-line recommender.
The supplied ranked_results were produced by a deterministic multi-omics model.
Never calculate, estimate, adjust, remove, or reorder a score or candidate.
Only explain supplied facts. Every claim must cite one or more supplied
evidence_ids. Every number in a claim must be present in its cited evidence.
A missing value means not measured, never low, absent, or zero. Follow the
output_language value in the supplied query when it is present. Return JSON
with one key, claims. Each claim contains text, claim_type, and evidence_ids."""


class Writer(Protocol):
    name: str

    def write(self, context: dict[str, Any]) -> dict[str, Any]: ...


def _key(name: str) -> str | None:
    if os.environ.get(name):
        return os.environ[name]
    if ENV_FILE.exists():
        for line in ENV_FILE.read_text(encoding="utf-8").splitlines():
            key, sep, value = line.partition("=")
            if sep and key.strip() == name:
                return value.strip().strip('"').strip("'") or None
    return None


class ScriptedWriter:
    """Deterministic offline writer; every sentence is copied from one card."""

    name = "scripted"

    def write(self, context: dict[str, Any]) -> dict[str, Any]:
        cards = context["cards"]
        selected: list[dict[str, Any]] = []
        first_of_type: dict[str, int] = {
            "RANKING": 1, "SCORE": 1, "MEASUREMENT": 3,
            "GAP": 10, "ALTERNATIVE": 3, "SUPPLEMENTARY": 1,
        }
        counts: dict[str, int] = {}
        for card in cards:
            kind = card["card_type"]
            counts[kind] = counts.get(kind, 0) + 1
            if counts[kind] > first_of_type.get(kind, 0):
                continue
            claim_type = {
                "GAP": "LIMITATION", "ALTERNATIVE": "ALTERNATIVE",
                "SUPPLEMENTARY": "CONTEXT",
            }.get(kind, "SUPPORT")
            selected.append({
                "text": card["text"],
                "claim_type": claim_type,
                "evidence_ids": [card["evidence_id"]],
            })
        return {"claims": selected}


class HttpWriter:
    endpoint = ""
    key_name = ""
    name = "http"

    def __init__(
        self,
        model: str,
        api_key: str | None = None,
        timeout: int = 90,
        endpoint: str | None = None,
        require_key: bool = True,
    ):
        self.model = model
        self.api_key = api_key or _key(self.key_name)
        self.timeout = timeout
        self.endpoint = endpoint or self.endpoint
        parsed = urlparse(self.endpoint)
        if parsed.scheme not in {"http", "https"} or not parsed.netloc:
            raise ValueError("Model API endpoint must be a valid http(s) URL")
        if parsed.username or parsed.password:
            raise ValueError("Do not put API credentials inside the endpoint URL")
        if require_key and not self.api_key:
            raise RuntimeError(f"{self.key_name} is not set in the environment or .env")

    @staticmethod
    def _json(text: str) -> dict[str, Any]:
        start, end = text.find("{"), text.rfind("}")
        if start < 0 or end <= start:
            raise ValueError("model output did not contain a JSON object")
        return json.loads(text[start:end + 1])


class OpenAIWriter(HttpWriter):
    endpoint = "https://api.openai.com/v1/chat/completions"
    key_name = "OPENAI_API_KEY"
    name = "openai"

    def __init__(
        self,
        model: str,
        api_key: str | None = None,
        timeout: int = 90,
        endpoint: str | None = None,
        require_key: bool = True,
    ) -> None:
        super().__init__(model, api_key, timeout, endpoint, require_key)

    def write(self, context: dict[str, Any]) -> dict[str, Any]:
        import requests
        body = {
            "model": self.model,
            "temperature": 0,
            "response_format": {"type": "json_object"},
            "messages": [
                {"role": "system", "content": SYSTEM_PROMPT},
                {"role": "user", "content": json.dumps(context, ensure_ascii=False)[:120000]},
            ],
        }
        headers = {"Content-Type": "application/json"}
        if self.api_key:
            headers["Authorization"] = f"Bearer {self.api_key}"
        response = requests.post(self.endpoint, json=body, headers=headers, timeout=self.timeout)
        response.raise_for_status()
        return self._json(response.json()["choices"][0]["message"]["content"])


class AnthropicWriter(HttpWriter):
    endpoint = "https://api.anthropic.com/v1/messages"
    key_name = "ANTHROPIC_API_KEY"
    name = "anthropic"

    def __init__(
        self,
        model: str,
        api_key: str | None = None,
        timeout: int = 90,
        endpoint: str | None = None,
    ) -> None:
        super().__init__(model, api_key, timeout, endpoint, require_key=True)

    def write(self, context: dict[str, Any]) -> dict[str, Any]:
        import requests
        body = {
            "model": self.model, "max_tokens": 5000, "temperature": 0,
            "system": SYSTEM_PROMPT,
            "messages": [{"role": "user", "content": json.dumps(context, ensure_ascii=False)[:120000]}],
        }
        response = requests.post(
            self.endpoint, json=body,
            headers={"x-api-key": self.api_key, "anthropic-version": "2023-06-01"},
            timeout=self.timeout,
        )
        response.raise_for_status()
        return self._json(response.json()["content"][0]["text"])


def get_writer(
    backend: str = "scripted",
    model: str | None = None,
    api_key: str | None = None,
    endpoint: str | None = None,
) -> Writer:
    backend = backend.lower()
    if backend == "scripted":
        return ScriptedWriter()
    if backend == "openai":
        return OpenAIWriter(model or "gpt-4o-mini", api_key=api_key, endpoint=endpoint)
    if backend == "openai_compatible":
        if not endpoint:
            raise ValueError("An API endpoint is required for an OpenAI-compatible model")
        return OpenAIWriter(
            model or "local-model",
            api_key=api_key,
            endpoint=endpoint,
            require_key=False,
        )
    if backend == "anthropic":
        return AnthropicWriter(
            model or "claude-3-5-haiku-latest",
            api_key=api_key,
            endpoint=endpoint,
        )
    raise ValueError(
        "agent backend must be scripted, openai, openai_compatible or anthropic"
    )
