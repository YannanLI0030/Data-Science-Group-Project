"""Local HTTP application for the canonical CellLineSelector backend.

The deterministic recommender owns filtering, scores and rank order.  The
optional language model receives an immutable saved run and can only produce
grounded, cited prose after ranking has finished.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import re
import threading
import time
import uuid
import webbrowser
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlparse

from dynamic_cellline_selector_gene_protein import (
    CONFIDENCE_WEIGHT,
    DEFAULT_CACHE_DIR,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_RAW_DATA_DIR,
    MAX_EXCLUSION_PENALTY,
    PROTEIN_ONLY_CONFIDENCE_WEIGHT,
    PROTEIN_ONLY_PROTEIN_WEIGHT,
    PROTEIN_WEIGHT,
    RNA_WEIGHT,
    DynamicMultiOmicsRecommender,
    build_disease_aliases,
    build_reason,
    determine_query_mode,
    save_ranked_results,
    score_candidates,
)
from src.agentic.output_agent import OutputAgent, ranking_fingerprint


ROOT = Path(__file__).resolve().parent
WEB_DIR = ROOT / "web"
RUN_DIR = ROOT / "runs" / "ui"
OUTPUT_DIR = Path(DEFAULT_OUTPUT_DIR)
MAX_BODY_BYTES = 1_000_000
MAX_TOP_N = 50
SYMBOL_RE = re.compile(r"^[A-Za-z0-9_.-]{1,64}$")
RUN_ID_RE = re.compile(r"^[a-f0-9]{32}$")
PROTEIN_ALIASES = {
    "HER2": "ERBB2",
    "HER-2": "ERBB2",
    "NEU": "ERBB2",
    "HER1": "EGFR",
    "ERBB1": "EGFR",
    "PD-L1": "CD274",
    "PDL1": "CD274",
    "PD-1": "PDCD1",
    "C-MET": "MET",
    "CMET": "MET",
    "C-KIT": "KIT",
    "CD117": "KIT",
    "P53": "TP53",
    "TRKA": "NTRK1",
}

_recommender: DynamicMultiOmicsRecommender | None = None
_data_lock = threading.RLock()


def _jsonable(value: Any) -> Any:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(v) for v in value]
    if hasattr(value, "item"):
        try:
            return _jsonable(value.item())
        except Exception:
            pass
    return str(value)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(
        json.dumps(_jsonable(payload), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    temp.replace(path)


def _clean_symbol(value: Any, label: str) -> str | None:
    text = str(value or "").strip().upper()
    if not text:
        return None
    if not SYMBOL_RE.fullmatch(text):
        raise ValueError(f"{label} may contain only letters, numbers, periods, underscores, or hyphens")
    return text


def _clean_query(body: dict[str, Any]) -> dict[str, Any]:
    gene = _clean_symbol(body.get("target_gene"), "Target gene")
    protein = _clean_symbol(body.get("target_protein"), "Target protein")
    protein = PROTEIN_ALIASES.get(protein, protein)
    exclusion = _clean_symbol(body.get("exclusion_gene"), "Exclusion gene")
    disease = str(body.get("disease") or "").strip()
    if not gene and not protein:
        raise ValueError("Enter at least one target gene or target protein")
    if not disease:
        raise ValueError("Disease or tissue is required")
    if len(disease) > 120:
        raise ValueError("The disease or tissue name is too long")
    try:
        top_n = int(body.get("top_n", 10))
    except (TypeError, ValueError):
        raise ValueError("The result count must be an integer") from None
    if not 1 <= top_n <= MAX_TOP_N:
        raise ValueError(f"The result count must be between 1 and {MAX_TOP_N}.")
    return {
        "target_gene": gene,
        "target_protein": protein,
        "exclusion_gene": exclusion,
        "disease": disease,
        "top_n": top_n,
    }


def _weights(query_mode: str) -> dict[str, Any]:
    if query_mode == "PROTEIN_ONLY":
        return {
            "rna": 0.0,
            "protein": PROTEIN_ONLY_PROTEIN_WEIGHT,
            "confidence": PROTEIN_ONLY_CONFIDENCE_WEIGHT,
            "exclusion_penalty_max": MAX_EXCLUSION_PENALTY,
            "note": "Protein is the only directly scored biological modality: 0.85× protein + 0.15× confidence. RNA contributes only supporting evidence to confidence.",
        }
    return {
        "rna": RNA_WEIGHT,
        "protein": PROTEIN_WEIGHT,
        "confidence": CONFIDENCE_WEIGHT,
        "exclusion_penalty_max": MAX_EXCLUSION_PENALTY,
        "note": "RNA and protein comprise the 0.85 biological-evidence weight. The weight of a missing modality is redistributed across the available biological modalities.",
    }


def _gaps(row: dict[str, Any], query_mode: str, target: str) -> list[str]:
    gaps: list[str] = []
    prefix = "supporting " if query_mode == "PROTEIN_ONLY" else ""
    if not row.get("hasDepMapRNA"):
        gaps.append(f"Missing {prefix}DepMap RNA-seq evidence.")
    if not row.get("hasHpaRNA"):
        gaps.append(f"Missing {prefix}HPA RNA-seq evidence.")
    if not row.get("hasGeoRNA"):
        gaps.append(f"Missing {prefix}GEO RNA-expression evidence.")
    if not row.get("hasProteomics"):
        gaps.append(f"Missing CCLE-Gygi proteomics evidence for {target}.")
    return gaps


def get_recommender() -> DynamicMultiOmicsRecommender:
    global _recommender
    with _data_lock:
        if _recommender is None:
            print("[backend] Initializing the dynamic multi-omics data layer...", flush=True)
            _recommender = DynamicMultiOmicsRecommender(
                Path(DEFAULT_RAW_DATA_DIR), Path(DEFAULT_CACHE_DIR), False
            )
            print("[backend] Dynamic multi-omics data layer is ready.", flush=True)
        return _recommender


def recommend(body: dict[str, Any]) -> dict[str, Any]:
    query = _clean_query(body)
    mode = determine_query_mode(query["target_gene"], query["target_protein"])
    aliases = build_disease_aliases(query["disease"])
    t0 = time.time()

    with _data_lock:
        data = get_recommender()
        candidates = data.fetch_candidate_evidence(
            target_gene=query["target_gene"],
            target_protein=query["target_protein"],
            disease_aliases=aliases,
            exclusion_gene=query["exclusion_gene"],
        )
        if not candidates:
            raise LookupError("No cell line satisfies both the disease hard filter and the target-evidence requirements")
        scored = score_candidates(candidates, query_mode=mode)
        if not scored:
            raise LookupError("Candidate cell lines lack the scoreable evidence required for this query mode")

        top_rows = scored[: query["top_n"]]
        top = top_rows[0]
        alternatives = data.fetch_similar_cell_lines(scored, 5)
        trace = data.fetch_evidence_trace(top)
        supplementary = data.get_supplementary_context(top["DepMap_ID"], 5)

    run_id = uuid.uuid4().hex
    request = {
        **query,
        "disease_aliases": aliases,
        "query_mode": mode,
        "output_language": "en",
    }
    target = query["target_protein"] or query["target_gene"] or "target"
    fingerprint = ranking_fingerprint(top_rows)
    payload = {
        "run_id": run_id,
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "query": request,
        "query_mode": mode,
        "weights": _weights(mode),
        "results": top_rows,
        "total_candidates": len(scored),
        "top_recommendation": top,
        "main_reasons": build_reason(
            top,
            query["target_gene"],
            query["target_protein"],
            query["exclusion_gene"],
            query["disease"],
            mode,
        ),
        "evidence_trace": trace,
        "data_gaps": _gaps(top, mode, target),
        "alternatives": alternatives,
        "supplementary": supplementary,
        "ranking_fingerprint": fingerprint,
        "ranking_authority": "dynamic_cellline_selector_gene_protein.score_candidates",
        "ranking_unchanged": True,
        "elapsed_seconds": round(time.time() - t0, 2),
    }
    payload = _jsonable(payload)
    _write_json(RUN_DIR / f"{run_id}.json", payload)
    save_ranked_results(scored, OUTPUT_DIR / f"{run_id}_ranking.csv")
    save_ranked_results(scored, OUTPUT_DIR / "dynamic_ranked_recommendations.csv")
    return payload


def _load_run(run_id: str) -> dict[str, Any]:
    if not RUN_ID_RE.fullmatch(run_id):
        raise ValueError("Invalid run ID")
    path = RUN_DIR / f"{run_id}.json"
    if not path.exists():
        raise FileNotFoundError("This recommendation run was not found; run the query again")
    return json.loads(path.read_text(encoding="utf-8"))


def candidate_detail(run_id: str, rank: int) -> dict[str, Any]:
    run = _load_run(run_id)
    row = next((x for x in run["results"] if int(x.get("rank", -1)) == rank), None)
    if row is None:
        raise LookupError("This rank is not included in the currently displayed results")
    with _data_lock:
        data = get_recommender()
        trace = data.fetch_evidence_trace(row)
        supplementary = data.get_supplementary_context(row["DepMap_ID"], 5)
    mode = run["query_mode"]
    target = run["query"].get("target_protein") or run["query"].get("target_gene") or "target"
    return _jsonable({
        "run_id": run_id,
        "candidate": row,
        "weights": run["weights"],
        "evidence_trace": trace,
        "data_gaps": _gaps(row, mode, target),
        "supplementary": supplementary,
        "ranking_fingerprint": run["ranking_fingerprint"],
        "ranking_unchanged": True,
    })


def _model_options(raw: Any) -> dict[str, Any]:
    raw = raw if isinstance(raw, dict) else {}
    provider = str(raw.get("provider") or "scripted").strip().lower()
    allowed = {"scripted", "openai", "openai_compatible", "anthropic"}
    if provider not in allowed:
        raise ValueError("Unsupported language-model provider")
    grounding = str(raw.get("grounding") or "l1l2").strip().lower()
    if grounding not in {"none", "l1l2", "full"}:
        raise ValueError("grounding must be none, l1l2, or full")
    model = str(raw.get("model") or "").strip() or None
    api_key = str(raw.get("api_key") or "").strip() or None
    endpoint = str(raw.get("endpoint") or "").strip() or None
    if endpoint:
        parsed = urlparse(endpoint)
        if parsed.scheme not in {"http", "https"} or not parsed.netloc:
            raise ValueError("The API endpoint must be a valid http(s) URL")
        if parsed.username or parsed.password:
            raise ValueError("Do not include the API key in the URL")
    if provider == "openai_compatible" and not endpoint:
        raise ValueError("An OpenAI-compatible provider requires a complete API endpoint")
    return {
        "provider": provider,
        "model": model,
        "api_key": api_key,
        "endpoint": endpoint,
        "grounding": grounding,
    }


def explain(body: dict[str, Any]) -> dict[str, Any]:
    run_id = str(body.get("run_id") or "")
    run = _load_run(run_id)
    current = ranking_fingerprint(run["results"])
    if current != run["ranking_fingerprint"]:
        raise RuntimeError("The saved ranking failed its integrity check, so the Agent stopped")
    opts = _model_options(body.get("model_config"))
    agent = OutputAgent(
        backend=opts["provider"],
        model=opts["model"],
        grounding_mode=opts["grounding"],
        api_key=opts["api_key"],
        endpoint=opts["endpoint"],
    )
    output = agent.run(
        request=run["query"],
        ranked_rows=run["results"],
        evidence_trace=run["evidence_trace"],
        alternatives=run["alternatives"],
        supplementary=run["supplementary"],
    ).to_dict()
    output["run_id"] = run_id
    output["ranking_unchanged"] = (
        output["ranking_fingerprint"] == run["ranking_fingerprint"]
    )
    output = _jsonable(output)
    _write_json(
        RUN_DIR / f"{run_id}_agent_{int(time.time())}.json",
        output,
    )
    return output


def ranking_csv(run_id: str) -> bytes:
    run = _load_run(run_id)
    rows = run["results"]
    if not rows:
        return b""
    keys = [
        "rank", "DepMap_ID", "cellLine", "lineage", "disease", "queryMode",
        "targetGene", "targetProtein", "finalScore", "confidenceScore",
        "recommendationLevel", "rnaScore", "proteinScore", "biologicalScore",
        "exclusionPenalty", "completenessScore", "sourceSupportScore",
    ]
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=keys, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    return ("\ufeff" + stream.getvalue()).encode("utf-8")


class CellLineHandler(SimpleHTTPRequestHandler):
    server_version = "CellLineSelector/1.0"

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, directory=str(WEB_DIR), **kwargs)

    def _send_json(self, payload: Any, status: int = 200) -> None:
        raw = json.dumps(_jsonable(payload), ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(raw)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(raw)

    def _error(self, status: int, message: str) -> None:
        self._send_json({"error": message}, status)

    def _body(self) -> dict[str, Any]:
        try:
            length = int(self.headers.get("Content-Length", "0"))
        except ValueError:
            raise ValueError("Invalid request length") from None
        if length <= 0 or length > MAX_BODY_BYTES:
            raise ValueError("The request body is empty or too large")
        raw = self.rfile.read(length)
        value = json.loads(raw.decode("utf-8"))
        if not isinstance(value, dict):
            raise ValueError("The request must be a JSON object")
        return value

    def do_POST(self) -> None:  # noqa: N802
        try:
            body = self._body()
            if self.path == "/api/recommend":
                self._send_json(recommend(body))
            elif self.path == "/api/explain":
                self._send_json(explain(body))
            else:
                self._error(HTTPStatus.NOT_FOUND, "Endpoint not found")
        except (ValueError, json.JSONDecodeError) as exc:
            self._error(HTTPStatus.BAD_REQUEST, str(exc))
        except FileNotFoundError as exc:
            self._error(HTTPStatus.NOT_FOUND, str(exc))
        except LookupError as exc:
            self._error(HTTPStatus.UNPROCESSABLE_ENTITY, str(exc))
        except Exception as exc:
            print(f"[api] {type(exc).__name__}: {exc}", flush=True)
            self._error(HTTPStatus.INTERNAL_SERVER_ERROR, f"{type(exc).__name__}: {exc}")

    def do_GET(self) -> None:  # noqa: N802
        try:
            path = unquote(urlparse(self.path).path)
            if path == "/api/health":
                self._send_json({"status": "ok", "ranking_authority": "score_candidates"})
                return
            if path == "/api/meta":
                self._send_json({
                    "app": "CellLineSelector Final",
                    "providers": ["scripted", "openai", "openai_compatible", "anthropic"],
                    "default_provider": "scripted",
                    "api_keys_persisted": False,
                    "protein_only_weights": {"protein": 0.85, "confidence": 0.15, "rna": 0.0},
                })
                return
            match = re.fullmatch(r"/api/runs/([a-f0-9]{32})/candidates/(\d+)", path)
            if match:
                self._send_json(candidate_detail(match.group(1), int(match.group(2))))
                return
            match = re.fullmatch(r"/api/runs/([a-f0-9]{32})/ranking\.csv", path)
            if match:
                raw = ranking_csv(match.group(1))
                self.send_response(HTTPStatus.OK)
                self.send_header("Content-Type", "text/csv; charset=utf-8")
                self.send_header(
                    "Content-Disposition",
                    f'attachment; filename="cellline_ranking_{match.group(1)[:8]}.csv"',
                )
                self.send_header("Content-Length", str(len(raw)))
                self.end_headers()
                self.wfile.write(raw)
                return
            if path.startswith("/api/"):
                self._error(HTTPStatus.NOT_FOUND, "Endpoint not found")
                return
            super().do_GET()
        except (ValueError, FileNotFoundError) as exc:
            self._error(HTTPStatus.NOT_FOUND, str(exc))
        except LookupError as exc:
            self._error(HTTPStatus.UNPROCESSABLE_ENTITY, str(exc))
        except Exception as exc:
            print(f"[api] {type(exc).__name__}: {exc}", flush=True)
            self._error(HTTPStatus.INTERNAL_SERVER_ERROR, f"{type(exc).__name__}: {exc}")

    def log_message(self, fmt: str, *args: Any) -> None:
        print(f"[http] {self.address_string()} - {fmt % args}", flush=True)


def run_server(host: str = "127.0.0.1", port: int = 8000, open_ui: bool = True) -> None:
    if not WEB_DIR.is_dir():
        raise FileNotFoundError(f"UI directory not found: {WEB_DIR}")
    server = ThreadingHTTPServer((host, port), CellLineHandler)
    server.daemon_threads = True
    url = f"http://{host}:{port}/"
    print("=" * 72)
    print("CellLineSelector Final has started")
    print(f"Open: {url}")
    print("score_candidates() determines the ranking; the Agent explains it but cannot reorder it.")
    print("Press Ctrl+C to stop.")
    print("=" * 72, flush=True)
    if open_ui:
        timer = threading.Timer(0.8, lambda: webbrowser.open(url))
        timer.daemon = True
        timer.start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nServer stopped.")
    finally:
        server.server_close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--no-browser", action="store_true")
    args = parser.parse_args()
    run_server(args.host, args.port, not args.no_browser)
