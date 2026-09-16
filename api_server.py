"""
api_server.py — 把 agentic 模块包成 HTTP 接口，供前端调用。

放到项目根目录（和 start.py 同级），运行：
    py -m pip install fastapi uvicorn
    py api_server.py                      # 默认 scripted 后端，完全免费
    py api_server.py --backend anthropic  # 付费，仅录 fixture / 演示时用

接口分两层，这是省钱的关键：
    POST /api/rank   纯统计排序，永远不花钱，前端 1→2 屏全靠它
    POST /api/ask    走 Agent（含 LLM 解释层），只有第 3 屏 AI 解读才调
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from dataclasses import asdict, is_dataclass
from datetime import date
from pathlib import Path
from typing import Any, Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from agentic import agent as agent_mod, llm, tools  # noqa: E402

# ─────────────────────────────────────────────────────────────
# 成本护栏
# ─────────────────────────────────────────────────────────────
MAX_LIVE_CALLS = 30                 # 付费后端每天上限
CACHE_DIR = ROOT / ".cache"
FIXTURE_DIR = ROOT / "fixtures"
COUNTER = CACHE_DIR / "usage.json"
CACHE_DIR.mkdir(exist_ok=True)

BACKEND = "scripted"
MODEL_ID: Optional[str] = None
GROUNDING = "full"
RECORD = False


def _calls_today() -> int:
    if not COUNTER.exists():
        return 0
    return json.loads(COUNTER.read_text()).get(str(date.today()), 0)


def _bump() -> int:
    d = json.loads(COUNTER.read_text()) if COUNTER.exists() else {}
    k = str(date.today())
    d[k] = d.get(k, 0) + 1
    COUNTER.write_text(json.dumps(d))
    return d[k]


def _key(*parts) -> str:
    raw = json.dumps(parts, sort_keys=True, ensure_ascii=False, default=str)
    return hashlib.sha256(raw.encode()).hexdigest()[:16]


def _plain(o: Any) -> Any:
    """dataclass / 对象 → 可 JSON 化的普通结构。"""
    if is_dataclass(o) and not isinstance(o, type):
        return asdict(o)
    if hasattr(o, "to_dict"):
        return o.to_dict()
    if isinstance(o, dict):
        return {k: _plain(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_plain(x) for x in o]
    return o


# ─────────────────────────────────────────────────────────────
# 单例：Toolbox 加载快照很贵，只做一次
# ─────────────────────────────────────────────────────────────
_toolbox = None
_agent = None


def toolbox():
    global _toolbox
    if _toolbox is None:
        print("加载快照中（首次约需十几秒）...", flush=True)
        _toolbox = tools.Toolbox()
        print("快照就绪", flush=True)
    return _toolbox


def agent():
    global _agent
    if _agent is None:
        _agent = agent_mod.Agent(
            toolbox=toolbox(),
            model=llm.get_model(BACKEND, MODEL_ID),
            grounding_mode=GROUNDING,
        )
    return _agent


# ─────────────────────────────────────────────────────────────
# 结构化条件 → 英文问句（Agent.run 只吃字符串）
# ─────────────────────────────────────────────────────────────
def build_question(q: "Query") -> str:
    joiner = " and " if (q.logic or "ALL").upper() == "ALL" else " or "
    s = "Which "
    if q.disease:
        s += f"{q.disease} "
    elif q.lineage:
        s += f"{q.lineage} "
    s += "cell lines express " + joiner.join(q.target_genes)
    if q.exclusion_genes:
        s += " but not " + " or ".join(q.exclusion_genes)
    return s + "?"


# ─────────────────────────────────────────────────────────────
# 归一化：转成前端契约
# ─────────────────────────────────────────────────────────────
def norm_candidate(c: dict) -> dict:
    per_gene = c.get("per_gene") or []
    has_protein = any(g.get("has_protein") for g in per_gene) if per_gene else None
    return {
        "rank": c.get("rank"),
        "cell_line": c.get("cell_line"),
        "depmap_id": c.get("depmap_id"),
        "meta": {
            "lineage": c.get("lineage"),
            "disease": c.get("primary_disease"),
        },
        "final_score": c.get("final_score"),
        # 分数分解：Candidate 上实际存在的三项
        "score_decomposition": {
            "biological": c.get("biological_score"),
            "confidence": c.get("confidence_score"),
            "penalty": c.get("exclusion_penalty") or 0.0,
        },
        "confidence_level": c.get("confidence_level"),
        "recommendation_level": c.get("recommendation_level"),
        "protein_measured": has_protein,
        "per_gene": [
            {
                "gene": g.get("gene"),
                "rna_percentile": g.get("rna_percentile"),
                "protein_percentile": g.get("protein_percentile"),
                "rna_sources": g.get("rna_sources"),
                "has_protein": g.get("has_protein"),
                "has_geo_support": g.get("has_geo_support"),
                "score": g.get("score"),
                "confidence": g.get("confidence"),
            }
            for g in per_gene
        ],
        "evidence_ids": c.get("evidence_ids") or [],
    }


def norm_ranking(payload: dict) -> dict:
    ranked = payload.get("ranked") or []
    return {
        "results": [norm_candidate(_plain(c)) for c in ranked],
        "excluded": payload.get("excluded") or [],
        "data_gaps": payload.get("data_gaps") or [],
        "discrimination": payload.get("discrimination") or {},
        "reference_populations": payload.get("reference_populations") or {},
        "query": payload.get("query") or {},
    }


def norm_answer(ans) -> dict:
    d = ans.to_dict() if hasattr(ans, "to_dict") else _plain(ans)
    recs = []
    for r in d.get("recommendations") or []:
        r = _plain(r)
        recs.append({
            "rank": r.get("rank"),
            "cell_line": r.get("cell_line"),
            "depmap_id": r.get("depmap_id"),
            "final_score": r.get("final_score"),
            "confidence_score": r.get("confidence_score"),
            "recommendation_level": r.get("recommendation_level"),
            "supporting_claims": [_plain(c) for c in (r.get("supporting_claims") or [])],
            "tradeoffs": [_plain(c) for c in (r.get("tradeoffs") or [])],
            "abstained_claims": r.get("abstained_claims") or [],
        })
    # grounding 统计藏在 trace 里 node == "ground" 的那条
    ground = [t for t in (d.get("trace") or []) if t.get("node") == "ground"]
    g = ground[-1] if ground else {}
    return {
        "query_interpretation": d.get("query_interpretation") or {},
        "recommendations": recs,
        "excluded_lines": d.get("excluded_lines") or [],
        "data_gaps": d.get("data_gaps") or [],
        "method_note": d.get("method_note") or "",
        "run_manifest": d.get("run_manifest") or {},
        "grounding": {
            "mode": g.get("mode"),
            "claims": g.get("claims"),
            "dropped": g.get("dropped"),
            "hallucination_rate": g.get("hallucination_rate"),
            "by_layer": g.get("by_layer"),
        },
    }


# 蛋白俗名 → 基因符号（蛋白矩阵按基因索引，必须先归一）
PROTEIN_ALIAS = {
    "HER2": "ERBB2", "HER-2": "ERBB2", "NEU": "ERBB2",
    "HER1": "EGFR", "ERBB1": "EGFR",
    "PD-L1": "CD274", "PDL1": "CD274", "PD-1": "PDCD1",
    "C-MET": "MET", "CMET": "MET",
    "C-KIT": "KIT", "CD117": "KIT",
    "P53": "TP53", "PTEN": "PTEN",
    "TRKA": "NTRK1", "ALK": "ALK",
}


def to_symbol(name: str) -> str:
    """蛋白名归一到基因符号。先查别名表，再试模块自带的解析。"""
    k = name.strip().upper()
    if k in PROTEIN_ALIAS:
        return PROTEIN_ALIAS[k]
    for fn in ("resolve_symbol", "normalise_symbol", "parse_symbol"):
        for mod in (llm, __import__("agentic.panel", fromlist=["*"])):
            if hasattr(mod, fn):
                try:
                    v = getattr(mod, fn)(k)
                    if v:
                        return str(v).upper()
                except Exception:
                    pass
    return k


# ─────────────────────────────────────────────────────────────
# HTTP
# ─────────────────────────────────────────────────────────────
app = FastAPI(title="CellLineSelector API")
app.add_middleware(CORSMiddleware, allow_origins=["*"],
                   allow_methods=["*"], allow_headers=["*"])


class Query(BaseModel):
    target_genes: list[str] = []
    require_protein: bool = False      # 目标蛋白模式：只保留有质谱测量的细胞系
    exclusion_genes: list[str] = []
    logic: str = "ALL"                 # ALL(几何平均) | ANY(noisy-OR)
    disease: Optional[str] = None
    lineage: Optional[str] = None
    top_n: int = 10


class AskQuery(Query):
    question: Optional[str] = None      # 给了就直接用，忽略上面的结构化字段


@app.get("/api/meta")
def meta():
    """顶栏数据版本 + 用量。只读文件，不花钱。"""
    from agentic import config
    out: dict = {"backend": BACKEND, "grounding": GROUNDING,
                 "live_calls_today": _calls_today(),
                 "live_calls_limit": None if BACKEND == "scripted" else MAX_LIVE_CALLS,
                 "cached": len([p for p in CACHE_DIR.glob("*.json") if p.name != "usage.json"])}
    mp = Path(config.MANIFEST) if hasattr(config, "MANIFEST") else ROOT / "snapshot/manifest.json"
    if Path(mp).exists():
        m = json.loads(Path(mp).read_text(encoding="utf-8"))
        out["snapshot"] = {
            "built_at": m.get("built_at"),
            "panel_genes": m.get("panel_genes"),
            "cell_lines": (m.get("cell_lines") or {}).get("rows"),
            "protein_lines": (m.get("protein") or {}).get("cell_lines"),
            "rna_lines": (m.get("rna_depmap") or {}).get("cell_lines"),
        }
    else:
        out["snapshot"] = None
    cp = ROOT / "snapshot/calibration.json"
    if cp.exists():
        c = json.loads(cp.read_text(encoding="utf-8"))
        out["calibration"] = {"scorer": c.get("scorer"), "alpha": c.get("alpha"),
                              "fitted_at": c.get("fitted_at")}
    return out


@app.post("/api/rank")
def rank(q: Query):
    """纯统计排序。不经过 LLM，永远免费。前端 1→2 屏用这个。"""
    if not q.target_genes:
        raise HTTPException(400, "至少需要一个 target gene")
    ck = CACHE_DIR / f"rank_{_key(q.model_dump())}.json"
    if ck.exists():
        out = json.loads(ck.read_text(encoding="utf-8"))
        out["_cached"] = True
        return out
    t0 = time.time()
    genes = [to_symbol(g) for g in q.target_genes]
    # 蛋白模式下要多取一些候选，因为随后要按"有蛋白测量"过滤
    want = q.top_n * 12 if q.require_protein else q.top_n
    try:
        res = toolbox().rank_cell_lines(
            target_genes=genes,
            logic=q.logic,
            exclusion_genes=[to_symbol(g) for g in q.exclusion_genes] or None,
            disease=q.disease,
            lineage=q.lineage,
            top_n=min(want, 500),
        )
        out = norm_ranking(_plain(res.payload))
    except Exception as e:
        raise HTTPException(500, f"{type(e).__name__}: {e}")

    if q.require_protein:
        total = len(out["results"])
        kept = [r for r in out["results"] if r.get("protein_measured")]
        out["results"] = kept[:q.top_n]
        # 保留统计模型给出的相对顺序，仅重新编号
        for i, r in enumerate(kept[:q.top_n], 1):
            r["display_rank"] = i
        out["protein_mode"] = {
            "requested": True,
            "candidates_before_filter": total,
            "candidates_with_protein": len(kept),
        }
        if not kept:
            out["data_gaps"].insert(0,
                f"蛋白层面查询：{total} 个候选中没有一个在 CCLE-Gygi 中测过该蛋白，"
                f"因此不作 RNA 代理推断。公开质谱覆盖 375 / 1479 个细胞系；"
                f"接入机构内部蛋白组数据可扩展覆盖范围。")
        else:
            out["data_gaps"].insert(0,
                f"蛋白层面查询：{total} 个候选中有 {len(kept)} 个经质谱实测，仅返回这些；"
                f"其余候选缺蛋白测量，未以 RNA 代理推断。"
                f"排序沿用统计模型原始顺序，前端不重新打分。")
    out.update({"_cached": False, "_elapsed_s": round(time.time() - t0, 2),
                "_cost": "free"})
    ck.write_text(json.dumps(out, ensure_ascii=False), encoding="utf-8")
    return out


@app.post("/api/ask")
def ask(q: AskQuery):
    """走完整 Agent，含解释断言与三层校验。付费后端下计入配额。"""
    question = q.question or build_question(q)
    ck = CACHE_DIR / f"ask_{_key(question, BACKEND, MODEL_ID, GROUNDING)}.json"
    if ck.exists():
        out = json.loads(ck.read_text(encoding="utf-8"))
        out["_cached"] = True
        return out

    if BACKEND != "scripted":
        used = _calls_today()
        if used >= MAX_LIVE_CALLS:
            raise HTTPException(429, f"今日实时调用已达上限 {MAX_LIVE_CALLS}（已用 {used}）。"
                                     f"改用 scripted 后端继续联调。")
        _bump()

    t0 = time.time()
    try:
        out = norm_answer(agent().run(question))
    except Exception as e:
        raise HTTPException(500, f"{type(e).__name__}: {e}")
    out.update({"_cached": False, "_question": question, "_backend": BACKEND,
                "_elapsed_s": round(time.time() - t0, 2)})
    ck.write_text(json.dumps(out, ensure_ascii=False), encoding="utf-8")
    if RECORD:
        FIXTURE_DIR.mkdir(exist_ok=True)
        tag = "_".join(q.target_genes) or "q"
        (FIXTURE_DIR / f"{tag}_{BACKEND}.json").write_text(
            json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    return out


def known_genes() -> set:
    """快照里认得的基因符号，parse_query 需要。"""
    tb = toolbox()
    for holder in (tb, getattr(tb, "store", None)):
        if holder is None:
            continue
        for attr in ("genes", "panel_genes", "gene_symbols", "panel"):
            if hasattr(holder, attr):
                v = getattr(holder, attr)
                if callable(v):
                    try:
                        v = v()
                    except Exception:
                        continue
                try:
                    return set(v)
                except TypeError:
                    continue
    return set()


class ParseBody(BaseModel):
    question: str


@app.post("/api/parse")
def parse(b: ParseBody):
    """自然语言 → 结构化条件。scripted 后端下是纯规则，免费。"""
    try:
        model = llm.get_model(BACKEND, MODEL_ID)
        parsed = model.parse_query(b.question, known_genes())
        disease, lineage = llm.disease_in(b.question)
        out = _plain(parsed) if isinstance(parsed, dict) else {"raw": str(parsed)}
        out.setdefault("disease", disease)
        out.setdefault("lineage", lineage)
        out["_backend"] = BACKEND
        return out
    except Exception as e:
        raise HTTPException(500, f"{type(e).__name__}: {e}")


@app.get("/api/availability")
def availability(genes: str):
    """某些基因有哪些模态的数据。免费，用于输入时提示。"""
    gl = [g.strip().upper() for g in genes.split(",") if g.strip()]
    try:
        return _plain(toolbox().check_data_availability(genes=gl).payload)
    except Exception as e:
        raise HTTPException(500, f"{type(e).__name__}: {e}")


@app.post("/api/evidence")
def evidence(body: dict):
    """按 evidence_id 取证据卡，第 3 屏详情用。免费。"""
    ids = body.get("evidence_ids") or []
    if not ids:
        raise HTTPException(400, "需要 evidence_ids")
    try:
        r = toolbox().get_evidence(evidence_ids=ids, gene=body.get("gene"))
        return {"payload": _plain(r.payload), "cards": [_plain(c) for c in r.cards]}
    except Exception as e:
        raise HTTPException(500, f"{type(e).__name__}: {e}")


@app.get("/api/usage")
def usage():
    return {"backend": BACKEND, "live_calls_today": _calls_today(),
            "limit": None if BACKEND == "scripted" else MAX_LIVE_CALLS,
            "cached_queries": len([p for p in CACHE_DIR.glob("*.json") if p.name != "usage.json"])}


if (ROOT / "web").is_dir():
    app.mount("/", StaticFiles(directory=str(ROOT / "web"), html=True), name="web")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", default="scripted",
                    choices=["scripted", "openai", "anthropic"])
    ap.add_argument("--model", default=None)
    ap.add_argument("--grounding", default="full", choices=["none", "l1l2", "full"])
    ap.add_argument("--record", action="store_true", help="结果存到 fixtures/")
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--warm", action="store_true", help="启动时就加载快照")
    a = ap.parse_args()
    BACKEND, MODEL_ID, GROUNDING, RECORD = a.backend, a.model, a.grounding, a.record

    if BACKEND != "scripted":
        key, where = llm.api_key_for(
            "ANTHROPIC_API_KEY" if BACKEND == "anthropic" else "OPENAI_API_KEY")
        print(f"\n[!] 付费后端 {BACKEND}，key 来自 {where}，"
              f"今日已用 {_calls_today()}/{MAX_LIVE_CALLS}。缓存命中不计费。\n")
    else:
        print("\n[ok] scripted 后端：离线、确定性、零成本。\n")

    if a.warm:
        toolbox()

    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=a.port)
