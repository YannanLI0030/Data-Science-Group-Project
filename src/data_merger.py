"""
data_merger.py  (v3)
====================
CellLineSelector —— 多数据源统一合并层 (DepMap + HPA + GEO)

v3 修复了 v2 在真实数据上暴露的问题
------------------------------------
1. [BUG] 列名匹配太死板 —— 你的 Cellosaurus 列叫 'Accession (CVCL_xxxx)' 而不是
   'Accession'。现在 `_find_col` 改用"规范化后模糊包含"匹配, 能吃掉括号后缀。
2. [BUG] 152,231 行 Cellosaurus 用 iterrows 会跑几分钟 —— 改为向量化, 秒级完成。
3. [严重] v2 的 build_long_table 会把整个 DepMap RNA 矩阵 melt 成 8000 万行,
   加上 HPA 的 2400 万行会直接 OOM。v3 把"目标基因集合"下推给每个数据源,
   源只加载需要的基因 -> 内存从 GB 级降到 MB 级。
4. [BUG] GEO File 3 实际是宽表 (Gene + 每个 GSM 一列), v2 假设长表。v3 自动识别。
5. [新增] 自动探测 DepMap sample_info 里的 RRID/Cellosaurus 列 —— 如果存在,
   它是比 Cellosaurus Cross-references 更直接的 CVCL<->ACH 桥。
6. [新增] 直接运行本文件会执行自检 (python data_merger.py <data_dir>)。

四条硬性要求的对应关系
------------------------------
1. 多数据源统一 ID 对齐  ←  CellLineIDResolver (CVCL 为枢纽)
2. 多组学统一结构化存储  ←  UNIFIED_COLUMNS 长表 + build_gene_table 宽表
3. 缺失值规范化处理      ←  外连接骨架 + has_<layer> 掩码 + redistribute_weights
4. 数据可扩展性预留      ←  OmicsSource 注册表模式
"""

from __future__ import annotations

__version__ = "3.1"

import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterable, Optional, Set

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# 0. 统一 schema
# ---------------------------------------------------------------------------

UNIFIED_COLUMNS = [
    "DepMap_ID", "gene_symbol", "ensembl_id",
    "omics_layer", "value_raw", "value_unit",
    "source", "source_version",
]

OMICS_LAYERS = ("rna", "protein", "mutation", "fusion", "mirna", "metabolite")


# ---------------------------------------------------------------------------
# 1. 通用工具：柔性列名匹配
# ---------------------------------------------------------------------------

def _norm_colname(c: str) -> str:
    """列名规范化: 小写、去掉非字母数字。'Accession (CVCL_xxxx)' -> 'accessioncvclxxxx'"""
    return re.sub(r"[^a-z0-9]", "", str(c).lower())


def _find_col(df: pd.DataFrame, keywords: list[str], what: str = "列") -> str:
    """按"规范化后包含关键词"找列, 容忍括号后缀/大小写/空格差异。

    先找完全相等的, 再找以关键词开头的, 最后找包含关键词的。
    """
    norm_map = {_norm_colname(c): c for c in df.columns}
    keys = [_norm_colname(k) for k in keywords]
    # 1) 精确
    for k in keys:
        if k in norm_map:
            return norm_map[k]
    # 2) 前缀
    for k in keys:
        for nc, orig in norm_map.items():
            if nc.startswith(k):
                return orig
    # 3) 包含
    for k in keys:
        for nc, orig in norm_map.items():
            if k in nc:
                return orig
    raise KeyError(
        f"找不到{what}: 期望列名含 {keywords} 之一。\n"
        f"实际列: {list(df.columns)[:15]}"
    )


def _find_col_optional(df: pd.DataFrame, keywords: list[str]) -> Optional[str]:
    try:
        return _find_col(df, keywords)
    except KeyError:
        return None


# ---------------------------------------------------------------------------
# 2. CellLineIDResolver —— 以 CVCL 为枢纽的 ID 对齐 (向量化)
# ---------------------------------------------------------------------------

class CellLineIDResolver:
    """把任意细胞系标识符 (CVCL / 细胞系名 / 别名 / CCLE_Name) 解析为 DepMap_ID。

    两张查找表:
      * cvcl_to_ach  : {'CVCL_0023': 'ACH-000681', ...}
      * name_to_cvcl : {规范化名字 -> 'CVCL_0023', ...}   (来自 Identifier + Synonyms)

    CVCL->ACH 的来源有两条, 按可靠性排序:
      (A) DepMap sample_info 自带的 RRID / Cellosaurus 列 —— 若存在, 这是最权威的;
      (B) Cellosaurus 'Cross-references' 里的 'DepMap; ACH-xxxxxx' 交叉引用。
    两条都会被使用, (A) 优先。
    """

    def __init__(self, cellosaurus: pd.DataFrame,
                 sample_info: Optional[pd.DataFrame] = None,
                 verbose: bool = True):
        acc_col = _find_col(cellosaurus, ["accession", "cvcl"], "Cellosaurus 的 CVCL 列")
        id_col = _find_col(cellosaurus, ["identifier", "cellline", "name"],
                           "Cellosaurus 的细胞系名列")
        xref_col = _find_col_optional(cellosaurus, ["crossreferences", "crossreference", "xref"])
        syn_col = _find_col_optional(cellosaurus, ["synonyms", "synonym"])

        if verbose:
            print(f"[resolver] Cellosaurus 列识别: CVCL={acc_col!r}, "
                  f"name={id_col!r}, xref={xref_col!r}, syn={syn_col!r}")

        cvcl = cellosaurus[acc_col].astype("string").str.strip()
        valid = cvcl.str.startswith("CVCL_", na=False)

        # ---------- 路径 A: sample_info 自带 CVCL 列 (最可靠) ----------
        self.cvcl_to_ach: dict[str, str] = {}
        n_from_sample_info = 0
        if sample_info is not None:
            si_cvcl_col = _find_col_optional(sample_info, ["rrid", "cellosaurus", "cvcl"])
            if si_cvcl_col is not None:
                si_cvcl = sample_info[si_cvcl_col].astype("string").str.strip()
                # RRID 常写成 'CVCL_0023' 或 'RRID:CVCL_0023'
                si_cvcl = si_cvcl.str.extract(r"(CVCL_[A-Za-z0-9]+)", expand=False)
                pairs = pd.DataFrame({"cvcl": si_cvcl,
                                      "ach": sample_info["DepMap_ID"].astype("string")}).dropna()
                self.cvcl_to_ach.update(dict(zip(pairs["cvcl"], pairs["ach"])))
                n_from_sample_info = len(self.cvcl_to_ach)
                if verbose:
                    print(f"[resolver] 路径 A: sample_info 的 {si_cvcl_col!r} 列提供 "
                          f"{n_from_sample_info} 条 CVCL->ACH (最权威)")
            elif verbose:
                print("[resolver] 路径 A: sample_info 里没找到 RRID/Cellosaurus 列, 跳过")

        # ---------- 路径 B: Cellosaurus Cross-references ----------
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
                print(f"[resolver] 路径 B: Cellosaurus Cross-references 提供 "
                      f"{n_from_xref} 条 CVCL->ACH "
                      f"(新增 {len(self.cvcl_to_ach) - before} 条)")

        if not self.cvcl_to_ach:
            print("[resolver] ⚠️  警告: 一条 CVCL->ACH 映射都没建立!\n"
                  "    请检查: (a) sample_info 是否有 RRID 列; "
                  "(b) Cellosaurus Cross-references 里是否真的有 'DepMap; ACH-xxx'。\n"
                  "    可以运行 diagnose_cellosaurus() 看看 xref 长什么样。")

        # ---------- 名字 -> CVCL 字典 (向量化) ----------
        names = cellosaurus.loc[valid, id_col].astype("string")
        cvcl_valid = cvcl[valid]
        self._name_to_cvcl: dict[str, str] = {}
        self._bulk_register(names, cvcl_valid)

        if syn_col is not None:
            syn = cellosaurus.loc[valid, syn_col].astype("string")
            # 同义词按 ';' 或 '||' 分隔 -> explode 成一行一别名
            exploded = (syn.str.split(r"\s*(?:\|\||;)\s*", regex=True)
                           .explode())
            exploded_cvcl = cvcl_valid.reindex(exploded.index)
            self._bulk_register(exploded, exploded_cvcl)

        # ---------- 把 sample_info 里的名字也挂到 CVCL 上 (供 GEO 用) ----------
        if sample_info is not None:
            ach_to_cvcl = {a: c for c, a in self.cvcl_to_ach.items()}
            si_cvcl_back = sample_info["DepMap_ID"].astype("string").map(ach_to_cvcl)
            for col in ("cell_line_name", "stripped_cell_line_name", "CCLE_Name"):
                if col in sample_info.columns:
                    self._bulk_register(sample_info[col].astype("string"), si_cvcl_back)

        if verbose:
            print(f"[resolver] 完成: CVCL->ACH {len(self.cvcl_to_ach)} 条, "
                  f"name->CVCL {len(self._name_to_cvcl)} 条")

    # ---- 批量注册名字 (向量化, 不用 iterrows) ----
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

    # ---- 公共 API ----
    def resolve_cvcl_to_ach(self, cvcl) -> Optional[str]:
        """HPA 主路径: File 11 已给 CVCL。"""
        if isinstance(cvcl, str) and cvcl.startswith("CVCL_"):
            return self.cvcl_to_ach.get(cvcl)
        return None

    def resolve_name(self, name) -> Optional[str]:
        """GEO 兜底路径: 只有名字 -> CVCL -> ACH。"""
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
    """打印 Cellosaurus Cross-references 里出现的资源名, 确认 DepMap 是否在列。"""
    xref_col = _find_col_optional(cellosaurus, ["crossreferences", "xref"])
    if xref_col is None:
        print("没有 Cross-references 列")
        return
    s = cellosaurus[xref_col].dropna().astype(str)
    print(f"Cross-references 非空行数: {len(s):,} / {len(cellosaurus):,}")
    # 抽取所有 'Resource;' 记号统计频次
    resources = (s.str.findall(r"([A-Za-z0-9_]+)\s*[;=]")
                  .explode().dropna())
    print(f"\n出现最多的 {n} 个交叉引用资源:")
    print(resources.value_counts().head(n).to_string())
    has_depmap = resources.str.lower().eq("depmap").any()
    print(f"\n是否包含 DepMap 引用: {has_depmap}")
    if has_depmap:
        sample = s[s.str.contains("depmap", case=False, na=False)].head(3)
        print("\nDepMap 引用样例:")
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
        """从 HPA File 1 建 symbol<->ENSG 字典。

        HPA File 1 是长表 (~2400 万行), 所以这里 **只读两列** 并去重,
        而不是 nrows=100000 (那样只能覆盖 ~83 个基因)。
        """
        head = pd.read_csv(path, sep="\t", nrows=5)
        ensg_col = _find_col(head, list(ensg_kw), "HPA 的 Gene(ENSG) 列")
        sym_col = _find_col(head, list(sym_kw), "HPA 的 Gene name 列")
        df = pd.read_csv(path, sep="\t", usecols=[ensg_col, sym_col],
                         dtype="string").drop_duplicates()
        self._bulk(df[ensg_col], df[sym_col])
        if verbose:
            print(f"[gene_res] 从 HPA 收录 {len(self._sym_to_ensg):,} 个 gene symbol")
        return self

    def to_symbol(self, ensg):
        return self._ensg_to_sym.get(ensg.split(".")[0]) if isinstance(ensg, str) else None

    def to_ensembl(self, symbol):
        return self._sym_to_ensg.get(symbol) if isinstance(symbol, str) else None


# ---------------------------------------------------------------------------
# 4. 数据源抽象接口  (v3: 支持基因过滤下推)
# ---------------------------------------------------------------------------

class OmicsSource(ABC):
    source: str
    layer: str
    version: str

    @abstractmethod
    def load_long(self, resolver: "CellLineIDResolver",
                  gene_resolver: "GeneIDResolver",
                  genes: Optional[Set[str]] = None) -> pd.DataFrame:
        """返回 UNIFIED_COLUMNS 长表。

        genes : 若给出, 数据源**只需**返回这些基因的数据 (内存关键)。
                实现时应尽量在读取阶段就过滤, 而不是读完再筛。
        """
        ...

    def _pack(self, df: pd.DataFrame) -> pd.DataFrame:
        out = df.reindex(columns=UNIFIED_COLUMNS)
        out["source"] = self.source
        out["source_version"] = self.version
        out["omics_layer"] = self.layer
        n_before = len(out)
        # 只丢主键缺失的行; 值为 NaN 的保留 (蛋白质谱未检测要留)
        out = out[out["DepMap_ID"].notna() & out["gene_symbol"].notna()]
        dropped = n_before - len(out)
        if dropped:
            print(f"  [{self.source}] {dropped:,} 行因 ID 未解析被跳过 "
                  f"({dropped / max(n_before, 1):.1%})")
        return out.reset_index(drop=True)


# ---------------------------------------------------------------------------
# 5. 具体适配器
# ---------------------------------------------------------------------------

def _split_symbol_ensg(columns: pd.Index) -> pd.DataFrame:
    """把 'EGFR (ENSG00000146648)' 这类列头拆成 symbol + ensembl_id。"""
    meta = pd.DataFrame({"orig": list(columns)})
    ext = meta["orig"].str.extract(r"^\s*(?P<sym>[^(]+?)\s*\((?P<ensg>[^)]+)\)\s*$")
    meta["gene_symbol"] = ext["sym"].fillna(meta["orig"].str.strip())
    meta["ensembl_id"] = ext["ensg"]
    return meta


class DepMapRNASource(OmicsSource):
    """DepMap RNA (复用 data_loader 的宽矩阵, 行=ACH, 列='SYMBOL (ENSG)')。

    v3: 只 melt 需要的基因列 —— 避免 8000 万行的全量 melt。
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
    """CCLE-Gygi 蛋白 (行=ACH, 列='UNIPROT (SYMBOL)')。NaN 保留。"""
    layer, source = "protein", "ccle_gygi_protein"

    def __init__(self, loader, version: str = "CCLE-Gygi-2020"):
        self.loader, self.version = loader, version
        self._colmeta: Optional[pd.DataFrame] = None

    def _meta(self, wide) -> pd.DataFrame:
        if self._colmeta is None:
            meta = pd.DataFrame({"orig": list(wide.columns)})
            # 'P00533 (EGFR)' -> symbol 在括号里
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
    """HPA 细胞系 RNA: File 1 (长表) + File 11 (描述表, 提供 CVCL)。

    ID 路径 (两跳精确查表, 不做字符串归一化):
        File 1 'Cell line' --File 11--> CVCL_xxxx --Cellosaurus--> ACH-xxxxxx

    v3: File 1 有约 2400 万行, 所以用 chunk 读取 + 按 genes 过滤,
        内存占用与目标基因数成正比, 而不是与文件大小成正比。
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
            name_col = _find_col(desc, ["cellline"], "File 11 的 Cell line 列")
            cvcl_col = _find_col(desc, ["cellosaurusid", "cellosaurus", "cvcl"],
                                 "File 11 的 Cellosaurus ID 列")
            d = pd.DataFrame({"n": desc[name_col].str.strip(),
                              "c": desc[cvcl_col].str.strip()}).dropna(subset=["n"])
            self._name_to_cvcl = dict(zip(d["n"], d["c"].fillna("")))
            if verbose:
                n_ok = sum(1 for v in self._name_to_cvcl.values()
                           if isinstance(v, str) and v.startswith("CVCL_"))
                print(f"  [hpa_rna] File 11: {len(self._name_to_cvcl)} 细胞系, "
                      f"{n_ok} 有 CVCL_id")
        return self._name_to_cvcl

    def load_long(self, resolver, gene_resolver, genes=None):
        name_to_cvcl = self._load_desc()

        head = pd.read_csv(self.expr_path, sep="\t", nrows=5)
        name_col = _find_col(head, ["cellline"], "File 1 的 Cell line 列")
        ensg_col = _find_col(head, ["gene"], "File 1 的 Gene(ENSG) 列")
        sym_col = _find_col(head, ["genename"], "File 1 的 Gene name 列")
        val_col = _find_col(head, [self.value_col], f"File 1 的 {self.value_col} 列")

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

        # 两跳: name -> CVCL -> ACH
        cvcl = expr[name_col].astype(str).str.strip().map(name_to_cvcl)
        ach = resolver.resolve_cvcl_series(cvcl)
        # 兜底: File 11 无 CVCL 的条目, 用名字直查 Cellosaurus 别名
        need_fb = ach.isna()
        if need_fb.any():
            fb = resolver.resolve_name_series(expr.loc[need_fb, name_col])
            ach = ach.copy()
            ach[need_fb] = fb
            if fb.notna().sum():
                print(f"  [hpa_rna] 兜底: {fb.notna().sum():,} 行经名字直查补回")

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
    """GEO 表达 (File 3, **宽表**: 行=基因/ENSG, 列=每个 GSM) + File 10 (样本元数据)。

    v3.1 针对真实数据的两处修正
    ----------------------------
    1. File 3 的基因列装的是 **ENSG**, 不是 gene symbol
       (`ENSG00000000003` 而不是 `TSPAN6`)。所以按 symbol 过滤前必须先把
       目标基因翻译成 ENSG, 拿到数据后再翻回 symbol。
    2. File 10 自带 **`Cellosaurus_ID` 列** (你导师已经用 Cellosaurus 的
       GEO 交叉引用预先对齐好了, 见 `Matching_Type` = 'Cello GEO GSM')。
       所以走 CVCL 精确查表, 与 HPA 同级可靠; 名字匹配仅作兜底。

    ID 路径 (精确):
        GSM --File 10 的 Cellosaurus_ID--> CVCL_xxxx --RRID/xref--> ACH-xxxxxx

    多个 GSM 常对应同一细胞系 (技术重复/不同 GSE), 合并时在 build_gene_table
    里按 DepMap_ID 取均值。
    """
    layer, source = "rna", "geo_rna"

    def __init__(self, expr_path, info_path,
                 version: str = "GEO-import",
                 log_transform: bool = True,
                 chunksize: int = 5000):
        """
        log_transform : File 3 的值是线性微阵列强度 (33.6 / 553.2 / 2182.3),
            而 DepMap 和 HPA 都是 log 尺度。默认做 log2(x+1) 让三个 RNA 源
            落在可比尺度上, `value_unit` 会记录为 'log2(intensity+1)' 以保持可追溯。
            设为 False 则保留原始线性值。
        """
        self.expr_path = Path(expr_path)
        self.info_path = Path(info_path)
        self.version = version
        self.log_transform = log_transform
        self.chunksize = chunksize
        self._gsm_to_ach: Optional[pd.Series] = None

    def _load_gsm_map(self, resolver, verbose=True) -> pd.Series:
        """从 File 10 建 GSM -> ACH。优先用 Cellosaurus_ID 列, 名字仅兜底。"""
        if self._gsm_to_ach is not None:
            return self._gsm_to_ach

        info = pd.read_csv(self.info_path, sep="\t", dtype="string")
        gsm_col = _find_col(info, ["geoaccession", "gsm", "sampleid", "accession"],
                            "File 10 的 GSM 列")
        cvcl_col = _find_col_optional(info, ["cellosaurusid", "cellosaurus", "cvcl"])
        cell_col = _find_col_optional(info, ["cellline", "sourcename", "title"])

        if verbose:
            print(f"  [geo_rna] File 10 列识别: GSM={gsm_col!r}, "
                  f"CVCL={cvcl_col!r}, name={cell_col!r}")

        gsm = info[gsm_col].astype(str).str.strip()
        ach = pd.Series([None] * len(info), index=info.index, dtype="object")

        # --- 主路径: Cellosaurus_ID 直查 (精确) ---
        n_cvcl = 0
        if cvcl_col is not None:
            cvcl = (info[cvcl_col].astype("string").str.strip()
                    .str.extract(r"(CVCL_[A-Za-z0-9]+)", expand=False))
            ach = resolver.resolve_cvcl_series(cvcl)
            n_cvcl = ach.notna().sum()
            if verbose:
                n_has_cvcl = cvcl.notna().sum()
                print(f"  [geo_rna] 主路径 CVCL: {n_has_cvcl:,} 个样本有 CVCL, "
                      f"其中 {n_cvcl:,} 个 CVCL 能对到 DepMap")

        # --- 兜底: 用细胞系名直查 Cellosaurus 别名 ---
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
                    print(f"  [geo_rna] 兜底路径: 名字直查再补回 "
                          f"{fb.notna().sum():,} 个样本")

        m = pd.Series(ach.values, index=gsm.values)
        m = m[m.notna()]
        m = m[~m.index.duplicated(keep="first")]
        self._gsm_to_ach = m
        if verbose:
            print(f"  [geo_rna] File 10 合计: {len(info):,} 样本 -> "
                  f"{len(m):,} 个 GSM 解析到 ACH ({len(m)/max(len(info),1):.1%}), "
                  f"覆盖 {m.nunique():,} 个不同细胞系")
        return self._gsm_to_ach

    def load_long(self, resolver, gene_resolver, genes=None):
        gsm_map = self._load_gsm_map(resolver)
        if not len(gsm_map):
            print("  [geo_rna] ⚠️  没有任何 GSM 能解析到 ACH, 返回空表")
            return self._pack(pd.DataFrame(columns=UNIFIED_COLUMNS))

        head = pd.read_csv(self.expr_path, sep="\t", nrows=5)
        gene_col = head.columns[0]        # File 3 第一列就是基因列

        # --- 判断基因列装的是 ENSG 还是 symbol ---
        sample_vals = head[gene_col].astype(str)
        is_ensg = sample_vals.str.startswith("ENSG").mean() > 0.5
        if is_ensg:
            print(f"  [geo_rna] 基因列 {gene_col!r} 装的是 ENSG, 将翻译为 symbol")
        else:
            print(f"  [geo_rna] 基因列 {gene_col!r} 装的是 gene symbol")

        # --- 把目标基因翻译到文件所用的 ID 空间 ---
        wanted = None
        if genes is not None:
            if is_ensg:
                wanted = {e for e in (gene_resolver.to_ensembl(g) for g in genes)
                          if e}
                missing = {g for g in genes if not gene_resolver.to_ensembl(g)}
                if missing:
                    print(f"  [geo_rna] ⚠️  这些基因没有 ENSG 映射, GEO 里查不到: "
                          f"{sorted(missing)}")
                if not wanted:
                    return self._pack(pd.DataFrame(columns=UNIFIED_COLUMNS))
            else:
                wanted = set(genes)

        # --- 只读 基因列 + 可解析的 GSM 列 (3267 列全读会很慢) ---
        keep_gsm = [c for c in head.columns[1:] if c in gsm_map.index]
        if not keep_gsm:
            print("  [geo_rna] ⚠️  File 3 的 GSM 列没有一个能在 File 10 里解析到 ACH")
            return self._pack(pd.DataFrame(columns=UNIFIED_COLUMNS))
        print(f"  [geo_rna] File 3: {len(head.columns)-1:,} 个 GSM 列, "
              f"其中 {len(keep_gsm):,} 个可解析 -> 只读这些")

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
            print("  [geo_rna] 目标基因在 File 3 中没有对应行")
            return self._pack(pd.DataFrame(columns=UNIFIED_COLUMNS))

        long = expr.melt(id_vars=gene_col, var_name="GSM", value_name="value_raw")
        long["DepMap_ID"] = long["GSM"].map(gsm_map)
        long = long[long["value_raw"].notna()]

        # --- 基因 ID 归位 ---
        gene_key = long[gene_col].astype(str).str.split(".").str[0]
        if is_ensg:
            long["ensembl_id"] = gene_key
            long["gene_symbol"] = gene_key.map(gene_resolver.to_symbol)
            n_unmapped = long["gene_symbol"].isna().sum()
            if n_unmapped:
                print(f"  [geo_rna] {n_unmapped:,} 行的 ENSG 无法翻译成 symbol, 跳过")
        else:
            long["gene_symbol"] = gene_key
            long["ensembl_id"] = gene_key.map(gene_resolver.to_ensembl)

        # --- 尺度对齐 ---
        if self.log_transform:
            v = pd.to_numeric(long["value_raw"], errors="coerce")
            long["value_raw"] = np.log2(v.clip(lower=0) + 1)
            long["value_unit"] = "log2(intensity+1)"
        else:
            long["value_unit"] = "linear_intensity"

        return self._pack(long)


# ---------------------------------------------------------------------------
# 6. 合并引擎
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
        """构建统一长表。

        genes : 强烈建议传入 (例如 {'EGFR','KRAS'}); 不传会尝试加载全部基因,
                在真实数据上可能需要数十 GB 内存。
        """
        key = frozenset(genes) if genes is not None else None
        if key in self._cache:
            return self._cache[key]
        if key is None:
            print("⚠️  未指定 genes —— 将尝试加载所有基因, 可能耗尽内存。\n"
                  "    建议: build_long_table(genes={'EGFR', 'KRAS', ...})")
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
        """单基因宽表: 一行一细胞系, 直接喂评分模型。"""
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


# ---------------------------------------------------------------------------
# 7. 缺失值 -> 权重再分配
# ---------------------------------------------------------------------------

def redistribute_weights(mask: dict[str, bool],
                         base_weights: dict[str, float]) -> dict[str, float]:
    """缺失层的权重按比例摊到有数据的层; 不清零、不删样本。"""
    available = {k: w for k, w in base_weights.items() if mask.get(k, False)}
    total = sum(available.values())
    if total == 0:
        return {k: 0.0 for k in base_weights}
    return {k: (available.get(k, 0.0) / total) for k in base_weights}


# ---------------------------------------------------------------------------
# 8. 自检入口: python data_merger.py <data_dir>
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    print("=" * 66)
    print(f"data_merger.py 自检   (__version__ = {__version__})")
    print("=" * 66)

    # 自动找 data 目录: 命令行参数 > 当前目录/data > 上级目录/data
    if len(sys.argv) > 1:
        data_dir = Path(sys.argv[1])
    else:
        here = Path(__file__).resolve().parent
        candidates = [Path.cwd() / "data", here / "data", here.parent / "data"]
        data_dir = next((c for c in candidates if c.exists()), candidates[0])

    print(f"data_dir = {data_dir.resolve()}")
    if not data_dir.exists():
        print("✗ 目录不存在。用法: python src/data_merger.py <你的 data 目录>")
        sys.exit(1)

    cellosaurus_path = data_dir / "nomenclature" / "7_cellosaurus.csv"
    print(f"\n[1/3] 读取 Cellosaurus: {cellosaurus_path.name}")
    cs = pd.read_csv(cellosaurus_path)
    print(f"      {cs.shape[0]:,} 行 x {cs.shape[1]} 列")

    print("\n[2/3] 诊断 Cross-references 是否含 DepMap 引用")
    diagnose_cellosaurus(cs)

    print("\n[3/3] 构建 CellLineIDResolver")
    try:
        from data_loader import CellLineDataLoader
        loader = CellLineDataLoader(data_dir)
        si = loader.sample_info
        print(f"      sample_info: {si.shape}")
        print(f"      sample_info 列: {list(si.columns)}")
    except Exception as e:
        print(f"      (跳过 data_loader: {e})")
        si = None

    res = CellLineIDResolver(cs, si)
    print(f"\n最终字典规模: {res.stats()}")

    print("\n抽查:")
    for probe in ["CVCL_0023", "A549", "A-549", "HeLa", "MCF7"]:
        if probe.startswith("CVCL_"):
            print(f"  CVCL {probe:12} -> {res.resolve_cvcl_to_ach(probe)}")
        else:
            print(f"  name {probe:12} -> {res.resolve_name(probe)}")
    print("\n✓ 自检完成")
