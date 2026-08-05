"""
00_diagnose.py —— 先跑这个, 再跑 notebook
==========================================
用法:
    cd D:\\User1\\UoB_Coursework\\FinalProject
    python src\\00_diagnose.py data

作用: 把接入 HPA/GEO 之前必须确认的几件事一次性打印出来, 免得在 notebook 里
一步一步撞 KeyError。它 **只读文件、不写任何东西**。
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).parent))


def sep(title):
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def main(data_dir: Path):
    print(f"data_dir = {data_dir.resolve()}")
    if not data_dir.exists():
        print("✗ 目录不存在"); return

    # ---------------- 1. DepMap sample_info ----------------
    sep("1. DepMap sample_info (File 9) —— 找有没有自带 CVCL/RRID 列")
    si_path = data_dir / "nomenclature" / "9_DepMap_sample_info.csv"
    si = pd.read_csv(si_path)
    print(f"形状: {si.shape}")
    print(f"全部 {len(si.columns)} 个列名:")
    for c in si.columns:
        print(f"    - {c}")
    cvcl_like = [c for c in si.columns
                 if any(k in c.lower() for k in ("rrid", "cellosaurus", "cvcl"))]
    if cvcl_like:
        print(f"\n✓ 找到疑似 CVCL 列: {cvcl_like}")
        for c in cvcl_like:
            print(f"    {c} 样例: {si[c].dropna().head(3).tolist()}")
        print("  → 这是最权威的 CVCL<->ACH 桥, resolver 会自动优先使用它")
    else:
        print("\n○ sample_info 没有 CVCL 列, 将只依赖 Cellosaurus Cross-references")

    # ---------------- 2. Cellosaurus ----------------
    sep("2. Cellosaurus (File 7) —— 确认 Cross-references 里有没有 DepMap")
    cs = pd.read_csv(data_dir / "nomenclature" / "7_cellosaurus.csv")
    print(f"形状: {cs.shape}")
    print(f"列名: {list(cs.columns)}")
    from data_merger import diagnose_cellosaurus
    diagnose_cellosaurus(cs)

    # ---------------- 3. HPA File 1 ----------------
    sep("3. HPA File 1 (1_4_hpa_rna_celline.tsv)")
    p1 = data_dir / "gene expression" / "1_4_hpa_rna_celline.tsv"
    h1 = pd.read_csv(p1, sep="\t", nrows=5)
    print(f"列名: {list(h1.columns)}")
    print(h1.to_string(index=False))
    # 文件有多大 / 多少行 (只数行, 不读进内存)
    with open(p1, "r", encoding="utf-8", errors="ignore") as f:
        n_lines = sum(1 for _ in f)
    print(f"\n总行数: {n_lines:,}  (长表: 每个基因 x 每个细胞系一行)")
    n_cells = pd.read_csv(p1, sep="\t", usecols=["Cell line"])["Cell line"].nunique()
    print(f"唯一细胞系数: {n_cells}")
    print(f"推算基因数: ~{n_lines // max(n_cells,1):,}")

    # ---------------- 4. HPA File 11 ----------------
    sep("4. HPA File 11 (11_hpa_rna_celline_description.tsv)")
    p11 = data_dir / "nomenclature" / "11_hpa_rna_celline_description.tsv"
    h11 = pd.read_csv(p11, sep="\t", dtype="string")
    print(f"形状: {h11.shape}")
    print(f"列名: {list(h11.columns)}")
    cvcl_col = [c for c in h11.columns if "cellosaurus" in c.lower()]
    if cvcl_col:
        c = cvcl_col[0]
        n_ok = h11[c].fillna("").str.startswith("CVCL_").sum()
        print(f"\n{c}: {n_ok} / {len(h11)} 有 CVCL_id "
              f"({len(h11)-n_ok} 个是 Uncategorized 等无 CVCL 条目)")

    # ---------------- 5. GEO File 3 ----------------
    sep("5. GEO File 3 (3_GEOexpression.txt) —— 宽表")
    p3 = data_dir / "gene expression" / "3_GEOexpression.txt"
    g3 = pd.read_csv(p3, sep="\t", nrows=3)
    print(f"前 6 列: {list(g3.columns[:6])}")
    print(f"总列数: {len(g3.columns)}  (1 个基因列 + {len(g3.columns)-1} 个 GSM 列)")
    print(g3.iloc[:, :5].to_string(index=False))

    # ---------------- 6. GEO File 10 ----------------
    sep("6. GEO File 10 (10_GEOInfo.txt) —— 这是 GEO 能不能接通的关键")
    p10 = data_dir / "nomenclature" / "10_GEOInfo.txt"
    g10 = pd.read_csv(p10, sep="\t", dtype="string", nrows=10)
    print(f"列名: {list(g10.columns)}")
    print(g10.to_string(index=False))
    print("\n→ 请确认: 哪一列是 GSM? 哪一列写着细胞系名?")
    print("  resolver 会尝试用关键词 gsm/sample/accession 和 "
          "cellline/cell/sourcename/title/characteristics 自动找。")

    sep("诊断完成")
    print("下一步: 打开 03_Data_Merge.ipynb 按顺序跑。")


if __name__ == "__main__":
    main(Path(sys.argv[1] if len(sys.argv) > 1 else "data"))
