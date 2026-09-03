import sys, os
from pathlib import Path

# 让脚本能找到 src 里的模块
SRC = Path(__file__).resolve().parent          # 本脚本在 src/ 里
sys.path.insert(0, str(SRC))
DATA_DIR = SRC.parent / 'data'                  # 假设 data/ 和 src/ 同级

import pandas as pd
import numpy as np
from data_loader import CellLineDataLoader
from data_merger import (CellLineIDResolver, GeneIDResolver, MultiOmicsMerger,
                         DepMapRNASource, HPARNASource, GEOWideSource,
                         MatrixProteinSource, MutationSource, FusionSource)

# ===== 1. 重建 loader / resolver / merger =====
print('初始化数据层...')
loader = CellLineDataLoader(DATA_DIR)
cell_res = CellLineIDResolver(
    pd.read_csv(DATA_DIR / 'nomenclature' / '7_cellosaurus.csv'),
    loader.sample_info)
gene_res = GeneIDResolver().build_from_hpa_file(
    DATA_DIR / 'gene expression' / '1_4_hpa_rna_celline.tsv')

merger = (MultiOmicsMerger(loader.sample_info['DepMap_ID'], cell_res, gene_res)
          .register(DepMapRNASource(loader))
          .register(HPARNASource(
              DATA_DIR / 'gene expression' / '1_4_hpa_rna_celline.tsv',
              DATA_DIR / 'nomenclature' / '11_hpa_rna_celline_description.tsv'))
          .register(GEOWideSource(
              DATA_DIR / 'gene expression' / '3_GEOexpression.txt',
              DATA_DIR / 'nomenclature' / '10_GEOInfo.txt'))
          .register(MatrixProteinSource(loader))
          .register(MutationSource(loader))
          .register(FusionSource(loader)))

# ===== 2. 选基因(排除22个,随机100个,固定种子) =====
EXCLUDE = {
    'AFP','ALB','AR','DES','ESR1','GFAP','MET','MITF','MSLN','PTPRC',
    'EGFR','KRAS','ERBB2','MYC','TP53','FGFR2','CD86','KLK4','GAPDH',
    'ASGR1','MUC1','CD3E'
}

rna_wide = loader._load_rna_matrix()
all_genes = [c.split(' (')[0] for c in rna_wide.columns]
candidates = [g for g in all_genes if g not in EXCLUDE]
print(f'候选基因数: {len(candidates)}')

np.random.seed(42)                              # 固定种子,保证可复现
selected = sorted(np.random.choice(candidates, size=100, replace=False))
assert not (set(selected) & EXCLUDE), '不应包含排除基因'
print(f'选中 100 个基因,已排除 {len(EXCLUDE)} 个')


# # planB:"四维完整"的基因,待定
#
# print(f'选中 {len(selected)} 个基因(已排除22个)')
# assert not (set(selected) & EXCLUDE), '不应包含排除基因'
# ===== 3. 一次性合并这100个基因(别循环) =====
print('合并中(会扫一遍 HPA,稍等)...')
long_100 = merger.build_long_table(genes=set(selected))

# ===== 4. 导出 =====
out_dir = DATA_DIR / 'merged'
out_dir.mkdir(parents=True, exist_ok=True)
long_100.to_parquet(out_dir / 'ablation_100genes.parquet')
pd.Series(selected, name='gene').to_csv(out_dir / 'ablation_100genes_list.csv', index=False)

print(f'\n完成:')
print(f'  数据: {out_dir / "ablation_100genes.parquet"}  ({len(long_100):,} 行)')
print(f'  基因清单: {out_dir / "ablation_100genes_list.csv"}')