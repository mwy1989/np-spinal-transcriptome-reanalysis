"""
v11 阶段 2 / 任务 34d：候选基因功能富集（GO-BP / KEGG）+ WGCNA 模块富集

相对 v10 的变化：
  - 输入基因集改为 v11 新候选名单（29 个）与新 WGCNA 模块
  - 使用 Enrichr 官方库（GO_Biological_Process / KEGG），非手工基因集
  - 明确标注：Enrichr 无 Rat 库，使用 Mouse 库（同源基因符号通用）
"""
import os
import numpy as np
import pandas as pd
import gseapy as gp

import os as _os
BASE = _os.environ.get("SCS_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
OUT = os.path.join(BASE, "output", "v11")

cand = pd.read_csv(os.path.join(OUT, "candidates_primary.csv"))
genes = [g for g in cand["gene"] if not g.startswith("AABR")]   # 剔除未注释的大鼠预测基因
print(f"候选基因: {len(cand)} 个，可用于富集: {len(genes)} 个")
print("列表:", ", ".join(genes))

gm = pd.read_csv(os.path.join(OUT, "wgcna_degs_gene_modules.csv"))
mods = {m: gm["gene"][gm["module"] == m].tolist() for m in ["turquoise", "blue", "brown"]}
print("\n模块基因数:", {k: len(v) for k, v in mods.items()})

GENE_SETS = ["GO_Biological_Process_2023", "KEGG_2019_Mouse",
             "MSigDB_Hallmark_2020", "Reactome_2022"]

def run_enrich(gl, tag, top=20):
    print(f"\n===== 富集分析: {tag} (n={len(gl)}) =====")
    try:
        res = gp.enrichr(gene_list=gl, gene_sets=GENE_SETS,
                         organism="mouse", outdir=None, no_plot=True)
    except Exception as e:
        print("  失败:", type(e).__name__, e)
        return None
    df = res.results.copy()
    df = df[df["Adjusted P-value"] < 0.25].sort_values("Adjusted P-value")
    if df.empty:
        print("  无 FDR<0.25 的条目")
        return df
    df.to_csv(os.path.join(OUT, f"enrichment_{tag}.csv"), index=False)
    show = df.head(top)
    for _, r in show.iterrows():
        print(f"  [{r['Gene_set'][:16]:16s}] {r['Term'][:58]:58s} FDR={r['Adjusted P-value']:.2e}  {r['Overlap']}")
    print(f"  ... 共 {len(df)} 条 FDR<0.25")
    return df

run_enrich(genes, "candidates")
run_enrich(mods["turquoise"][:1500], "module_turquoise")
run_enrich(mods["blue"][:1500], "module_blue")

print("\n完成。")
