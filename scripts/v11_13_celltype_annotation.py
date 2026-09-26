"""
v11 阶段 2 / 任务 36：候选基因细胞类型表达注释（纯描述性）

相对 v10 的关键变化：
  1. v10 强制把 10 种细胞类型二分（glial vs neuronal），并用 Oligodendrocyte
     作"神经元 proxy"——但 GSE189070 原始注释中**根本不含神经元**
     -> v11 取消双分区标签，只报告各真实细胞类型中的表达与检测率
  2. 不写"空间隔离""胶质抵抗"等超出单细胞表达数据支持范围的表述
  3. 如实报告：该参考不含神经元；部分候选基因在该参考中缺失
"""
import os
import numpy as np
import pandas as pd
import anndata as ad

import os as _os
BASE = _os.environ.get("SCS_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
OUT = os.path.join(BASE, "output", "v11")
os.makedirs(OUT, exist_ok=True)

# ---------- 1. 读单细胞参考 ----------
print("读取 GSE189070 注释数据...")
adata = ad.read_h5ad(os.path.join(BASE, "output", "scRNA_GSE189070", "GSE189070_annotated.h5ad"))
print(f"  {adata.shape[0]} 细胞 x {adata.shape[1]} 基因")
cts = adata.obs["celltype"].astype(str)
print("  细胞类型:", cts.value_counts().to_dict())
print("  !!! 该参考不含神经元细胞类型（原始注释 10 类均为胶质/免疫/血管/室管膜）")

# ---------- 2. 候选基因 ----------
cand = pd.read_csv(os.path.join(OUT, "candidates_primary.csv"))
genes = cand["gene"].tolist()
var_upper = {v.upper(): v for v in adata.var_names}
present, missing = [], []
for g in genes:
    if g.upper() in var_upper:
        present.append((g, var_upper[g.upper()]))
    else:
        missing.append(g)
print(f"\n候选基因 {len(genes)} 个；在该参考中可检测 {len(present)} 个；缺失 {len(missing)} 个: {missing}")

# ---------- 3. 每细胞类型的平均表达 + 检测率 ----------
CT_ORDER = ["Microglia", "Macrophage", "Neutrophil", "Astrocyte", "Oligodendrocyte",
            "Ependymal", "Endothelial", "Pericyte", "B_cell", "T_cell"]

mean_expr = pd.DataFrame(index=[g for g, _ in present], columns=CT_ORDER, dtype=float)
det_rate = mean_expr.copy()

for ct in CT_ORDER:
    mask = (cts == ct).values
    if mask.sum() == 0:
        continue
    sub = adata.X[mask, :]
    for g, vname in present:
        j = adata.var_names.get_loc(vname)
        col = np.asarray(sub[:, j].todense()).ravel() if hasattr(sub[:, j], "todense") \
              else np.asarray(sub[:, j]).ravel()
        mean_expr.loc[g, ct] = col.mean()
        det_rate.loc[g, ct] = (col > 0).mean()

mean_expr.round(4).to_csv(os.path.join(OUT, "candidate_celltype_mean_expression.csv"))
det_rate.round(4).to_csv(os.path.join(OUT, "candidate_celltype_detection_rate.csv"))

# ---------- 4. 优势细胞类型 + tau 特异性指数 ----------
def tau(v):
    v = np.asarray(v, dtype=float)
    mx = v.max()
    if mx <= 0 or len(v) < 2:
        return np.nan
    return float(np.sum(1 - v / mx) / (len(v) - 1))

summary = []
for g in mean_expr.index:
    row = mean_expr.loc[g].astype(float)
    dom = row.idxmax()
    dr = float(det_rate.loc[g].max())
    # 可评估性：最高检测率 > 10% 且该类型平均表达 > 0.1
    assessable = (dr > 0.10) and (row.max() > 0.1)
    summary.append({
        "gene": g,
        "dominant_celltype": dom if assessable else "not_assessable",
        "max_mean_expr": round(row.max(), 4),
        "second_celltype": row.drop(dom).idxmax() if len(row) > 1 else None,
        "tau": round(tau(row.values), 3),
        "max_detection_rate": round(dr, 3),
        "n_celltypes_detected": int((row > 0.05).sum()),
        "assessable": assessable,
    })
summ = pd.DataFrame(summary).sort_values("gene")
summ.to_csv(os.path.join(OUT, "candidate_dominant_celltype.csv"), index=False)

print("\n===== 候选基因细胞类型注释（纯描述性）=====")
print(summ.to_string(index=False))

n_ass = int(summ["assessable"].sum())
print(f"\n在该参考中可可靠评估: {n_ass} / {len(summ)}")
print(f"不可评估（检测率<10%）: {summ.loc[~summ['assessable'],'gene'].tolist()}")

# ---------- 5. 汇总：优势细胞类型分布（仅可评估者）----------
print("\n===== 优势细胞类型计数（仅可评估基因）=====")
print(summ.loc[summ["assessable"], "dominant_celltype"].value_counts().to_string())

print("\n===== 必须写入文章的事实说明 =====")
print(f"  1. 缺失基因（该参考无注释）: {missing}")
print(f"  2. 该参考不含神经元细胞类型 → 对神经元来源基因无法评估，不做任何推断")
print(f"  3. 检测率<10% 的 {int((~summ['assessable']).sum())} 个基因不做细胞类型归属"
      f"（低表达下 idxmax 属噪声，v10 将其强制归入某一'区室'是过度解读）")

# 供 Fig 5 使用的宽表（按 tau 排序，便于展示）
fig5 = mean_expr.loc[summ.sort_values("tau", ascending=False)["gene"]]
fig5.to_csv(os.path.join(OUT, "fig5_celltype_heatmap_matrix.csv"))
print("\n输出:")
for f in ["candidate_celltype_mean_expression.csv", "candidate_celltype_detection_rate.csv",
          "candidate_dominant_celltype.csv", "fig5_celltype_heatmap_matrix.csv"]:
    print("  ", f)
