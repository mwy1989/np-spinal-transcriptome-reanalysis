"""
v11 阶段 2 / 任务 34b：PROGENy 通路活性重跑（阈值无关，作 GSVA 的正交验证）

关键修正（相对 v10）：
  1. v10 的 progeny_dorothea_step1.py 用「列位置」分配样本分组
     （for i,g: for j in range(3): idx = i*3+j）——本次实测恰好正确，
     但属运气；v11 改为从 metadata 显式取分组
  2. 统一口径为 log2(FPKM + 0.1)，与 limma / GSVA 主分析一致
  3. 输出到 output/v11/，只做打分；统计统一在 R 中做 limma
"""
import os
import numpy as np
import pandas as pd
import decoupler as dc

import os as _os
BASE = _os.environ.get("SCS_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
OUT = os.path.join(BASE, "output", "v11")
os.makedirs(OUT, exist_ok=True)

# ---------- 1. 表达矩阵（显式对齐 metadata）----------
expr_lin = pd.read_csv(os.path.join(BASE, "output", "GSE175760_expression_matrix.csv"),
                       index_col=0)
meta = pd.read_csv(os.path.join(BASE, "output", "GSE175760_metadata.csv"))
meta = meta.set_index("sample_id")

# 显式按 metadata 顺序取列（不再按位置猜）
common = [s for s in meta.index if s in expr_lin.columns]
assert len(common) == len(meta), "样本对齐失败"
expr_lin = expr_lin[common].dropna()          # 与 limma 主分析一致：剔除含缺失的行
print(f"剔除含 NaN 的行后: {expr_lin.shape[0]} 基因")
mat = np.log2(expr_lin + 0.1)
print(f"表达矩阵: {mat.shape[0]} 基因 x {mat.shape[1]} 样本")
print("分组（来自 metadata）:", meta.loc[common, "group"].value_counts().to_dict())
print("口径: log2(FPKM + 0.1)")

# ---------- 2. PROGENy 网络，人 -> 大鼠 ----------
net = dc.op.progeny(organism="human")
print(f"\nPROGENy 网络（人）: {net.shape[0]} 条 边, {net['source'].nunique()} 条通路")
rat_upper = {g.upper(): g for g in mat.index}
net["target"] = net["target"].str.upper().map(rat_upper)
n_before = len(net)
net = net.dropna(subset=["target"])
print(f"映射到大鼠后: {len(net)} / {n_before} 条边（覆盖率 {len(net)/n_before:.1%}）")
cov = net.groupby("source")["target"].nunique().sort_values()
print("\n各通路可用基因数:")
print(cov.to_string())

# ---------- 3. PROGENy 通路活性（multivariate linear model）----------
print("\n运行 PROGENy MLM...")
# decoupler 2.x 期望 samples x features（与 1.x 相反）
mat_t = mat.T.copy()
mat_t.index.name = "sample_id"
scores, pvals = dc.mt.mlm(mat_t, net, tmin=5)
assert scores.shape == (18, 14), f"意外形状: {scores.shape}"   # samples x pathways
scores.index.name = "sample_id"
print("通路活性矩阵:", scores.shape)

# 附分组
scores_out = scores.copy()
scores_out.insert(0, "group", meta.loc[scores.index, "group"])
scores_out.to_csv(os.path.join(OUT, "progeny_scores.csv"))
pd.DataFrame(pvals, index=scores.index).to_csv(os.path.join(OUT, "progeny_pvals_raw.csv"))
print("输出: progeny_scores.csv, progeny_pvals_raw.csv")

# ---------- 4. 快速查看各时间点均值 ----------
GROUPS = ["Sham", "CCI_0.5d", "CCI_1d", "CCI_3d", "CCI_7d", "CCI_14d"]
means = scores.groupby(meta.loc[scores.index, "group"]).mean().reindex(GROUPS)
delta = means.subtract(means.loc["Sham"], axis=1)
print("\n===== 各时间点 vs Sham 的通路活性变化（z 分数差）=====")
print(delta.round(3).T.to_string())
delta.T.to_csv(os.path.join(OUT, "progeny_delta_quicklook.csv"))

# ---------- 5. 与 Hallmark GSVA 的方向一致性（正/负号）----------
gsva_h = pd.read_csv(os.path.join(OUT, "gsva_hallmark_delta.csv"))
KEY = {
    "JAK-STAT": "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "NFkB": "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "TNFa": "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "Hypoxia": "HALLMARK_HYPOXIA",
    "p53": "HALLMARK_P53_PATHWAY",
    "PI3K": "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
    "TGFb": "HALLMARK_TGF_BETA_SIGNALING",
    "MAPK": "HALLMARK_KRAS_SIGNALING_UP",
}
rows = []
gsva_h = gsva_h.set_index("pathway")
for pw, hm in KEY.items():
    if pw not in delta.columns or hm not in gsva_h.index:
        continue
    a = delta.loc["CCI_7d", pw]
    b = gsva_h.loc[hm, "7d"]
    rows.append({"progeny": pw, "PROGENy_delta_7d": round(a, 3),
                 "hallmark": hm, "GSVA_delta_7d": round(b, 3),
                 "sign_agree": (np.sign(a) == np.sign(b))})
if rows:
    cons = pd.DataFrame(rows)
    print("\n===== PROGENy vs GSVA(Hallmark) 7d 方向一致性 =====")
    print(cons.to_string(index=False))
    print(f"符号一致: {cons.sign_agree.sum()}/{len(cons)}")
    cons.to_csv(os.path.join(OUT, "progeny_gsva_consistency.csv"), index=False)

print("\n完成。")
