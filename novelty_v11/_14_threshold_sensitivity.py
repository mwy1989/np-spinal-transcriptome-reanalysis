"""
方案 B：阈值 / 效应量敏感性分析

核心目的
  1. 量化"同一数据集在不同分析口径下能得到多少个 DEG"—— 这个谱从未有人系统报告
  2. 精确重现 He 2025 的口径（p<=0.05 & FC>=1.5 & within-group FPKM mean>=0.5）
     看能否复现其报告的 7d=529 / 14d=352 —— 若接近，说明差异纯粹来自阈值而非数据
  3. 检验本稿 29 个候选基因是否"阈值脆弱"（换口径就掉出去 = 危险信号）

输入： output/v11/limma_all_contrasts.csv（13,996 基因 x 5 时点，含 P_raw / P_BH / log2FC / AveExpr）
      output/v11/candidates_primary.csv（29 候选）
输出： B_threshold_spectrum.csv / B_jaccard_matrix.csv / B_candidate_robustness.csv
"""
import os
import numpy as np
import pandas as pd

# --- 复现包路径（可移植化改写；原始脚本使用本机绝对路径）---
_PKG = os.environ.get("SCS_ROOT") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(_PKG, "results")               # 原始: F:\scs research\output\v11
OUT = os.path.join(_PKG, "results")   # 原始: F:\scs research\outputs\novelty_v11

d = pd.read_csv(os.path.join(BASE, "limma_all_contrasts.csv"))
cand = pd.read_csv(os.path.join(BASE, "candidates_primary.csv"))["gene"].tolist()
print(f"limma 结果 {d.shape}，候选基因 {len(cand)} 个")

TP_ORDER = ["0.5d", "1d", "3d", "7d", "14d"]
# 本稿口径：|log2FC| > 0.58（手稿 Methods 逐字："> 0.58 (a 1.5-fold change)"）
FC15 = 0.58
# He 2025 口径：原文写 "a fold change of at least 1.5" => |log2FC| >= log2(1.5)
FC15_HE = np.log2(1.5)   # 0.58496
LOG2_2   = np.log2(2.0)   # 1.0
LOG2_3   = np.log2(3.0)   # 1.585

# AveExpr 是 log2(FPKM+0.1) 的均值 => FPKM>=0.5 等价于 AveExpr >= log2(0.6)
FPKM05 = np.log2(0.6)     # -0.737

CRITERIA = [
    ("① BH<0.05 & |FC|>0.58（本稿严格集）",   lambda x: (x.P_BH < 0.05) & (x.log2FC.abs() > FC15)),
    ("② P<0.05 & |FC|>0.58（本稿探索集）",    lambda x: (x.P_raw < 0.05) & (x.log2FC.abs() > FC15)),
    ("③ BH<0.05（不限 FC）",                 lambda x: x.P_BH < 0.05),
    ("④ P<0.05（不限 FC）",                  lambda x: x.P_raw < 0.05),
    ("⑤ P<0.01 & |FC|>0.58",                lambda x: (x.P_raw < 0.01) & (x.log2FC.abs() > FC15)),
    ("⑥ ★He 2025 口径: P<=0.05 & FC>=1.5 & FPKM>=0.5",
                                            lambda x: (x.P_raw <= 0.05) & (x.log2FC.abs() >= FC15_HE) & (x.AveExpr >= FPKM05)),
    ("⑦ BH<0.05 & |FC|>=2",                 lambda x: (x.P_BH < 0.05) & (x.log2FC.abs() >= LOG2_2)),
    ("⑧ P<0.05 & |FC|>=2",                  lambda x: (x.P_raw < 0.05) & (x.log2FC.abs() >= LOG2_2)),
    ("⑨ P<0.01 & |FC|>=2（最严）",            lambda x: (x.P_raw < 0.01) & (x.log2FC.abs() >= LOG2_2)),
    ("⑩ BH<0.05 & |FC|>0.58 & FPKM>=0.5",   lambda x: (x.P_BH < 0.05) & (x.log2FC.abs() > FC15) & (x.AveExpr >= FPKM05)),
]

# ---------- 1. DEG 数量谱 ----------
rows = []
sets = {}   # criterion -> set(genes unioned across timepoints)
for name, fn in CRITERIA:
    mask = fn(d)
    row = {"criterion": name, "total_rows": int(mask.sum())}
    genes_union = set()
    for tp in TP_ORDER:
        m = mask & (d.tp == tp)
        row[tp] = int(m.sum())
        genes_union |= set(d.loc[m, "gene"])
    row["union"] = len(genes_union)
    sets[name] = genes_union
    rows.append(row)

spec = pd.DataFrame(rows)
spec.to_csv(os.path.join(OUT, "B_threshold_spectrum.csv"), index=False)

print("\n" + "=" * 100)
print("===== B1. DEG 数量谱（每个时点分别计数 + 跨时点去重并集）=====")
print(spec.to_string(index=False))

# ---------- 2. He 2025 对照 ----------
he = spec[spec.criterion.str.contains("He 2025")]
if len(he):
    h = he.iloc[0]
    print("\n" + "=" * 100)
    print("===== B2. 与 He 2025 的对照（同一数据集 GSE175760）=====")
    print(f"{'时点':<8}{'He 2025 报告':>14}{'本稿同口径':>12}{'本稿 BH<0.05':>14}{'本稿名义P<0.05':>16}")
    he_rep = {"0.5d": None, "1d": None, "3d": None, "7d": 529, "14d": 352}
    for tp in TP_ORDER:
        bh = spec.loc[spec.criterion.str.contains("① "), tp].iloc[0]
        nm = spec.loc[spec.criterion.str.contains("② "), tp].iloc[0]
        rep = he_rep[tp] if he_rep[tp] else "n/a"
        print(f"{tp:<8}{str(rep):>14}{h[tp]:>12}{bh:>14}{nm:>16}")
    print(f"\n并集：He 7d+14d = 529+352 = {529+352}；本稿同口径 7d+14d = {h['7d']+h['14d']}")

# ---------- 3. Jaccard 稳定性矩阵 ----------
names = [n for n, _ in CRITERIA]
J = pd.DataFrame(index=names, columns=names, dtype=float)
for a in names:
    for b in names:
        A, B = sets[a], sets[b]
        J.loc[a, b] = len(A & B) / len(A | B) if (A | B) else np.nan
J.to_csv(os.path.join(OUT, "B_jaccard_matrix.csv"))

print("\n" + "=" * 100)
print("===== B3. 口径间 Jaccard 稳定性（基因集两两重叠）=====")
print(J.round(2).iloc[:, :5].to_string())

# ---------- 4. 候选基因稳健性 ----------
print("\n" + "=" * 100)
print("===== B4. 29 个候选基因在各口径下是否仍被选中 =====")
rows = []
for g in cand:
    sub = d[d.gene == g]
    if sub.empty:
        rows.append({"gene": g, "in_any": 0, "note": "未在 limma 结果中"})
        continue
    rec = {"gene": g}
    for name, fn in CRITERIA:
        m = fn(sub)
        # 至少在 1 个时点满足
        rec[name.split("（")[0].strip()] = int(m.any())
    # 该基因最宽松口径（名义 P<0.05）下的最小 P 与最大 FC
    rec["min_P_raw"] = sub.P_raw.min()
    rec["min_P_BH"] = sub.P_BH.min()
    rec["max_abs_log2FC"] = sub.log2FC.abs().max()
    rec["min_AveExpr"] = sub.AveExpr.min()
    rec["n_criteria_pass"] = sum(rec[k] for k in rec if k.startswith(("①", "②", "③", "④", "⑤", "⑥", "⑦", "⑧", "⑨", "⑩")))
    rows.append(rec)

rob = pd.DataFrame(rows)
rob.to_csv(os.path.join(OUT, "B_candidate_robustness.csv"), index=False)

crit_keys = [c for c in rob.columns if c[0] in "①②③④⑤⑥⑦⑧⑨⑩"]
print(f"\n{'基因':<18}{'通过口径数':>10}  " + " ".join(f"{k:>5}" for k in crit_keys))
for _, r in rob.iterrows():
    if "n_criteria_pass" not in r or pd.isna(r.get("n_criteria_pass")):
        print(f"{r['gene']:<18}{'—':>10}")
        continue
    print(f"{r['gene']:<18}{int(r['n_criteria_pass']):>10}  " + " ".join(f"{int(r[k]):>5}" for k in crit_keys))

n_robust = int((rob["n_criteria_pass"] == len(CRITERIA)).sum()) if "n_criteria_pass" in rob else 0
n_bh_only = int(rob[crit_keys[0]].sum())
print(f"\n【关键】在全部 {len(CRITERIA)} 种口径下都被选中的候选基因：{n_robust} / {len(cand)}")
print(f"【关键】在严格集 BH<0.05 口径下被选中的候选基因：{n_bh_only} / {len(cand)}")
print(f"\n输出：B_threshold_spectrum.csv / B_jaccard_matrix.csv / B_candidate_robustness.csv")
