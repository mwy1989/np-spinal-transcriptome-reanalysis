"""
Fig 2 — Temporal expression structure and pathway programme
A) 四个时间软聚类轨迹（fuzzy c-means on the exploratory set）
B) WGCNA 模块 × 时间聚类的交叉分布（两种独立框架的收敛）
C) GSVA Hallmark 通路时程热图（7 d 达峰）
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v11_figstyle import (apply_style, save, despine, panel_tag, V11,
                          TEXT, GRID, SPINE, NEUTRAL, TP_ORDER, TP_LABEL,
                          CLUSTER_COLOR, CLUSTER_NAME, MODULE_COLOR,
                          div_cmap, SOFT_BLUE)

apply_style()
print("Fig 2 ...")

prof = pd.read_csv(os.path.join(V11, "cluster_k4_profiles.csv"))
asg = pd.read_csv(os.path.join(V11, "cluster_k4_assignments.csv"))
gm = pd.read_csv(os.path.join(V11, "wgcna_degs_gene_modules.csv"))
long = pd.read_csv(os.path.join(V11, "gsva_hallmark_delta_long.csv"))
wide = pd.read_csv(os.path.join(V11, "gsva_hallmark_delta.csv"))

fig = plt.figure(figsize=(7.2, 8.0))
gs = GridSpec(2, 5, figure=fig, height_ratios=[0.72, 1.28],
              hspace=0.40, wspace=0.55)

# ================= A. Cluster trajectories =================
axA = fig.add_subplot(gs[0, :3])
xs = np.arange(len(TP_ORDER))
for c in sorted(prof["cluster"].unique()):
    sub = prof[prof["cluster"] == c].set_index("tp").loc[TP_ORDER]
    n = int(sub["n_genes"].iloc[0])
    axA.plot(xs, sub["z_profile"].values, "-o", color=CLUSTER_COLOR[c],
             lw=1.6, ms=4.0, mec="white", mew=0.6, zorder=3,
             label=f"{CLUSTER_NAME[c]}  (n = {n})")
axA.axhline(0, color=SPINE, lw=0.6, ls=(0, (3, 3)), zorder=1)
axA.set_xticks(xs)
axA.set_xticklabels([TP_LABEL[t] for t in TP_ORDER])
axA.set_xlabel("Time after CCI")
axA.set_ylabel("Mean z-scored expression")
axA.set_title("Temporal expression programmes", loc="left")
axA.legend(loc="lower left", fontsize=5.9, handlelength=1.4,
           bbox_to_anchor=(-0.005, 0.0), labelspacing=0.35)
despine(axA)
axA.grid(True, axis="y", lw=0.5, alpha=0.5)
axA.set_axisbelow(True)
axA.set_ylim(-2.35, 1.80)
axA.set_yticks([-2, -1, 0, 1])
panel_tag(axA, "A")

# ================= B. module x cluster =================
axB = fig.add_subplot(gs[0, 3:])
m = asg.merge(gm[["gene", "module"]], on="gene", how="left")
m["module"] = m["module"].fillna("not_assigned")
mods = ["turquoise", "blue", "brown", "grey"]
ct = pd.crosstab(m["module"], m["cluster"]).reindex(mods).fillna(0)
pct = ct.div(ct.sum(axis=1).replace(0, np.nan), axis=0) * 100

ys = np.arange(len(mods))[::-1]
left = np.zeros(len(mods))
for c in sorted(prof["cluster"].unique()):
    v = pct[c].values if c in pct.columns else np.zeros(len(mods))
    axB.barh(ys, v, left=left, height=0.62, color=CLUSTER_COLOR[c],
             edgecolor="white", linewidth=0.5, zorder=3,
             label=f"C{c}" if False else None)
    for yi, (l, w) in enumerate(zip(left, v)):
        if w >= 12:
            axB.text(l + w / 2, ys[yi], f"{int(ct[c].values[yi])}", ha="center",
                     va="center", fontsize=6.0, color=TEXT)
    left = left + v

axB.set_yticks(ys)
axB.set_yticklabels([f"{mm}\n(n = {int(ct.loc[mm].sum())})" for mm in mods], fontsize=6.4)
axB.set_xlim(0, 100)
axB.set_xlabel("Genes in module (%)")
axB.set_title("WGCNA modules vs temporal clusters", loc="left", fontsize=8.4)
despine(axB, keep=("bottom",))
axB.grid(True, axis="x", lw=0.5, alpha=0.5)
axB.set_axisbelow(True)

handles = [plt.Rectangle((0, 0), 1, 1, fc=CLUSTER_COLOR[c]) for c in sorted(prof["cluster"].unique())]
axB.legend(handles, [f"C{c}" for c in sorted(prof["cluster"].unique())],
           ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.34),
           fontsize=6.2, handlelength=0.8, columnspacing=0.8)
panel_tag(axB, "B", dx=-0.24)

# ================= C. GSVA hallmark heatmap =================
axC = fig.add_subplot(gs[1, :])
sig = long[long["FDR"] < 0.05]["pathway"].unique()
H = wide.set_index("pathway")
sel = H.loc[[p for p in H.index if p in set(sig)]].copy()
sel = sel.sort_values("7d", ascending=False)
mat = sel[["0.5d", "1d", "3d", "7d", "14d"]].values

vmax = np.nanmax(np.abs(mat))
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
im = axC.imshow(mat, cmap=div_cmap(), norm=norm, aspect="auto")

# FDR 标注
fdr = long.pivot(index="pathway", columns="timepoint", values="FDR").reindex(sel.index)
COLS = ["0.5d", "1d", "3d", "7d", "14d"]
fdr = fdr[COLS]
star_log = {}
for i in range(mat.shape[0]):
    hits = []
    for j in range(mat.shape[1]):
        v = fdr.values[i, j]
        if not np.isfinite(v):
            continue
        if v < 0.001:
            axC.text(j, i, "***", ha="center", va="center", fontsize=5.6,
                     color="#FFFFFF" if abs(mat[i, j]) > vmax * 0.55 else TEXT)
        elif v < 0.01:
            axC.text(j, i, "**", ha="center", va="center", fontsize=5.6,
                     color="#FFFFFF" if abs(mat[i, j]) > vmax * 0.55 else TEXT)
        elif v < 0.05:
            axC.text(j, i, "*", ha="center", va="center", fontsize=6.4,
                     color="#FFFFFF" if abs(mat[i, j]) > vmax * 0.55 else TEXT)
        else:
            continue
        hits.append(COLS[j])
    if hits:
        star_log[sel.index[i]] = hits

axC.set_xticks(range(5))
axC.set_xticklabels(["0.5 d", "1 d", "3 d", "7 d", "14 d"])
axC.set_yticks(range(len(sel)))
axC.set_yticklabels(sel.index, fontsize=6.3)
axC.set_xlabel("Time after CCI")
axC.set_title("Hallmark pathway programmes (GSVA, FDR < 0.05 at >= 1 time point)", loc="left")
for s in axC.spines.values():
    s.set_visible(False)
axC.tick_params(length=0)
cb = fig.colorbar(im, ax=axC, fraction=0.022, pad=0.015)
cb.set_label("\u0394 GSVA score vs Sham", fontsize=6.6)
cb.ax.tick_params(labelsize=6.2)
cb.outline.set_linewidth(0.5)
cb.outline.set_edgecolor(SPINE)
panel_tag(axC, "C", dx=-0.075)

save(fig, "Figure_2_Temporal_Structure_and_Pathways")
print(f"  pathways plotted: {len(sel)}")
print("  --- 显著性星号位置核对 ---")
for k, v in star_log.items():
    print(f"    {k:42s} {v}")
print("done.")
