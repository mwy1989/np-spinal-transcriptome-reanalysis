"""
Fig 3 — Candidate genes and cross-dataset concordance
A) 29 个候选基因 × 6 时间点表达热图（行 z-score），左侧标注 WGCNA 模块
B) 发现队列 vs GSE5296 效应量散点（7 d）
C) 各时间点配对的方向一致率与效应量相关
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
                          MODULE_COLOR, SOFT_RED, div_cmap)

apply_style()
print("Fig 3 ...")

BASE = "F:/scs research"
cand = pd.read_csv(os.path.join(V11, "candidates_primary.csv"))
expr = pd.read_csv(os.path.join(V11, "expr_log2_pseudo0.1.csv"), index_col=0)
meta = pd.read_csv(os.path.join(BASE, "output", "GSE175760_metadata.csv"))
asg = pd.read_csv(os.path.join(V11, "cluster_k4_assignments.csv"))
conc = pd.read_csv(os.path.join(V11, "gse5296_candidate_concordance.csv"))
det = pd.read_csv(os.path.join(V11, "gse5296_7d_candidate_detail.csv"))

# ---------- 行排序：模块 -> 时间簇 -> 效应量 ----------
MOD_ORDER = ["turquoise", "blue", "brown", "grey", "not_assigned"]
cand["mod_rank"] = cand["WGCNA_module"].map({m: i for i, m in enumerate(MOD_ORDER)}).fillna(9)
cand = cand.merge(asg[["gene", "cluster"]], on="gene", how="left")
cand["cluster"] = cand["cluster"].fillna(0)
cand = cand.sort_values(["mod_rank", "cluster", "max_abs_log2FC"],
                        ascending=[True, True, False]).reset_index(drop=True)

genes = [g for g in cand["gene"] if g in expr.index]
cand = cand[cand["gene"].isin(genes)].reset_index(drop=True)

# ---------- 表达矩阵 ----------
gm = pd.DataFrame({t: expr[meta[meta["group"] == t]["sample_id"]].mean(axis=1)
                   for t in TP_ORDER}).loc[cand["gene"]]
z = gm.sub(gm.mean(axis=1), axis=0).div(gm.std(axis=1).replace(0, np.nan), axis=0)

fig = plt.figure(figsize=(7.2, 8.6))
gs = GridSpec(2, 3, figure=fig, width_ratios=[1.45, 1.0, 1.0],
              height_ratios=[1.0, 0.78], hspace=0.46, wspace=0.62)

# ================= A. heatmap =================
axA = fig.add_subplot(gs[:, 0])
vmax = np.nanmax(np.abs(z.values))
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
im = axA.imshow(z.values, cmap=div_cmap(), norm=norm, aspect="auto")
axA.set_xticks(range(6))
axA.set_xticklabels(["Sham", "0.5", "1", "3", "7", "14"], rotation=45, ha="right")
axA.set_yticks(range(len(cand)))
axA.set_yticklabels(cand["gene"], fontsize=5.6)
axA.set_xlabel("Time after CCI (d)")
axA.set_title("Candidate gene expression", loc="left")
for s in axA.spines.values():
    s.set_visible(False)
axA.tick_params(length=0)

# 基因名按 WGCNA 模块着色
for i, m in enumerate(cand["WGCNA_module"]):
    t = axA.get_yticklabels()[i]
    t.set_color(MODULE_COLOR.get(m, NEUTRAL))
    t.set_fontweight("bold")

# BH 显著列（独立窄轴，避免遮挡热图）
axH = axA.inset_axes([1.035, 0, 0.055, 1], transform=axA.transAxes)
axH.set_xlim(-0.5, 0.5)
axH.set_ylim(len(cand) - 0.5, -0.5)
axH.set_xticks([])
axH.set_yticks([])
for s in axH.spines.values():
    s.set_visible(False)
for i in range(len(cand)):
    nbh = int(cand["n_tp_bh"].iloc[i])
    if nbh > 0:
        axH.scatter([0], [i], s=13, c="#8F3F3C", marker="o", lw=0)
        axH.text(0.42, i, str(nbh), ha="left", va="center", fontsize=5.0,
                 color="#8F3F3C")
axH.text(0, -1.35, "BH", ha="center", va="bottom", fontsize=5.6,
         color="#8F3F3C", fontweight="bold")

cb = fig.colorbar(im, ax=axA, fraction=0.032, pad=0.215)
cb.set_label("Row z-score", fontsize=6.4)
cb.ax.tick_params(labelsize=6.0)
cb.outline.set_linewidth(0.5)
cb.outline.set_edgecolor(SPINE)

handles = [plt.Rectangle((0, 0), 1, 1, fc=MODULE_COLOR[m]) for m in MOD_ORDER[:4]]
axA.legend(handles, [m.capitalize() for m in MOD_ORDER[:4]], ncol=2,
           loc="upper left", bbox_to_anchor=(0.0, -0.125), fontsize=5.8,
           handlelength=0.8, columnspacing=0.8, title="WGCNA module",
           title_fontsize=5.8)
panel_tag(axA, "A", dx=-0.155)

# ================= B. effect-size scatter =================
axB = fig.add_subplot(gs[0, 1:])
x = det["log2FC_rat175760"].values
y = det["log2FC_mouse5296"].values
ok = det["concordant"].values.astype(bool)

axB.axhline(0, color=SPINE, lw=0.6, ls=(0, (3, 3)), zorder=1)
axB.axvline(0, color=SPINE, lw=0.6, ls=(0, (3, 3)), zorder=1)
lim = 6.4
axB.plot([-lim, lim], [-lim, lim], color="#CFC9C2", lw=0.7,
         ls=(0, (4, 3)), zorder=1)
axB.scatter(x[ok], y[ok], s=30, c=SOFT_RED, edgecolors="white", lw=0.6,
            zorder=3, label=f"Concordant ({int(ok.sum())})")
axB.scatter(x[~ok], y[~ok], s=30, c="#FFFFFF", edgecolors=NEUTRAL, lw=0.8,
            zorder=3, label=f"Discordant ({int((~ok).sum())})")
axB.set_xlim(-lim, lim)
axB.set_ylim(-2.2, 5.2)
axB.set_xlabel("log$_2$FC, GSE175760 (rat, 7 d)")
axB.set_ylabel("log$_2$FC, GSE5296\n(mouse, 7 d)")
axB.set_title("Cross-species effect-size agreement", loc="left")
axB.legend(loc="upper left", fontsize=5.8, handletextpad=0.3)
despine(axB)
axB.grid(True, lw=0.5, alpha=0.5)
axB.set_axisbelow(True)
axB.text(0.975, 0.06, "Spearman \u03c1 = 0.649\n20 / 23 concordant",
         transform=axB.transAxes, ha="right", va="bottom", fontsize=6.2,
         color=TEXT, linespacing=1.5,
         bbox=dict(boxstyle="round,pad=0.30", fc="#FBF7F4", ec=GRID, lw=0.6))
panel_tag(axB, "B")

# ================= C. concordance by time point =================
axC = fig.add_subplot(gs[1, 1:])
lab = [f"{a}\nvs {b}" for a, b in zip(conc["gse5296_tp"], conc["gse175760_tp"])]
xs = np.arange(len(conc))
rate = conc["concordance_rate"].values * 100
bars = axC.bar(xs, rate, 0.56, color=SOFT_RED, edgecolor="white", lw=0.5, zorder=3)
for xi, v in zip(xs, rate):
    axC.text(xi, v + 2.5, f"{v:.0f}%", ha="center", va="bottom", fontsize=6.0, color=TEXT)
axC.set_xticks(xs)
axC.set_xticklabels(lab, fontsize=5.8)
axC.set_xlabel("GSE5296 time point  vs  discovery time point")
axC.set_ylabel("Direction concordance (%)")
axC.set_ylim(0, 118)
axC.set_yticks([0, 25, 50, 75, 100])
axC.set_title(f"Concordance across the injury time course  (n = "
              f"{int(conc['n_candidates_detected'].iloc[0])} detectable candidates)", loc="left")
despine(axC)
axC.grid(True, axis="y", lw=0.5, alpha=0.5)
axC.set_axisbelow(True)

ax2 = axC.twinx()
ax2.plot(xs, conc["spearman_rho"].values, "-o", color="#8F3F3C", lw=1.2, ms=3.4,
         mec="white", mew=0.5, zorder=4)
ax2.set_ylabel("Spearman \u03c1", fontsize=7.5, color="#8F3F3C")
ax2.set_ylim(0, 0.95)
ax2.tick_params(axis="y", colors="#8F3F3C")
despine(ax2, keep=("right",))
ax2.spines["right"].set_color("#8F3F3C")
ax2.spines["right"].set_linewidth(0.7)
axC.text(0.985, 0.10, "\u2014 \u25cf \u2014  Spearman \u03c1", transform=axC.transAxes,
         ha="right", va="bottom", fontsize=5.8, color="#8F3F3C")
panel_tag(axC, "C")

save(fig, "Figure_3_Candidates_and_CrossDataset")
print(f"  candidates in heatmap: {len(cand)}")
print("done.")
