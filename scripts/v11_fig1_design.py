"""
Fig 1 — Study design and data overview
A) 四队列构成与各自可回答的问题
B) 发现队列 18 样本 PCA（含分组质心轨迹）
C) limma 双轨 DEG 计数（探索集 raw P<0.05 / 严格集 BH<0.05）
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.gridspec import GridSpec
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v11_figstyle import (apply_style, save, despine, panel_tag, V11,
                          TEXT, GRID, SPINE, NEUTRAL, TP_ORDER, TP_LABEL, TP_COLOR,
                          STRICT_COLOR, EXPLOR_COLOR, SOFT_RED, AMBER, SOFT_BLUE)

import os as _os
BASE = _os.environ.get("SCS_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
apply_style()
print("Fig 1 ...")

expr = pd.read_csv(os.path.join(V11, "expr_log2_pseudo0.1.csv"), index_col=0)
meta = pd.read_csv(os.path.join(BASE, "output", "GSE175760_metadata.csv"))
deg = pd.read_csv(os.path.join(V11, "DEG_summary_by_timepoint.csv"))

fig = plt.figure(figsize=(7.2, 6.6))
gs = GridSpec(2, 2, figure=fig, height_ratios=[1.06, 0.94],
              hspace=0.42, wspace=0.30)

# ================= A. Cohort design =================
axA = fig.add_subplot(gs[0, :])
axA.set_axis_off()
axA.set_xlim(0, 11.2)
axA.set_ylim(0, 5.0)

def box(x, y, w, h, title, lines, fc="#FFFFFF", ec=SPINE, lw=0.9, tcolor=TEXT,
        fs_t=7.8, fs_b=6.6):
    axA.add_patch(FancyBboxPatch((x, y), w, h,
                                 boxstyle="round,pad=0.06,rounding_size=0.12",
                                 fc=fc, ec=ec, lw=lw))
    axA.text(x + w / 2, y + h - 0.28, title, ha="center", va="top",
             fontsize=fs_t, fontweight="bold", color=tcolor)
    axA.text(x + w / 2, y + h - 0.60, "\n".join(lines), ha="center", va="top",
             fontsize=fs_b, color=TEXT, linespacing=1.5)

def arrow(p1, p2, color=SPINE):
    axA.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=8,
                                  lw=0.8, color=color, shrinkA=1, shrinkB=1))

# 发现队列
box(0.10, 1.55, 3.30, 3.05,
    "Discovery cohort  |  GSE175760",
    ["Rat, chronic constriction injury",
     "Lumbar spinal cord, bulk RNA-seq",
     "6 groups \u00d7 n = 3  (18 samples)",
     "Sham + 0.5 / 1 / 3 / 7 / 14 d",
     "",
     "\u2192 time-course differential expression",
     "\u2192 co-expression modules",
     "\u2192 candidate gene selection"],
    fc="#FDF4F3", ec=SOFT_RED, lw=1.1)

# 三个外部数据集（同列纵向排列，箭头不交叉）
box(4.30, 3.50, 3.20, 1.30,
    "Cross-platform  |  GSE5296",
    ["Mouse, spinal cord injury, microarray",
     "96 samples; impact region",
     "\u2192 direction concordance of candidates"],
    fc="#FFFFFF", ec=SPINE, fs_b=6.4)

box(4.30, 2.05, 3.20, 1.30,
    "Stimulation  |  GSE243038",
    ["Mouse, spinal cord motoneurons, Smart-seq2",
     "6 groups \u00d7 3 stimulation frequencies",
     "\u2192 stimulation-associated trends"],
    fc="#FFFFFF", ec=SPINE, fs_b=6.4)

box(4.30, 0.60, 3.20, 1.30,
    "Cell-type reference  |  GSE189070",
    ["Mouse, spinal cord scRNA-seq",
     "10 annotated cell types; no neuronal cluster",
     "\u2192 descriptive expression background"],
    fc="#FFFFFF", ec=SPINE, fs_b=6.4)

arrow((3.42, 4.05), (4.27, 4.15), color=SOFT_RED)
arrow((3.42, 3.05), (4.27, 2.70), color=SOFT_RED)
arrow((3.42, 2.05), (4.27, 1.25), color=SOFT_RED)

axA.text(0.10, 0.18,
         "All analyses re-run on a single unified expression scale  [ log$_2$(FPKM + 0.1) ]  with limma-trend;  "
         "GSE5296 re-run after restoring the log$_2$ scale of the RMA matrix.",
         fontsize=6.4, color=NEUTRAL, ha="left", va="center")
panel_tag(axA, "A", dx=-0.015, dy=1.0)

# ================= B. PCA =================
axB = fig.add_subplot(gs[1, 0])
samples = [s for s in meta["sample_id"] if s in expr.columns]
X = expr[samples].T.values
X = StandardScaler().fit_transform(X)
pca = PCA(n_components=2).fit(X)
P = pca.transform(X)
gmap = dict(zip(meta["sample_id"], meta["group"]))

cent = {}
for g in TP_ORDER:
    idx = [i for i, s in enumerate(samples) if gmap[s] == g]
    cent[g] = P[idx].mean(axis=0)
    axB.scatter(P[idx, 0], P[idx, 1], s=26, c=TP_COLOR[g],
                edgecolors="white", linewidths=0.6, zorder=3,
                label=TP_LABEL[g])

for i in range(len(TP_ORDER) - 1):
    a, b = cent[TP_ORDER[i]], cent[TP_ORDER[i + 1]]
    axB.annotate("", xy=b, xytext=a,
                 arrowprops=dict(arrowstyle="-|>", color="#9A948E", lw=0.9,
                                 alpha=0.75, mutation_scale=7,
                                 connectionstyle="arc3,rad=0.12"), zorder=2)

axB.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f} %)")
axB.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f} %)")
axB.set_title("Global transcriptome structure", loc="left")
axB.legend(ncol=2, loc="lower right", fontsize=6.4, handletextpad=0.2,
           columnspacing=0.8, borderpad=0.2)
despine(axB)
axB.grid(True, lw=0.5, alpha=0.55)
axB.set_axisbelow(True)
panel_tag(axB, "B")

# ================= C. DEG counts =================
axC = fig.add_subplot(gs[1, 1])
lab = ["0.5 d", "1 d", "3 d", "7 d", "14 d"]
x = np.arange(len(lab))
w = 0.38
expl = deg["n_raw_total"].values
strict = deg["n_bh_total"].values

b1 = axC.bar(x - w / 2, expl, w, color=EXPLOR_COLOR, edgecolor="white",
             linewidth=0.5, label="Exploratory set\n(nominal P < 0.05)", zorder=3)
b2 = axC.bar(x + w / 2, strict, w, color=STRICT_COLOR, edgecolor="white",
             linewidth=0.5, label="Strict set\n(BH FDR < 0.05)", zorder=3)

for xi, v in zip(x - w / 2, expl):
    axC.text(xi, v + 26, f"{int(v)}", ha="center", va="bottom", fontsize=6.2, color=TEXT)
for xi, v in zip(x + w / 2, strict):
    if v == 0:
        axC.text(xi, 26, "0", ha="center", va="bottom", fontsize=6.2, color=STRICT_COLOR,
                 fontweight="bold")
    else:
        axC.text(xi, v + 26, f"{int(v)}", ha="center", va="bottom", fontsize=6.2,
                 color=STRICT_COLOR, fontweight="bold")

axC.set_xticks(x)
axC.set_xticklabels(lab)
axC.set_xlabel("Time after CCI")
axC.set_ylabel("Differentially expressed genes")
axC.set_ylim(0, 1520)
axC.set_yticks([0, 250, 500, 750, 1000, 1250])
axC.set_title("DEG yield before and after multiplicity control", loc="left")
axC.legend(loc="upper left", fontsize=6.2, handlelength=0.9, borderpad=0.25,
           labelspacing=0.7, bbox_to_anchor=(0.0, 1.0))
despine(axC)
axC.grid(True, axis="y", lw=0.5, alpha=0.55)
axC.set_axisbelow(True)

# 注释框
axC.text(0.985, 0.975,
         "union across time points\n"
         f"{int(deg['union_explor'].iloc[0])} exploratory   |   "
         f"{int(deg['union_strict'].iloc[0])} strict\n"
         "strict set empty at 0.5 and 1 d",
         transform=axC.transAxes, ha="right", va="top", fontsize=6.2,
         color=TEXT, linespacing=1.55,
         bbox=dict(boxstyle="round,pad=0.32", fc="#FBF7F4", ec=GRID, lw=0.6))
panel_tag(axC, "C")

save(fig, "Figure_1_Study_Design_and_Data_Overview")
print("done.")
