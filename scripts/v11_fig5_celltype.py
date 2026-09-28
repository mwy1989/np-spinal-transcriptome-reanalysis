"""
Fig 5 — Cell-type expression background (GSE189070 scRNA-seq reference)
A) Relative expression of each detected candidate across 10 annotated spinal cell types
B) Dominant cell type per candidate, after pooling the three myeloid subtypes
C) Detection rate per candidate, with the assessability threshold
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v11_figstyle import (apply_style, save, panel_tag, V11,
                          TEXT, GRID, SPINE, NEUTRAL, SOFT_RED, seq_cmap)

apply_style()
print("Fig 5 ...")

mat = pd.read_csv(os.path.join(V11, "fig5_celltype_heatmap_matrix.csv"), index_col=0)
det = pd.read_csv(os.path.join(V11, "candidate_celltype_detection_rate.csv"), index_col=0)
dom = pd.read_csv(os.path.join(V11, "candidate_dominant_celltype.csv")).set_index("gene")
assess = dom["assessable"].astype(bool)

# 合并髓系后的归属（脚本 scripts/v11_14_myeloid_merge.py 产出）
mrg = pd.read_csv(os.path.join(V11, "myeloid_merge_dominant_celltype.csv")).set_index("gene")

CELLS = ["Microglia", "Macrophage", "Neutrophil", "B_cell", "T_cell",
         "Astrocyte", "Oligodendrocyte", "Ependymal", "Endothelial", "Pericyte"]
DISP = {"Microglia": "Microglia", "Macrophage": "Macrophage", "Neutrophil": "Neutrophil",
        "B_cell": "B cell", "T_cell": "T cell", "Astrocyte": "Astrocyte",
        "Oligodendrocyte": "Oligodendrocyte", "Ependymal": "Ependymal",
        "Endothelial": "Endothelial", "Pericyte": "Pericyte",
        "Myeloid": "Myeloid"}

AMBER = "#D99A2B"   # Gapt: passes the 10 % detection threshold but fails the 0.1 expression one
TEAL = "#3E8C72"

# 合并后失去髓系优势、且原始标签为小胶质者：加星号
FLIP = {"Cfh", "Nrp1", "Lmo2", "Arhgap25"}


def merged_label(g):
    lab = DISP[mrg.loc[g, "dominant_merged_cells"]]
    return lab + "*" if g in FLIP else lab


ord_a = dom[assess].assign(_lab=[DISP[mrg.loc[g, "dominant_merged_cells"]] for g in dom[assess].index]) \
    .sort_values(["_lab", "max_mean_expr"], ascending=[True, False]).index.tolist()
ord_n = dom[~assess].sort_values("max_mean_expr", ascending=False).index.tolist()
order = ord_a + ord_n
N = len(order)
print(f"  assessable {len(ord_a)} / not assessable {len(ord_n)}")

M = mat.loc[order, CELLS]
rel = M.div(M.max(axis=1).replace(0, np.nan), axis=0)

fig = plt.figure(figsize=(7.2, 4.55))
outer = GridSpec(1, 3, figure=fig, width_ratios=[1.48, 0.36, 0.86],
                 wspace=0.30, left=0.132, right=0.975, top=0.885, bottom=0.155)

# ================= A. heatmap =================
axA = fig.add_subplot(outer[0, 0])
im = axA.imshow(rel.values, cmap=seq_cmap(SOFT_RED, "rel", light="#FFFFFF"),
                vmin=0, vmax=1, aspect="auto")
axA.set_xticks(range(len(CELLS)))
axA.set_xticklabels([DISP[c] for c in CELLS], rotation=45, ha="right", fontsize=5.8)
axA.set_yticks(range(N))
axA.set_yticklabels(order, fontsize=5.6)
for i, g in enumerate(order):
    axA.get_yticklabels()[i].set_color(TEXT if assess[g] else NEUTRAL)
    axA.get_yticklabels()[i].set_fontweight("bold" if assess[g] else "normal")
axA.axhline(len(ord_a) - 0.5, color=TEXT, lw=0.9)
for s in axA.spines.values():
    s.set_visible(False)
axA.tick_params(length=0)
axA.set_title("Relative expression across spinal cell types", loc="left", pad=16)

cax = axA.inset_axes([0.63, 1.035, 0.37, 0.022])
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label("Fraction of gene maximum", fontsize=5.8, labelpad=2)
cb.ax.tick_params(labelsize=5.4, length=1.8, pad=1)
cb.outline.set_linewidth(0.5)
cb.outline.set_edgecolor(SPINE)
panel_tag(axA, "A", dx=-0.155)

# ================= B. dominant cell type (myeloid pooled) =================
axD = fig.add_subplot(outer[0, 1])
axD.set_xlim(0, 1)
axD.set_ylim(N - 0.5, -0.5)
axD.set_xticks([])
axD.set_yticks(range(N))
labels = [merged_label(g) if assess[g] else "" for g in order]
axD.set_yticklabels(labels, fontsize=5.5)
for i, g in enumerate(order):
    t = axD.get_yticklabels()[i]
    is_my = assess[g] and mrg.loc[g, "dominant_merged_cells"] == "Myeloid"
    t.set_color(TEAL if is_my else TEXT)
    t.set_fontweight("bold" if is_my else "normal")
for s in axD.spines.values():
    s.set_visible(False)
axD.tick_params(length=0)
axD.axhline(len(ord_a) - 0.5, color=TEXT, lw=0.9)
axD.text(0.0, len(ord_a) + 0.6, "not assessed", fontsize=5.2, color=NEUTRAL)
axD.text(-0.02, 1.028, "B   Dominant cell type\n      (myeloid pooled)", transform=axD.transAxes,
         fontsize=7.4, fontweight="bold", color=TEXT, ha="left", va="bottom", linespacing=1.35)

# ================= C. detection rate =================
axB = fig.add_subplot(outer[0, 2])
ys = np.arange(N)
rates = [det.loc[g].max() * 100 for g in order]
cols = ["#5FBFA8" if assess[g] else (AMBER if g == "Gapt" else NEUTRAL) for g in order]
axB.barh(ys, rates, height=0.64, color=cols, edgecolor="white", lw=0.4, zorder=3)
axB.axvline(10, color="#8F3F3C", lw=0.8, ls=(0, (4, 3)), zorder=4)
axB.axhline(len(ord_a) - 0.5, color=TEXT, lw=0.9)
_gapt_y = order.index("Gapt")
axB.annotate("Gapt: detection 11 %\nbut max mean 0.09", xy=(11.4, _gapt_y),
             xytext=(26, _gapt_y + 2.1), fontsize=5.0, color="#A87213",
             ha="left", va="center",
             arrowprops=dict(arrowstyle="-", lw=0.6, color="#A87213",
                             shrinkA=1, shrinkB=1))
axB.set_yticks(ys)
axB.set_yticklabels([])
axB.set_ylim(N - 0.5, -0.5)
axB.set_xlim(0, 108)
axB.set_xticks([0, 50, 100])
axB.set_xlabel("Max. detection rate (%)\ndashed line = 10 % detection threshold", fontsize=6.4)
axB.text(-0.03, 1.028, "C   Detectability in the reference", transform=axB.transAxes,
         fontsize=7.4, fontweight="bold", color=TEXT, ha="left", va="bottom")
for s in ("top", "right", "left"):
    axB.spines[s].set_visible(False)
axB.grid(True, axis="x", lw=0.5, alpha=0.5)
axB.set_axisbelow(True)

axB.text(0.0, 1.20,
         "* attribution differs from the original subtype label (panel B)",
         transform=axB.transAxes, fontsize=5.2, color=NEUTRAL, ha="left", va="bottom")

save(fig, "Figure_5_CellType_Expression_Background")
print("done.")
