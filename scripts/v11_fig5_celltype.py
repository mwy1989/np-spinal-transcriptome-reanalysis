"""
Fig 5 — Cell-type expression background (GSE189070 scRNA-seq reference)
A) Relative expression of each detected candidate across 10 annotated spinal cell types
B) Dominant cell type per candidate
C) Detection rate per candidate, with the assessability threshold
D) Evidence-boundary statement
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
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

CELLS = ["Microglia", "Macrophage", "Neutrophil", "B_cell", "T_cell",
         "Astrocyte", "Oligodendrocyte", "Ependymal", "Endothelial", "Pericyte"]
DISP = {"Microglia": "Microglia", "Macrophage": "Macrophage", "Neutrophil": "Neutrophil",
        "B_cell": "B cell", "T_cell": "T cell", "Astrocyte": "Astrocyte",
        "Oligodendrocyte": "Oligodendrocyte", "Ependymal": "Ependymal",
        "Endothelial": "Endothelial", "Pericyte": "Pericyte"}

ord_a = dom[assess].sort_values(["dominant_celltype", "max_mean_expr"],
                                ascending=[True, False]).index.tolist()
ord_n = dom[~assess].sort_values("max_mean_expr", ascending=False).index.tolist()
order = ord_a + ord_n
N = len(order)
print(f"  assessable {len(ord_a)} / not assessable {len(ord_n)}")

M = mat.loc[order, CELLS]
rel = M.div(M.max(axis=1).replace(0, np.nan), axis=0)

fig = plt.figure(figsize=(7.2, 7.2))
outer = GridSpec(2, 3, figure=fig, width_ratios=[1.48, 0.34, 0.86],
                 height_ratios=[1.0, 0.40], hspace=0.44, wspace=0.30,
                 left=0.115, right=0.975, top=0.905, bottom=0.055)

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

# ================= B. dominant cell type =================
axD = fig.add_subplot(outer[0, 1])
axD.set_xlim(0, 1)
axD.set_ylim(N - 0.5, -0.5)
axD.set_xticks([])
axD.set_yticks(range(N))
labels = []
for g in order:
    if assess[g]:
        labels.append(DISP[dom.loc[g, "dominant_celltype"]])
    else:
        labels.append("")
axD.set_yticklabels(labels, fontsize=5.5)
for i, g in enumerate(order):
    t = axD.get_yticklabels()[i]
    t.set_color("#3E8C72" if (assess[g] and
                              dom.loc[g, "dominant_celltype"] == "Microglia") else TEXT)
    t.set_fontweight("bold" if assess[g] and
                     dom.loc[g, "dominant_celltype"] == "Microglia" else "normal")
for s in axD.spines.values():
    s.set_visible(False)
axD.tick_params(length=0)
axD.axhline(len(ord_a) - 0.5, color=TEXT, lw=0.9)
axD.text(0.0, len(ord_a) + 0.6, "not assessed", fontsize=5.2, color=NEUTRAL)
axD.text(-0.02, 1.030, "B   Dominant cell type", transform=axD.transAxes,
         fontsize=7.4, fontweight="bold", color=TEXT, ha="left", va="bottom")

# ================= C. detection rate =================
axB = fig.add_subplot(outer[0, 2])
ys = np.arange(N)
rates = [det.loc[g].max() * 100 for g in order]
cols = ["#5FBFA8" if assess[g] else NEUTRAL for g in order]
axB.barh(ys, rates, height=0.64, color=cols, edgecolor="white", lw=0.4, zorder=3)
axB.axvline(10, color="#8F3F3C", lw=0.8, ls=(0, (4, 3)), zorder=4)
axB.axhline(len(ord_a) - 0.5, color=TEXT, lw=0.9)
axB.set_yticks(ys)
axB.set_yticklabels([])
axB.set_ylim(N - 0.5, -0.5)
axB.set_xlim(0, 108)
axB.set_xticks([0, 50, 100])
axB.set_xlabel("Max. detection rate (%)\ndashed line = 10 % threshold", fontsize=6.4)
axB.text(-0.03, 1.030, "C   Detectability in the reference", transform=axB.transAxes,
         fontsize=7.4, fontweight="bold", color=TEXT, ha="left", va="bottom")
for s in ("top", "right", "left"):
    axB.spines[s].set_visible(False)
axB.grid(True, axis="x", lw=0.5, alpha=0.5)
axB.set_axisbelow(True)

# ================= D. evidence boundary =================
axC = fig.add_subplot(outer[1, :])
axC.set_axis_off()
axC.set_xlim(0, 10)
axC.set_ylim(0, 2.35)

axC.add_patch(FancyBboxPatch((0.05, 0.15), 4.32, 2.00,
                             boxstyle="round,pad=0.06,rounding_size=0.10",
                             fc="#F6FBF8", ec="#5FBFA8", lw=1.0))
axC.text(2.21, 1.92, "Supported by the present data", ha="center", va="top",
         fontsize=7.3, fontweight="bold", color="#3E8C72")
axC.text(0.26, 1.52,
         "\u2022  Relative expression of each candidate across the ten\n"
         "    annotated cell types of an adult mouse spinal cord reference\n"
         "\u2022  Dominant cell type and detection rate per candidate\n"
         "\u2022  16 of 29 candidates are assessable; complement and\n"
         "    microglial activation genes are microglia-dominant",
         fontsize=6.1, color=TEXT, va="top", linespacing=1.55)

axC.add_patch(FancyBboxPatch((5.55, 0.15), 4.40, 2.00,
                             boxstyle="round,pad=0.06,rounding_size=0.10",
                             fc="#FDF7F6", ec=SOFT_RED, lw=1.0))
axC.text(7.75, 1.92, "Not addressed by the present data", ha="center", va="top",
         fontsize=7.3, fontweight="bold", color="#B4635C")
axC.text(5.76, 1.52,
         "\u2022  The reference contains no neuronal cluster, so neuronal\n"
         "    expression cannot be assessed for any candidate\n"
         "\u2022  Single-cell expression is not spatial information\n"
         "\u2022  No glial response to stimulation was measured\n"
         "    (the stimulation dataset contains motoneurons only)",
         fontsize=6.1, color=TEXT, va="top", linespacing=1.55)

axC.add_patch(FancyArrowPatch((4.44, 1.12), (5.50, 1.12), arrowstyle="-|>",
                              mutation_scale=9, lw=0.9, color=NEUTRAL))
axC.text(4.97, 1.30, "evidence\nboundary", fontsize=5.5, color=NEUTRAL,
         ha="center", va="bottom", linespacing=1.3)
panel_tag(axC, "D", dx=-0.006, dy=1.10)

save(fig, "Figure_5_CellType_Expression_Background")
print("done.")
