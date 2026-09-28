"""
Fig 6 — Analysis-criterion spectrum, framework dependence and candidate-gene robustness
A) DEG yield across ten analysis criteria (union), sorted; 13.4-fold span
B) Head-to-head with the published criterion of the data originators (He 2025)
C) DEG counts per criterion x time point (criteria grouped by multiple-testing strategy)
D) Pairwise Jaccard overlap between criteria
E) Pass/fail grid of the 29 candidates across the ten criteria
F) Which criteria exclude the nine exceptions (effect size vs criteria passed)

Input : outputs/novelty_v11/B_threshold_spectrum.csv
        outputs/novelty_v11/B_jaccard_matrix.csv
        outputs/novelty_v11/B_candidate_robustness.csv
Output: Figures_v11/Figure_6_Analysis_Criterion_Spectrum.{png,pdf}
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle, FancyBboxPatch
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v11_figstyle import (apply_style, save, panel_tag, V11,
                          TEXT, GRID, SPINE, NEUTRAL, SOFT_RED, AMBER, seq_cmap)

import os as _os
BASE = _os.environ.get("SCS_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
apply_style()
print("Fig 6 ...")

OUT = os.path.join(BASE, "outputs", "novelty_v11")
spec = pd.read_csv(os.path.join(OUT, "B_threshold_spectrum.csv"))
jac = pd.read_csv(os.path.join(OUT, "B_jaccard_matrix.csv"), index_col=0)
rob = pd.read_csv(os.path.join(OUT, "B_candidate_robustness.csv"))

TPC = ["0.5d", "1d", "3d", "7d", "14d"]
TPD = ["0.5 d", "1 d", "3 d", "7 d", "14 d"]

# ---- 口径定义（row index in spec -> 展示标签） -------------------------------
# 分为三组：FDR 校正组 / 已发表口径 / 名义 P 组
DEFS = {
    "C1":  ("BH FDR < 0.05 & |log2FC| > 0.58",                          "fdr"),
    "C3":  ("BH FDR < 0.05, no effect-size filter",                     "fdr"),
    "C7":  ("BH FDR < 0.05 & |log2FC| >= 1 (two-fold)",                 "fdr"),
    "C10": ("BH FDR < 0.05 & |log2FC| > 0.58 & FPKM >= 0.5",            "fdr"),
    "C6":  ("P <= 0.05 & FC >= 1.5 & FPKM >= 0.5 (published, He 2025)", "pub"),
    "C2":  ("nominal P < 0.05 & |log2FC| > 0.58",                       "nom"),
    "C4":  ("nominal P < 0.05, no effect-size filter",                  "nom"),
    "C5":  ("nominal P < 0.01 & |log2FC| > 0.58",                       "nom"),
    "C8":  ("nominal P < 0.05 & |log2FC| >= 1 (two-fold)",              "nom"),
    "C9":  ("nominal P < 0.01 & |log2FC| >= 1 (two-fold)",              "nom"),
}
# spec.csv 行序：1,2,3,4,5,6,7,8,9,10 -> 映射到 C 编号
rowmap = {"C1": 0, "C2": 1, "C3": 2, "C4": 3, "C5": 4,
          "C6": 5, "C7": 6, "C8": 7, "C9": 8, "C10": 9}
ORDER = list(DEFS.keys())                      # 展示顺序（按组）
GCOL = {"fdr": "#8F3F3C", "pub": AMBER, "nom": SOFT_RED}

uni = {c: int(spec.iloc[rowmap[c]]["union"]) for c in ORDER}
per_tp = {c: [int(spec.iloc[rowmap[c]][t]) for t in TPC] for c in ORDER}

fig = plt.figure(figsize=(7.2, 9.2))
gs = GridSpec(4, 2, figure=fig, height_ratios=[1.00, 0.92, 1.02, 0.74],
              hspace=0.62, wspace=0.42,
              left=0.115, right=0.975, top=0.945, bottom=0.055)

# ============================================================ A 谱
axA = fig.add_subplot(gs[0, 0])
srt = sorted(ORDER, key=lambda c: uni[c])
y = np.arange(len(srt))
axA.barh(y, [uni[c] for c in srt], height=0.66,
         color=[GCOL[DEFS[c][1]] for c in srt],
         edgecolor="white", lw=0.5, zorder=3)
for i, c in enumerate(srt):
    axA.text(uni[c] * 1.06, i, f"{uni[c]:,}", va="center", ha="left",
             fontsize=5.8, color=TEXT, zorder=4)
axA.set_yticks(y)
axA.set_yticklabels(srt, fontsize=6.2)
for i, c in enumerate(srt):
    axA.get_yticklabels()[i].set_fontweight("bold" if DEFS[c][1] == "pub" else "normal")
    axA.get_yticklabels()[i].set_color(GCOL[DEFS[c][1]])
axA.set_xscale("log")
axA.set_xlim(230, 12000)
axA.set_xticks([300, 1000, 3000, 10000])
axA.set_xticklabels(["300", "1,000", "3,000", "10,000"], fontsize=6)
axA.set_xlabel("Differentially expressed genes (union across time points)",
               fontsize=6.6)
axA.grid(True, axis="x", lw=0.5, alpha=0.45)
axA.set_axisbelow(True)
for s in ("top", "right", "left"):
    axA.spines[s].set_visible(False)
axA.tick_params(length=0)
axA.set_title("Yield of one dataset under ten analysis criteria", loc="left", pad=16)
axA.text(0.985, 0.055, "13.4-fold span", transform=axA.transAxes, fontsize=6.4,
         color="#8F3F3C", fontweight="bold", ha="right", va="bottom")
_h = [plt.Line2D([], [], marker="s", ls="", ms=4.2, mfc=c, mec="white", label=l)
      for c, l in [("#8F3F3C", "FDR-controlled"), (AMBER, "published criterion"),
                   (SOFT_RED, "nominal P")]]
axA.legend(handles=_h, loc="lower right", fontsize=5.3, handlelength=0.7,
           handletextpad=0.35, labelspacing=0.22, borderaxespad=0.15,
           bbox_to_anchor=(1.0, 0.10))
panel_tag(axA, "A", dx=-0.10, dy=1.02)

# ============================================================ B 与已发表口径对照
axB = fig.add_subplot(gs[0, 1])
x = np.arange(len(TPC))
vals = per_tp["C6"]
axB.bar(x, vals, width=0.60, color=AMBER, edgecolor="white", lw=0.5, zorder=3,
        label="This study, published criterion applied")
axB.plot([3, 4], [529, 352], ls="none", marker="D", ms=4.4, mfc="white",
         mec="#8F3F3C", mew=1.0, zorder=5, label="Reported by the data originators")
for xi, v in zip(x, vals):
    axB.text(xi, v + 26, f"{v}", ha="center", va="bottom", fontsize=5.6, color=TEXT)
for xi, v in [(3, 529), (4, 352)]:
    axB.text(xi - 0.12, v, f"{v}", ha="right", va="center", fontsize=5.6,
             color="#8F3F3C", fontweight="bold")
axB.set_xticks(x)
axB.set_xticklabels(TPD, fontsize=6.4)
axB.set_ylim(0, 1040)
axB.set_ylabel("Genes passing the criterion", fontsize=6.6)
axB.set_xlabel("Time after injury", fontsize=6.6)
axB.grid(True, axis="y", lw=0.5, alpha=0.45)
axB.set_axisbelow(True)
for s in ("top", "right"):
    axB.spines[s].set_visible(False)
axB.tick_params(length=2.2)
axB.set_title("Same threshold, different implementation", loc="left", pad=16)
axB.legend(loc="upper left", fontsize=5.5, handlelength=1.1,
           borderaxespad=0.15, labelspacing=0.25)
axB.text(0.985, 0.055, "1.6-2.0x more genes", transform=axB.transAxes,
         fontsize=6.4, color="#8F3F3C", fontweight="bold", ha="right", va="bottom")
panel_tag(axB, "B", dx=-0.13, dy=1.02)

# ============================================================ C 逐时点计数
axC = fig.add_subplot(gs[1, :])
M = np.array([per_tp[c] + [uni[c]] for c in ORDER], dtype=float)
disp = M / np.log1p(M.max())          # log 压缩后归一化到 0-1
im = axC.imshow(disp, cmap=seq_cmap("#8F3F3C", "cnt", light="#FFFFFF"),
                vmin=0, vmax=1, aspect="auto")
axC.set_xticks(range(6))
axC.set_xticklabels(TPD + ["union"], fontsize=6.6)
axC.set_yticks(range(len(ORDER)))
axC.set_yticklabels([f"{c}   {DEFS[c][0]}" for c in ORDER], fontsize=5.8)
for i, c in enumerate(ORDER):
    axC.get_yticklabels()[i].set_color(GCOL[DEFS[c][1]])
    axC.get_yticklabels()[i].set_fontweight("bold" if DEFS[c][1] == "pub" else "normal")
for i in range(len(ORDER)):
    for j in range(6):
        v = int(M[i, j])
        axC.text(j, i, f"{v:,}" if v else "0", ha="center", va="center",
                 fontsize=5.5,
                 color="#FFFFFF" if disp[i, j] > 0.55 else (NEUTRAL if v == 0 else TEXT),
                 fontweight="bold" if v == 0 else "normal")
for i in (3.5, 4.5):                  # 组间分隔线
    axC.axhline(i, color=TEXT, lw=1.0)
axC.axvline(4.5, color=TEXT, lw=0.9)
# 圈出「FDR 校正口径下 0.5/1 d 为空」区块
axC.add_patch(Rectangle((-0.5, -0.5), 2.0, 4.0, fill=False,
                        ec="#8F3F3C", lw=1.1, ls=(0, (3, 2)), zorder=6))
axC.text(0.5, -0.62, "empty under every FDR-controlled criterion",
         ha="center", va="bottom", fontsize=5.8, color="#8F3F3C", fontweight="bold")
axC.text(6.10, 0.5, "FDR-controlled", rotation=270, fontsize=5.8, color="#8F3F3C",
         ha="left", va="center", fontweight="bold")
axC.text(6.10, 4.0, "published", rotation=270, fontsize=5.8, color="#B07A22",
         ha="left", va="center", fontweight="bold")
axC.text(6.10, 7.0, "nominal P", rotation=270, fontsize=5.8, color="#B4635C",
         ha="left", va="center", fontweight="bold")
for s in axC.spines.values():
    s.set_visible(False)
axC.tick_params(length=0)
axC.set_title("Differential-expression yield by criterion and time point",
              loc="left", pad=18)
panel_tag(axC, "C", dx=-0.145, dy=1.035)

# ============================================================ D Jaccard
axD = fig.add_subplot(gs[2, 0])
J = jac.loc[[spec.iloc[rowmap[c]]["criterion"] for c in ORDER],
            [spec.iloc[rowmap[c]]["criterion"] for c in ORDER]].values
imD = axD.imshow(J, cmap=seq_cmap("#5B8FC0", "jac", light="#FFFFFF"),
                 vmin=0, vmax=1, aspect="equal")
axD.set_xticks(range(len(ORDER)))
axD.set_yticks(range(len(ORDER)))
axD.set_xticklabels(ORDER, fontsize=5.8, rotation=90)
axD.set_yticklabels(ORDER, fontsize=5.8)
for i in range(len(ORDER)):
    for j in range(len(ORDER)):
        axD.text(j, i, f"{J[i, j]:.2f}", ha="center", va="center", fontsize=4.7,
                 color="#FFFFFF" if J[i, j] > 0.62 else TEXT)
    axD.get_yticklabels()[i].set_color(GCOL[DEFS[ORDER[i]][1]])
for p in ((0, 5), (5, 0)):            # C1 <-> C2
    axD.add_patch(Rectangle((p[1] - 0.5, p[0] - 0.5), 1, 1, fill=False,
                            ec="#8F3F3C", lw=1.2, zorder=6))
for s in axD.spines.values():
    s.set_visible(False)
axD.tick_params(length=0)
axD.set_title("Overlap between criteria", loc="left", pad=16)
axD.text(1.04, 0.30,
         "The strict set is a subset of the\nexploratory set: 386 of 1,940 genes\n(19.9 %) survive correction\n"
         "(C1 vs C2: Jaccard 0.20)",
         transform=axD.transAxes, fontsize=5.7, color="#8F3F3C", ha="left",
         va="center", linespacing=1.5)
axD.text(1.04, 0.86, "Jaccard\nindex", transform=axD.transAxes, fontsize=5.7,
         color=TEXT, fontweight="bold", ha="left", va="center", linespacing=1.4)
panel_tag(axD, "D", dx=-0.145, dy=1.045)

# ============================================================ F 脆弱成因
axF = fig.add_subplot(gs[2, 1])
colmap = {c: f"{c} {DEFS[c][0]}" for c in ORDER}
rob["npass"] = rob["n_criteria_pass"]
xs = rob["max_abs_log2FC"].values
ys = rob["npass"].values
cols = ["#C4685F" if n >= 10 else (AMBER if n == 8 else SOFT_RED) for n in ys]
axF.axvline(1.0, color="#8F3F3C", lw=0.9, ls=(0, (4, 3)), zorder=2)
axF.scatter(xs, ys, s=17, c=cols, edgecolor="white", lw=0.5, zorder=4)
for _, r in rob.iterrows():                 # 仅单独标注低表达导致的掉队者
    if r["npass"] == 8:
        axF.annotate(r["gene"], (r["max_abs_log2FC"], r["npass"]),
                     textcoords="offset points", xytext=(4.5, 1.0),
                     fontsize=5.4, color=TEXT, fontweight="bold")
axF.set_xscale("log")
axF.set_xlim(0.4, 13)
axF.set_ylim(6.4, 10.8)
axF.set_yticks([7, 8, 10])
axF.set_xticks([0.5, 1, 2, 5, 10])
axF.set_xticklabels(["0.5", "1", "2", "5", "10"], fontsize=6)
axF.set_xlabel("Largest |log$_2$ fold change| of the gene", fontsize=6.2)
axF.set_ylabel("Number of criteria passed", fontsize=6.2)
for s in ("top", "right"):
    axF.spines[s].set_visible(False)
axF.grid(True, axis="y", lw=0.5, alpha=0.4)
axF.set_axisbelow(True)
axF.tick_params(length=2.2)
axF.set_title("Effect size accounts for the nine exceptions", loc="left", pad=16)
axF.annotate("8 genes below the two-fold\nthreshold lose only the\ncriteria that require |log2FC| >= 1",
             xy=(1.30, 7.05), xytext=(3.05, 9.35), fontsize=5.4, color="#8F3F3C",
             ha="left", va="center", linespacing=1.5,
             arrowprops=dict(arrowstyle="-", lw=0.7, color="#8F3F3C",
                             shrinkA=2, shrinkB=2))
axF.text(2.55, 7.72, "1 gene (Chrna1) is lost by the\nexpression filter (low abundance)",
         fontsize=5.4, color=TEXT, ha="left", va="center", linespacing=1.5)
panel_tag(axF, "F", dx=-0.145, dy=1.045)

# ============================================================ E 候选稳健性网格
axE = fig.add_subplot(gs[3, :])
robs = rob.sort_values(["npass", "max_abs_log2FC"], ascending=[False, False])
genes = robs["gene"].tolist()
# rob 的前 10 个数据列即按口径 ①…⑩ 排列（列名与 spec 措辞略有差异，故按位置映射）
crit_cols = list(rob.columns[1:11])
assert len(crit_cols) == 10, rob.columns.tolist()
key = {c: crit_cols[rowmap[c]] for c in ORDER}
G = np.array([[int(robs.iloc[i][key[c]]) for c in ORDER] for i in range(len(genes))])
assert G.sum() == robs["npass"].sum(), (G.sum(), robs["npass"].sum())
axE.imshow(G.T, cmap=plt.matplotlib.colors.ListedColormap(["#F2EFEC", "#C4685F"]),
           vmin=0, vmax=1, aspect="auto")
axE.set_xticks(range(len(genes)))
axE.set_xticklabels(genes, rotation=90, fontsize=5.1)
for i, g in enumerate(genes):
    axE.get_xticklabels()[i].set_color(TEXT if robs.iloc[i]["npass"] == 10 else "#B4635C")
axE.set_yticks(range(len(ORDER)))
axE.set_yticklabels(ORDER, fontsize=6.0)
for i, c in enumerate(ORDER):
    axE.get_yticklabels()[i].set_color(GCOL[DEFS[c][1]])
for i in (3.5, 4.5):
    axE.axhline(i, color=TEXT, lw=1.0)
axE.axvline(19.5, color=TEXT, lw=0.9)
axE.axvline(20.5, color=TEXT, lw=0.9)
axE.text(9.5, -0.85, "pass all ten criteria (20)", ha="center", va="bottom",
         fontsize=5.4, color="#8F3F3C", fontweight="bold")
axE.text(20.5, -0.85, "1 gene passes 8 (Chrna1)", ha="center", va="bottom",
         fontsize=5.4, color=TEXT)
axE.text(25.0, -0.85, "8 genes pass seven", ha="center", va="bottom",
         fontsize=5.4, color="#B4635C")
for s in axE.spines.values():
    s.set_visible(False)
axE.tick_params(length=0)
axE.set_title("Pass / fail of the 29 candidate genes across the ten criteria",
              loc="left", pad=30)
panel_tag(axE, "E", dx=-0.145, dy=1.16)

save(fig, "Figure_6_Analysis_Criterion_Spectrum")
print("done.")
