"""
Fig 4 — Stimulation-associated expression trends (GSE243038)
A) Slope plot (broken y-axis): Sham -> SCI -> SCS group means, all evaluable candidates
B) Per-sample expression for candidates trending back toward Sham
   (purely descriptive; no P values are reported)
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v11_figstyle import (apply_style, save, despine, panel_tag, V11,
                          TEXT, GRID, SPINE, NEUTRAL, SOFT_GREEN, SOFT_BLUE)

apply_style()
print("Fig 4 ...")

d = pd.read_csv(os.path.join(V11, "scs_candidate_expression.csv")).set_index("gene")
FLOOR = -3.321928094887362
KEYS = ["Sham", "SCI", "SCS"]
GROUPS = {
    "Sham": ("Sham_mean", ["Sham_Sham1", "Sham_Sham2", "Sham_Sham3"]),
    "SCI": ("SCI_untrained_mean", ["SCI_untrained_Untrained1",
                                   "SCI_untrained_Untrained2",
                                   "SCI_untrained_Untrained3"]),
    "SCS": ("SCS_10_15_20Hz_mean", ["SCS_10_15_20Hz_SCS1",
                                    "SCS_10_15_20Hz_SCS2",
                                    "SCS_10_15_20Hz_SCS3"]),
}

ev = d[d["SCI_shift"].abs() > 0.1].copy()
ev = ev.sort_values("SCI_shift", ascending=False)
tr = ev[ev["SCS_10_15_20Hz_toward_sham"]].index.tolist()
GREEN = "#4E9E7E"
print(f"  evaluable: {len(ev)}   toward-Sham trend: {len(tr)}")

fig = plt.figure(figsize=(7.2, 7.2))
outer = GridSpec(1, 2, figure=fig, width_ratios=[1.0, 1.20], wspace=0.30,
                 left=0.095, right=0.985, top=0.875, bottom=0.085)

lcol = outer[0].subgridspec(2, 1, height_ratios=[0.14, 0.86], hspace=0.10)
axAT = fig.add_subplot(lcol[0])
axAB = fig.add_subplot(lcol[1])
xs = [0, 1, 2]

for ax in (axAT, axAB):
    for g in ev.index:
        vals = [d.loc[g, GROUPS[k][0]] for k in KEYS]
        is_tr = bool(d.loc[g, "SCS_10_15_20Hz_toward_sham"])
        ax.plot(xs, vals, "-o", ms=2.6, lw=1.1,
                color=SOFT_GREEN if is_tr else "#D5CFC8",
                mec="white", mew=0.4, zorder=3 if is_tr else 2)
    despine(ax, keep=("left",))
    ax.grid(True, axis="y", lw=0.5, alpha=0.5)
    ax.set_axisbelow(True)
    ax.set_xlim(-0.22, 3.30)

axAT.set_ylim(7.9, 10.1)
axAT.set_yticks([8, 9])
axAT.set_xticks([])
axAB.set_ylim(-4.35, 7.3)
axAB.set_yticks([-3, 0, 3, 6])
axAB.set_xticks(xs)
axAB.set_xticklabels(["Sham", "SCI\n(untrained)", "SCI +\nSCS"])
axAB.set_ylabel("Expression  [ log$_2$(FPKM + 0.1) ]")

kw = dict(color=SPINE, lw=0.7, clip_on=False)
axAT.plot([-0.012, 0.012], [0, 0], transform=axAT.transAxes, **kw)
axAB.plot([-0.012, 0.012], [1, 1], transform=axAB.transAxes, **kw)

lab_y = {"top": {}, "bot": {}}
for g in tr:
    v = d.loc[g, GROUPS["SCS"][0]]
    lab_y["top" if v > 7.3 else "bot"][g] = v

for key, ax, gap in [("top", axAT, 0.55), ("bot", axAB, 0.62)]:
    order = sorted(lab_y[key], key=lambda k: lab_y[key][k])
    ys = [lab_y[key][k] for k in order]
    for i in range(1, len(ys)):
        if ys[i] - ys[i - 1] < gap:
            ys[i] = ys[i - 1] + gap
    for g, yl in zip(order, ys):
        v = d.loc[g, GROUPS["SCS"][0]]
        ax.annotate(g, xy=(2, v), xytext=(2.17, yl), fontsize=5.5, va="center",
                    ha="left", color=GREEN, fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color=SOFT_GREEN, lw=0.5,
                                    shrinkA=0.5, shrinkB=1.5)
                    if abs(v - yl) > 0.14 else None)

axAB.axhline(FLOOR, color=SPINE, lw=0.6, ls=(0, (3, 3)))
axAB.text(0.10, FLOOR - 0.42, "zero expression  [ log$_2$(0 + 0.1) ]", fontsize=5.3, color=NEUTRAL)
h = [plt.Line2D([], [], color=SOFT_GREEN, marker="o", ms=2.6, lw=1.1, mec="white"),
     plt.Line2D([], [], color="#D5CFC8", marker="o", ms=2.6, lw=1.1, mec="white")]
axAB.legend(h, ["trend toward Sham", "no trend"], loc="upper left", fontsize=5.8,
            handlelength=1.5)
axAT.set_title("Group-level expression,\nall evaluable candidates", loc="left",
               fontsize=8.4, pad=6)
panel_tag(axAT, "A", dx=-0.10, dy=1.20)

inner = outer[1].subgridspec(4, 2, hspace=0.68, wspace=0.34)
for k, g in enumerate(tr):
    ax = fig.add_subplot(inner[k // 2, k % 2])
    for xi, key in enumerate(KEYS):
        _, cols = GROUPS[key]
        vals = [d.loc[g, c] for c in cols]
        ax.scatter(np.full(3, xi) + np.array([-0.11, 0, 0.11]), vals,
                   s=12, c=SOFT_BLUE if key == "SCS" else "#CFC9C2",
                   edgecolors="white", lw=0.4, zorder=3)
        ax.plot([xi - 0.24, xi + 0.24], [np.mean(vals)] * 2, color=TEXT,
                lw=1.1, zorder=4)
    ax.plot([0, 1, 2], [d.loc[g, GROUPS[k2][0]] for k2 in KEYS],
            color="#D5CFC8", lw=0.8, zorder=2)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(KEYS, fontsize=5.4)
    ax.set_title(g, fontsize=6.8, pad=2.5, color=GREEN)
    ax.set_ylim(FLOOR - 0.5, 6.8)
    ax.set_yticks([-3, 0, 3, 6])
    ax.tick_params(labelsize=5.2)
    despine(ax)
    if k % 2 == 0:
        ax.set_ylabel("log$_2$(FPKM+0.1)", fontsize=5.4)
    else:
        ax.set_yticklabels([])

fig.text(0.455, 0.935, "Per-sample expression, candidates trending back toward Sham",
         fontsize=8.4, fontweight="bold", color=TEXT)
fig.text(0.985, 0.028,
         "SCS samples (left to right) received 10, 15 and 20 Hz; horizontal bars denote group means. "
         "Descriptive only \u2014 no P values are reported.",
         fontsize=5.6, color=NEUTRAL, ha="right")
fig.text(0.415, 0.935, "B", fontsize=10, fontweight="bold", color=TEXT)

save(fig, "Figure_4_Stimulation_Associated_Trends")
print("done.")
