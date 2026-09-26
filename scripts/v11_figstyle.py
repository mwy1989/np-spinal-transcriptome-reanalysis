"""
v11 统一绘图样式基线
--------------------
约定（与项目既有 SCI 图规范一致）：
  - 白底、无阴影、无发光、无 emoji
  - 字体 Arial；正文色 #4A4540；网格 #E5E2DE
  - Mfuzz 对齐柔色系：软红 #E8A9A4 / 软绿 #7DB89A / 软蓝 #5BC0DE / 琥珀 #F0AD4E
输出目录：<SCS_ROOT>/Figures_v11/
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt

import os as _os
BASE = _os.environ.get("SCS_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
V11 = os.path.join(BASE, "output", "v11")
FIGDIR = os.path.join(BASE, "Figures_v11")
os.makedirs(FIGDIR, exist_ok=True)

# ---------------- 色板 ----------------
TEXT = "#4A4540"
GRID = "#E5E2DE"
PANEL = "#FFFFFF"
SPINE = "#C9C3BC"

SOFT_RED = "#E8A9A4"     # 上调 / 早期爆发 / 免疫上调
SOFT_GREEN = "#7DB89A"    # 渐进上调
SOFT_BLUE = "#5BC0DE"     # 持续下调
AMBER = "#F0AD4E"         # 中晚期达峰
NEUTRAL = "#B8B2AC"       # Sham / 未分配

# 时间点梯度（Sham -> 14d）：暖色递进，晚期加深
TP_ORDER = ["Sham", "CCI_0.5d", "CCI_1d", "CCI_3d", "CCI_7d", "CCI_14d"]
TP_LABEL = {"Sham": "Sham", "CCI_0.5d": "0.5 d", "CCI_1d": "1 d",
            "CCI_3d": "3 d", "CCI_7d": "7 d", "CCI_14d": "14 d"}
TP_COLOR = {"Sham": NEUTRAL, "CCI_0.5d": "#F6D9B8", "CCI_1d": AMBER,
            "CCI_3d": SOFT_RED, "CCI_7d": "#C4685F", "CCI_14d": "#8F3F3C"}

# 四个时间软聚类
CLUSTER_COLOR = {1: SOFT_RED, 2: SOFT_GREEN, 3: AMBER, 4: SOFT_BLUE}
CLUSTER_NAME = {1: "Early burst (0.5 d)", 2: "Progressive up (14 d)",
                3: "Mid-late peak (7 d)", 4: "Persistent down (Sham)"}

# WGCNA 模块
MODULE_COLOR = {"turquoise": "#5FBFA8", "blue": "#5B8FC0",
                "brown": "#B08050", "grey": NEUTRAL,
                "not_assigned": NEUTRAL}

# 双轨 DEG
STRICT_COLOR = "#8F3F3C"
EXPLOR_COLOR = "#E8A9A4"


def apply_style():
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.titleweight": "bold",
        "axes.labelsize": 8,
        "axes.labelcolor": TEXT,
        "axes.edgecolor": SPINE,
        "axes.linewidth": 0.7,
        "axes.facecolor": PANEL,
        "figure.facecolor": PANEL,
        "savefig.facecolor": PANEL,
        "text.color": TEXT,
        "xtick.color": TEXT,
        "ytick.color": TEXT,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "legend.fontsize": 7,
        "legend.frameon": False,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.grid": False,
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "legend.handletextpad": 0.35,
        "legend.columnspacing": 0.9,
    })


def despine(ax, keep=("left", "bottom")):
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in keep)


def panel_tag(ax, tag, dx=-0.085, dy=1.045):
    ax.text(dx, dy, tag, transform=ax.transAxes, fontsize=10,
            fontweight="bold", color=TEXT, ha="left", va="bottom")


def save(fig, name):
    """同时保存 PNG(tight) 与 PDF(矢量)。"""
    png = os.path.join(FIGDIR, name + ".png")
    pdf = os.path.join(FIGDIR, name + ".pdf")
    fig.savefig(png)
    fig.savefig(pdf)
    plt.close(fig)
    print(f"  saved -> {png}")
    return png


def seq_cmap(color, name=None, light="#FFFFFF"):
    """由单色生成浅->深顺序色图，用于热图。"""
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(name or "seq", [light, color])


def div_cmap(neg=SOFT_BLUE, pos="#8F3F3C", mid="#FFFFFF", name="div"):
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(name, [neg, mid, pos])
