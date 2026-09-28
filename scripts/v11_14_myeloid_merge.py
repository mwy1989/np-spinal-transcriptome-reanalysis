"""
v11 阶段 2 / 任务 37：髓系亚型合并后的优势细胞类型重算

背景：
  GSE189070 原始注释把髓系细胞拆成 Microglia / Macrophage / Neutrophil 三类，
  但该参考自身 26 个簇中有 4 个与其被赋予的髓系亚类标签不完全吻合
  （见 Table 1 局限栏与 Methods 2.11），因此"小胶质 vs 巨噬细胞"这一层的
  归属可疑，而"是否属于髓系"这一层是稳的。

本脚本回答一个明确的问题：
  把 Microglia + Macrophage + Neutrophil 三类细胞**合并后重新计算**平均表达
  （而不是仅把标签改称髓系），各候选基因的优势细胞类型是否改变？

两种合并方式都算，以检验结论对合并权重的敏感性：
  A. 细胞级合并（cells）：把三类细胞的所有细胞合成一个集合后重算均值，
     等价于按各类细胞数加权 —— 这是"重新计算"的严格形式。
  B. 亚型均值等权（equal-weight）：三类各自均值的算术平均。

实现说明：
  该 h5ad 由较新版本 anndata 写出（obs/_index 与 var/GeneSymbol 为
  nullable-string 编码），anndata 0.10.8 无法直接读入，故用 h5py 直读
  CSR 三件套并做列切片。数值与 v11_13_celltype_annotation.py 的产出做了
  交叉校验（原始 10 类口径下优势细胞类型必须逐基因一致）。

输出：
  output/v11/myeloid_merge_dominant_celltype.csv
  output/v11/myeloid_merge_summary.csv
"""
import os
import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csr_matrix

BASE = os.environ.get("SCS_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "output", "v11")
H5 = os.path.join(BASE, "output", "scRNA_GSE189070", "GSE189070_annotated.h5ad")

CT_ORDER = ["Microglia", "Macrophage", "Neutrophil", "Astrocyte", "Oligodendrocyte",
            "Ependymal", "Endothelial", "Pericyte", "B_cell", "T_cell"]
MYELOID = ["Microglia", "Macrophage", "Neutrophil"]
NON_MYELOID = [c for c in CT_ORDER if c not in MYELOID]


def _dec(x):
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)


# ---------- 1. 读参考（h5py 直读）----------
print("读取 GSE189070 注释数据...")
f = h5py.File(H5, "r")
shape = tuple(int(v) for v in f["X"].attrs["shape"])
X = csr_matrix((f["X/data"][:], f["X/indices"][:], f["X/indptr"][:]), shape=shape)
X = X.tocsc()                                   # 按列切片 O(nnz_col)
cats = [_dec(c) for c in f["obs/celltype/categories"][:]]
codes = np.asarray(f["obs/celltype/codes"][:])
names = [_dec(x) for x in f["var/GeneSymbol/values"][:]]
f.close()
print("  细胞 x 基因:", shape)
print("  细胞类型:", dict(zip(cats, [int((codes == i).sum()) for i in range(len(cats))])))

celltype = np.array([cats[c] for c in codes])
mye_mask = np.isin(celltype, MYELOID)
print("  髓系细胞合计:", int(mye_mask.sum()))

upper = {}
for i, n in enumerate(names):
    if n and n.upper() not in ("NAN", ""):
        upper.setdefault(n.upper(), i)

cand = pd.read_csv(os.path.join(OUT, "candidates_primary.csv"))
genes = cand["gene"].tolist()
present = [(g, upper[g.upper()]) for g in genes if g.upper() in upper]
absent = [g for g in genes if g.upper() not in upper]
print(f"  候选 {len(genes)} 个；参考中可检测 {len(present)}；缺失 {len(absent)}: {absent}")

orig = pd.read_csv(os.path.join(OUT, "candidate_dominant_celltype.csv")).set_index("gene")


def col_mean(j, mask):
    v = np.asarray(X[:, j].todense()).ravel()
    return float(v[mask].mean()) if mask.sum() else np.nan


def col_det(j, mask):
    v = np.asarray(X[:, j].todense()).ravel()
    return float((v[mask] > 0).mean()) if mask.sum() else np.nan


# ---------- 2. 逐基因计算 ----------
masks = {ct: (celltype == ct) for ct in CT_ORDER}
rows = []
for g, j in present:
    v = np.asarray(X[:, j].todense()).ravel()
    m10 = {ct: float(v[masks[ct]].mean()) if masks[ct].sum() else np.nan for ct in CT_ORDER}
    d10 = {ct: float((v[masks[ct]] > 0).mean()) if masks[ct].sum() else np.nan for ct in CT_ORDER}
    mye_v = v[mye_mask]
    mye_cells_mean = float(mye_v.mean())
    mye_equal = float(np.mean([m10[c] for c in MYELOID]))

    m8_cells = {"Myeloid": mye_cells_mean, **{c: m10[c] for c in NON_MYELOID}}
    m8_equal = {"Myeloid": mye_equal, **{c: m10[c] for c in NON_MYELOID}}
    dom10 = max(m10, key=m10.get)
    best_nonmye = max(NON_MYELOID, key=lambda c: m10[c])

    rows.append({
        "gene": g,
        "dominant_original_10ct": dom10,
        "dominant_merged_cells": max(m8_cells, key=m8_cells.get),
        "dominant_merged_equalweight": max(m8_equal, key=m8_equal.get),
        "max_mean_expr_10ct": round(m10[dom10], 4),
        "max_detection_rate_10ct": round(d10[dom10], 3),
        "myeloid_mean_cells": round(mye_cells_mean, 4),
        "myeloid_mean_equalweight": round(mye_equal, 4),
        "best_nonmyeloid_type": best_nonmye,
        "best_nonmyeloid_mean": round(m10[best_nonmye], 4),
        "is_myeloid_dominant_cells": max(m8_cells, key=m8_cells.get) == "Myeloid",
        "is_myeloid_dominant_equalweight": max(m8_equal, key=m8_equal.get) == "Myeloid",
    })

df = pd.DataFrame(rows).sort_values("gene").reset_index(drop=True)
df.to_csv(os.path.join(OUT, "myeloid_merge_dominant_celltype.csv"), index=False)

# ---------- 3. 与既有注释表交叉校验 ----------
chk = df.set_index("gene").join(orig[["dominant_celltype", "assessable"]], how="left")
mismatch = chk[chk["dominant_original_10ct"] != chk["dominant_celltype"]]
print("\n交叉校验（原始 10 类口径 vs 既有表）:",
      "全部一致 ✓" if mismatch.empty else "不一致 %d 个:\n%s" % (len(mismatch), mismatch.to_string()))

# ---------- 4. 只在可评估基因上统计（与正文口径一致）----------
sub = chk[chk["assessable"].astype(bool)].copy()
sub["changed_designation"] = sub["dominant_original_10ct"] != sub["dominant_merged_cells"]
n_ass = len(sub)
n_micro = int((sub["dominant_original_10ct"] == "Microglia").sum())
n_mye = int(sub["is_myeloid_dominant_cells"].sum())
n_mye_eq = int(sub["is_myeloid_dominant_equalweight"].sum())
changed = sub.loc[sub["changed_designation"], ["dominant_original_10ct", "dominant_merged_cells"]]

print(f"\n===== 可评估候选 {n_ass} 个 =====")
print("  原始 10 类口径下优势细胞类型分布:")
print("   " + sub["dominant_original_10ct"].value_counts().to_string().replace("\n", "\n   "))

# --- 三口径对照：仅改标签 / 细胞级重算 / 等权重算 ---
label_only = int(sub["dominant_original_10ct"].isin(MYELOID).sum())
table = pd.DataFrame([
    {"口径": "仅改标签（原亚型归属改称髓系，不重算）", "髓系优势": label_only},
    {"口径": "合并后细胞级重算（按细胞数加权）", "髓系优势": n_mye},
    {"口径": "合并后亚型均值等权重算", "髓系优势": n_mye_eq},
])
print(f"\n===== 三种口径对照（分母 {n_ass}）=====")
print(table.to_string(index=False))
detail = sub[["dominant_original_10ct", "myeloid_mean_cells", "myeloid_mean_equalweight",
              "best_nonmyeloid_type", "best_nonmyeloid_mean",
              "is_myeloid_dominant_cells", "is_myeloid_dominant_equalweight"]]
print("\n  逐基因明细（可评估者）:")
print(detail.sort_values("myeloid_mean_cells", ascending=False).to_string())
n_flip = int((sub["is_myeloid_dominant_cells"] != sub["is_myeloid_dominant_equalweight"]).sum())
print(f"\n  两种重算口径的结论差异基因数: {n_flip}")
print(f"\n  原始小胶质优势                        : {n_micro} / {n_ass}")
print(f"  合并髓系（细胞级重算）后髓系优势      : {n_mye} / {n_ass}")
print(f"  合并髓系（亚型等权均值）后髓系优势    : {n_mye_eq} / {n_ass}")
print(f"\n  归属发生改变的基因（{len(changed)} 个）:")
print(changed.to_string())

summary = pd.DataFrame([
    {"item": "assessable_candidates", "value": n_ass},
    {"item": "microglia_dominant_original", "value": n_micro},
    {"item": "myeloid_dominant_after_merge_cells", "value": n_mye},
    {"item": "myeloid_dominant_after_merge_equalweight", "value": n_mye_eq},
    {"item": "genes_changing_designation", "value": int(sub["changed_designation"].sum())},
    {"item": "genes_changing_designation_list", "value": ";".join(changed.index.tolist())},
])
summary.to_csv(os.path.join(OUT, "myeloid_merge_summary.csv"), index=False)
print("\n输出：myeloid_merge_dominant_celltype.csv, myeloid_merge_summary.csv")
