"""
v11 阶段 2：下载标准基因集（MSigDB Hallmark 2020 + KEGG）并检查大鼠映射覆盖率。
产出 GMT 供 R GSVA 使用。
"""
import os
import numpy as np
import pandas as pd
import gseapy as gp

import os as _os
BASE = _os.environ.get("SCS_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
OUT = os.path.join(BASE, "output", "v11")
os.makedirs(OUT, exist_ok=True)

expr = pd.read_csv(os.path.join(BASE, "output", "GSE175760_expression_matrix.csv"),
                   index_col=0)
rat_upper = {g.upper(): g for g in expr.index}
print(f"大鼠表达矩阵基因数: {len(expr.index)}")

for lib, tag in [("MSigDB_Hallmark_2020", "hallmark"),
                 ("KEGG_2019_Mouse", "kegg")]:
    print(f"\n===== {lib} =====")
    try:
        gmt = gp.get_library(name=lib, organism="Mouse")
    except Exception as e:
        print("  下载失败:", type(e).__name__, e)
        continue

    if gmt is None:
        print("  未返回数据")
        continue

    rows = []
    hits = {}
    for name, genes in gmt.items():
        genes = list(dict.fromkeys(genes))
        found = [rat_upper[g.upper()] for g in genes if g.upper() in rat_upper]
        hits[name] = found
        rows.append({"geneset": name, "n_genes": len(genes),
                     "n_matched_rat": len(found),
                     "coverage": round(len(found) / max(len(genes), 1), 3)})
    cov = pd.DataFrame(rows).sort_values("coverage", ascending=False)
    print(f"  基因集数: {len(cov)}")
    print(f"  覆盖率: median={cov.coverage.median():.3f}  min={cov.coverage.min():.3f}  max={cov.coverage.max():.3f}")
    print(f"  覆盖率>=0.5 的集数: {(cov.coverage>=0.5).sum()}")

    # 只保留覆盖 >= 30%（GSVA 需要足够基因）
    keep = {k: v for k, v in hits.items() if len(v) >= 10}
    print(f"  保留（>=10 个匹配基因）: {len(keep)} 个")

    gmt_path = os.path.join(OUT, f"v11_{tag}_rat.gmt")
    with open(gmt_path, "w", encoding="utf-8") as f:
        for name, genes in keep.items():
            f.write(f"{name}\t{lib}\t" + "\t".join(genes) + "\n")
    print("  写出:", gmt_path)
    cov.to_csv(os.path.join(OUT, f"v11_{tag}_coverage.csv"), index=False)
    print("\n  覆盖率最高的 10 个:")
    print(cov.head(10).to_string(index=False))
