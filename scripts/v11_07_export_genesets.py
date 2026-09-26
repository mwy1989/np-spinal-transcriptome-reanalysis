"""
v11 阶段 2：从 gsva_pilot_v2.py 提取 curated 通路基因集，导出为 GMT。
避免手工复制 28 个通路的大字典造成错误。
"""
import ast
import os

BASE = r"F:\scs research"
SRC = os.path.join(BASE, "SCS-PCR", "scripts", "gsva_pilot_v2.py")
OUT = os.path.join(BASE, "output", "v11")
os.makedirs(OUT, exist_ok=True)

tree = ast.parse(open(SRC, encoding="utf-8").read())

PATHWAYS = None
for node in tree.body:
    if isinstance(node, ast.Assign):
        for t in node.targets:
            if isinstance(t, ast.Name) and t.id == "PATHWAYS":
                PATHWAYS = ast.literal_eval(node.value)

assert PATHWAYS, "未找到 PATHWAYS 定义"
print(f"提取到 {len(PATHWAYS)} 个通路基因集")

gmt_path = os.path.join(OUT, "v11_pathways_custom.gmt")
with open(gmt_path, "w", encoding="utf-8") as f:
    for name, genes in PATHWAYS.items():
        # GMT: 名称 \t 描述 \t 基因...
        f.write(f"{name}\tcurated_v11\t" + "\t".join(dict.fromkeys(genes)) + "\n")
print("写出:", gmt_path)

# 同时输出一个通路 -> 基因数汇总，便于 Methods 描述
import csv
with open(os.path.join(OUT, "v11_pathway_gene_counts.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["pathway", "n_genes_defined"])
    for name, genes in PATHWAYS.items():
        w.writerow([name, len(set(genes))])
print("写出: v11_pathway_gene_counts.csv")
for name, genes in PATHWAYS.items():
    print(f"  {name:36s} {len(set(genes)):4d} genes")
