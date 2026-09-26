"""
方案 C：29 个候选基因 × 6 篇竞品全文 的基因级交叉比对

输出：
  1. 每篇文献提及了哪些候选基因（= 独立复现候选）
  2. 每条提及的上下文句子（判断是"作为 DEG 报告"还是"泛泛讨论"）
  3. 交叉表矩阵
"""
import os
import re
import json
import pandas as pd

# --- 复现包路径（可移植化改写；原始脚本使用本机绝对路径）---
# 注意：本脚本需把竞品论文全文的纯文本缓存放在本目录下（因版权不随包分发）；
#       论文清单与提及情况见 results/C_gene_crosswalk_summary.csv。
BASE = os.path.dirname(os.path.abspath(__file__))
OUT = BASE

CANDIDATES = [
    "Reg3b", "C1qa", "Cfh", "Gapt", "Hexb", "Nrp1", "Laptm5", "Nlrc4", "Mmp3",
    "Ecel1", "Clec7a", "Arhgap25", "Mrpl10", "Fcrl2", "Cd68", "Lmo2", "Rfx7",
    "Zc3h12a", "C4a", "Cd84", "Man2b2", "Klhl6", "Fbxw4", "Fbxl7", "Dnase2",
    "AABR07034445.1", "Chrna1", "Svop", "Rassf3",
]

PAPERS = {
    "41008561": "Li 2025 Biomolecules 15:1254（同刊同题，含qRT-PCR）",
    "41184936": "He 2025 J Transl Med 23:1209（GSE175760 原始数据作者）",
    "41733593": "Cao 2026 Brain Behav（GSE175760 作验证集）",
    "42039142": "未标注 PMID 42039142",
    "42606631": "Metab Brain Dis 2026（概念撞车）",
    "42706380": "Sci Rep 2026 16:27944（小胶质表型时程）",
}


def split_sentences(text):
    """把长文本切成句子（保留位置）"""
    # 按句末标点 + 空格切分，保留较完整的语义单元
    parts = re.split(r"(?<=[.!?])\s+(?=[A-Z(])", text)
    return parts


def find_mentions(text, gene):
    """找出基因在文本中所有出现位置的上下文句子"""
    # 词边界：前后不能是字母或数字（避免 C4a 匹配到 C4ap 之类）
    pat = re.compile(r"(?<![A-Za-z0-9])" + re.escape(gene) + r"(?![A-Za-z0-9])", re.IGNORECASE)
    hits = []
    for m in pat.finditer(text):
        s = max(0, m.start() - 350)
        e = min(len(text), m.end() + 350)
        ctx = text[s:e].replace("\n", " ")
        hits.append(ctx)
    return hits


rows = []
detail = []

for pmid, label in PAPERS.items():
    fp = os.path.join(BASE, f"_fulltext_{pmid}.txt")
    if not os.path.exists(fp):
        print(f"[跳过] {pmid} 无全文文件")
        continue
    text = open(fp, encoding="utf-8", errors="ignore").read()
    print(f"\n{'='*70}\n{pmid}  {label}   全文 {len(text)} 字符")

    mentioned = []
    for g in CANDIDATES:
        hits = find_mentions(text, g)
        if hits:
            mentioned.append(g)
            for h in hits:
                detail.append({"paper": pmid, "gene": g, "context": h})
    print(f"  提及候选基因 {len(mentioned)}/{len(CANDIDATES)}: {mentioned}")
    rows.append({"paper": pmid, "label": label, "n_mentioned": len(mentioned),
                 "mentioned": "; ".join(mentioned)})

cross = pd.DataFrame(rows)
cross.to_csv(os.path.join(OUT, "C_gene_crosswalk_summary.csv"), index=False)

det = pd.DataFrame(detail)
det.to_csv(os.path.join(OUT, "C_gene_crosswalk_contexts.csv"), index=False)

# 基因 × 文献 矩阵
mat = pd.DataFrame(0, index=CANDIDATES, columns=list(PAPERS.keys()))
for d in detail:
    mat.loc[d["gene"], d["paper"]] = 1
mat["n_papers"] = mat.sum(axis=1)
mat = mat.sort_values("n_papers", ascending=False)
mat.to_csv(os.path.join(OUT, "C_gene_crosswalk_matrix.csv"))

print(f"\n{'='*70}\n===== 交叉矩阵（1 = 该文献提及该基因）=====")
print(mat.to_string())
print(f"\n出现上下文共 {len(det)} 条 -> C_gene_crosswalk_contexts.csv")
print("汇总 -> C_gene_crosswalk_summary.csv")
