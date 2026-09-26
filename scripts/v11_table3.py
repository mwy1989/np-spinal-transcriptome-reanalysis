"""
Table 3 — Gene-level cross-reference of the 29 candidate genes against six competing studies
输出：
  Figures_v11/Table_3_Gene_CrossReference.xlsx
  Figures_v11/Table_3_Gene_CrossReference.md          （独立版，表头用 PMID）
  Figures_v11/Table_3_Gene_CrossReference_embed.md    （稿件嵌入版，表头用 [@key] 引用）
输入：output/v11/candidates_primary.csv
      output/v11/cluster_k4_assignments.csv
      outputs/novelty_v11/C_gene_crosswalk_matrix.csv
"""
import os, sys
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v11_figstyle import V11, FIGDIR

import os as _os
BASE = _os.environ.get("SCS_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
NOV = os.path.join(BASE, "outputs", "novelty_v11")
print("Table 3 ...")

cand = pd.read_csv(os.path.join(V11, "candidates_primary.csv"))
asg = pd.read_csv(os.path.join(V11, "cluster_k4_assignments.csv"))
mat = pd.read_csv(os.path.join(NOV, "C_gene_crosswalk_matrix.csv"), index_col=0)

CL = {1: "Early burst (0.5 d)", 2: "Progressive up (14 d)",
      3: "Mid-late peak (7 d)", 4: "Persistent down (sham)"}
clu = asg.set_index("gene")["cluster"].to_dict()

# ---- 六篇竞品（列顺序 = 相关度） --------------------------------------------
# pmid -> (稿件引用 key, 短标签, 独立版表头)
PAPERS = [
    ("41184936", "he2025",   "He 2025 (GSE175760 originators)"),
    ("41008561", "li2025",   "Li 2025 (same journal)"),
    ("41733593", "cao2026",  "Cao 2026"),
    ("42706380", "qian2026", "Sci Rep 2026"),
    ("42039142", "zhou2026", "Zhou 2026 (review)"),
    ("42606631", "park2026", "Metab Brain Dis 2026"),
]

# ---- 命中语境（依据 C_gene_crosswalk_contexts.csv 逐条人工归纳） --------------
HIT = {
    ("41184936", "C1qa"):  "day-7 top-20 up",
    ("41184936", "C4a"):   "day-7 top-20 up",
    ("41184936", "Cd68"):  "day-14 lysosomal gene",
    ("41184936", "Reg3b"): "PPI network (20 genes)",
    ("41008561", "C1qa"):  "top-10 PPI hub",
    ("42706380", "C1qa"):  "up at PID3",
    ("42039142", "C4a"):   "concept (C4A-C3 axis)",
    ("42039142", "Cd68"):  "concept (CD68+ lysosome)",
}

rows = []
order = list(cand["gene"])
hit_genes = sorted({g for (_, g) in HIT}, key=lambda g: (-sum(g == g2 for (_, g2) in HIT), g))
ordered = hit_genes + [g for g in order if g not in hit_genes]
for g in ordered:
    r = {"Candidate gene": g,
         "Temporal cluster": CL.get(clu.get(g), "unassigned")}
    n = 0
    for pmid, key, _ in PAPERS:
        code = HIT.get((pmid, g), "\u2014")
        if code != "\u2014":
            n += 1
        r[pmid] = code
    r["Studies reporting"] = n
    rows.append(r)
t3 = pd.DataFrame(rows)

NOTES = {
    "41184936": "He 2025 (PMID 41184936) is the group that generated GSE175760.",
    "41008561": "Li 2025 (PMID 41008561) is a same-journal study on an independent SNI cohort.",
    "41733593": "Cao 2026 (PMID 41733593) used GSE175760 as a validation set.",
    "42706380": "Sci Rep 2026 (PMID 42706380) reports microglial temporal dynamics in five unrelated datasets.",
    "42039142": "Zhou 2026 (PMID 42039142) is a review; \u201cconcept\u201d marks a gene that appears in its conceptual framework without new data.",
    "42606631": "Park 2026 (PMID 42606631) is a metabolomic study; RNA-seq datasets do not overlap with the present work.",
}

t3_note = (
    "Cells give the context in which the gene was reported; \u201c\u2014\u201d means the gene was not mentioned in the "
    "retrieved full text of that study and is not evidence of absence. Matching was case-insensitive and "
    "word-bounded against the gene symbol, covering rat, mouse and human nomenclature. Only the body text of the six "
    "studies was searched; their supplementary gene tables were not retrieved, so a gene reported exclusively in a "
    "supplementary file would be missed. " + " ".join(NOTES[p] for p, _, _ in PAPERS)
)

# ---- 输出 -------------------------------------------------------------------
with pd.ExcelWriter(os.path.join(FIGDIR, "Table_3_Gene_CrossReference.xlsx"),
                    engine="openpyxl") as w:
    out = t3.rename(columns={p: l for p, _, l in PAPERS})
    out.to_excel(w, sheet_name="Table 3", index=False)
    pd.DataFrame({"Note": [t3_note]}).to_excel(w, sheet_name="Note", index=False)


def md(df, note):
    s = ["| " + " | ".join(df.columns) + " |",
         "|" + "|".join(["---"] * len(df.columns)) + "|"]
    for _, r in df.iterrows():
        s.append("| " + " | ".join(str(r[c]).replace("|", "\\|") for c in df.columns) + " |")
    return "\n".join(s) + "\n\n**Note.** " + note + "\n"


with open(os.path.join(FIGDIR, "Table_3_Gene_CrossReference.md"), "w", encoding="utf-8") as f:
    f.write("# Table 3. Gene-level cross-reference of the 29 candidate genes against six competing studies\n\n")
    f.write(md(t3.rename(columns={p: l for p, _, l in PAPERS}), t3_note))

# 稿件嵌入版：表头用 [@key] 引用，便于 build_v11_manuscript.py 自动编号
emb = t3.rename(columns={p: f"{l.split(' (')[0]} [@{k}]" + (f" ({l.split(' (')[1]}"
                    if " (" in l else "")
                    for p, k, l in PAPERS})
with open(os.path.join(FIGDIR, "Table_3_Gene_CrossReference_embed.md"), "w", encoding="utf-8") as f:
    f.write(md(emb, t3_note))

print(t3.to_string(index=False))
print(f"\n  rows={len(t3)}  hit genes={len(hit_genes)}  total hits={int(t3['Studies reporting'].sum())}")
print("done.")
