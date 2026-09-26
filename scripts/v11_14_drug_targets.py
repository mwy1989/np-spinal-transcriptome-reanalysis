"""
v11 阶段 2 / 任务 37：候选基因与现有药物靶点的交叉比对（补充材料级）

相对 v10 的变化：
  1. 输入改为 v11 新候选基因（29 个）与新 WGCNA 模块
  2. 明确删除"协同治疗/联合用药"推断——本分析仅为描述性靶点比对
  3. 增加补体/免疫调节药物一组（因候选基因富集于补体级联）
  4. 输出到 output/v11/
"""
import os
import pandas as pd

BASE = r"F:\scs research"
OUT = os.path.join(BASE, "output", "v11")

cand = pd.read_csv(os.path.join(OUT, "candidates_primary.csv"))
genes = cand["gene"].tolist()
gm = pd.read_csv(os.path.join(OUT, "wgcna_degs_gene_modules.csv"))
mods = {m: set(gm["gene"][gm["module"] == m]) for m in ["turquoise", "blue", "brown"]}

# ---------- 药物靶点库 ----------
NP_DRUGS = {
    "Gabapentin": ["Cacna2d1", "Cacna2d2"],
    "Pregabalin": ["Cacna2d1"],
    "Duloxetine": ["Slc6a4", "Slc6a2"],
    "Amitriptyline": ["Slc6a4", "Hrh1", "Chrm1", "Adra1a"],
    "Tramadol": ["Oprm1", "Slc6a4"],
    "Morphine": ["Oprm1", "Oprk1", "Oprd1"],
    "Oxycodone": ["Oprm1", "Oprk1"],
    "Carbamazepine": ["Scn2a", "Scn3a"],
    "Baclofen": ["Gabbr1", "Gabbr2"],
    "NSAIDs": ["Ptgs1", "Ptgs2"],
    "Ketamine": ["Grin1", "Grin2a", "Grin2b"],
    "Ziconotide": ["Cacna1b"],
    "Capsaicin": ["Trpv1"],
    "Lidocaine": ["Scn9a", "Scn10a"],
}
IMMUNE_DRUGS = {
    "Tofacitinib": ["Jak1", "Jak3"],
    "Ruxolitinib": ["Jak1", "Jak2"],
    "Baricitinib": ["Jak1", "Jak2"],
    "Eculizumab": ["C5"],
    "Pegcetacoplan": ["C3"],
    "Avacopan": ["C5ar1"],
    "Anakinra": ["Il1r1"],
    "Canakinumab": ["Il1b"],
    "Tocilizumab": ["Il6r"],
    "Sarilumab": ["Il6r"],
    "Minocycline": ["Casp1", "Mmp9"],
    "Colchicine": ["Nlrp3", "Tubb1"],
    "Dapansutrile": ["Nlrp3"],
    "Maraviroc": ["Ccr5"],
}

def cross(drugs, label):
    rows = []
    for drug, targets in drugs.items():
        for t in targets:
            hit = [g for g in genes if g.upper() == t.upper()]
            rows.append({"drug_class": label, "drug": drug, "target": t,
                         "candidate_hit": hit[0] if hit else None,
                         "match": bool(hit)})
    return pd.DataFrame(rows)

df = pd.concat([cross(NP_DRUGS, "conventional_NP"), cross(IMMUNE_DRUGS, "immune_modulator")],
               ignore_index=True)
df.to_csv(os.path.join(OUT, "drug_target_cross.csv"), index=False)

all_targets = set(t.upper() for t in list(NP_DRUGS.values()) + list(IMMUNE_DRUGS.values())[0]
                  for t in ([t] if isinstance(t, str) else t))
all_targets = set()
for d in (NP_DRUGS, IMMUNE_DRUGS):
    for tlist in d.values():
        all_targets.update(x.upper() for x in tlist)

print("=" * 66)
print("候选基因 × 现有药物靶点 交叉比对（描述性）")
print("=" * 66)
print(f"候选基因数: {len(genes)}")
print(f"药物靶点数: {len(all_targets)}（{len(NP_DRUGS)} 种常规 NP 药 + {len(IMMUNE_DRUGS)} 种免疫调节药）")

hits = df[df["match"]]
print(f"\n直接命中: {len(hits)} 条")
if len(hits):
    print(hits[["drug_class", "drug", "target", "candidate_hit"]].to_string(index=False))
else:
    print("  无直接重叠")

# 模块层面（放宽到模块，因模块含数百基因）
print("\n===== WGCNA 模块层面的靶点命中数 =====")
for m, gs in mods.items():
    h = sorted(all_targets & set(x.upper() for x in gs))
    print(f"  {m:10s} (n={len(gs):4d}): 命中 {len(h)} 个靶点 -> {h if h else '无'}")

# 结论声明
print("\n" + "=" * 66)
print("【必须写入文章的限制声明】")
print("  本分析仅为描述性靶点比对。候选基因与现有药物靶点的重叠（或缺乏重叠）")
print("  不能用于推断联合治疗、协同效应或临床转化潜力；无药理学实验验证。")
print("=" * 66)

df.to_csv(os.path.join(OUT, "drug_target_cross.csv"), index=False)
print("\n输出: drug_target_cross.csv")
