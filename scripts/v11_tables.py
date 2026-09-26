"""
v11 — Table 1 (cohorts) and Table 2 (candidate gene evidence)
输出：Figures_v11/Table_1_Cohorts.xlsx|.md, Table_2_Candidate_Evidence.xlsx|.md
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v11_figstyle import V11, FIGDIR

import os as _os
BASE = _os.environ.get("SCS_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
print("Tables ...")

# ============================================================
# Table 1 — Cohorts and the questions they can answer
# ============================================================
t1 = pd.DataFrame([
    {
        "Cohort (accession)": "Discovery\nGSE175760",
        "Species": "Rat",
        "Model / tissue": "Chronic constriction injury (CCI); lumbar spinal cord; bulk RNA-seq",
        "Groups and n": "6 groups x n = 3 (18 samples)\nSham, 0.5, 1, 3, 7, 14 d after CCI",
        "Stimulation": "None",
        "Role in this study": "Time-course differential expression; co-expression modules; candidate gene selection",
        "Principal limitation": "Single injury model; whole-tissue bulk data cannot resolve cell types; n = 3 per time point",
    },
    {
        "Cohort (accession)": "Cross-platform\nGSE5296",
        "Species": "Mouse (A/J, C57BL/6, ICR)",
        "Model / tissue": "Spinal cord injury; impact region; microarray (RMA, supplied on a linearised scale)",
        "Groups and n": "96 samples total\n0.5 h - 28 d; injured vs sham",
        "Stimulation": "None",
        "Role in this study": "Direction concordance of candidate genes across platform and species",
        "Principal limitation": "Three strains pooled (strain fitted as a covariate); only the impact region is available in the downloaded series; array platform differs from the discovery cohort",
    },
    {
        "Cohort (accession)": "Stimulation\nGSE243038",
        "Species": "Mouse",
        "Model / tissue": "Spinal cord motoneurons; Smart-seq2",
        "Groups and n": "6 groups x 3 samples\nSham, untrained, MS, SCS, ineffective DES, effective DES",
        "Stimulation": "Dual electrical stimulation at defined frequencies (each group: 10, 15 and 20 Hz; ineffective DES: 1, 30 and 40 Hz)",
        "Role in this study": "Descriptive stimulation-associated expression trends for the candidate genes",
        "Principal limitation": "Motoneuron-enriched rather than whole tissue; the three samples per group are three different stimulation frequencies, not biological replicates; the original study assessed volitional motor output, not pain behaviour",
    },
    {
        "Cohort (accession)": "Cell-type reference\nGSE189070",
        "Species": "Mouse",
        "Model / tissue": "Spinal cord; single-cell RNA-seq",
        "Groups and n": "10 annotated cell types\nMicroglia, macrophage, neutrophil, astrocyte, oligodendrocyte, ependymal, endothelial, pericyte, B cell, T cell",
        "Stimulation": "None",
        "Role in this study": "Descriptive expression background and dominant cell type for each candidate",
        "Principal limitation": "No neuronal cluster is annotated, so neuronal expression cannot be assessed; the reference contains no stimulation arm; single-cell data are not spatial",
    },
])

# 脚注：被排除的队列
t1_note = (
    "GSE155610 (rat, motor cortex stimulation after spinal cord injury, n = 2 vs 2) was reviewed and excluded: "
    "the stimulation target is the motor cortex rather than the spinal cord, and the group size does not support "
    "any meaningful comparison. All analyses were performed on the unified expression scale log2(FPKM + 0.1) with "
    "limma-trend; the GSE5296 matrix was analysed after restoring its log2 scale."
)

# ============================================================
# Table 2 — Candidate gene evidence summary
# ============================================================
cand = pd.read_csv(os.path.join(V11, "candidates_primary.csv"))
asg = pd.read_csv(os.path.join(V11, "cluster_k4_assignments.csv"))
det = pd.read_csv(os.path.join(V11, "gse5296_7d_candidate_detail.csv"))
scs = pd.read_csv(os.path.join(V11, "scs_candidate_expression.csv")).set_index("gene")
dom = pd.read_csv(os.path.join(V11, "candidate_dominant_celltype.csv")).set_index("gene")

CL = {1: "Early burst (0.5 d)", 2: "Progressive up (14 d)",
      3: "Mid-late peak (7 d)", 4: "Persistent down (Sham)"}
MOD = {"turquoise": "turquoise", "blue": "blue", "brown": "brown",
       "grey": "grey", "not_assigned": "unassigned"}

t2 = cand.merge(asg[["gene", "cluster"]], on="gene", how="left")
t2["Temporal cluster"] = t2["cluster"].map(CL).fillna("unassigned")

dt = det.set_index("gene")
sc = scs
rows = []
for _, r in t2.iterrows():
    g = r["gene"]
    src = []
    if bool(r["in_LASSO"]):
        src.append("LASSO")
    if bool(r["in_SVM_RFE"]):
        src.append("SVM-RFE")
    ext = "not detected"
    if g in dt.index:
        d = dt.loc[g]
        ext = "concordant" if bool(d["concordant"]) else "discordant"
    if g in sc.index:
        s = sc.loc[g]
        sa = "yes \u2014 trend toward Sham" if bool(s["SCS_10_15_20Hz_toward_sham"]) else "yes \u2014 no trend"
        if abs(s["SCI_shift"]) <= 0.1:
            sa = "yes \u2014 shift too small to evaluate"
    else:
        sa = "not detected"
    if g in dom.index and bool(dom.loc[g, "assessable"]):
        ct = dom.loc[g, "dominant_celltype"]
        rate = f"{dom.loc[g, 'max_detection_rate']*100:.0f} %"
    else:
        ct, rate = "not assessable", "< 10 %"
    rows.append({
        "Gene": g,
        "Selected by": " + ".join(src),
        "WGCNA module": MOD.get(r["WGCNA_module"], r["WGCNA_module"]),
        "Temporal cluster": r["Temporal cluster"],
        "Max abs. log2FC (discovery)": float(r["max_abs_log2FC"]),
        "Min raw P (discovery)": float(r["min_P_raw"]),
        "Min BH P (discovery)": float(r["min_P_BH"]),
        "Time points BH-significant": int(r["n_tp_bh"]),
        "GSE5296 direction (7 d)": ext,
        "GSE243038 evaluable": sa,
        "Dominant cell type (reference)": ct,
        "Detection rate": rate,
    })
t2 = pd.DataFrame(rows)
t2["Max abs. log2FC (discovery)"] = t2["Max abs. log2FC (discovery)"].round(2)

t2_note = (
    "BH P values are Benjamini-Hochberg adjusted across all genes at each time point (limma-trend). "
    "\"Time points BH-significant\" counts the number of time points at which a gene passed BH FDR < 0.05 with "
    "abs(log2FC) > 0.58. GSE243038 trends are descriptive only and carry no P values: the three samples per group "
    "are three different stimulation frequencies rather than biological replicates, and the dataset contains "
    "motoneurons only. \"Not assessable\" in the cell-type columns means the gene was not detected reliably in the "
    "single-cell reference (detection rate below 10 %); it does not mean that the gene is absent or unresponsive. "
    "AABR07034445.1 is an unannotated Rattus norvegicus transcript."
)

# ============================================================
# 输出
# ============================================================
def md_table(df, note=None, align_left=None):
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |",
           "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                if 0 < v < 0.001:
                    v = f"{v:.2e}"
                else:
                    v = f"{v:g}"
            cells.append(str(v).replace("\n", "<br>").replace("|", "\\|"))
        out.append("| " + " | ".join(cells) + " |")
    s = "\n".join(out)
    if note:
        s += "\n\n" + note
    return s


with pd.ExcelWriter(os.path.join(FIGDIR, "Table_1_Cohorts.xlsx"), engine="openpyxl") as w:
    t1.to_excel(w, sheet_name="Table 1", index=False)
    pd.DataFrame({"Note": [t1_note]}).to_excel(w, sheet_name="Note", index=False)

with pd.ExcelWriter(os.path.join(FIGDIR, "Table_2_Candidate_Evidence.xlsx"), engine="openpyxl") as w:
    t2.to_excel(w, sheet_name="Table 2", index=False)
    pd.DataFrame({"Note": [t2_note]}).to_excel(w, sheet_name="Note", index=False)

with open(os.path.join(FIGDIR, "Table_1_Cohorts.md"), "w", encoding="utf-8") as f:
    f.write("# Table 1. Cohorts and the questions they can answer\n\n")
    f.write(md_table(t1))
    f.write("\n\n**Note.** " + t1_note + "\n")

with open(os.path.join(FIGDIR, "Table_2_Candidate_Evidence.md"), "w", encoding="utf-8") as f:
    f.write("# Table 2. Evidence summary for the 29 candidate genes\n\n")
    f.write(md_table(t2))
    f.write("\n\n**Note.** " + t2_note + "\n")

print(f"  Table 1: {t1.shape}")
print(f"  Table 2: {t2.shape}")
print(f"  BH-significant >=1 time point: {(t2['Time points BH-significant']>0).sum()}/{len(t2)}")
print("done.")
