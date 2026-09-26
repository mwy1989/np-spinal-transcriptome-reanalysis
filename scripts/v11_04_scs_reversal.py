# ============================================================
# v11 Step 4: SCS 部分重建（GSE243038）
#   事实认定（据 SOFT 元数据）:
#     - 组织: 脊髓运动神经元 (Smart-seq2)，非全脊髓、非胶质
#     - 设计: Sham / Untrained(SCI) / MS / SCS / IneffectiveDES / EffectiveDES
#     - 每"组"3 个样本实为 3 个不同刺激频率 (10/15/20Hz 或 1/30/40Hz)
#     - 原研究目的: 运动功能恢复，非疼痛
#   ⇒ 本分析仅作"跨模型电刺激相关表达趋势"的描述性报告，不报 P 值、
#      不推断镇痛机制、不合并不同频率为生物学重复
# 用法: python scripts/v11_04_scs_reversal.py
# ============================================================
import gzip, numpy as np, pandas as pd

BASE = r"F:\scs research"
OUT  = BASE + r"\output"
V11  = OUT + r"\v11"

CAND = pd.read_csv(f"{V11}/candidates_genes.csv")["gene"].tolist()
print(f"候选基因数: {len(CAND)}")

# ---- 读取 GSE243038 FPKM ----
with gzip.open(f"{BASE}/GSE243038_DES_fpkm.anno.xls.gz", 'rt',
               encoding='utf-8', errors='ignore') as f:
    df = pd.read_csv(f, sep='\t', low_memory=False)

meta_cols = ['GeneName','Biotype','Position']
sample_cols = [c for c in df.columns if c not in meta_cols
               and not any(c.startswith(p) for p in
               ('NR:','NT:','Uniprot:','COG:','Pfam:','GO:','KEGG:'))]
print(f"样本列 ({len(sample_cols)}): {sample_cols}")

expr = df[sample_cols].apply(pd.to_numeric, errors='coerce')
expr.index = df['GeneName'].astype(str)
# 同名基因取均值
expr = expr.groupby(level=0).mean()
expr_log = np.log2(expr + 0.1)
print(f"矩阵: {expr.shape[0]} 基因 x {expr.shape[1]} 样本\n")

# ---- 跨物种匹配 (rat -> mouse) ----
rat2mouse = {g.upper(): g for g in expr.index}
found, missing = {}, []
for g in CAND:
    key = g.upper()
    if key in rat2mouse:
        found[g] = rat2mouse[key]
    else:
        # 常见格式差异: 首字母大写其余小写
        alt = g[0].upper() + g[1:].lower()
        if alt in expr.index:
            found[g] = alt
        else:
            missing.append(g)
print(f"跨物种匹配成功: {len(found)}/{len(CAND)}")
print(f"未匹配: {missing}\n")

GROUPS = {
    'Sham':            ['Sham1','Sham2','Sham3'],
    'SCI_untrained':   ['Untrained1','Untrained2','Untrained3'],
    'MS_10_15_20Hz':   ['MS1','MS2','MS3'],
    'SCS_10_15_20Hz':  ['SCS1','SCS2','SCS3'],
    'DES_ineff_1_30_40Hz': ['IneffectiveDES1','IneffectiveDES2','IneffectiveDES3'],
    'DES_eff_10_15_20Hz':  ['EffectiveDES1','EffectiveDES2','EffectiveDES3'],
}

rows = []
for g, mg in found.items():
    r = {'gene': g, 'mouse_symbol': mg}
    for gname, cols in GROUPS.items():
        vals = expr_log.loc[mg, cols].values.astype(float)
        r[f'{gname}_mean'] = np.nanmean(vals)
        for i, c in enumerate(cols, 1):
            r[f'{gname}_{c}'] = vals[i-1]
    sham = r['Sham_mean']; sci = r['SCI_untrained_mean']
    r['SCI_shift'] = sci - sham
    for gname in ['SCS_10_15_20Hz','DES_eff_10_15_20Hz','MS_10_15_20Hz']:
        stim = r[f'{gname}_mean']
        r[f'{gname}_delta_vs_SCI'] = stim - sci
        # 回归判据: 与 SCI 偏移同号方向相反(即朝 Sham 移动) 且 幅度更小
        d_sci = r['SCI_shift']; d_stim = stim - sham
        same_side = np.sign(d_stim) == np.sign(d_sci)
        smaller = abs(d_stim) < abs(d_sci)
        r[f'{gname}_toward_sham'] = bool(same_side and smaller and abs(d_sci) > 0.1)
    rows.append(r)

res = pd.DataFrame(rows)
res.to_csv(f"{V11}/scs_candidate_expression.csv", index=False)

# ============================================================
show = res[['gene','SCI_shift','SCS_10_15_20Hz_delta_vs_SCI',
            'SCS_10_15_20Hz_toward_sham','DES_eff_10_15_20Hz_toward_sham']].copy()
show = show.sort_values('SCI_shift', key=abs, ascending=False)
print("=== 候选基因在 GSE243038 中的表达变化（log2 FPKM 尺度）===")
print(f"{'gene':>10} {'SCI-Sham':>9} {'SCS-SCI':>9} {'朝Sham(SCS)':>12} {'朝Sham(DES_eff)':>16}")
for _, r in show.iterrows():
    print(f"{r['gene']:>10} {r['SCI_shift']:>+9.3f} "
          f"{r['SCS_10_15_20Hz_delta_vs_SCI']:>+9.3f} "
          f"{('是' if r['SCS_10_15_20Hz_toward_sham'] else '-'):>12} "
          f"{('是' if r['DES_eff_10_15_20Hz_toward_sham'] else '-'):>16}")

n_scs = int(res['SCS_10_15_20Hz_toward_sham'].sum())
n_des = int(res['DES_eff_10_15_20Hz_toward_sham'].sum())
n_eval = int((res['SCI_shift'].abs() > 0.1).sum())
print(f"\n可评估基因(|SCI-Sham|>0.1): {n_eval}/{len(res)}")
print(f"SCS 组呈向 Sham 回归趋势: {n_scs}")
print(f"DES(最优) 组呈向 Sham 回归趋势: {n_des}")
print(f"\n输出: {V11}\\scs_candidate_expression.csv")
