# ============================================================
# v11 Step 3: 候选基因筛选（LASSO + SVM-RFE）
#   输入: output/v11/wgcna_degs_chronic_module_genes.csv (turquoise 渐进上调模块)
#   方法1 LASSO: 18 样本, 时间 0-5 为连续标签
#   方法2 SVM-RFE: 14d vs Sham (n=6)
#   输出: output/v11/candidates_*
# 用法: python scripts/v11_03_candidates.py
# ============================================================
import pandas as pd, numpy as np
from sklearn.linear_model import LassoCV
from sklearn.svm import SVC
from sklearn.feature_selection import RFECV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, KFold
import warnings; warnings.filterwarnings('ignore')

BASE = r"F:\scs research"
OUT  = BASE + r"\output"
V11  = OUT + r"\v11"

expr = pd.read_csv(f"{OUT}/GSE175760_expression_matrix.csv", index_col=0).dropna()
meta = pd.read_csv(f"{OUT}/GSE175760_metadata.csv")
log2expr = np.log2(expr + 0.1)

mod_genes = pd.read_csv(f"{V11}/DEG_exploratory_genes.csv")["gene"].tolist()
mod_genes = [g for g in mod_genes if g in expr.index]
print(f"候选筛选输入基因（探索集，不限定模块）: {len(mod_genes)}")

# 模块归属作为注释（非筛选条件）
gene_mod = pd.read_csv(f"{V11}/wgcna_degs_gene_modules.csv").set_index("gene")["module"]

time_map = {'Sham':0,'CCI_0.5d':1,'CCI_1d':2,'CCI_3d':3,'CCI_7d':4,'CCI_14d':5}
grp = meta.set_index('sample_id')['group']
X_all = log2expr.loc[mod_genes, grp.index].T.values
y_time = grp.map(time_map).values

# ---------- 方法 1: LASSO (连续时间) ----------
sc = StandardScaler(); Xs = sc.fit_transform(X_all)
lasso = LassoCV(cv=KFold(5, shuffle=True, random_state=42),
                random_state=42, max_iter=20000,
                alphas=np.logspace(-4, 1, 200)).fit(Xs, y_time)
coef = pd.Series(lasso.coef_, index=mod_genes)
nz = coef[coef.abs() > 1e-6].sort_values(key=abs, ascending=False)
print(f"\nLASSO: alpha={lasso.alpha_:.4f}, 非零系数基因 = {len(nz)}")

# ---------- 方法 2: SVM-RFE (14d vs Sham) ----------
sham = meta[meta['group']=='Sham']['sample_id'].tolist()
c14  = meta[meta['group']=='CCI_14d']['sample_id'].tolist()
X_cs = log2expr.loc[mod_genes, sham + c14].T.values
y_cs = np.array([0]*len(sham) + [1]*len(c14))
X_cs = StandardScaler().fit_transform(X_cs)
rfecv = RFECV(estimator=SVC(kernel='linear', random_state=42), step=1,
              cv=StratifiedKFold(n_splits=3, shuffle=True, random_state=42),
              scoring='accuracy', min_features_to_select=5).fit(X_cs, y_cs)
svm_sel = [mod_genes[i] for i in np.where(rfecv.support_)[0]]
print(f"SVM-RFE: 选中 = {len(svm_sel)}")

# ---------- 并集 ----------
hub = sorted(set(nz.index) | set(svm_sel))
print(f"\n并集候选基因 = {len(hub)}")

# ---------- 附上 limma 统计 ----------
lim = pd.read_csv(f"{V11}/limma_all_contrasts.csv")
lim = lim[lim['gene'].isin(hub)]
agg = lim.groupby('gene').agg(
    min_P_raw=('P_raw','min'), min_P_BH=('P_BH','min'),
    max_abs_log2FC=('log2FC', lambda s: s.abs().max()),
    n_tp_raw=('dir_explor', lambda s:(s!='ns').sum()),
    n_tp_bh=('dir_strict', lambda s:(s!='ns').sum())).reset_index()

cand = pd.DataFrame({'gene': hub})
cand['in_LASSO'] = cand['gene'].isin(nz.index)
cand['in_SVM_RFE'] = cand['gene'].isin(svm_sel)
cand['LASSO_coef'] = cand['gene'].map(coef)
cand['WGCNA_module'] = cand['gene'].map(gene_mod).fillna('not_assigned')
cand = cand.merge(agg, on='gene', how='left')

# ---- 第三重证据: 与时间的相关性（全 18 样本）----
from scipy import stats as sps
tc = []
for g in hub:
    v = log2expr.loc[g, grp.index].values
    r, p = sps.spearmanr(y_time, v)
    tc.append({'gene': g, 'rho_time': r, 'P_time_corr': p})
tc = pd.DataFrame(tc)
cand = cand.merge(tc, on='gene', how='left')

# ---- 主候选名单: ML 选中 且 至少1个时间点 BH<0.05（双重要求）----
primary = cand[cand['n_tp_bh'] > 0].copy()
primary = primary.sort_values('min_P_BH')

# 三重证据计数
primary['n_evidence'] = (primary['in_LASSO'].astype(int)
                         + primary['in_SVM_RFE'].astype(int)
                         + (primary['P_time_corr'] < 0.05).astype(int))
cand = cand.sort_values(['in_LASSO','in_SVM_RFE','min_P_BH'],
                        ascending=[False,False,True])
cand.to_csv(f"{V11}/candidates_ml_union_all.csv", index=False)
primary.to_csv(f"{V11}/candidates_primary.csv", index=False)
pd.DataFrame({'gene': primary['gene']}).to_csv(f"{V11}/candidates_genes.csv", index=False)
pd.DataFrame({'gene': nz.index, 'coef': nz.values}).to_csv(
    f"{V11}/candidates_lasso.csv", index=False)
pd.DataFrame({'gene': svm_sel}).to_csv(f"{V11}/candidates_svmrfe.csv", index=False)

print("\n=== 主候选基因名单 (ML 选中 且 至少1时间点 BH<0.05) ===")
show = primary.copy()
show['min_P_raw'] = show['min_P_raw'].map(lambda x: f"{x:.2e}")
show['min_P_BH']  = show['min_P_BH'].map(lambda x: f"{x:.2e}")
show['LASSO_coef'] = show['LASSO_coef'].round(4)
show['max_abs_log2FC'] = show['max_abs_log2FC'].round(3)
print(show[['gene','in_LASSO','in_SVM_RFE','WGCNA_module','LASSO_coef',
            'max_abs_log2FC','min_P_raw','min_P_BH','n_tp_raw','n_tp_bh',
            'rho_time','P_time_corr','n_evidence']].to_string(index=False))

print(f"\n  ML 并集(全部): {len(cand)}")
print(f"  主候选名单(含 BH 显著): {len(primary)}")
print(f"  其中 LASSO 与 SVM-RFE 共同选中: {int((primary['in_LASSO']&primary['in_SVM_RFE']).sum())}")
print(f"  其中与时间显著相关 (Spearman P<0.05): {int((primary['P_time_corr']<0.05).sum())}")
print(f"  三重证据齐备: {int((primary['n_evidence']==3).sum())}; 双重: {int((primary['n_evidence']==2).sum())}")
print(f"  主名单模块归属: {primary['WGCNA_module'].value_counts().to_dict()}")
print(f"\n输出: {V11}\\candidates_*.csv")
