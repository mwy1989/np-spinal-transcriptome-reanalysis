# Reproducibility package — Temporal spinal transcriptomic reprogramming and stimulation-associated candidate genes in neuropathic pain

This package accompanies the manuscript:

> **Temporal spinal transcriptomic reprogramming and stimulation-associated candidate genes in neuropathic pain: a criterion-explicit re-analysis framework across four public datasets**
> Wen-Yi Mao, MD. Department of Neurosurgery, Nantong Second People's Hospital, Nantong, Jiangsu, China.

It contains every analysis script and every intermediate result table from which each number, figure and table in the manuscript is derived. All analyses use **only public data**; no new data were generated and no ethics approval was required.

---

## 1. Data sources

| Accession | Species | Assay / tissue | Role in this study |
|---|---|---|---|
| **GSE175760** | Rat | Bulk RNA-seq, lumbar spinal cord, CCI; 6 groups × n = 3 (Sham, 0.5, 1, 3, 7, 14 d) | **Discovery cohort**: time-course differential expression, co-expression network, candidate gene selection |
| **GSE5296** | Mouse (A/J, C57BL/6, ICR) | Microarray (Affymetrix Mouse Genome 430 2.0), spinal cord, contusion; 96 samples, 0.5 h – 28 d | Cross-platform, cross-species **direction concordance** |
| **GSE243038** | Mouse | Smart-seq2, spinal motoneurons; 6 conditions × 3 samples | **Descriptive** stimulation-associated expression trends |
| **GSE189070** | Mouse | Single-cell RNA-seq, spinal cord; 10 annotated cell types | **Cell-type expression reference** |

All four series were downloaded from the NCBI Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/) and placed under `data/` (not redistributed here). Download them with:

```bash
mkdir -p data && cd data
for g in GSE175760 GSE5296 GSE243038 GSE189070; do
  wget -r -np -nH --cut-dirs=4 -R "index.html*" \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/${g%???}nnn/$g/"
done
```

GSE155610 was reviewed and **excluded** from the study: its stimulation target is the motor cortex rather than the spinal cord, and the group size (n = 2 vs 2) does not support any meaningful comparison.

---

## 2. Environment

### R (analyses, clustering, network, GSVA)
- R 4.3.3
- limma 3.58.1, WGCNA 1.74, GSVA 1.50.5, e1071, clusterProfiler dependencies
- `requirements_R.txt`

```r
install.packages(c("BiocManager", "WGCNA", "e1071"))
BiocManager::install(c("limma", "GSVA", "GSEABase"))
```

### Python (machine learning, enrichment, plotting)
- Python 3.11
- scikit-learn, pandas, numpy, matplotlib, gseapy, anndata, decoupler 2.2.0
- `requirements_python.txt`

Figures 1–5 in this package were rendered with matplotlib 3.10.9 and Figure 6 with 3.11.0. Because figure geometry is computed from a tight bounding box, another matplotlib version may shift the canvas by a pixel or two; no plotted value, colour or label is affected. Pin `matplotlib==3.10.9` if byte-identical output is required.

```bash
pip install -r requirements_python.txt
```

### Path configuration
Scripts resolve their working root from the environment variable `SCS_ROOT`. If it is not set, the root defaults to the **parent directory of `scripts/`**, so **no path editing is required**:

```bash
export SCS_ROOT=/path/to/your/copy        # optional; overrides the default
```

Under that root the pipeline expects the layout below. The provided package already contains every intermediate table, so the layout can be staged in three commands:

| Path under `$SCS_ROOT` | Contents | Where it comes from |
|---|---|---|
| `output/v11/` | Intermediate analysis tables | the shipped `results/*.csv` |
| `outputs/novelty_v11/` | Criterion-sensitivity and cross-reference tables | the shipped `results/B_*.csv` and `results/C_*.csv` |
| `Figures_v11/` | Rendered figures | the shipped `figures/` |
| `output/GSE175760_expression_matrix.csv`, `output/GSE175760_metadata.csv` | Discovery-cohort expression matrix and sample sheet | **not redistributed** — obtain from GEO (section 1) |
| Single-cell object for `v11_13_celltype_annotation.py` | GSE189070 annotated cell-by-gene matrix | **not redistributed** — obtain from GEO (section 1) |

```bash
mkdir -p output/v11 outputs/novelty_v11 Figures_v11
cp results/*.csv                    output/v11/
cp results/B_*.csv results/C_*.csv  outputs/novelty_v11/
cp figures/*                        Figures_v11/
```

Because every intermediate table ships with the package, **each reported number can be verified without re-running the pipeline**; re-execution is only needed to regenerate a table from scratch. The mapping from each reported value to its source file is given in section 5, and a machine-readable version of the script-to-output mapping is in `MANIFEST.tsv`.

---

## 3. Analysis pipeline — execution order

### Stage A — Differential expression (discovery cohort)
| # | Script | Produces |
|---|---|---|
| A1 | `limma_standard_v11.R` | `limma_standard_all.csv`, `limma_all_contrasts.csv` — limma-trend (`eBayes(trend=TRUE)`) per time point on the unified scale log2(FPKM + 0.1) |
| A2 | `v11_01_deg_finalize.R` | `DEG_strict_genes.csv` (386), `DEG_exploratory_genes.csv` (1,940), `DEG_summary_by_timepoint.csv` |

**Unified scale.** Every analysis in the manuscript uses `log2(FPKM + 0.1)`. This single choice is applied to differential expression, clustering, network construction, pathway scoring and every table and figure. No other scale, pseudo-count or transformation is used anywhere.

**Dual thresholds.** Strict set = BH FDR < 0.05 and |log2FC| > 0.58; exploratory set = nominal P < 0.05 and |log2FC| > 0.58.

### Stage B — Temporal structure and co-expression
| # | Script | Produces |
|---|---|---|
| B1 | `v11_05_temporal_clusters.R` | `cluster_k4_*.csv` (K = 4), `cluster_k6_*.csv` (K = 6) — fuzzy c-means on z-standardised group means |
| B2 | `v11_02_wgcna.R`, `v11_02b_wgcna_degs.R` | `wgcna_degs_*` — signed network, power = 14 (R² = 0.805), `blockwiseModules` |

Note: the module–cluster cross-tabulation reported in Figure 2B compares frameworks built on **different information** (trajectory shape vs. co-expression topology). It is reported as **convergence**, not as independent replication.

### Stage C — Candidate genes
| # | Script | Produces |
|---|---|---|
| C1 | `v11_03_candidates.py` | `candidates_lasso.csv`, `candidates_svmrfe.csv`, `candidates_ml_union_all.csv` (62), `candidates_primary.csv` (29) |

Primary list = selected by ≥1 machine-learning method **and** passing BH FDR < 0.05 with |log2FC| > 0.58 at ≥1 time point.

### Stage D — Pathway activity
| # | Script | Produces |
|---|---|---|
| D1 | `v11_07_export_genesets.py`, `v11_08_fetch_genesets.py` | `v11_hallmark_rat.gmt` (orthologue-mapped), coverage tables |
| D2 | `v11_09_gsva.R` | `gsva_hallmark_scores.csv`, `gsva_hallmark_delta.csv` — GSVA 1.50.5 with `gsvaParam`, Gaussian kernel. **Pathways are always indexed by name, never by position.** |
| D3 | `v11_10_progeny.py` | `progeny_scores.csv`, `progeny_delta_limma.csv` — PROGENy via decoupleR 2.2.0 (`dc.mt.mlm`, tmin = 5) |
| D4 | `v11_11_pathway_consolidate.R` | `pathway_summary_hallmark.csv`, `pathway_cross_method_consistency.csv` (9/12 agree at 7 d) |
| D5 | `v11_12_enrichment.py` | `enrichment_candidates.csv`, `enrichment_module_turquoise.csv`, `enrichment_module_blue.csv` |

### Stage E — External datasets
| # | Script | Produces |
|---|---|---|
| E1 | `v11_06_gse5296.R` | `gse5296_limma_all.csv`, `gse5296_candidate_concordance.csv`, `gse5296_module_overrepresentation.csv`. The distributed matrix is a **linearised back-transformation of RMA log2 values** and is log2-transformed before any computation (median after log2 = 6.08, range −0.69 to 11.87); model is `~ strain + condition` |
| E2 | `v11_04_scs_reversal.py` | `scs_candidate_expression.csv`. **No P value, confidence interval or restoration percentage is computed** — the three samples per condition are three stimulation frequencies, not replicates |
| E3 | `v11_13_celltype_annotation.py` | `candidate_dominant_celltype.csv`, `candidate_celltype_detection_rate.csv`, `fig5_celltype_heatmap_matrix.csv` |
| E4 | `v11_14_drug_targets.py` | `drug_target_cross.csv` (supplementary only) |

### Stage F — Figures and tables
| # | Script | Produces |
|---|---|---|
| F1 | `v11_figstyle.py` | Shared style (soft red `#E8A9A4`, soft green `#7DB89A`, soft blue `#5BC0DE`, amber `#F0AD4E`; text `#4A4540`) |
| F2 | `v11_fig1_design.py` | Figure 1 |
| F3 | `v11_fig2_temporal.py` | Figure 2 |
| F4 | `v11_fig3_candidates.py` | Figure 3 |
| F5 | `v11_fig4_scs_trends.py` | Figure 4 |
| F6 | `v11_fig5_celltype.py` | Figure 5 |
| F7 | `v11_fig6_threshold.py` | Figure 6 |
| F8 | `v11_tables.py` | Table 1, Table 2 |
| F9 | `v11_table3.py` | Table 3 |

---

## 4. Figure / table → script manifest

| Manuscript item | Producing script | Primary result table |
|---|---|---|
| Figure 1 | `v11_fig1_design.py` | `DEG_summary_by_timepoint.csv` |
| Figure 2 | `v11_fig2_temporal.py` | `cluster_k4_profiles.csv`, `wgcna_degs_module_trait.csv`, `gsva_hallmark_delta.csv` |
| Figure 3 | `v11_fig3_candidates.py` | `candidates_primary.csv`, `gse5296_candidate_concordance.csv` |
| Figure 4 | `v11_fig4_scs_trends.py` | `scs_candidate_expression.csv` |
| Figure 5 | `v11_fig5_celltype.py` | `candidate_dominant_celltype.csv` |
| Figure 6 | `v11_fig6_threshold.py` | `B_threshold_spectrum.csv`, `B_jaccard_matrix.csv`, `B_candidate_robustness.csv` |
| Table 1 | `v11_tables.py` | (cohort metadata) |
| Table 2 | `v11_tables.py` | `candidates_primary.csv`, `candidate_dominant_celltype.csv` |
| Table 3 | `v11_table3.py` | `C_gene_crosswalk_matrix.csv`, `C_gene_crosswalk_summary.csv` |

A machine-readable version of this mapping is in `MANIFEST.tsv`.

---

## 5. Key reported numbers and where they come from

| Reported value | Source file |
|---|---|
| 386 strict-set DEGs; 239 / 273 / 183 at 3 / 7 / 14 d; 0 at 0.5 and 1 d | `DEG_summary_by_timepoint.csv`, `DEG_strict_genes.csv` |
| 1,940 exploratory-set DEGs | `DEG_exploratory_genes.csv` |
| WGCNA turquoise n = 569, r = 0.852; blue n = 238, r = −0.815; power = 14, R² = 0.805 | `wgcna_degs_module_trait.csv`, `wgcna_degs_softthreshold.csv` |
| Cluster sizes 307 / 299 / 828 / 506 | `cluster_k4_assignments.csv` |
| 62 machine-learning union → 29 candidates | `candidates_ml_union_all.csv`, `candidates_primary.csv` |
| GSVA significant pathways 0 / 0 / 5 / 15 / 6 | `pathway_summary_hallmark.csv` |
| PROGENy JAK-STAT +6.746 at 7 d | `progeny_delta_limma.csv` |
| GSE5296 concordance 20/23 at 7 d, ρ = 0.649 | `gse5296_candidate_concordance.csv` |
| 8 of 22 candidates with stimulation-associated trends | `scs_candidate_expression.csv` |
| 10 of 16 assessable candidates microglia-dominant | `candidate_dominant_celltype.csv` |
| DEG yield spans 301–4,047 genes (13.4-fold) across ten analysis criteria | `B_threshold_spectrum.csv` |
| Published criterion applied to identical contrasts: 861 / 708 vs 529 / 352 reported | `B_threshold_spectrum.csv` |
| Jaccard 0.20 between the strict and the exploratory set (about 80 % of the list replaced) | `B_jaccard_matrix.csv` |
| 20 of 29 candidates pass all ten criteria; 29 of 29 pass the strict criterion | `B_candidate_robustness.csv` |
| 25 of 29 candidates not mentioned in the body text of six competing studies | `C_gene_crosswalk_matrix.csv` |

---

## 6. Scope of inference (reproduced here for clarity)

- The discovery cohort has **three animals per time point**; although limma's variance moderation improves power, no gene passes strict correction at 0.5 or 1 d.
- **FPKM with a pseudo-count of 0.1** is a legitimate and internally consistent choice, but it is not a count-based model; a different pseudo-count would change fold-change magnitudes for lowly expressed genes.
- The stimulation dataset (**GSE243038**) contains motoneurons only and its three samples per condition are three different stimulation frequencies, not biological replicates. It supports **no inferential claim** of any kind.
- The single-cell reference (**GSE189070**) contains **no neuronal cluster** and is not spatial. Neuronal expression could not be assessed for any gene, and no spatial inference is made.
- **No wet-laboratory validation** was performed. Every statement in the manuscript is a hypothesis or a dataset description, not a demonstration of mechanism.

---

## 7. Licence and citation

- **Code**: MIT (see `LICENSE`).
- **Derived result tables and figures**: CC BY 4.0.
- The underlying expression data remain subject to the terms of their original GEO submitters.

Citation metadata is provided in `CITATION.cff`; Zenodo deposit metadata is in `.zenodo.json`.

---

## 8. Contact

Wen-Yi Mao, MD — maowenyi89@hotmail.com
Department of Neurosurgery, Nantong Second People's Hospital, Nantong, Jiangsu 226002, China
