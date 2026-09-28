# Supplementary Table S1. Sensitivity of the results to the pseudo-count, to the omission of a single animal, and to the pooling of the myeloid subtypes

Bilingual companion table to Sections 2.11, 2.15, 3.5 and 3.6 of the manuscript. Panels A–C were produced by `scripts/v11b_pcnt_loo_sensitivity.R`; Panel D was produced by `scripts/v11_14_myeloid_merge.py`. No data beyond `output/GSE175760_expression_matrix.csv`, the single-cell reference `GSE189070_annotated.h5ad` and the 29 candidates in `output/v11/candidates_primary.csv` were used. The expression pipeline is identical to the main chain (log2(FPKM + pseudo-count), limma-trend, contrasts CCI_t − Sham at 0.5, 1, 3, 7 and 14 d).

## Panel A. Pseudo-count sensitivity

| Pseudo-count | Strict set (BH FDR < 0.05, \|log2FC\| > 0.58) | Exploratory set (P < 0.05, \|log2FC\| > 0.58) | Candidates FDR-significant (of 29) | Jaccard of the strict set against pseudo-count 0.1 |
|---|---|---|---|---|
| 0.01 | 357 | 1977 | 29 | 0.905 |
| 0.10 | 386 | 1940 | 29 | — (reference) |
| 1.00 | 352 | 1355 | 29 | 0.84 |

## Panel B. Peak log2 fold change of each candidate under three pseudo-counts

Values are the log2 fold change at the time point of largest absolute effect for that candidate (as defined in Section 2.15). The peak time point and the direction of change were identical under all three pseudo-counts for all 29 candidates.

| Candidate gene | Pseudo-count 0.01 | Pseudo-count 0.10 | Pseudo-count 1.00 |
|---|---|---|---|
| Reg3b | 6.05 | 5.62 | 3.85 |
| C1qa | 2.73 | 2.73 | 2.71 |
| Cfh | 1.25 | 1.25 | 1.17 |
| Gapt | 3.16 | 2.89 | 1.69 |
| Hexb | 0.67 | 0.67 | 0.66 |
| Nrp1 | 1.14 | 1.12 | 0.92 |
| Laptm5 | 1.58 | 1.58 | 1.51 |
| Nlrc4 | 2.16 | 1.91 | 0.92 |
| Mmp3 | 3.37 | 2.80 | 1.24 |
| Ecel1 | 0.81 | 0.80 | 0.69 |
| Clec7a | 4.05 | 3.46 | 1.70 |
| Arhgap25 | 1.51 | 1.45 | 1.09 |
| Mrpl10 | 1.64 | 1.63 | 1.53 |
| Fcrl2 | 4.61 | 4.12 | 2.45 |
| Cd68 | 2.67 | 2.59 | 2.01 |
| Lmo2 | 0.85 | 0.84 | 0.75 |
| Rfx7 | -3.07 | -2.91 | -2.00 |
| Zc3h12a | 1.42 | 1.30 | 0.72 |
| C4a | 3.34 | 3.30 | 3.00 |
| Cd84 | 1.19 | 1.17 | 1.00 |
| Man2b2 | 0.95 | 0.95 | 0.88 |
| Klhl6 | 1.59 | 1.52 | 1.06 |
| Fbxw4 | 0.93 | 0.92 | 0.83 |
| Fbxl7 | -2.43 | -1.81 | -0.60 |
| Dnase2 | 0.91 | 0.90 | 0.84 |
| AABR07034445.1 | -12.57 | -9.25 | -5.95 |
| Chrna1 | 1.36 | 1.15 | 0.46 |
| Svop | -0.86 | -0.86 | -0.80 |
| Rassf3 | 0.61 | 0.60 | 0.51 |

## Panel C. Leave-one-animal-out stability of each candidate

The differential-expression model was refitted 18 times, each time omitting one of the 18 animals, and the adjusted significance and peak time point of each of the 29 already-selected candidates were recorded in every replicate. Clustering, network construction and machine-learning selection were not re-run, so this panel describes the stability of the differential-expression step for a fixed candidate list and not the stability of the selection procedure. "Resamples with adjusted significance" counts how many of the 18 replicates retained BH FDR < 0.05 at at least one time point. "Peak time points observed" lists the time points of largest absolute effect among those significant replicates; an entry with two time points means the assignment changed when an animal was omitted.

| Candidate gene | Resamples with adjusted significance (of 18) | Fraction | Peak time points observed |
|---|---|---|---|
| AABR07034445.1 | 13 | 72 % | 14d |
| Chrna1 | 13 | 72 % | 14d |
| Svop | 13 | 72 % | 7d |
| Dnase2 | 14 | 78 % | 3d/7d |
| Rassf3 | 14 | 78 % | 14d/3d |
| Fbxl7 | 15 | 83 % | 14d |
| Fbxw4 | 16 | 89 % | 14d |
| Man2b2 | 16 | 89 % | 7d |
| Arhgap25 | 18 | 100 % | 14d/7d |
| C1qa | 18 | 100 % | 7d |
| C4a | 18 | 100 % | 7d |
| Cd68 | 18 | 100 % | 14d/7d |
| Cd84 | 18 | 100 % | 7d |
| Cfh | 18 | 100 % | 7d |
| Clec7a | 18 | 100 % | 14d/7d |
| Ecel1 | 18 | 100 % | 14d |
| Fcrl2 | 18 | 100 % | 3d/7d |
| Gapt | 18 | 100 % | 7d |
| Hexb | 18 | 100 % | 7d |
| Klhl6 | 18 | 100 % | 3d/7d |
| Laptm5 | 18 | 100 % | 14d/7d |
| Lmo2 | 18 | 100 % | 14d/7d |
| Mmp3 | 18 | 100 % | 14d/7d |
| Mrpl10 | 18 | 100 % | 14d |
| Nlrc4 | 18 | 100 % | 7d |
| Nrp1 | 18 | 100 % | 14d/7d |
| Reg3b | 18 | 100 % | 14d/7d |
| Rfx7 | 18 | 100 % | 14d |
| Zc3h12a | 18 | 100 % | 14d/7d |

## Panel D. Cell-type attribution of the assessable candidates before and after pooling the myeloid subtypes

The single-cell reference annotates three separate myeloid subtypes (microglia, macrophages, neutrophils). Panel D reports the dominant type of each candidate under those original labels (left) and after pooling the three subtypes into one myeloid category and recomputing the mean expression over the pooled cells (right), using two weightings: the mean over all pooled cells, and the mean of the three subtype means with equal weight. Both weightings gave the same dominant type for every assessable candidate; for one of the eight that failed the criterion (Gapt) the two weightings fall on opposite sides of the myeloid boundary. Only the 16 assessable candidates are shown; genes that could not be assigned a dominant type are listed in the notes. Under the original labels 10 of the 16 assessable candidates were microglia-dominant; after pooling and recomputation 9 were myeloid-dominant. Counting the original subtype assignments as myeloid without recomputing them would give 13 of 16, but that measures the labels rather than the expression and is not used in the manuscript.

| Candidate gene | Original dominant type (10 labels) | Pooled dominant type (cell-weighted) | Pooled dominant type (equal weight) | Myeloid mean (cell-weighted) | Best non-myeloid type (mean) | Myeloid-dominant after pooling |
|---|---|---|---|---|---|---|
| C1qa | Microglia | Myeloid | Myeloid | 3.11 | Endothelial (1.12) | yes |
| Hexb | Microglia | Myeloid | Myeloid | 2.62 | Endothelial (0.93) | yes |
| Laptm5 | Microglia | Myeloid | Myeloid | 2.05 | T_cell (1.69) | yes |
| Cd68 | Macrophage | Myeloid | Myeloid | 1.67 | Pericyte (0.40) | yes |
| Cd84 | Microglia | Myeloid | Myeloid | 0.89 | T_cell (0.57) | yes |
| Cfh | Microglia | Pericyte | Pericyte | 0.56 | Pericyte (0.59) | no |
| Nrp1 | Microglia | Endothelial | Endothelial | 0.54 | Endothelial (0.62) | no |
| Lmo2 | Microglia | Endothelial | Endothelial | 0.54 | Endothelial (0.63) | no |
| Klhl6 | B_cell | B_cell | B_cell | 0.28 | B_cell (0.54) | no |
| Fbxw4 | Microglia | Myeloid | Myeloid | 0.24 | T_cell (0.12) | yes |
| Arhgap25 | Microglia | B_cell | B_cell | 0.24 | B_cell (0.27) | no |
| Rassf3 | Neutrophil | Myeloid | Myeloid | 0.20 | Endothelial (0.15) | yes |
| Zc3h12a | Neutrophil | Myeloid | Myeloid | 0.13 | B_cell (0.10) | yes |
| Man2b2 | Microglia | Myeloid | Myeloid | 0.10 | T_cell (0.05) | yes |
| Mrpl10 | Ependymal | Ependymal | Ependymal | 0.09 | Ependymal (0.18) | no |
| Rfx7 | B_cell | B_cell | B_cell | 0.06 | B_cell (0.20) | no |

**Notes.** (i) Pseudo-count 0.10 is the value used throughout the manuscript. (ii) The strict set is the union across the five post-injury time points; the 29 candidates were defined before this analysis and were not re-selected here. (iii) Across the 18 leave-one-animal-out replicates the number of candidates retaining adjusted significance ranged from 22 to 29. (iv) No candidate in any replicate had its peak at 0.5 d or 1 d. (v) A count-based re-analysis could not be performed on these data: the processed matrix used here contains no original read counts, so a count model could not be fitted to it; such an analysis would require re-quantifying the deposited reads. (vi) Five of the 29 candidates (Clec7a, Fcrl2, C4a, Dnase2, AABR07034445.1) were not found in the annotation matrix of the single-cell reference, and eight further candidates failed the assessability criterion, which requires both a maximum detection rate above 10 % and a maximum mean expression above 0.1: seven had maximum detection rates below 10 % (Fbxl7 (0.1 %); Reg3b (0.1 %); Svop (0.3 %); Ecel1 (0.5 %); Mmp3 (1.1 %); Chrna1 (2.1 %); Nlrc4 (2.5 %)) and one (Gapt) reached a maximum detection rate of 11.0 % with a maximum mean expression of 0.09, below the 0.1 required. Neither group is shown in Panel D. (vii) Pooling the three myeloid subtypes changes the label of the 13 assessable candidates that were attributed to microglia, macrophages or neutrophils, but the dominant type after recomputing the mean over the pooled cells is myeloid for only 9 of the 16 assessable candidates; the other four (Cfh, Nrp1, Lmo2, Arhgap25) fall just behind a non-myeloid type. The values for all detected candidates, including those that are not assessable, are provided in `results/myeloid_merge_dominant_celltype.csv`. (viii) The two pooled weightings gave the same dominant type for all 16 assessable candidates; among the eight that failed the criterion they disagree for Gapt only (cell-weighted B cell, equal weight myeloid), which is one of the reasons that gene is not interpreted.
