# Supplementary Table S1. Sensitivity of the results to the pseudo-count and to the omission of a single animal

Bilingual companion table to Section 2.15 and Section 3.6 of the manuscript. Both panels were produced by `scripts/v11b_pcnt_loo_sensitivity.R` from the deposited expression matrix; no data beyond `output/GSE175760_expression_matrix.csv` and the 29 candidates in `output/v11/candidates_primary.csv` were used. The analysis pipeline is identical to the main chain (log2(FPKM + pseudo-count), limma-trend, contrasts CCI_t − Sham at 0.5, 1, 3, 7 and 14 d).

## Panel A. Pseudo-count sensitivity

| Pseudo-count | Strict set (BH FDR < 0.05, \|log2FC\| > 0.58) | Exploratory set (P < 0.05, \|log2FC\| > 0.58) | Candidates FDR-significant (of 29) | Jaccard of the strict set against pseudo-count 0.1 |
|---|---|---|---|---|
| 0.01 | 357 | 1977 | 29 | 0.905 |
| 0.10 | 386 | 1940 | 29 | — (reference) |
| 1.00 | 352 | 1355 | 29 | 0.84 |

## Panel B. Peak log2 fold change of each candidate under three pseudo-counts

Values are the log2 fold change at the time point of largest absolute effect for that candidate. The peak time point and the direction of change were identical under all three pseudo-counts for all 29 candidates.

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

The whole pipeline was repeated 18 times, each time omitting one of the 18 animals. "Resamples with adjusted significance" counts how many of the 18 replicates retained BH FDR < 0.05 at at least one time point. "Peak time points observed" lists the time points of largest absolute effect among those significant replicates; an entry with two time points means the assignment changed when an animal was omitted.

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

**Notes.** (i) Pseudo-count 0.10 is the value used throughout the manuscript. (ii) The strict set is the union across the five post-injury time points; the 29 candidates were defined before this analysis and were not re-selected here. (iii) Across the 18 leave-one-animal-out replicates the number of candidates retaining adjusted significance ranged from 22 to 29. (iv) No candidate in any replicate had its peak at 0.5 d or 1 d. (v) A count-based re-analysis could not be performed on these data: the sequencing reads of GSE175760 were deposited as FPKM values only, so such an analysis would require re-quantifying the original reads.
