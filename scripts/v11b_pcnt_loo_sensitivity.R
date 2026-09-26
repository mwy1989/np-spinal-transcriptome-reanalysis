# ============================================================
# v11b ③a + ③b：伪计数敏感性 + 留一动物 jackknife
#
# 目的（回应第二份 GPT 审稿意见"计数型稳健性检查"）：
#   GSE175760 官方只发布 FPKM（StringTie + Ballgown），GEO 中不存在 count 矩阵，
#   因此真正的 count-based 重分析须自行从 SRA 重新定量。本脚本改为检验同一组
#   顾虑的两个可低成本回答的部分：
#     ③a 结果是否依赖 log 变换的伪计数选择（0.01 / 0.1 / 1.0）
#     ③b 结果是否由单只动物驱动（18 次留一重采样）
#   两者都只用已存档的表达矩阵，不改动任何主链输出。
#
# 输入: output/GSE175760_expression_matrix.csv （线性 FPKM）
#       output/GSE175760_metadata.csv
#       output/v11/candidates_primary.csv    （29 候选）
# 输出: output/v11/v11b_pseudocount_sensitivity.csv  （③a，3 行）
#       output/v11/v11b_pseudocount_gene_fc.csv      （③a，29 行 × 3 口径峰值 logFC）
#       output/v11/v11b_loo_per_resample.csv         （③b，18 行）
#       output/v11/v11b_loo_candidate_stability.csv  （③b，29 行）
#
# 与主链完全一致的设定：log2(FPKM + P) + limma-trend，对比 CCI_t − Sham，
# 时间点 0.5/1/3/7/14 d，候选集取 candidates_primary.csv。
# 运行：Rscript scripts/v11b_pcnt_loo_sensitivity.R
# ============================================================
suppressPackageStartupMessages(library(limma))

BASE <- Sys.getenv("SCS_ROOT", unset = {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(.f)) normalizePath(file.path(dirname(.f[1]), "..")) else normalizePath(".")
})
OUT <- file.path(BASE, "output")
V11 <- file.path(OUT, "v11")

expr <- read.csv(file.path(OUT, "GSE175760_expression_matrix.csv"), row.names = 1, check.names = FALSE)
meta <- read.csv(file.path(OUT, "GSE175760_metadata.csv"), stringsAsFactors = FALSE)
expr <- expr[complete.cases(expr), , drop = FALSE]
meta$group <- factor(meta$group, levels = c("Sham", "CCI_0.5d", "CCI_1d", "CCI_3d", "CCI_7d", "CCI_14d"))
rownames(meta) <- meta$sample_id
meta <- meta[colnames(expr), ]

cand <- read.csv(file.path(V11, "candidates_primary.csv"))$gene
tps  <- c("0.5d", "1d", "3d", "7d", "14d")
FC_FLOOR <- 0.58   # 1.5 倍：主链探索集/严格集的效应量下限

cat("矩阵:", dim(expr)[1], "基因 x", dim(expr)[2], "样本 | 候选:", length(cand), "\n\n")

# 对给定表达矩阵与分组跑同一套 limma-trend，返回 logFC / BH / raw P
fit_block <- function(e, m) {
  d  <- model.matrix(~ 0 + group, data = m); colnames(d) <- levels(m$group)
  cm <- makeContrasts(CCI_0.5d - Sham, CCI_1d - Sham, CCI_3d - Sham, CCI_7d - Sham, CCI_14d - Sham, levels = d)
  f  <- eBayes(contrasts.fit(lmFit(e, d), cm), trend = TRUE)
  tt <- lapply(1:5, function(i) topTable(f, coef = i, number = Inf, sort.by = "none"))
  fc <- sapply(tt, function(x) x$logFC)
  bh <- sapply(tt, function(x) x$adj.P.Val)
  pv <- sapply(tt, function(x) x$P.Value)
  dimnames(fc) <- dimnames(bh) <- dimnames(pv) <- list(rownames(e), tps)
  list(fc = fc, bh = bh, pv = pv)
}

# ---------------- ③a 伪计数敏感性 ----------------
PSEUDO <- c(0.01, 0.1, 1.0)
a_rows <- list(); gene_fc <- NULL; strict_sets <- list()
for (P in PSEUDO) {
  r <- fit_block(log2(as.matrix(expr) + P), meta)
  strict <- rownames(r$fc)[apply(abs(r$fc) > FC_FLOOR & r$bh < 0.05, 1, any)]
  expl   <- rownames(r$fc)[apply(abs(r$fc) > FC_FLOOR & r$pv < 0.05, 1, any)]
  peak   <- tps[apply(abs(r$fc[cand, , drop = FALSE]), 1, which.max)]
  dirn   <- sapply(seq_along(cand), function(i) sign(r$fc[cand[i], peak[i]]))
  a_rows[[length(a_rows) + 1]] <- data.frame(
    pseudocount = P,
    strict_n = length(strict),
    exploratory_n = length(expl),
    candidates_FDR_significant_of_29 = sum(apply(r$bh[cand, , drop = FALSE] < 0.05, 1, any))
  )
  gene_fc <- cbind(gene_fc, r$fc[cand, ][cbind(seq_along(cand), match(peak, tps))])
  strict_sets[[as.character(P)]] <- strict
  assign(paste0("peak_", P), peak); assign(paste0("dir_", P), dirn)
}
acc <- do.call(rbind, a_rows)
gene_fc <- data.frame(gene = cand, gene_fc, check.names = FALSE)
colnames(gene_fc) <- c("gene", paste0("peak_log2FC_pcnt", PSEUDO))
acc$strict_Jaccard_vs_0.1 <- sapply(PSEUDO, function(P) {
  a <- strict_sets[[as.character(P)]]; b <- strict_sets[["0.1"]]
  round(length(intersect(a, b)) / length(union(a, b)), 3)
})
acc$peak_timepoint_identical_to_0.1 <- c(NA, all(peak_0.1 == peak_0.01), all(peak_0.1 == peak_1))
acc$direction_identical_to_0.1      <- c(NA, all(dir_0.1 == dir_0.01),  all(dir_0.1 == dir_1))
write.csv(acc, file.path(V11, "v11b_pseudocount_sensitivity.csv"), row.names = FALSE)
write.csv(gene_fc, file.path(V11, "v11b_pseudocount_gene_fc.csv"), row.names = FALSE)

cat("=== ③a 伪计数敏感性 ===\n"); print(acc, row.names = FALSE)
cat("\n示例候选峰值 log2FC：\n")
for (g in c("Reg3b", "Gapt", "Nlrc4", "C1qa", "Cfh", "Hexb"))
  cat(sprintf("  %-8s %s\n", g, paste(sprintf("%.2f", unlist(gene_fc[gene_fc$gene == g, -1])), collapse = " -> ")))
cat("全部 29 候选峰值时点与 0.1 一致:", all(peak_0.1 == peak_0.01) && all(peak_0.1 == peak_1), "\n\n")

# ---------------- ③b 留一动物 jackknife ----------------
loo <- NULL; per <- NULL
for (s in colnames(expr)) {
  cols <- setdiff(colnames(expr), s)
  m2 <- droplevels(meta[cols, ])
  if (!all(levels(meta$group) %in% levels(m2$group))) next
  r <- fit_block(log2(as.matrix(expr[, cols, drop = FALSE]) + 0.1), m2)
  sig  <- apply(r$bh[cand, , drop = FALSE] < 0.05, 1, any)
  peak <- tps[apply(abs(r$fc[cand, , drop = FALSE]), 1, which.max)]
  loo <- rbind(loo, data.frame(dropped_sample = s, gene = cand,
                               FDR_significant = as.integer(sig), peak_timepoint = peak))
  per <- rbind(per, data.frame(dropped_sample = s, n_candidates_FDR_significant = sum(sig)))
}
write.csv(loo, file.path(V11, "v11b_loo_per_candidate.csv"), row.names = FALSE)
write.csv(per, file.path(V11, "v11b_loo_per_resample.csv"), row.names = FALSE)
n_rep <- length(unique(loo$dropped_sample))
stab <- do.call(rbind, lapply(cand, function(g) {
  sub <- loo[loo$gene == g, ]
  data.frame(gene = g,
             n_resamples_FDR_significant = sum(sub$FDR_significant),
             n_resamples = n_rep,
             fraction = round(sum(sub$FDR_significant) / n_rep, 3),
             peak_timepoints_observed = paste(sort(unique(sub$peak_timepoint[sub$FDR_significant == 1])), collapse = "/"))
}))
stab <- stab[order(stab$n_resamples_FDR_significant), ]
write.csv(stab, file.path(V11, "v11b_loo_candidate_stability.csv"), row.names = FALSE)

cat("=== ③b 留一动物（", n_rep, " 次重采样）===\n", sep = "")
cat("在全部重采样中均保持 FDR 显著的候选:", sum(stab$n_resamples_FDR_significant == n_rep), "/", length(cand), "\n")
cat("最低保留次数:", min(stab$n_resamples_FDR_significant), "/", n_rep,
    " (", sprintf("%.0f%%", 100 * min(stab$fraction)), ")\n", sep = "")
cat("每次重采样仍显著的候选数范围:", min(per$n_candidates_FDR_significant), "-",
    max(per$n_candidates_FDR_significant), "\n")
flip <- stab[grepl("/", stab$peak_timepoints_observed), ]
cat("峰值时点不唯一（随删样本改变）的候选:", nrow(flip), "/", length(cand), "\n")
cat("其中峰值时点集合均为 7d 与 14d 者:", sum(flip$peak_timepoints_observed == "14d/7d"), "\n")
cat("峰值出现在 0.5d 或 1d 者:",
    sum(grepl("0.5d|1d", stab$peak_timepoints_observed)), "（注：'1d' 子串需人工核对）\n")
early <- grepl("(^|/)0\\.5d|(^|/)1d", stab$peak_timepoints_observed)
cat("严格判定——峰值集合含 0.5d 或 1d 的候选:", sum(early), "\n")
cat("\n最脆弱的 8 个候选:\n"); print(head(stab[, c("gene", "n_resamples_FDR_significant", "peak_timepoints_observed")], 8), row.names = FALSE)
cat("\n输出已写入:", V11, "\n")
