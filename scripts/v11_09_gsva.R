# ============================================================
# v11 阶段 2 / 任务 34：通路分析重跑（真 GSVA + 标准基因集）
# ------------------------------------------------------------
# 关键修正（相对 v10）：
#   1. 方法名坐实：v10 脚本自称 "GSVA" 实则用 mean z-score
#      → 本脚本用 R GSVA 1.50.5 官方实现（gsvaParam + gsva）
#   2. 基因集标准化：v10 用手写 curated 基因集（来源不明）
#      → 主分析改用 MSigDB Hallmark 2020（50 集，大鼠映射覆盖率中位 86.5%）
#      → curated 基因集仅作神经/代谢补充
#   3. 标签错位修复：v10 的 gsva_pilot_v2.py 在 sort_values 后用
#      list(pathway_info.index).index(n) 取位置（恒为 0..n），导致
#      "排序后的通路名" 配 "未排序的数值"
#      → 本脚本一律用名称索引（result[pathway_name, timepoint]）
# ============================================================

suppressMessages({
  library(GSVA)
  library(limma)
})

BASE <- "F:/scs research"
OUT  <- file.path(BASE, "output", "v11")

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  sets <- lapply(lines, function(l) {
    p <- strsplit(l, "\t")[[1]]
    unique(p[-c(1, 2)])
  })
  names(sets) <- vapply(lines, function(l) strsplit(l, "\t")[[1]][1], character(1))
  sets
}

# ---------- 1. 表达矩阵（与主链同口径）----------
expr_lin <- read.csv(file.path(BASE, "output/GSE175760_expression_matrix.csv"),
                     row.names = 1, check.names = FALSE)
meta <- read.csv(file.path(BASE, "output/GSE175760_metadata.csv"), stringsAsFactors = FALSE)
expr <- as.matrix(log2(expr_lin + 0.1))
cat("表达矩阵:", nrow(expr), "基因 x", ncol(expr), "样本\n")
cat("口径: log2(FPKM + 0.1)  [与 limma 主分析一致]\n")

GROUPS <- c("Sham", "CCI_0.5d", "CCI_1d", "CCI_3d", "CCI_7d", "CCI_14d")
TP_LAB <- c("0.5d", "1d", "3d", "7d", "14d")
grp <- factor(meta$group, levels = GROUPS)
names(grp) <- meta$sample_id
expr <- expr[, names(grp)]
cat("分组:", paste(levels(grp), as.integer(table(grp)), sep = "=", collapse = "  "), "\n\n")

# ---------- 2. GSVA（真实现）----------
run_gsva <- function(gmt_path, tag) {
  gs <- read_gmt(gmt_path)
  cat(sprintf("===== GSVA: %s（%d 个基因集）=====\n", tag, length(gs)))
  keep <- vapply(gs, function(g) sum(g %in% rownames(expr)) >= 10, logical(1))
  gs <- gs[keep]
  cat("  可用基因集（>=10 基因）:", length(gs), "\n")

  param <- gsvaParam(expr, gs, kcdf = "Gaussian")
  sc <- gsva(param, verbose = FALSE)
  cat("  分数矩阵:", nrow(sc), "通路 x", ncol(sc), "样本\n")
  write.csv(data.frame(pathway = rownames(sc), sc, check.names = FALSE),
            file.path(OUT, sprintf("gsva_%s_scores.csv", tag)), row.names = FALSE)

  # ---------- 3. limma 检验（Sham vs 各时间点）----------
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- levels(grp)
  fit <- eBayes(lmFit(sc, design))
  ctrs <- paste0("CCI_", TP_LAB, " - Sham")
  cm <- makeContrasts(contrasts = ctrs, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cm))

  delta <- matrix(NA_real_, nrow(sc), length(TP_LAB),
                  dimnames = list(rownames(sc), TP_LAB))
  pval <- padj <- delta
  for (j in seq_along(TP_LAB)) {
    tt <- topTable(fit2, coef = j, number = Inf, sort.by = "none")
    idx <- match(rownames(sc), rownames(tt))          # <- 按名称对齐，不按位置
    delta[, j] <- tt$logFC[idx]
    pval[, j]  <- tt$P.Value[idx]
    padj[, j]  <- tt$adj.P.Val[idx]
  }

  out <- data.frame(
    pathway = rep(rownames(sc), times = length(TP_LAB)),
    timepoint = rep(TP_LAB, each = nrow(sc)),
    delta = as.vector(delta),
    P = as.vector(pval),
    FDR = as.vector(padj),
    stringsAsFactors = FALSE)
  out <- out[order(out$FDR), ]
  write.csv(out, file.path(OUT, sprintf("gsva_%s_delta_long.csv", tag)), row.names = FALSE)

  # 宽表（按名称索引生成，杜绝错位）
  dw <- as.data.frame(delta)
  dw$pathway <- rownames(sc)
  write.csv(dw[, c("pathway", TP_LAB)],
            file.path(OUT, sprintf("gsva_%s_delta.csv", tag)), row.names = FALSE)

  cat("\n  各时间点 FDR<0.05 的通路数: ",
      paste(sprintf("%s=%d", TP_LAB, colSums(padj < 0.05, na.rm = TRUE)), collapse = "  "), "\n")
  cat("  各时间点 FDR<0.25 的通路数: ",
      paste(sprintf("%s=%d", TP_LAB, colSums(padj < 0.25, na.rm = TRUE)), collapse = "  "), "\n")

  # 7d / 14d 最强通路
  cat("\n  7d 变化最强的 10 条通路:\n")
  o <- order(-abs(delta[, "7d"]))
  for (i in head(o, 10)) {
    cat(sprintf("    %-46s delta=%+.3f  FDR=%.2e\n",
                rownames(sc)[i], delta[i, "7d"], padj[i, "7d"]))
  }
  cat("\n  14d 变化最强的 10 条通路:\n")
  o <- order(-abs(delta[, "14d"]))
  for (i in head(o, 10)) {
    cat(sprintf("    %-46s delta=%+.3f  FDR=%.2e\n",
                rownames(sc)[i], delta[i, "14d"], padj[i, "14d"]))
  }
  list(scores = sc, delta = delta, padj = padj)
}

hall <- run_gsva(file.path(OUT, "v11_hallmark_rat.gmt"), "hallmark")
cat("\n")
cur  <- run_gsva(file.path(OUT, "v11_pathways_custom.gmt"), "custom")

# ---------- 4. 交叉一致性检查：Hallmark 与 curated 重叠通路方向是否一致 ----------
cat("\n===== 交叉检查：同义通路在两套基因集中的方向一致性 =====\n")
PAIRS <- list(c("HALLMARK_INFLAMMATORY_RESPONSE", "Inflammatory response"),
              c("HALLMARK_COMPLEMENT", "Complement & coagulation"),
              c("HALLMARK_IL6_JAK_STAT3_SIGNALING", "JAK-STAT signaling"),
              c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "TNF-a via NF-kB"),
              c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "Oxidative phosphorylation"),
              c("HALLMARK_GLYCOLYSIS", "Glycolysis"),
              c("HALLMARK_APOPTOSIS", "Apoptosis"),
              c("HALLMARK_HYPOXIA", "Hypoxia"),
              c("HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "ER stress / UPR"),
              c("HALLMARK_P53_PATHWAY", "p53 pathway"))
cons <- data.frame()
for (pr in PAIRS) {
  if (!(pr[1] %in% rownames(hall$delta)) || !(pr[2] %in% rownames(cur$delta))) next
  a <- hall$delta[pr[1], ]; b <- cur$delta[pr[2], ]
  agree <- sum(sign(a) == sign(b), na.rm = TRUE)
  cons <- rbind(cons, data.frame(hallmark = pr[1], custom = pr[2],
                                 sign_agree = sprintf("%d/5", agree),
                                 rho = round(cor(a, b, method = "spearman"), 3)))
}
print(cons, row.names = FALSE)
write.csv(cons, file.path(OUT, "gsva_cross_geneset_consistency.csv"), row.names = FALSE)

cat("\n完成。输出:\n")
cat("  gsva_hallmark_scores.csv / _delta.csv / _delta_long.csv\n")
cat("  gsva_custom_scores.csv   / _delta.csv / _delta_long.csv\n")
cat("  gsva_cross_geneset_consistency.csv\n")
