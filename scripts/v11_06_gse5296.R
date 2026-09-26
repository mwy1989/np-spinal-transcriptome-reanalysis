# ============================================================
# v11 阶段 2 / 任务 35：GSE5296 跨队列验证重跑
# ------------------------------------------------------------
# 关键修正（相对 v10）：
#   1. 矩阵尺度：GSE5296_gene_expression.csv 是 RMA log2 值的 2^x 线性还原
#      （实测 log2 后 median=6.08, 范围 -0.69~11.87，正是 RMA 区间）
#      → 必须先取 log2；v10 直接用线性值算 log2FC 属尺度错误
#   2. 三品系池化问题：用 limma 以 strain 作协变量（v10 直接合并均值）
#   3. 删除 Fisher 合并 P 值（存在选择依赖，不是独立验证）
#   4. region 实测只有 impact（96 样本），Methods 不得写"三个脊髓区域"
# ============================================================

suppressMessages({
  library(limma)
})

BASE <- "F:/scs research"
OUT  <- file.path(BASE, "output", "v11")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---------- 1. 读入并 log2 还原 ----------
cat("===== 1. 数据读入与尺度处理 =====\n")
expr_lin <- read.csv(file.path(BASE, "output/GSE5296/GSE5296_gene_expression.csv"),
                     row.names = 1, check.names = FALSE)
meta <- read.csv(file.path(BASE, "output/GSE5296/GSE5296_sample_metadata.csv"),
                 row.names = 1, stringsAsFactors = FALSE)

cat("线性矩阵:", dim(expr_lin)[1], "基因 x", dim(expr_lin)[2], "样本\n")
cat(sprintf("线性尺度: min=%.3f max=%.3f median=%.3f\n",
            min(expr_lin), max(expr_lin), median(as.matrix(expr_lin))))

expr <- log2(as.matrix(expr_lin))
cat(sprintf("log2 尺度: min=%.3f max=%.3f median=%.3f  <- RMA 典型区间\n",
            min(expr), max(expr), median(expr)))

# 样本对齐
common <- intersect(colnames(expr), rownames(meta))
expr <- expr[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]
cat("对齐后样本数:", length(common), "\n")
cat("region 分布:", paste(names(table(meta$region)), table(meta$region), sep = "=", collapse = ", "), "\n")
cat("区域数:", length(unique(meta$region)),
    "(v10 稿件声称 3 个区域 impact/rostral/caudal -> 需修正)\n\n")

# ---------- 2. 逐时间点 limma（strain 作协变量）----------
cat("===== 2. 逐时间点 limma（~ strain + condition）=====\n")
TPS <- c("0.5h", "4h", "24h", "72h", "7d", "28d")

allres <- list()
for (tp in TPS) {
  sub <- meta[meta$timepoint == tp, , drop = FALSE]
  if (nrow(sub) == 0) next
  e <- expr[, rownames(sub), drop = FALSE]
  grp <- factor(sub$condition, levels = c("Sham", "Injured"))
  st  <- factor(sub$strain)
  cat(sprintf("%-5s: n=%d (Injured=%d, Sham=%d)  品系=%s\n",
              tp, nrow(sub), sum(grp == "Injured"), sum(grp == "Sham"),
              paste(levels(st), collapse = "+")))

  design <- model.matrix(~ st + grp)
  fit <- eBayes(lmFit(e, design))
  tt <- topTable(fit, coef = ncol(design), number = Inf, sort.by = "none")

  res <- data.frame(gene = rownames(tt), tp = tp,
                    log2FC = tt$logFC, AveExpr = tt$AveExpr,
                    t = tt$t, P_raw = tt$P.Value, P_BH = tt$adj.P.Val,
                    row.names = NULL)
  allres[[tp]] <- res
}
allres <- do.call(rbind, allres)
write.csv(allres, file.path(OUT, "gse5296_limma_all.csv"), row.names = FALSE)
cat("\n输出: gse5296_limma_all.csv\n\n")

# ---------- 3. 各时间点 DEG 计数 ----------
cat("===== 3. GSE5296 各时间点 DEG 计数 =====\n")
cat(sprintf("%-6s %10s %14s %14s\n", "tp", "仅|FC|>0.58", "P<0.05+FC", "BH<0.05+FC"))
for (tp in TPS) {
  d <- allres[allres$tp == tp, ]
  cat(sprintf("%-6s %10d %14d %14d\n", tp,
              sum(abs(d$log2FC) > 0.58),
              sum(abs(d$log2FC) > 0.58 & d$P_raw  < 0.05),
              sum(abs(d$log2FC) > 0.58 & d$P_BH   < 0.05)))
}

# ---------- 4. 新候选基因跨队列方向一致性 ----------
cat("\n===== 4. 候选基因跨队列方向一致性 =====\n")
cand <- read.csv(file.path(OUT, "candidates_primary.csv"), stringsAsFactors = FALSE)$gene
cat("候选基因数:", length(cand), "\n")

# GSE175760 limma 结果
g175 <- read.csv(file.path(OUT, "limma_all_contrasts.csv"), stringsAsFactors = FALSE)

# 时间点配对（最接近的对应关系）
PAIRS <- list(c("7d", "7d"), c("7d", "14d_rat"), c("28d", "14d_rat"))
# 明确的配对： GSE5296 tp  <->  GSE175760 tp
MATCH <- data.frame(
  gse5296 = c("0.5h", "4h", "24h", "72h", "7d", "28d"),
  gse175760 = c("0.5d", "0.5d", "1d", "3d", "7d", "14d"),
  stringsAsFactors = FALSE
)

cons_all <- list()
for (i in seq_len(nrow(MATCH))) {
  t5296 <- MATCH$gse5296[i]; t175 <- MATCH$gse175760[i]
  a <- allres[allres$tp == t5296, c("gene", "log2FC", "P_raw", "P_BH")]
  b <- g175[g175$tp == t175, c("gene", "log2FC", "P_raw", "P_BH")]
  m <- merge(a, b, by = "gene", suffixes = c("_5296", "_175760"))
  m <- m[m$gene %in% cand, ]
  if (nrow(m) < 3) next

  ok <- is.finite(m$log2FC_5296) & is.finite(m$log2FC_175760)
  rho <- cor(m$log2FC_5296[ok], m$log2FC_175760[ok], method = "spearman")
  conc <- sum(sign(m$log2FC_5296[ok]) == sign(m$log2FC_175760[ok]))
  n <- sum(ok)

  cons_all[[length(cons_all) + 1]] <- data.frame(
    gse5296_tp = t5296, gse175760_tp = t175,
    n_candidates_detected = n,
    direction_concordant = conc,
    concordance_rate = round(conc / n, 3),
    spearman_rho = round(rho, 3),
    stringsAsFactors = FALSE)
}
cons <- do.call(rbind, cons_all)
cat("\n候选基因跨队列方向一致性（描述性，不报合并 P）:\n")
print(cons, row.names = FALSE)
write.csv(cons, file.path(OUT, "gse5296_candidate_concordance.csv"), row.names = FALSE)

# ---------- 5. 7d 详细对照（主结论时间点）----------
cat("\n===== 5. 7d 配对详细对照 =====\n")
a <- allres[allres$tp == "7d", c("gene", "log2FC", "P_raw", "P_BH")]
b <- g175[g175$tp == "7d", c("gene", "log2FC", "P_raw", "P_BH")]
m <- merge(a, b, by = "gene", suffixes = c("_mouse5296", "_rat175760"))
m <- m[m$gene %in% cand, ]
m <- m[order(-abs(m$log2FC_rat175760)), ]
m$dir_5296  <- ifelse(m$log2FC_mouse5296 > 0, "up", "down")
m$dir_175760 <- ifelse(m$log2FC_rat175760 > 0, "up", "down")
m$concordant <- m$dir_5296 == m$dir_175760
write.csv(m, file.path(OUT, "gse5296_7d_candidate_detail.csv"), row.names = FALSE)
cat(sprintf("7d: 可评估候选基因 %d 个，方向一致 %d 个（%.0f%%）\n",
            nrow(m), sum(m$concordant), 100 * mean(m$concordant)))
print(m[, c("gene", "log2FC_mouse5296", "log2FC_rat175760", "concordant")], row.names = FALSE)

# ---------- 6. WGCNA 模块过表达检验（替代 M_yellow）----------
cat("\n===== 6. 新 WGCNA 模块在 GSE5296 DEG 中的过表达 =====\n")
gm <- read.csv(file.path(OUT, "wgcna_degs_gene_modules.csv"), stringsAsFactors = FALSE)
MODS <- c("turquoise", "blue", "brown")

hyp <- function(mod_genes, deg_genes, bg_genes) {
  mod_genes <- intersect(mod_genes, bg_genes)
  deg_genes <- intersect(deg_genes, bg_genes)
  N <- length(bg_genes); K <- length(mod_genes); n <- length(deg_genes)
  k <- length(intersect(mod_genes, deg_genes))
  if (K == 0 || n == 0) return(c(k = 0, expected = 0, fold = NA, p = NA))
  p <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
  c(k = k, expected = K * n / N, fold = k / (K * n / N), p = p)
}

rows <- list()
for (tp in TPS) {
  d <- allres[allres$tp == tp, ]
  deg <- d$gene[abs(d$log2FC) > 0.58 & d$P_BH < 0.05]
  bg <- rownames(expr)
  for (mo in MODS) {
    mg <- gm$gene[gm$module == mo]
    h <- hyp(mg, deg, bg)
    rows[[length(rows) + 1]] <- data.frame(
      tp = tp, module = mo, n_module = length(intersect(mg, bg)),
      n_deg = length(intersect(deg, bg)), overlap = h["k"],
      expected = round(h["expected"], 1), fold_enrichment = round(h["fold"], 2),
      hypergeom_P = h["p"], stringsAsFactors = FALSE)
  }
}
ovr <- do.call(rbind, rows)
ovr <- ovr[!is.na(ovr$hypergeom_P), ]
ovr <- ovr[order(ovr$hypergeom_P), ]
print(ovr, row.names = FALSE)
write.csv(ovr, file.path(OUT, "gse5296_module_overrepresentation.csv"), row.names = FALSE)

cat("\n完成。\n")
cat("  - gse5296_limma_all.csv              逐时间点 limma 结果（strain 校正）\n")
cat("  - gse5296_candidate_concordance.csv  候选基因跨队列方向一致性\n")
cat("  - gse5296_7d_candidate_detail.csv    7d 配对明细\n")
cat("  - gse5296_module_overrepresentation.csv  模块过表达\n")
