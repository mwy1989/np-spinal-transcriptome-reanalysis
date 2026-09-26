# ============================================================
# v11 标准口径重跑：真 limma (limma-trend) 替代自制近似
# 输入: output/GSE175760_expression_matrix.csv (线性 FPKM)
# 输出: output/v11/limma_standard_all.csv  +  控制台诊断
# ============================================================
suppressPackageStartupMessages(library(limma))

BASE <- "F:/scs research"
OUT  <- file.path(BASE, "output")
V11  <- file.path(OUT, "v11")
dir.create(V11, showWarnings = FALSE, recursive = TRUE)

expr <- read.csv(file.path(OUT, "GSE175760_expression_matrix.csv"), row.names = 1, check.names = FALSE)
meta <- read.csv(file.path(OUT, "GSE175760_metadata.csv"), stringsAsFactors = FALSE)
expr <- expr[complete.cases(expr), , drop = FALSE]

cat("矩阵:", dim(expr)[1], "基因 x", dim(expr)[2], "样本\n")
cat("表达范围(线性): min=", round(min(expr), 3), " max=", round(max(expr), 3),
    " median=", round(median(as.matrix(expr)), 3), "\n\n", sep = "")

# ---- 样本分组 ----
meta$group <- factor(meta$group, levels = c("Sham", "CCI_0.5d", "CCI_1d", "CCI_3d", "CCI_7d", "CCI_14d"))
rownames(meta) <- meta$sample_id
meta <- meta[colnames(expr), ]
cat("分组:\n"); print(table(meta$group)); cat("\n")

# ---- 对数转换（两种伪计数，主用 0.1 与现主链一致）----
PSEUDO <- 0.1
log2expr <- log2(as.matrix(expr) + PSEUDO)

# ============================================================
# limma-trend
# ============================================================
design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)

fit  <- lmFit(log2expr, design)
cm   <- makeContrasts(
  CCI_0.5d - Sham,
  CCI_1d   - Sham,
  CCI_3d   - Sham,
  CCI_7d   - Sham,
  CCI_14d  - Sham,
  levels = design
)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2, trend = TRUE)

tps <- c("0.5d", "1d", "3d", "7d", "14d")
res <- list()
for (i in seq_along(tps)) {
  tt <- topTable(fit2, coef = i, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt$tp   <- tps[i]
  res[[i]] <- tt[, c("gene", "tp", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
}
allres <- do.call(rbind, res)
write.csv(allres, file.path(V11, "limma_standard_all.csv"), row.names = FALSE)

cat("============================================================\n")
cat("limma-trend 结果  (pseudo-count =", PSEUDO, ")\n")
cat("============================================================\n\n")

FC <- 0.58
cat(sprintf("%-8s %10s %12s %12s %12s\n", "时间点", "|FC|>0.58", "P<.05+FC", "BH<.05+FC", "BH<.25+FC"))
for (tp in tps) {
  d <- allres[allres$tp == tp, ]
  n_onlyfc <- sum(abs(d$logFC) > FC)
  n_raw    <- sum(abs(d$logFC) > FC & d$P.Value   < 0.05)
  n_bh05   <- sum(abs(d$logFC) > FC & d$adj.P.Val < 0.05)
  n_bh25   <- sum(abs(d$logFC) > FC & d$adj.P.Val < 0.25)
  cat(sprintf("%-8s %10d %12d %12d %12d\n", tp, n_onlyfc, n_raw, n_bh05, n_bh25))
}

u <- function(mask_fn) {
  s <- character(0)
  for (tp in tps) { d <- allres[allres$tp == tp, ]; s <- union(s, d$gene[mask_fn(d)]) }
  length(s)
}
cat("\nunion 基因数:\n")
cat("  仅 |FC|>0.58        :", u(function(d) abs(d$logFC) > FC), "\n")
cat("  P<0.05 & |FC|>0.58  :", u(function(d) abs(d$logFC) > FC & d$P.Value   < 0.05), "\n")
cat("  BH<0.05 & |FC|>0.58 :", u(function(d) abs(d$logFC) > FC & d$adj.P.Val < 0.05), "\n")
cat("  BH<0.25 & |FC|>0.58 :", u(function(d) abs(d$logFC) > FC & d$adj.P.Val < 0.25), "\n")

# ---- 22 hub 基因存活 ----
HUB <- c("Il10rb","Jak3","Cd37","Per1","Dera","Dbr1","Lpxn","Nat9","Pbx2","Fbxw4",
         "Edem1","Nrp1","Itga7","Casp3","Svop","Mrpl10","Mtif3","Smc2","Akap13",
         "Emilin2","Apool","Abhd1")
HUB <- HUB[HUB %in% rownames(expr)]
cat("\n22 hub 基因在 limma 下的显著性（P<0.05 & |FC|>0.58）:\n")
cat(sprintf("%-9s", ""), sprintf("%7s", tps), "\n", sep = "")
surv_raw <- c(); surv_bh <- c()
for (g in HUB) {
  d <- allres[allres$gene == g, ]; d <- d[match(tps, d$tp), ]
  mk <- ifelse(abs(d$logFC) > FC & d$P.Value < 0.05, "  Y", "  .")
  cat(sprintf("%-9s", g), sprintf("%7s", mk), "\n", sep = "")
  if (any(abs(d$logFC) > FC & d$P.Value   < 0.05)) surv_raw <- c(surv_raw, g)
  if (any(abs(d$logFC) > FC & d$adj.P.Val < 0.05)) surv_bh  <- c(surv_bh,  g)
}
cat("\n  至少1个时间点 P<0.05 :", length(surv_raw), "/22\n")
cat("  至少1个时间点 BH<0.05:", length(surv_bh),  "/22\n")

cat("\n完成。结果已写入: output/v11/limma_standard_all.csv\n")
