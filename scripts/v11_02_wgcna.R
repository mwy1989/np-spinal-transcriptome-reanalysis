# ============================================================
# v11 Step 2: 真 WGCNA（替代原 step4_modules.py 的 KMeans 改名做法）
#   - 输入: output/GSE175760_expression_matrix.csv (线性 FPKM)
#   - 流程: goodSamplesGenes -> pickSoftThreshold -> blockwiseModules
#           -> module-trait correlation
#   - 输出: output/v11/wgcna_*
# 用法: Rscript scripts/v11_02_wgcna.R
# ============================================================
suppressPackageStartupMessages({
  library(WGCNA)
})
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

BASE <- Sys.getenv("SCS_ROOT", unset = {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(.f)) normalizePath(file.path(dirname(.f[1]), "..")) else normalizePath(".")
})
OUT  <- file.path(BASE, "output")
V11  <- file.path(OUT, "v11")
dir.create(V11, showWarnings = FALSE, recursive = TRUE)

expr <- read.csv(file.path(OUT, "GSE175760_expression_matrix.csv"),
                 row.names = 1, check.names = FALSE)
meta <- read.csv(file.path(OUT, "GSE175760_metadata.csv"), stringsAsFactors = FALSE)
expr <- expr[complete.cases(expr), , drop = FALSE]

meta$group <- factor(meta$group,
  levels = c("Sham","CCI_0.5d","CCI_1d","CCI_3d","CCI_7d","CCI_14d"))
rownames(meta) <- meta$sample_id
meta <- meta[colnames(expr), , drop = FALSE]

log2expr <- log2(as.matrix(expr) + 0.1)

# ---- 基因过滤：去掉低表达/低变异基因 ----
# 保留至少在 3 个样本中 log2 表达 > 1 的基因（即 FPKM > 1）
keep <- rowSums(log2expr > 1) >= 3
datExpr <- t(log2expr[keep, ])
cat("基因过滤: ", nrow(log2expr), " -> ", ncol(datExpr), "\n", sep = "")
cat("样本数: ", nrow(datExpr), "\n\n", sep = "")

# ---- 样本聚类检查离群 ----
gsg <- goodSamplesGenes(datExpr, verbose = 0)
cat("goodSamplesGenes 通过: ", gsg$allOK, "\n", sep = "")
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  cat("  剔除后: ", nrow(datExpr), " 样本 x ", ncol(datExpr), " 基因\n", sep = "")
}

# ---- 软阈值选择 ----
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                         networkType = "signed", verbose = 0)
sftTab <- data.frame(power = sft$fitIndices$Power,
                     SFT_R2 = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
                     slope = sft$fitIndices$slope,
                     mean_k = sft$fitIndices$mean.k.)
write.csv(sftTab, file.path(V11, "wgcna_softthreshold.csv"), row.names = FALSE)
cat("\n=== 软阈值拟合 ===\n")
print(round(sftTab, 3), row.names = FALSE)

# 选第一个 R2 >= 0.8 的 power（若无则取 R2 最大者）
above <- sftTab$power[sftTab$SFT_R2 >= 0.8]
if (length(above) > 0) {
  power <- min(above)
} else {
  power <- sftTab$power[which.max(sftTab$SFT_R2)]
  cat("  注意: 无 power 达 R2>=0.8, 取 R2 最大者\n")
}
cat("\n选定 power =", power, "\n\n")

# ---- 模块识别 ----
net <- blockwiseModules(datExpr, power = power,
                        networkType = "signed",
                        TOMType = "signed",
                        minModuleSize = 30,
                        deepSplit = 2,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        maxBlockSize = 20000,
                        verbose = 0)

moduleColors <- labels2colors(net$colors)
names(moduleColors) <- colnames(datExpr)
tab <- table(moduleColors)
cat("=== 模块大小 ===\n"); print(tab)
cat("总模块数(含 grey):", length(tab), "\n\n")

# ---- 模块 eigengene 与性状相关 ----
MEs <- net$MEs
time_map <- c(Sham = 0, CCI_0.5d = 1, CCI_1d = 2, CCI_3d = 3, CCI_7d = 4, CCI_14d = 5)
traits <- data.frame(
  Time      = time_map[as.character(meta$group)],
  CCI_any   = as.integer(meta$group != "Sham"),
  Chronic14 = as.integer(meta$group == "CCI_14d"),
  CCI_7d    = as.integer(meta$group == "CCI_7d"),
  CCI_3d    = as.integer(meta$group == "CCI_3d"),
  row.names = rownames(meta)
)
modTrait <- cor(MEs, traits, use = "p")
modTraitP <- corPvalueStudent(modTrait, nrow(datExpr))
write.csv(round(modTrait, 4), file.path(V11, "wgcna_module_trait_cor.csv"))
write.csv(signif(modTraitP, 4), file.path(V11, "wgcna_module_trait_p.csv"))

cat("=== 模块-性状相关 (r) ===\n")
print(round(modTrait, 3))
cat("\n=== 模块-性状 P ===\n")
print(signif(modTraitP, 3))

# ---- 输出模块基因 ----
geneModule <- data.frame(gene = colnames(datExpr), module = moduleColors)
write.csv(geneModule, file.path(V11, "wgcna_gene_modules.csv"), row.names = FALSE)
for (m in setdiff(unique(moduleColors), "grey")) {
  g <- geneModule$gene[geneModule$module == m]
  write.csv(data.frame(gene = g),
            file.path(V11, paste0("wgcna_module_", m, "_genes.csv")), row.names = FALSE)
}
write.csv(data.frame(sample = rownames(MEs), MEs),
          file.path(V11, "wgcna_module_eigengenes.csv"), row.names = FALSE)

# ---- 与既有模块的对照 ----
old_yellow <- tryCatch(
  read.csv(file.path(OUT, "WGCNA_M_yellow_genes.csv"))$gene,
  error = function(e) character(0))
if (length(old_yellow) > 0) {
  # 找与 Chronic14 相关最强的模块
  best <- names(which.max(abs(modTrait[, "Chronic14"])))
  new_top <- geneModule$gene[geneModule$module == sub("ME", "", best)]
  ov <- length(intersect(old_yellow, new_top))
  cat("\n=== 与旧 M_yellow 对照 ===\n")
  cat("旧 M_yellow:", length(old_yellow), "基因\n")
  cat("新 Chronic14 最强模块 (", best, "):", length(new_top), "基因\n", sep = "")
  cat("重叠:", ov, " (Jaccard =",
      round(ov / length(union(old_yellow, new_top)), 3), ")\n", sep = "")
}

saveRDS(list(net = net, moduleColors = moduleColors, MEs = MEs,
             modTrait = modTrait, power = power),
        file.path(V11, "wgcna_result.rds"))
cat("\n完成。输出目录: output/v11/\n")
