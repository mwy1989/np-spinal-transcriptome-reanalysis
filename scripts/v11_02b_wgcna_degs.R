# ============================================================
# v11 Step 2b: WGCNA on 探索集基因（1940）
#   理由: 全基因(11640) x n=18 会产生 27 个模块的过度切分；
#         改用探索集基因使基因数/样本数比例合理，模块更稳，
#         且与"从差异基因出发识别共表达模块"的叙事连贯。
#   用法: Rscript scripts/v11_02b_wgcna_degs.R
# ============================================================
suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

BASE <- Sys.getenv("SCS_ROOT", unset = {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(.f)) normalizePath(file.path(dirname(.f[1]), "..")) else normalizePath(".")
})
OUT  <- file.path(BASE, "output")
V11  <- file.path(OUT, "v11")

expr <- read.csv(file.path(OUT, "GSE175760_expression_matrix.csv"),
                 row.names = 1, check.names = FALSE)
meta <- read.csv(file.path(OUT, "GSE175760_metadata.csv"), stringsAsFactors = FALSE)
expr <- expr[complete.cases(expr), , drop = FALSE]
meta$group <- factor(meta$group,
  levels = c("Sham","CCI_0.5d","CCI_1d","CCI_3d","CCI_7d","CCI_14d"))
rownames(meta) <- meta$sample_id
meta <- meta[colnames(expr), , drop = FALSE]

# ---- 基因集：探索集 ∩ 通过表达过滤 ----
degs <- read.csv(file.path(V11, "DEG_exploratory_genes.csv"))$gene
log2expr_all <- log2(as.matrix(expr) + 0.1)
keep_expr <- rowSums(log2expr_all > 1) >= 3
genes <- intersect(degs, rownames(expr)[keep_expr])
cat("探索集基因: ", length(degs), " -> 表达过滤后参与 WGCNA: ", length(genes), "\n\n", sep = "")

datExpr <- t(log2expr_all[genes, ])
gsg <- goodSamplesGenes(datExpr, verbose = 0)
cat("goodSamplesGenes:", gsg$allOK, "\n")
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# ---- 软阈值 ----
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                         networkType = "signed", verbose = 0)
sftTab <- data.frame(power = sft$fitIndices$Power,
                     SFT_R2 = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
                     slope = sft$fitIndices$slope,
                     mean_k = sft$fitIndices$mean.k.)
write.csv(sftTab, file.path(V11, "wgcna_degs_softthreshold.csv"), row.names = FALSE)
cat("\n=== 软阈值 ===\n"); print(round(sftTab, 3), row.names = FALSE)

above <- sftTab$power[sftTab$SFT_R2 >= 0.8]
power <- if (length(above) > 0) min(above) else sftTab$power[which.max(sftTab$SFT_R2)]
cat("\n选定 power =", power, "\n\n")

# ---- 模块 ----
net <- blockwiseModules(datExpr, power = power,
                        networkType = "signed", TOMType = "signed",
                        minModuleSize = 40, deepSplit = 2,
                        mergeCutHeight = 0.25, numericLabels = TRUE,
                        pamRespectsDendro = FALSE, maxBlockSize = 5000,
                        verbose = 0)
moduleColors <- labels2colors(net$colors)
names(moduleColors) <- colnames(datExpr)
cat("=== 模块大小 ===\n"); print(table(moduleColors)); cat("\n")

# ---- 模块-性状 ----
MEs <- net$MEs
time_map <- c(Sham = 0, CCI_0.5d = 1, CCI_1d = 2, CCI_3d = 3, CCI_7d = 4, CCI_14d = 5)
traits <- data.frame(
  Time       = time_map[as.character(meta$group)],
  CCI_any    = as.integer(meta$group != "Sham"),
  Chronic14  = as.integer(meta$group == "CCI_14d"),
  CCI_7d     = as.integer(meta$group == "CCI_7d"),
  CCI_3d     = as.integer(meta$group == "CCI_3d"),
  row.names  = rownames(meta))
modTrait  <- cor(MEs, traits, use = "p")
modTraitP <- corPvalueStudent(modTrait, nrow(datExpr))

geneModule <- data.frame(gene = colnames(datExpr),
                         module_label = as.integer(net$colors),
                         module = moduleColors)
sizes <- table(net$colors)                      # names are "0","1",...
lab <- sub("^ME", "", rownames(modTrait))       # "1","2",...
out <- data.frame(module = rownames(modTrait),
                  color = labels2colors(as.integer(lab)),
                  n_genes = as.integer(sizes[lab]),
                  modTrait, P_Time = modTraitP[, "Time"],
                  P_Chronic14 = modTraitP[, "Chronic14"], row.names = NULL)
out <- out[order(-abs(out$Chronic14)), ]
write.csv(out, file.path(V11, "wgcna_degs_module_trait.csv"), row.names = FALSE)
cat("=== 模块-性状（按 |Chronic14| 排序）===\n")
show <- out[, c("module","color","n_genes","Time","Chronic14","P_Time","P_Chronic14")]
for (cc in c("Time","Chronic14")) show[[cc]] <- round(show[[cc]], 4)
for (cc in c("P_Time","P_Chronic14")) show[[cc]] <- signif(show[[cc]], 4)
print(show, row.names = FALSE)

# ---- 输出 ----
write.csv(geneModule, file.path(V11, "wgcna_degs_gene_modules.csv"), row.names = FALSE)
write.csv(data.frame(sample = rownames(MEs), MEs),
          file.path(V11, "wgcna_degs_eigengenes.csv"), row.names = FALSE)

# ---- 候选模块: 与时间正相关最强（渐进上调）----
chronic_lab <- as.integer(sub("^ME", "", out$module[which.max(out$Chronic14)]))
prog_lab    <- as.integer(sub("^ME", "", out$module[which.max(out$Time)]))
cand <- geneModule$gene[geneModule$module_label == prog_lab]
write.csv(data.frame(gene = cand),
          file.path(V11, "wgcna_degs_chronic_module_genes.csv"), row.names = FALSE)
cat("\n与 Chronic14 最相关模块:", out$module[which.max(out$Chronic14)],
    " (", out$color[which.max(out$Chronic14)], ", n=",
    out$n_genes[which.max(out$Chronic14)], ")\n", sep = "")
cat("与 Time 最相关(渐进上调)模块:", out$module[which.max(out$Time)],
    " (", out$color[which.max(out$Time)], ", n=", out$n_genes[which.max(out$Time)],
    ", r=", round(max(out$Time), 3), ", P=", signif(out$P_Time[which.max(out$Time)], 3), ")\n", sep = "")
cat("已输出该模块基因列表 (wgcna_degs_chronic_module_genes.csv)\n")

# ---- 与旧 M_yellow 对照 ----
old_yellow <- tryCatch(read.csv(file.path(OUT, "WGCNA_M_yellow_genes.csv"))$gene,
                       error = function(e) character(0))
if (length(old_yellow) > 0) {
  ov <- length(intersect(old_yellow, cand))
  cat("\n=== 与旧 M_yellow(265) 对照 ===\n")
  cat("新模块基因数:", length(cand), " 重叠:", ov,
      " Jaccard =", round(ov / length(union(old_yellow, cand)), 3), "\n")
}

saveRDS(list(net = net, moduleColors = moduleColors, MEs = MEs,
             modTrait = modTrait, modTraitP = modTraitP,
             power = power, genes = genes),
        file.path(V11, "wgcna_degs_result.rds"))
cat("\n完成。\n")
