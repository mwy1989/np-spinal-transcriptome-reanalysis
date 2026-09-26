# ============================================================
# v11 Step 5: 时间轨迹软聚类（fuzzy c-means）
#   说明: Mfuzz 的核心算法即 e1071::cmeans + 基因级标准化；
#         此处直接调用同样算法，方法学描述与实现一致。
#   输入: 探索集基因 x 6 时间点组均值
#   输出: output/v11/cluster_*
# 用法: Rscript scripts/v11_05_temporal_clusters.R
# ============================================================
suppressPackageStartupMessages(library(e1071))

BASE <- Sys.getenv("SCS_ROOT", unset = {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(.f)) normalizePath(file.path(dirname(.f[1]), "..")) else normalizePath(".")
}); OUT <- file.path(BASE,"output"); V11 <- file.path(OUT,"v11")

expr <- read.csv(file.path(OUT,"GSE175760_expression_matrix.csv"),
                 row.names=1, check.names=FALSE)
meta <- read.csv(file.path(OUT,"GSE175760_metadata.csv"), stringsAsFactors=FALSE)
expr <- expr[complete.cases(expr), , drop=FALSE]
genes <- read.csv(file.path(V11,"DEG_exploratory_genes.csv"))$gene
genes <- intersect(genes, rownames(expr))

log2expr <- log2(as.matrix(expr) + 0.1)
TPL <- c("Sham","CCI_0.5d","CCI_1d","CCI_3d","CCI_7d","CCI_14d")
tm <- sapply(TPL, function(g) rowMeans(log2expr[, meta$sample_id[meta$group==g], drop=FALSE]))
tm <- tm[genes, ]
cat("输入基因:", nrow(tm), " 时间点:", ncol(tm), "\n")

# ---- 基因级 z-score ----
z <- t(scale(t(tm)))
ok <- rowSums(is.na(z)) == 0 & apply(z, 1, sd) > 0
z <- z[ok, ]
cat("标准化后可用基因:", nrow(z), "\n\n")

results <- list()
for (K in c(4, 6)) {
  set.seed(42)
  cm <- cmeans(z, centers = K, m = 1.25, iter.max = 200, method = "cmeans")
  results[[as.character(K)]] <- cm
  cat(strrep("=",64), "\n")
  cat("c-means  K =", K, "   iter =", cm$iter, "   within.error =", round(cm$withinerror,2), "\n")
  cat(strrep("=",64), "\n")
  for (k in seq_len(K)) {
    idx <- cm$cluster == k
    prof <- colMeans(tm[names(cm$cluster)[idx], , drop=FALSE])  # 用原始 log2 均值画轨迹
    peak_tp <- TPL[which.max(prof)]
    cat(sprintf("  簇 %d  n=%4d  峰时间=%-8s  轨迹: %s\n", k, sum(idx), peak_tp,
                paste(sprintf("%+.2f", prof), collapse=" -> ")))
  }
  cat("\n")
  # 保存
  asg <- data.frame(gene = names(cm$cluster), cluster = cm$cluster,
                    membership = apply(cm$membership, 1, max))
  write.csv(asg, file.path(V11, sprintf("cluster_k%d_assignments.csv", K)), row.names=FALSE)
  prof_df <- do.call(rbind, lapply(seq_len(K), function(k) {
    idx <- cm$cluster == k
    data.frame(cluster = k, tp = TPL,
               mean_log2 = colMeans(tm[names(cm$cluster)[idx], , drop=FALSE]),
               z_profile = colMeans(z[names(cm$cluster)[idx], , drop=FALSE]),
               n_genes = sum(idx))
  }))
  write.csv(prof_df, file.path(V11, sprintf("cluster_k%d_profiles.csv", K)), row.names=FALSE)
  # 每簇代表基因（隶属度最高）
  reps <- do.call(rbind, lapply(seq_len(K), function(k) {
    idx <- names(cm$cluster)[cm$cluster == k]
    m <- cm$membership[idx, k]
    top <- names(sort(m, decreasing = TRUE))[1:10]
    data.frame(cluster = k, rank = 1:10, gene = top,
               membership = round(as.numeric(m[top]), 3))
  }))
  write.csv(reps, file.path(V11, sprintf("cluster_k%d_representatives.csv", K)), row.names=FALSE)
}

# ---- 与 WGCNA 模块对照 (K=4) ----
asg <- read.csv(file.path(V11,"cluster_k4_assignments.csv"))
gm  <- read.csv(file.path(V11,"wgcna_degs_gene_modules.csv"))
mg  <- merge(asg, gm, by="gene")
cat("=== 软聚类簇 vs WGCNA 模块 交叉表 ===\n")
print(table(mg$cluster, mg$module))
cat("\n完成。输出目录: output/v11/\n")
