# ============================================================
# v11 Step 1: 固化唯一口径的 DEG 结果（双轨）
#   - 主口径: log2(FPKM + 0.1), limma-trend
#   - 严格集: BH FDR < 0.05  & |log2FC| > 0.58
#   - 探索集: raw P < 0.05   & |log2FC| > 0.58
# 所有下游分析、正文数字、图注、附表均须从本脚本输出生成
# 用法: Rscript scripts/v11_01_deg_finalize.R
# ============================================================
suppressPackageStartupMessages(library(limma))

BASE <- Sys.getenv("SCS_ROOT", unset = {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(.f)) normalizePath(file.path(dirname(.f[1]), "..")) else normalizePath(".")
})
OUT  <- file.path(BASE, "output")
V11  <- file.path(OUT, "v11")
dir.create(V11, showWarnings = FALSE, recursive = TRUE)

FC_CUT <- 0.58
PSEUDO <- 0.1

expr <- read.csv(file.path(OUT, "GSE175760_expression_matrix.csv"),
                 row.names = 1, check.names = FALSE)
meta <- read.csv(file.path(OUT, "GSE175760_metadata.csv"), stringsAsFactors = FALSE)
expr <- expr[complete.cases(expr), , drop = FALSE]

meta$group <- factor(meta$group,
  levels = c("Sham","CCI_0.5d","CCI_1d","CCI_3d","CCI_7d","CCI_14d"))
rownames(meta) <- meta$sample_id
meta <- meta[colnames(expr), , drop = FALSE]

log2expr <- log2(as.matrix(expr) + PSEUDO)
write.csv(round(log2expr, 4), file.path(V11, "expr_log2_pseudo0.1.csv"))

tps <- c("0.5d","1d","3d","7d","14d")

# ---------- limma-trend ----------
design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)
fit  <- lmFit(log2expr, design)
cm   <- makeContrasts(CCI_0.5d - Sham, CCI_1d - Sham, CCI_3d - Sham,
                      CCI_7d - Sham, CCI_14d - Sham, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE)

res <- do.call(rbind, lapply(seq_along(tps), function(i) {
  tt <- topTable(fit2, coef = i, number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), tp = tps[i],
             log2FC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t,
             P_raw = tt$P.Value, P_BH = tt$adj.P.Val, B = tt$B,
             stringsAsFactors = FALSE)
}))
res$dir_strict <- ifelse(res$P_BH < 0.05 & abs(res$log2FC) > FC_CUT,
                         ifelse(res$log2FC > 0, "up", "down"), "ns")
res$dir_explor <- ifelse(res$P_raw < 0.05 & abs(res$log2FC) > FC_CUT,
                         ifelse(res$log2FC > 0, "up", "down"), "ns")
write.csv(res, file.path(V11, "limma_all_contrasts.csv"), row.names = FALSE)

# ---------- 各时间点计数汇总 ----------
summ <- do.call(rbind, lapply(tps, function(tp) {
  d <- res[res$tp == tp, ]
  data.frame(tp = tp,
    n_fc_only       = sum(abs(d$log2FC) > FC_CUT),
    n_raw_up        = sum(d$dir_explor == "up"),
    n_raw_down      = sum(d$dir_explor == "down"),
    n_raw_total     = sum(d$dir_explor != "ns"),
    n_bh_up         = sum(d$dir_strict == "up"),
    n_bh_down       = sum(d$dir_strict == "down"),
    n_bh_total      = sum(d$dir_strict != "ns"),
    stringsAsFactors = FALSE)
}))
uni <- function(col) {
  s <- unique(res$gene[res[[col]] != "ns"]); s
}
summ$union_explor <- length(uni("dir_explor"))
summ$union_strict <- length(uni("dir_strict"))
write.csv(summ, file.path(V11, "DEG_summary_by_timepoint.csv"), row.names = FALSE)

# ---------- 双轨基因列表 ----------
write.csv(data.frame(gene = uni("dir_strict")),
          file.path(V11, "DEG_strict_genes.csv"), row.names = FALSE)
write.csv(data.frame(gene = uni("dir_explor")),
          file.path(V11, "DEG_exploratory_genes.csv"), row.names = FALSE)
write.csv(res[res$dir_explor != "ns", ],
          file.path(V11, "DEG_exploratory_long.csv"), row.names = FALSE)
write.csv(res[res$dir_strict != "ns", ],
          file.path(V11, "DEG_strict_long.csv"), row.names = FALSE)

# ---------- 重点基因（广义）统计: 旧22 + 免疫核心 + SCS候选 ----------
FOCUS <- c("Il10rb","Jak3","Cd37","Per1","Dera","Dbr1","Lpxn","Nat9","Pbx2","Fbxw4",
           "Edem1","Nrp1","Itga7","Casp3","Svop","Mrpl10","Mtif3","Smc2","Akap13",
           "Emilin2","Apool","Abhd1")
FOCUS <- FOCUS[FOCUS %in% rownames(expr)]
fs <- res[res$gene %in% FOCUS, ]
fs <- fs[order(fs$gene, match(fs$tp, tps)), ]
write.csv(fs, file.path(V11, "focus_genes_limma_stats.csv"), row.names = FALSE)

nsig_raw <- sapply(FOCUS, function(g) sum(fs$gene == g & fs$dir_explor != "ns"))
nsig_bh  <- sapply(FOCUS, function(g) sum(fs$gene == g & fs$dir_strict != "ns"))
summary_focus <- data.frame(
  gene = FOCUS,
  n_tp_sig_raw = nsig_raw,
  n_tp_sig_bh  = nsig_bh,
  max_abs_log2FC = sapply(FOCUS, function(g) round(max(abs(fs$log2FC[fs$gene == g])), 3)),
  min_P_raw   = sapply(FOCUS, function(g) signif(min(fs$P_raw[fs$gene == g]), 3)),
  min_P_BH    = sapply(FOCUS, function(g) signif(min(fs$P_BH[fs$gene == g]), 3)),
  row.names = NULL
)
write.csv(summary_focus, file.path(V11, "focus_genes_summary.csv"), row.names = FALSE)

# ============================================================
cat(strrep("=", 62), "\n")
cat("v11 双轨 DEG 结果已固化\n")
cat(strrep("=", 62), "\n\n")
cat("口径: log2(FPKM +", PSEUDO, ") + limma-trend;  |log2FC| >", FC_CUT, "\n\n")
print(summ[, c("tp","n_fc_only","n_raw_total","n_bh_total")], row.names = FALSE)
cat("\nunion 探索集 (raw P<0.05 & FC): ", length(uni("dir_explor")), "\n", sep = "")
cat("union 严格集 (BH<0.05  & FC): ", length(uni("dir_strict")), "\n", sep = "")
cat("\n重点基因显著性计数:\n")
print(summary_focus[, c("gene","n_tp_sig_raw","n_tp_sig_bh","max_abs_log2FC","min_P_BH")],
      row.names = FALSE)
cat("\n至少1时间点 raw 显著: ", sum(summary_focus$n_tp_sig_raw > 0), "/", nrow(summary_focus), "\n", sep = "")
cat("至少1时间点 BH  显著: ", sum(summary_focus$n_tp_sig_bh  > 0), "/", nrow(summary_focus), "\n", sep = "")
cat("\n输出目录: output/v11/\n")
