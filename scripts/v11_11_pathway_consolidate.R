# ============================================================
# v11 阶段 2 / 任务 34c：通路分析统一统计与交叉验证
#   1. PROGENy 分数的 limma 统计（Sham vs 各时间点）
#   2. GSVA(Hallmark) / GSVA(curated) / PROGENy 三方法方向一致性
#   3. 输出统一的通路结果表（供 Fig 2 使用）
# ============================================================

suppressMessages(library(limma))

BASE <- "F:/scs research"
OUT  <- file.path(BASE, "output", "v11")
TP_LAB <- c("0.5d", "1d", "3d", "7d", "14d")
GROUPS <- c("Sham", paste0("CCI_", TP_LAB))

# ---------- 1. PROGENy limma 统计 ----------
ps <- read.csv(file.path(OUT, "progeny_scores.csv"), row.names = 1, check.names = FALSE)
grp <- factor(ps$group, levels = GROUPS)
sc <- t(as.matrix(ps[, setdiff(colnames(ps), "group"), drop = FALSE]))  # pathways x samples
cat("PROGENy 通路活性矩阵:", nrow(sc), "通路 x", ncol(sc), "样本\n")

design <- model.matrix(~ 0 + grp)
colnames(design) <- levels(grp)
fit <- eBayes(lmFit(sc, design))
cm <- makeContrasts(contrasts = paste0("CCI_", TP_LAB, " - Sham"), levels = design)
fit2 <- eBayes(contrasts.fit(fit, cm))

prod <- do.call(rbind, lapply(seq_along(TP_LAB), function(j) {
  tt <- topTable(fit2, coef = j, number = Inf, sort.by = "none")
  idx <- match(rownames(sc), rownames(tt))
  data.frame(pathway = rownames(sc), timepoint = TP_LAB[j],
             delta = tt$logFC[idx], P = tt$P.Value[idx], FDR = tt$adj.P.Val[idx],
             stringsAsFactors = FALSE)
}))
prod <- prod[order(prod$FDR), ]
write.csv(prod, file.path(OUT, "progeny_delta_limma.csv"), row.names = FALSE)

dw <- reshape(prod[, c("pathway", "timepoint", "delta")],
              idvar = "pathway", timevar = "timepoint", direction = "wide")
colnames(dw) <- sub("^delta\\.", "", colnames(dw))
write.csv(dw[, c("pathway", TP_LAB)], file.path(OUT, "progeny_delta.csv"), row.names = FALSE)

cat("\nPROGENy 各时间点 FDR<0.05 通路数: ",
    paste(sprintf("%s=%d", TP_LAB, sapply(TP_LAB, function(t)
      sum(prod$FDR[prod$timepoint == t] < 0.05))), collapse = "  "), "\n")

cat("\nPROGENy 7d 变化最强:\n")
d7 <- prod[prod$timepoint == "7d", ]; d7 <- d7[order(-abs(d7$delta)), ]
for (i in head(seq_len(nrow(d7)), 8))
  cat(sprintf("  %-10s delta=%+7.3f  FDR=%.2e\n", d7$pathway[i], d7$delta[i], d7$FDR[i]))

# ---------- 2. 三方法交叉验证 ----------
gh <- read.csv(file.path(OUT, "gsva_hallmark_delta.csv"), check.names = FALSE)
gc <- read.csv(file.path(OUT, "gsva_custom_delta.csv"), check.names = FALSE)
rownames(gh) <- gh$pathway; rownames(gc) <- gc$pathway
ro <- dw$pathway
rownames(dw) <- dw$pathway

MAP <- data.frame(
  hallmark = c("IL-6/JAK/STAT3 Signaling", "TNF-alpha Signaling via NF-kB",
               "Inflammatory Response", "Complement", "Oxidative Phosphorylation",
               "Hypoxia", "Apoptosis", "p53 Pathway", "Unfolded Protein Response",
               "Glycolysis", "TGF-beta Signaling", "Wnt-beta Catenin Signaling"),
  custom   = c("JAK-STAT signaling", "TNF-a via NF-kB", "Inflammatory response",
               "Complement & coagulation", "Oxidative phosphorylation", "Hypoxia",
               "Apoptosis", "p53 pathway", "ER stress / UPR", "Glycolysis",
               NA, NA),
  progeny  = c("JAK-STAT", "TNFa", NA, NA, NA, "Hypoxia", NA, "p53", NA, NA,
               "TGFb", "WNT"),
  stringsAsFactors = FALSE)

rows <- list()
for (i in seq_len(nrow(MAP))) {
  a <- if (!is.na(MAP$hallmark[i]) && MAP$hallmark[i] %in% rownames(gh)) gh[MAP$hallmark[i], TP_LAB] else NULL
  b <- if (!is.na(MAP$custom[i])   && MAP$custom[i]   %in% rownames(gc)) gc[MAP$custom[i], TP_LAB]   else NULL
  c3 <- if (!is.na(MAP$progeny[i]) && MAP$progeny[i] %in% rownames(dw))  dw[MAP$progeny[i], TP_LAB]   else NULL
  rec <- list(hallmark = MAP$hallmark[i], custom = MAP$custom[i], progeny = MAP$progeny[i])
  for (nm in c("hallmark", "custom", "progeny")) rec[[paste0("d7_", nm)]] <- NA
  if (!is.null(a)) rec$d7_hallmark <- round(as.numeric(a["7d"]), 3)
  if (!is.null(b)) rec$d7_custom   <- round(as.numeric(b["7d"]), 3)
  if (!is.null(c3)) rec$d7_progeny <- round(as.numeric(c3["7d"]), 3)
  # 7d 符号一致性
  s <- c(sign(rec$d7_hallmark), sign(rec$d7_custom), sign(rec$d7_progeny))
  s <- s[!is.na(s)]
  rec$n_methods <- length(s)
  rec$sign_consistent_7d <- length(unique(s)) == 1
  rows[[length(rows) + 1]] <- as.data.frame(rec, stringsAsFactors = FALSE)
}
cv <- do.call(rbind, rows)
cat("\n===== 三方法 7d 方向交叉验证 =====\n")
print(cv, row.names = FALSE)
write.csv(cv, file.path(OUT, "pathway_cross_method_consistency.csv"), row.names = FALSE)

multi <- cv[cv$n_methods > 1, ]
cat(sprintf("\n可比较通路: %d 条，7d 方向完全一致: %d 条 (%.0f%%)\n",
            nrow(multi), sum(multi$sign_consistent_7d),
            100 * mean(multi$sign_consistent_7d)))

# ---------- 3. 汇总表（Fig 2 用）----------
cat("\n===== 通路结果汇总（7d, 三方法）=====\n")
gw <- gh[, c("pathway", TP_LAB)]
colnames(gw)[-1] <- paste0("hallmark_", TP_LAB)
cw <- gc[, c("pathway", TP_LAB)]
colnames(cw)[-1] <- paste0("custom_", TP_LAB)
pw <- dw[, c("pathway", TP_LAB)]
colnames(pw)[-1] <- paste0("progeny_", TP_LAB)
write.csv(gw, file.path(OUT, "pathway_summary_hallmark.csv"), row.names = FALSE)
write.csv(cw, file.path(OUT, "pathway_summary_custom.csv"), row.names = FALSE)
write.csv(pw, file.path(OUT, "pathway_summary_progeny.csv"), row.names = FALSE)

cat("\n完成。\n")
