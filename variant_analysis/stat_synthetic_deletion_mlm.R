#!/usr/bin/env Rscript

library(ggplot2)

df <- read.delim("synthetic_deletion_core_scored_mlm.tsv", stringsAsFactors = FALSE)

score_col <- "delta_pll_mut_minus_ref"
#df$score_drop_mut_minus_ref <- -df$score_drop_ref_minus_mut

df$sv_type <- factor(
  df$sv_type,
  levels = c(
    "Deep_intron_deletion",
    "TIS_upstream_deletion",
    "CDS_inframe_deletion",
    "CDS_frameshift_deletion"
  )
)

# -----------------------------
# Summary statistics
# -----------------------------
sv_types <- levels(df$sv_type)
out <- data.frame()

for (sv in sv_types) {
  x <- df[df$sv_type == sv, score_col]
  x <- x[!is.na(x)]

  n <- length(x)
  n_positive <- sum(x > 0)
  frac_positive <- n_positive / n

  wil <- wilcox.test(
    x,
    mu = 0,
    alternative = "greater",
    exact = FALSE
  )

  bin <- binom.test(
    n_positive,
    n,
    p = 0.5,
    alternative = "greater"
  )

  out <- rbind(
    out,
    data.frame(
      sv_type = sv,
      n = n,
      mean = mean(x),
      median = median(x),
      sd = sd(x),
      min = min(x),
      q25 = quantile(x, 0.25),
      q75 = quantile(x, 0.75),
      max = max(x),
      n_positive = n_positive,
      frac_positive = frac_positive,
      wilcoxon_p = wil$p.value,
      sign_test_p = bin$p.value
    )
  )
}

out$wilcoxon_p_adj_BH <- p.adjust(out$wilcoxon_p, method = "BH")
out$sign_test_p_adj_BH <- p.adjust(out$sign_test_p, method = "BH")

write.table(
  out,
  file = "synthetic_deletion_core_score_statistics_mlm.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(out)

# -----------------------------
# Violin + boxplot
# -----------------------------
# -----------------------------
# Violin plot only
# -----------------------------
p <- ggplot(df, aes(x = sv_type, y = delta_pll_mut_minus_ref, fill = sv_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_violin(trim = FALSE, alpha = 0.75, color = "black", linewidth = 0.3) +
  geom_jitter(
    width = 0.12,
    size = 0.18,
    alpha = 0.06,
    color = "black"
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 2.8,
    fill = "white",
    color = "black"
  ) +
  labs(
    title = "Synthetic deletion score drop",
    x = NULL,
    y = "Score drop: ref junction log-probability − mutant junction log-probability"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(
  "synthetic_deletion_score_violin_mlm.pdf",
  p,
  width = 7,
  height = 4.5
)

