# Figure 3 / Supporting Figures 1–2 statistics

The author confirms one-sided Welch's t-tests with Benjamini–Hochberg correction for group comparisons in Figure 3 and Supporting Figures 1–2. Do not replace these with two-sided tests or other procedures.

This confirmation does not specify every alternative direction or the grouping of p-values into BH families. Preserve those details from the exact final plotting code rather than guessing from star annotations. The earlier bundled plot_recovery.R defines directional exposure contrasts, with separate primary/secondary correction families; its exact match to every final panel is not yet established.

Correlation panels (e.g. unique-gene fraction versus recovery in S1) are distinct from group comparisons; the figure displays Spearman rho. The statement about Welch tests is not interpreted as replacing a correlation test with a t-test.
