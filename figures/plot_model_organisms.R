#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2); library(readr)})
d <- read_tsv('results/model_organisms/exact_recovery_long.tsv', show_col_types=FALSE)
p <- ggplot(d,aes(x=model,y=species,fill=exact_recovery)) + geom_tile() +
  facet_wrap(~task) + scale_fill_gradient(limits=c(0,1),low='white',high='#246847') +
  labs(x=NULL,y=NULL,fill='Exact recovery') + theme_minimal(base_size=10) +
  theme(axis.text.x=element_text(angle=60,hjust=1),panel.grid=element_blank())
ggsave('results/model_organisms/exact_recovery_overview.pdf',p,width=12,height=7)
