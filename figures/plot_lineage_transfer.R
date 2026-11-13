#!/usr/bin/env Rscript
# Run from repository root. Plots observed values only; no significance tests.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readr)})
d <- read_csv('figures/source_data/lineage_transfer/combined_exact_recovery.csv', show_col_types=FALSE)
stopifnot(all(is.finite(d$Exact_recovery)), all(d$Exact_recovery>=0 & d$Exact_recovery<=1))
dir.create('results/lineage_transfer', recursive=TRUE, showWarnings=FALSE)
p <- ggplot(d,aes(x=Species,y=Exact_recovery,colour=Model)) +
  geom_point(position=position_dodge(width=.5),size=1.8) + facet_wrap(~Task,ncol=2) +
  scale_y_continuous(limits=c(0,1)) + labs(x=NULL,y='Exact recovery') +
  theme_bw(base_size=10) + theme(axis.text.x=element_text(angle=65,hjust=1),legend.position='bottom')
ggsave('results/lineage_transfer/22_species_overview.pdf',p,width=15,height=9)
# The 18-species analysis excludes the four model angiosperms only.
excluded <- c('A_thaliana_Araport11','S_lycopersicum_ITAG4','O_sativa_MSU7','Z_mays_APGv4',
              'Arabidopsis_thaliana','Solanum_lycopersicum','Oryza_sativa','Zea_mays')
extension <- filter(d,!Species %in% excluded)
if(n_distinct(extension$Species)!=18) stop('Species identifiers did not yield the intended 18-species subset; inspect labels.')
write_tsv(extension,'results/lineage_transfer/18_species_values.tsv')
ggsave('results/lineage_transfer/18_species_overview.pdf',p %+% extension,width=14,height=9)
# Current main Figure 6 excludes all seven model organisms, leaving 15 species.
bryophytes <- c('Anthoceros_agrestis','Marchantia_polymorpha','Physcomitrium_patens',
                'A_agrestis_v10','M_polymorpha_v61','P_patens_v33')
main_extension <- filter(extension,!Species %in% bryophytes)
if(n_distinct(main_extension$Species)!=15) stop('Expected 15 additional species for current Figure 6')
write_tsv(main_extension,'results/lineage_transfer/15_species_Figure6_values.tsv')
ggsave('results/lineage_transfer/15_species_Figure6_overview.pdf',p %+% main_extension,width=13,height=9)
