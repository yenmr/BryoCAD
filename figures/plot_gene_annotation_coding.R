#!/usr/bin/env Rscript
# Coding-only metrics for Figures 5A / 5D; no UTR-sensitive full_tx substitution.
suppressPackageStartupMessages(library(ggplot2))
out <- 'results/gene_annotation_coding'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
select_coding <- function(path) {
  d <- read.delim(path,check.names=FALSE,na.strings='NA')
  d <- d[d$section=='CDS' & d$metric %in% c('cds_ex','tx') & d$status=='available',]
  d$Metric <- factor(d$metric,levels=c('cds_ex','tx'),labels=c('CDS segments','Coding transcripts (UTRs excluded)'))
  if(!nrow(d)) stop('No available coding metrics: ',path)
  d
}
a <- select_coding('figures/source_data/gene_annotation/comparison.metrics.tsv')
write.table(a,file.path(out,'Figure5A_coding_source.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
p <- ggplot(a,aes(recall,precision,colour=Model,shape=Head)) + geom_point(size=2.5) +
  facet_grid(Metric~Species) + coord_equal(xlim=c(0,1),ylim=c(0,1)) + theme_bw() +
  labs(x='Recall',y='Precision') + theme(legend.position='bottom')
ggsave(file.path(out,'Figure5A_coding_overview.pdf'),p,width=11,height=6)
d <- select_coding('figures/source_data/pangenome_candidates/pangenome.metrics.tsv')
if(length(unique(d$short_name))!=123) stop('Expected 123 pangenome assemblies')
if(anyDuplicated(d[c('short_name','Model','Head','metric')])) stop('Duplicate assembly/model/head/metric')
write.table(d,file.path(out,'Figure5D_coding_source.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
p <- ggplot(d,aes(division,f1,colour=division)) + geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=.12,height=0,alpha=.4,size=.8) + facet_wrap(~Metric) +
  scale_y_continuous(limits=c(0,1)) + theme_bw() +
  labs(x=NULL,y='F1',caption='123 assemblies; selected whole chromosomes/contigs; approximately 50 Mb target') +
  theme(legend.position='none')
set.seed(42)
ggsave(file.path(out,'Figure5D_coding_overview.pdf'),p,width=10,height=5)
