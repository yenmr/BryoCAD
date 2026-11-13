#!/usr/bin/env Rscript
# v2: PDF-only, species-level tests, count-pooled gene subsets, source TSVs.
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(readr)})
a <- commandArgs(TRUE)
if(!length(a) %in% c(2,3)) stop('Usage: Rscript plot_recovery.R processed_directory figure_directory [gene_category_counts.tsv]')
input <- a[1]; output <- a[2]; dir.create(output,recursive=TRUE,showWarnings=FALSE)
rd <- function(n) read_tsv(file.path(input,n),show_col_types=FALSE,na=c('NA',''))
s <- rd('species_motif_summary.tsv'); g <- rd('species_gene_category_summary.tsv'); meta <- rd('species_metadata.tsv'); t <- rd('tts_codon_confusion.tsv')
set.seed(20261008)
# Display-only cleanup; original IDs remain intact in data and source TSVs.
species_label <- function(x) sub('_[vV][0-9]+([.][0-9]+)*$', '', x)
groups <- c('trained','unseen_species_seen_order','unseen_order')
subsets <- c('all','core_accessory','core')
basekeys <- c('short_name','species','division','order','training_group','motif')
pool <- function(cats,label) g %>% filter(gene_category %in% cats) %>% group_by(across(all_of(basekeys))) %>%
  summarise(n_sites=sum(n_sites),n_exact=sum(n_exact),.groups='drop') %>% mutate(exact_recovery_rate=n_exact/n_sites,subset=label)
x <- bind_rows(s %>% mutate(subset='all'),pool(c('core','accessory'),'core_accessory'),pool('core','core')) %>%
  mutate(subset=factor(subset,levels=subsets),training_group=factor(training_group,levels=groups),motif=factor(motif,levels=c('acceptor','donor','tis','tts')))
if(!all(c('core','accessory') %in% g$gene_category)) stop('Expected exact category labels core and accessory')
write_tsv(x,file.path(output,'all_analysis_subsets.tsv'))
theme_set(theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),strip.background=element_rect(fill='grey95')))
export <- function(p,name,data,w=12,h=6,stats=NULL) {
  write_tsv(data,file.path(output,paste0(name,'.tsv')))
  if(!is.null(stats)) write_tsv(stats,file.path(output,paste0(name,'_statistics.tsv')))
  ggsave(file.path(output,paste0(name,'.pdf')),p,width=w,height=h,limitsize=FALSE)
  # ggsave(file.path(output,paste0(name,'.png')),p,width=w,height=h,dpi=300,limitsize=FALSE)
}
boot <- function(z) {
  z<-z[is.finite(z)]; n<-length(z)
  if(n<2) return(tibble(n_species=n,mean=mean(z),lower=NA_real_,upper=NA_real_))
  b<-replicate(2000,mean(sample(z,n,TRUE)))
  tibble(n_species=n,mean=mean(z),lower=unname(quantile(b,.025)),upper=unname(quantile(b,.975)))
}
# Directional hypothesis: trained > unseen. Primary contrast pools unseen species.
comparisons <- list(c('trained','unseen_all'),c('trained','unseen_species_seen_order'),c('trained','unseen_order'))
pairtest <- function(d) bind_rows(lapply(comparisons,function(pair) {
  u<-d$exact_recovery_rate[d$training_group=='trained']
  v<-d$exact_recovery_rate[if(pair[2]=='unseen_all') d$training_group!='trained' else d$training_group==pair[2]]
  fit<-if(length(u)>=2 && length(v)>=2) tryCatch(t.test(u,v,alternative='greater',var.equal=FALSE),error=function(e) NULL) else NULL
  tibble(group1=pair[1],group2=pair[2],n1=length(u),n2=length(v),mean_difference_pp=100*(mean(u)-mean(v)),p=if(is.null(fit)) NA_real_ else fit$p.value)
}))
tests <- x %>% group_by(subset,motif) %>% group_modify(~pairtest(.x)) %>% ungroup() %>% mutate(family=if_else(group2=='unseen_all','primary_pooled','secondary_stratified')) %>% group_by(family) %>% mutate(q=p.adjust(p,'BH'),test='one-sided Welch t-test: trained > unseen; species unit') %>% ungroup()
# One BH family across all subset/motif/group contrasts.
write_tsv(tests,file.path(output,'training_group_tests.tsv'))
ci <- x %>% group_by(subset,motif,training_group) %>% group_modify(~boot(.x$exact_recovery_rate)) %>% ungroup()
write_tsv(ci,file.path(output,'species_macro_mean_CI.tsv'))
short <- c(trained='Trained',unseen_species_seen_order='Unseen species\nseen order',unseen_order='Unseen order')
for(ss in subsets) {
  d<-filter(x,subset==ss); st<-filter(tests,subset==ss)
  ann<-st %>% group_by(motif) %>% mutate(y=1.04+.065*(row_number()-1),label=paste0(ifelse(group2=='unseen_all','Trained > all unseen',ifelse(group2=='unseen_order','Trained > unseen order','Trained > seen-order unseen')),': q=',ifelse(is.na(q),'NA',format.pval(q,digits=2,eps=.001)))) %>% ungroup()
  p<-ggplot(d,aes(training_group,exact_recovery_rate))+geom_boxplot(outlier.shape=NA,width=.55,fill='grey95')+geom_point(aes(color=division),position=position_jitter(width=.14,height=0,seed=3),alpha=.7)+
    geom_text(data=ann,aes(x=2,y=y,label=label),inherit.aes=FALSE,size=2.6)+facet_wrap(~motif,nrow=1)+
    scale_x_discrete(labels=short)+scale_y_continuous(breaks=seq(0,1,.2))+coord_cartesian(ylim=c(0,1.21))+
    labs(x=NULL,y='Exact recovery',title=paste('Gene subset:',ss),subtitle='One-sided Welch tests: trained > unseen; BH q. Primary: pooled unseen; secondary: each unseen group.')
  export(p,paste0('02_training_groups_',ss),d,14,5,st)

}
# Gene categories: absolute rates, fixed biological order; matched species tests.
catdata<-g %>% filter(gene_category %in% c('core','accessory','unique')) %>% mutate(gene_category=factor(gene_category,levels=c('core','accessory','unique')))
export(ggplot(catdata,aes(gene_category,exact_recovery_rate))+geom_boxplot(outlier.shape=NA,width=.55,fill='grey95')+
  geom_point(aes(color=division),position=position_jitter(width=.12,height=0,seed=3),alpha=.6,size=1.3)+
  facet_wrap(~motif,nrow=1)+labs(x=NULL,y='Exact recovery',title='Core / accessory / unique'),
  '04_gene_categories',catdata,14,5)
# Unique fraction among evaluated genes, matched to each overall motif result.
uf<-g %>% group_by(short_name,motif) %>% summarise(
  n_unique_genes=sum(n_genes[gene_category=='unique']),
  n_classified_genes=sum(n_genes[gene_category %in% c('core','accessory','unique')]),
  n_total_genes=sum(n_genes),.groups='drop') %>%
  mutate(unique_gene_fraction=n_unique_genes/n_total_genes,
         classified_gene_fraction=n_classified_genes/n_total_genes) %>%
  inner_join(s,by=c('short_name','motif')) %>% filter(n_total_genes>0)
if(any(uf$n_total_genes!=uf$n_genes)) stop('Gene counts disagree between category and overall summary')
ufstats<-uf %>% group_by(motif) %>% group_modify(~{
  fit<-if(nrow(.x)>=4 && n_distinct(.x$unique_gene_fraction)>1 && n_distinct(.x$exact_recovery_rate)>1)
    suppressWarnings(cor.test(.x$unique_gene_fraction,.x$exact_recovery_rate,method='spearman',exact=FALSE)) else NULL
  tibble(n_species=nrow(.x),rho=if(is.null(fit)) NA_real_ else unname(fit$estimate),p=if(is.null(fit)) NA_real_ else fit$p.value)
}) %>% ungroup() %>% mutate(q=p.adjust(p,'BH'),label=sprintf('n=%d; rho=%.2f; q=%.3g',n_species,rho,q))
export(ggplot(uf,aes(unique_gene_fraction,exact_recovery_rate))+
  geom_point(aes(color=division),alpha=.75,size=1.8)+
  geom_smooth(method='lm',formula=y~x,color='black',linewidth=.5,se=TRUE)+
  geom_text(data=ufstats,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.03,vjust=1.3,size=3)+
  facet_wrap(~motif,nrow=1)+coord_cartesian(ylim=c(0,1))+
  labs(x='Unique genes / all evaluated genes (within species and motif)',y='Overall exact recovery',
       subtitle='One point per species; Spearman rho and BH q; OLS trend with 95% CI'),
  '04d_unique_gene_fraction_vs_recovery',uf,14,5,ufstats)
# Supplementary composition of evaluated genes, separately by motif.
# Counts cannot be added across motifs; unknown mappings remain in denominators.
composition_genes<-g %>% group_by(short_name,species,division,order,training_group,motif) %>%
  mutate(total_evaluated_genes=sum(n_genes),gene_percentage=100*n_genes/total_evaluated_genes) %>% ungroup()
write_tsv(composition_genes,file.path(output,'04b_evaluated_gene_composition_all_categories.tsv'))
# Fill absent known categories with zero only for species/motif pairs with data.
comp<-composition_genes %>% select(short_name,species,division,order,training_group,motif,total_evaluated_genes,gene_category,n_genes,gene_percentage) %>%
  group_by(short_name,species,division,order,training_group,motif,total_evaluated_genes) %>%
  complete(gene_category=union(as.character(gene_category),c('core','accessory','unique')),fill=list(n_genes=0,gene_percentage=0)) %>% ungroup()
known<-comp %>% filter(gene_category %in% c('core','accessory','unique')) %>%
  mutate(gene_category=factor(gene_category,levels=c('core','accessory','unique')))
for(measure in c('n_genes','gene_percentage')) {
  export(ggplot(known,aes(gene_category,.data[[measure]]))+geom_boxplot(outlier.shape=NA,fill='grey95')+
    geom_point(aes(color=division),position=position_jitter(width=.13,height=0,seed=8),alpha=.55,size=1)+
    facet_grid(motif~division)+labs(x=NULL,y=if(measure=='n_genes') 'Number of evaluated genes' else 'Percent of all evaluated genes',
      subtitle='Per species and motif; unmapped/ambiguous genes remain in percentage denominator'),
    paste0('04b_gene_composition_',measure),known,12,10)
}
comp$gene_category<-factor(comp$gene_category,levels=c('core','accessory','unique',sort(setdiff(unique(comp$gene_category),c('core','accessory','unique')))))
sporder<-meta %>% arrange(division,order,short_name) %>% pull(short_name)
comp$short_name<-factor(comp$short_name,levels=rev(sporder))
export(ggplot(comp,aes(gene_percentage,short_name,fill=gene_category))+scale_y_discrete(labels=species_label)+geom_col(width=.85)+facet_wrap(~motif,nrow=1)+
  labs(x='Percent of evaluated genes',y=NULL,fill='Gene category')+theme(axis.text.y=element_text(size=5)),
  '04c_species_gene_composition',comp,14,max(10,length(sporder)*.12))
# Whole orthogroup-index gene counts; do not sum motif-specific gene counts.
gene_count_file<-if(length(a)==3) a[3] else file.path(input,'gene_category_counts.tsv')
if(!file.exists(gene_count_file)) stop('Missing gene_category_counts.tsv. Run export_gene_counts.py on orthogroup_index.sqlite first.')
gc<-read_tsv(gene_count_file,show_col_types=FALSE)
if(anyDuplicated(gc[c('short_name','gene_category')])) stop('Duplicated species/category gene count rows')
gc<-gc %>% filter(gene_category %in% c('core','accessory','unique')) %>%
  group_by(short_name) %>% complete(gene_category=c('core','accessory','unique'),fill=list(n_genes=0)) %>% ungroup() %>%
  left_join(meta %>% select(short_name,species,division,order,training_group),by='short_name') %>%
  mutate(gene_category=factor(gene_category,levels=c('core','accessory','unique')))
if(anyNA(gc$division)) stop('Gene-count species missing metadata')
if(any(!meta$short_name %in% gc$short_name)) stop('Some metadata species have no gene counts; check OG mapping coverage')
cols<-c(core='#4477AA',accessory='#EEAA33',unique='#AA3377')
export(ggplot(gc,aes(gene_category,n_genes,fill=gene_category))+geom_boxplot(outlier.shape=NA,width=.6)+
  geom_point(position=position_jitter(width=.12,height=0,seed=6),size=1.2,alpha=.6)+facet_wrap(~division)+
  scale_fill_manual(values=cols)+labs(x=NULL,y='Number of classified genes',subtitle='Orthogroup-index counts; one point per species; absolute counts')+
  theme(legend.position='none'), '04b_gene_counts_by_division',gc,11,5)
gc<-gc %>% mutate(short_name=factor(short_name,levels=meta %>% arrange(division,order,short_name) %>% pull(short_name)))
export(ggplot(gc,aes(short_name,n_genes,fill=gene_category))+scale_x_discrete(labels=species_label)+geom_col(position=position_stack(reverse=TRUE),width=.85)+
  facet_grid(~division,scales='free_x',space='free_x')+scale_fill_manual(values=cols)+
  labs(x=NULL,y='Number of classified genes',fill='Gene category',subtitle='Stacked counts: core at base, then accessory, then unique; taxonomic order')+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=5.5)), '04c_gene_counts_by_species',gc,22,7)
count_summary<-gc %>% group_by(division,gene_category) %>% summarise(n_species=n(),median=median(n_genes),q25=quantile(n_genes,.25),q75=quantile(n_genes,.75),min=min(n_genes),max=max(n_genes),.groups='drop')
write_tsv(count_summary,file.path(output,'04b_gene_counts_by_division_summary.tsv'))
# Independent lineage comparison, two-sided Welch tests (no directional hypothesis).
divtests<-x %>% group_by(subset,motif) %>% group_modify(~{
  lev<-sort(unique(.x$division))
  bind_rows(lapply(combn(lev,2,simplify=FALSE),function(pair){
    u<-.x$exact_recovery_rate[.x$division==pair[1]];v<-.x$exact_recovery_rate[.x$division==pair[2]]
    fit<-if(length(u)>=2 && length(v)>=2) tryCatch(t.test(u,v),error=function(e) NULL) else NULL
    tibble(group1=pair[1],group2=pair[2],n1=length(u),n2=length(v),mean_difference_pp=100*(mean(u)-mean(v)),p=if(is.null(fit)) NA_real_ else fit$p.value)
  }))
}) %>% ungroup() %>% mutate(q=p.adjust(p,'BH'),test='two-sided Welch t-test; species unit')
for(ss in subsets) {
  d<-filter(x,subset==ss); st<-filter(divtests,subset==ss)
  ann<-st %>% group_by(motif) %>% mutate(y=1.04+.065*(row_number()-1),label=paste0(group1,' vs ',group2,': q=',ifelse(is.na(q),'NA',format.pval(q,digits=2,eps=.001)))) %>% ungroup()
  export(ggplot(d,aes(division,exact_recovery_rate))+geom_boxplot(outlier.shape=NA,fill='grey95')+
    geom_point(aes(color=training_group),position=position_jitter(width=.12,height=0,seed=4),alpha=.7)+
    geom_text(data=ann,aes(x=2,y=y,label=label),inherit.aes=FALSE,size=2.4)+facet_wrap(~motif,nrow=1)+
    coord_cartesian(ylim=c(0,1.21))+labs(x=NULL,y='Exact recovery',title=paste('Lineage comparison:',ss)),
    paste0('05_division_',ss),d,14,5,st)
}
# Supplementary species overview with metadata in source table.
ord<-meta %>% arrange(division,order,short_name) %>% pull(short_name)
hm<-s %>% mutate(short_name=factor(short_name,levels=rev(ord)))
export(ggplot(hm,aes(motif,short_name,fill=exact_recovery_rate))+scale_y_discrete(labels=species_label)+geom_tile()+facet_grid(division~.,scales='free_y',space='free_y')+
  scale_fill_viridis_c(limits=c(0,1))+theme(axis.text.y=element_text(size=6))+labs(x=NULL,y=NULL),
  '01_species_supplementary',hm,7,max(8,length(ord)*.13))
# Genomic associations. OLS curve/95% CI is descriptive; Spearman evaluates monotonic association.
z<-x %>% left_join(meta %>% select(short_name,starts_with('genome_')),by='short_name')
vars<-c('genome_assembly_span_bp','genome_gc_fraction','genome_coding_genes_per_Mb')
assoc<-list(); datasets<-list()
for(v in intersect(vars,names(z))) for(ss in subsets) for(mode in 'all_species') {
  if(mode!='all_species' && v!='genome_coding_genes_per_Mb') next
  d<-filter(z,subset==ss);d$genome_value<-suppressWarnings(as.numeric(d[[v]]))
  d<-filter(d,is.finite(genome_value),is.finite(exact_recovery_rate))
  if(mode!='all_species') d<-filter(d,as.numeric(genome_gff_target_missing_sequence)==0)
  if(v=='genome_assembly_span_bp') d<-filter(d,genome_value>0)
  d$plot_x<-if(v=='genome_assembly_span_bp') log10(d$genome_value) else d$genome_value
  if(!nrow(d)) next
  id<-paste(v,ss,mode,sep='__');datasets[[id]]<-d
  assoc[[id]]<-d %>% group_by(motif) %>% group_modify(~{
    ok<-nrow(.x)>=3 && n_distinct(.x$plot_x)>1 && n_distinct(.x$exact_recovery_rate)>1
    fit<-if(ok) suppressWarnings(cor.test(.x$plot_x,.x$exact_recovery_rate,method='spearman',exact=FALSE)) else NULL
    tibble(n_species=nrow(.x),rho=if(ok) unname(fit$estimate) else NA_real_,p=if(ok) fit$p.value else NA_real_)
  }) %>% ungroup() %>% mutate(id=id,variable=v,subset=ss,sensitivity=mode)
}
at<-bind_rows(assoc) %>% mutate(q=p.adjust(p,'BH'))
write_tsv(at,file.path(output,'genome_association_tests.tsv'))
for(id in names(datasets)) {
  d<-datasets[[id]];st<-filter(at,.data$id==!!id)
  st$label<-sprintf('n=%d; rho=%.2f; q=%s',st$n_species,st$rho,ifelse(is.na(st$q),'NA',format.pval(st$q,digits=2,eps=.001)))
  p<-ggplot(d,aes(plot_x,exact_recovery_rate))+geom_point(aes(color=division),alpha=.7)+
    geom_smooth(method='lm',formula=y~x,se=TRUE,color='black',linewidth=.5)+
    geom_text(data=st,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.03,vjust=1.3,size=3)+
    facet_wrap(~motif)+labs(x=if(st$variable[1]=='genome_assembly_span_bp') 'log10(assembly span bp)' else st$variable[1],y='Exact recovery',title=paste(st$subset[1],st$sensitivity[1]),subtitle='OLS line with 95% CI; Spearman correlation, BH-adjusted q')
  export(p,paste0(if(st$sensitivity[1]=='all_species') '09_' else '10_',id),d,11,7,st)
}
# TTS retains all-gene scope; existing confusion table has no OG stratification.
if(nrow(t)) {
  tt<-s %>% filter(motif=='tts') %>% pivot_longer(c(exact_recovery_rate,canonical_recovery_rate),names_to='metric',values_to='rate')
  export(ggplot(tt,aes(metric,rate,group=short_name))+geom_line(alpha=.2)+geom_point(aes(color=division))+facet_wrap(~training_group)+labs(x=NULL,y='TTS recovery'), '06_tts_recovery',tt,13,5)
  usage<-t %>% group_by(short_name,division,order,reference_motif) %>% summarise(n=sum(n_sites),.groups='drop') %>%
    group_by(short_name,division,order) %>% complete(reference_motif=c('TAA','TAG','TGA'),fill=list(n=0)) %>% mutate(proportion=n/sum(n)) %>% ungroup()
  write_tsv(usage,file.path(output,'07_tts_order_usage_species.tsv'))
  ou<-usage %>% group_by(division,order,reference_motif) %>% summarise(proportion=mean(proportion),n_species=n_distinct(short_name),.groups='drop')
  export(ggplot(ou,aes(order,proportion,fill=reference_motif))+geom_col()+coord_flip()+labs(x=NULL,y='Species mean reference proportion'),'07_tts_order_usage',ou,8,max(5,n_distinct(t$order)*.23))
  export(ggplot(usage,aes(reference_motif,proportion))+geom_boxplot(outlier.shape=NA)+
    geom_point(aes(color=division),position=position_jitter(width=.12,height=0,seed=5),alpha=.65)+
    facet_wrap(~division)+labs(x='Reference stop codon',y='Within-species reference frequency'),
    '07b_reference_codon_distribution',usage,11,5)
  # Denominator: all evaluated sites with this reference codon in each species.
  codon_recovery<-t %>% group_by(short_name,division,order,training_group,reference_motif) %>%
    summarise(n_exact=sum(n_sites[predicted_motif==reference_motif]),
      n_canonical=sum(n_sites[predicted_motif %in% c('TAA','TAG','TGA')]),n_sites=sum(n_sites),.groups='drop') %>%
    mutate(exact_recovery=n_exact/n_sites,canonical_recovery=n_canonical/n_sites,
      reference_motif=factor(reference_motif,levels=c('TAA','TAG','TGA'))) %>%
    pivot_longer(c(exact_recovery,canonical_recovery),names_to='metric',values_to='recovery_rate') %>%
    mutate(metric=factor(metric,levels=c('exact_recovery','canonical_recovery'),labels=c('Exact recovery','Canonical recovery')))
  export(ggplot(codon_recovery,aes(reference_motif,recovery_rate))+
    geom_boxplot(outlier.shape=NA,fill='grey95',width=.6)+
    geom_point(aes(color=division),position=position_jitter(width=.12,height=0,seed=5),alpha=.65)+
    facet_grid(metric~division)+scale_y_continuous(limits=c(0,1))+
    labs(x='Reference stop codon',y='Recovery rate',subtitle='One point per species; denominator includes all predictions for the reference codon'),
    '07c_reference_codon_recovery',codon_recovery,12,8)
  cf<-t %>% mutate(predicted_class=if_else(predicted_motif %in% c('TAA','TAG','TGA'),predicted_motif,'Other')) %>%
    group_by(short_name,division,order,reference_motif,predicted_class) %>% summarise(n=sum(n_sites),.groups='drop') %>%
    group_by(short_name,division,order,reference_motif) %>% complete(predicted_class=c('TAA','TAG','TGA','Other'),fill=list(n=0)) %>% mutate(proportion=n/sum(n)) %>% ungroup()
  write_tsv(cf,file.path(output,'08_tts_conditional_order_species.tsv'))
  # Main matrix: rows reference, columns prediction; species-equal conditional rates.
  cm<-cf %>% group_by(reference_motif,predicted_class) %>% summarise(proportion=mean(proportion),n_species=n_distinct(short_name),.groups='drop') %>%
    mutate(reference_motif=factor(reference_motif,levels=c('TAA','TAG','TGA')),predicted_class=factor(predicted_class,levels=c('TAA','TAG','TGA','Other')))
  export(ggplot(cm,aes(predicted_class,reference_motif,fill=proportion))+geom_tile()+geom_text(aes(label=sprintf('%.1f%%',100*proportion)))+
    scale_fill_viridis_c(limits=c(0,1))+labs(x='Predicted TTS',y='Reference TTS',title='TTS conditional predictions (species-equal means)'),
    '08a_tts_confusion',cm,7,4)
  # Condition on predicted canonical, NOT on exact recovery.
  canon<-cf %>% filter(predicted_class!='Other') %>% group_by(short_name,division,order,reference_motif) %>%
    mutate(n_canonical=sum(n)) %>% filter(n_canonical>0) %>% mutate(proportion=n/n_canonical) %>% ungroup()
  write_tsv(canon,file.path(output,'08a_canonical_only_species.tsv'))
  cc<-canon %>% group_by(reference_motif,predicted_class) %>% summarise(proportion=mean(proportion),n_species=n_distinct(short_name),.groups='drop')
  export(ggplot(cc,aes(predicted_class,reference_motif,fill=proportion))+geom_tile()+geom_text(aes(label=sprintf('%.1f%%',100*proportion)))+
    scale_fill_viridis_c(limits=c(0,1))+labs(x='Predicted stop codon',y='Reference stop codon',title='Conditional on a canonical prediction; species-equal means'),
    '08a_canonical_only',cc,7,4)
  co<-cf %>% group_by(division,order,reference_motif,predicted_class) %>% summarise(proportion=mean(proportion),n_species=n_distinct(short_name),.groups='drop')
  export(ggplot(co,aes(predicted_class,order,fill=proportion))+geom_tile()+facet_wrap(~reference_motif,nrow=1)+scale_fill_viridis_c(limits=c(0,1)),
    '08_tts_conditional_order',co,11,max(5,n_distinct(t$order)*.23))
}

# Species-stratified fixed-marginal null: shuffle predictions among sites within species.
# r2dtable samples the exact contingency-table distribution induced by label permutation.
B<-2000L; nullmeans<-numeric(B); observed<-numeric(); rows<-list(); cells<-list()
for(sp in unique(t$short_name)) {
  d<-filter(t,short_name==sp) %>% mutate(predicted_class=if_else(predicted_motif %in% c('TAA','TAG','TGA'),predicted_motif,'Other'))
  mat<-matrix(0,3,4,dimnames=list(c('TAA','TAG','TGA'),c('TAA','TAG','TGA','Other')))
  for(i in seq_len(nrow(d))) mat[d$reference_motif[i],d$predicted_class[i]]<-mat[d$reference_motif[i],d$predicted_class[i]]+d$n_sites[i]
  n<-sum(mat); rr<-rowSums(mat);cl<-colSums(mat);expected<-outer(rr,cl)/n
  obs<-sum(diag(mat[,1:3]))/n; chance<-sum(diag(expected[,1:3]))/n
  sims<-vapply(r2dtable(B,as.integer(rr),as.integer(cl)),function(m) sum(diag(m[,1:3]))/n,numeric(1))
  nullmeans<-nullmeans+sims;observed<-c(observed,obs)
  rows[[sp]]<-tibble(short_name=sp,division=d$division[1],n_sites=n,observed_exact=obs,expected_exact=chance,
    excess_pp=100*(obs-chance),null_lower=unname(quantile(sims,.025)),null_upper=unname(quantile(sims,.975)),p=(1+sum(sims>=obs))/(B+1),n_permutations=B)
  cd<-as.data.frame(as.table(mat));names(cd)<-c('reference_motif','predicted_class','observed_n')
  cd$short_name<-sp;cd$expected_n<-as.vector(expected);cd$excess_conditional_pp<-100*(as.vector(mat-expected)/rep(rr,4))
  cells[[sp]]<-cd
}
perm<-bind_rows(rows) %>% mutate(q=p.adjust(p,'BH'))
nullmeans<-nullmeans/length(observed)
global<-tibble(n_species=length(observed),observed_macro_exact=mean(observed),expected_macro_exact=mean(perm$expected_exact),
  p=(1+sum(nullmeans>=mean(observed)))/(B+1),n_permutations=B)
write_tsv(global,file.path(output,'08c_shuffle_global.tsv'))
write_tsv(tibble(permutation=seq_len(B),null_macro_exact=nullmeans),file.path(output,'08c_shuffle_global_null.tsv'))
r<-range(c(perm$expected_exact,perm$observed_exact),na.rm=TRUE)
padding<-diff(r)*0.05
if(padding==0) padding<-0.01
r<-c(r[1]-padding,r[2]+padding)
export(ggplot(perm,aes(expected_exact,observed_exact,color=division))+geom_abline(slope=1,intercept=0,linetype=2)+geom_point()+
  coord_equal(xlim=r,ylim=r)+labs(x='Expected exact recovery under within-species shuffle',y='Observed exact recovery',subtitle='Fixed reference and prediction frequencies; dashed line = chance expectation'),
  '08c_shuffle_exact',perm,6,6,perm)
cell<-bind_rows(cells);write_tsv(cell,file.path(output,'08d_excess_over_shuffle_species.tsv'))
ce<-cell %>% group_by(reference_motif,predicted_class) %>% summarise(excess_pp=mean(excess_conditional_pp,na.rm=TRUE),.groups='drop')
export(ggplot(ce,aes(predicted_class,reference_motif,fill=excess_pp))+geom_tile()+geom_text(aes(label=sprintf('%+.1f pp',excess_pp)))+
  scale_fill_gradient2(low='#2166ac',mid='white',high='#b2182b')+labs(x='Prediction',y='Reference',title='Observed minus shuffle expectation; species-equal means'),
  '08d_excess_over_shuffle',ce,8,5)

# Lineage-stratified associations and partial associations controlling division/training.
adjusted<-list()
for(v in intersect(vars,names(z))) for(ss in subsets) {
  d<-filter(z,subset==ss); d$genome_value<-suppressWarnings(as.numeric(d[[v]]))
  d<-filter(d,is.finite(genome_value),is.finite(exact_recovery_rate))
  if(v=='genome_assembly_span_bp') d<-filter(d,genome_value>0)
  d$plot_x<-if(v=='genome_assembly_span_bp') log10(d$genome_value) else d$genome_value
  st<-d %>% group_by(motif,division) %>% group_modify(~{
    fit<-if(nrow(.x)>=4 && n_distinct(.x$plot_x)>1 && n_distinct(.x$exact_recovery_rate)>1) suppressWarnings(cor.test(.x$plot_x,.x$exact_recovery_rate,method='spearman',exact=FALSE)) else NULL
    tibble(n=nrow(.x),rho=if(is.null(fit)) NA_real_ else unname(fit$estimate),p=if(is.null(fit)) NA_real_ else fit$p.value)
  }) %>% ungroup() %>% mutate(q=p.adjust(p,'BH'),label=sprintf('n=%d; rho=%.2f; q=%.3g',n,rho,q))
  export(ggplot(d,aes(plot_x,exact_recovery_rate))+geom_point(aes(color=training_group),alpha=.65)+
    geom_smooth(method='lm',formula=y~x,color='black',linewidth=.4)+facet_grid(division~motif,scales='free_x')+
    geom_text(data=st,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.02,vjust=1.3,size=2.3)+
    labs(x=if(v=='genome_assembly_span_bp') 'log10(assembly span bp)' else v,y='Exact recovery',title=paste('Within-division associations:',ss)),
    paste0('09b_within_division_',v,'_',ss),d,15,9,st)
  for(m in unique(d$motif)) {
    dd<-filter(d,motif==m);dd$standardized_x<-as.numeric(scale(dd$plot_x))
    fit<-tryCatch(lm(exact_recovery_rate~standardized_x+division+training_group,data=dd),error=function(e) NULL)
    if(is.null(fit)) next
    co<-summary(fit)$coefficients
    if(!'standardized_x' %in% rownames(co)) next
    ci2<-confint(fit,'standardized_x')
    adjusted[[paste(v,ss,m)]]<-tibble(variable=v,subset=ss,motif=as.character(m),n_species=nrow(dd),estimate_pp=100*co['standardized_x',1],
      lower_pp=100*ci2[1],upper_pp=100*ci2[2],p=co['standardized_x',4],rank_deficient=fit$rank<length(coef(fit)))
  }
}
ad<-bind_rows(adjusted) %>% mutate(q=p.adjust(p,'BH'))
export(ggplot(ad,aes(estimate_pp,variable,color=subset))+geom_vline(xintercept=0,linetype=2)+
  geom_point(position=position_dodge(width=.5))+geom_errorbar(aes(xmin=lower_pp,xmax=upper_pp),orientation='y',width=.2,position=position_dodge(width=.5))+
  facet_wrap(~motif)+labs(x='Recovery difference (pp) per 1 SD predictor; adjusted for division and training group',y=NULL),
  '09c_adjusted_associations',ad,14,7,ad)
# Species counts show sampling/training imbalance, not actual training token exposure.
composition<-meta %>% count(division,training_group,name='n_species')
export(ggplot(composition,aes(division,n_species,fill=training_group))+geom_col(position='dodge')+labs(x=NULL,y='Number of species in metadata'),
  '05b_training_composition',composition,9,5)
writeLines(capture.output(sessionInfo()),file.path(output,'R_sessionInfo.txt'))
message('Finished: ',output)
