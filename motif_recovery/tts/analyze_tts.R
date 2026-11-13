#!/usr/bin/env Rscript
# Independent TTS analysis. Rscript analyze_tts.R processed_dir output_dir [model_label]
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(readr)})
a<-commandArgs(TRUE)
if(!length(a) %in% c(2,3)) stop('Usage: Rscript analyze_tts.R processed_dir output_dir [model_label]')
input<-a[1];output<-a[2];model_label<-if(length(a)==3) a[3] else basename(normalizePath(input))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
rd<-function(n) read_tsv(file.path(input,n),show_col_types=FALSE,na=c('NA',''))
t<-rd('tts_codon_confusion.tsv');meta<-rd('species_metadata.tsv')
req<-c('short_name','division','order','training_group','reference_motif','predicted_motif','n_sites')
if(!all(req %in% names(t))) stop('Missing required confusion columns: ',paste(setdiff(req,names(t)),collapse=', '))
if(!nrow(t) || anyNA(t[req]) || any(t$n_sites<0 | t$n_sites!=floor(t$n_sites))) stop('Invalid/empty confusion data')
if(any(!t$reference_motif %in% c('TAA','TAG','TGA'))) stop('Confusion input must contain canonical references only')
if(any(!grepl('^[ACGT]{3}$',t$predicted_motif))) stop('Invalid predicted motif')
if(anyDuplicated(meta$short_name)) stop('Duplicate metadata species')
if(anyDuplicated(t[c('short_name','reference_motif','predicted_motif')])) stop('Duplicate confusion cells')
t<-filter(t,n_sites>0)
if(!all(unique(t$short_name) %in% meta$short_name)) stop('Species missing from metadata')
if(!'genome_gc_fraction' %in% names(meta)) stop('species_metadata.tsv needs genome_gc_fraction (0..1)')
meta$genome_gc_fraction<-suppressWarnings(as.numeric(meta$genome_gc_fraction))
if(any(meta$genome_gc_fraction<0 | meta$genome_gc_fraction>1,na.rm=TRUE)) stop('GC must be a fraction, not percent')
s<-t %>% group_by(short_name,division,order,training_group) %>%
 summarise(n_exact=sum(n_sites[reference_motif==predicted_motif]),n_canonical=sum(n_sites[predicted_motif %in% c('TAA','TAG','TGA')]),n_sites=sum(n_sites),.groups='drop') %>%
 mutate(motif='tts',exact_recovery_rate=n_exact/n_sites,canonical_recovery_rate=n_canonical/n_sites)
set.seed(20261008)
theme_set(theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),strip.background=element_rect(fill='grey95')))
export<-function(p,name,data,w=12,h=6,stats=NULL) {
 write_tsv(data,file.path(output,paste0(name,'.tsv')))
 if(!is.null(stats)) write_tsv(stats,file.path(output,paste0(name,'_statistics.tsv')))
 ggsave(file.path(output,paste0(name,'.pdf')),p+labs(caption=model_label),width=w,height=h,limitsize=FALSE)
 # PNG deliberately disabled.
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


# GC analyses: canonical composition and unconditional probabilities are distinct.
gcmeta<-meta %>% select(short_name,genome_gc_fraction) %>% mutate(gc_percent=100*genome_gc_fraction)
gcqc<-gcmeta %>% filter(short_name %in% t$short_name) %>% mutate(gc_available=is.finite(gc_percent))
write_tsv(gcqc,file.path(output,'GC_metadata_QC.tsv'))
if(any(!gcqc$gc_available)) warning(sum(!gcqc$gc_available),' species omitted from GC analyses (missing GC)')
# Both reference-stratified and pooled counts; pooled changes can reflect reference composition.
byref<-cf %>% transmute(short_name,division,reference_group=as.character(reference_motif),predicted_class,n)
pooled<-byref %>% group_by(short_name,division,predicted_class) %>% summarise(n=sum(n),.groups='drop') %>% mutate(reference_group='All references')
prob<-bind_rows(byref,pooled) %>% group_by(short_name,division,reference_group) %>%
 mutate(n_total=sum(n),n_canonical=sum(n[predicted_class!='Other'])) %>% ungroup() %>% filter(predicted_class!='Other') %>%
 mutate(all_sites=n/n_total,canonical_only=if_else(n_canonical>0,n/n_canonical,NA_real_)) %>%
 pivot_longer(c(all_sites,canonical_only),names_to='denominator',values_to='probability') %>%
 left_join(gcmeta,by='short_name') %>% filter(is.finite(gc_percent),is.finite(probability)) %>%
 mutate(predicted_class=factor(predicted_class,levels=c('TAA','TAG','TGA')),
 reference_group=factor(reference_group,levels=c('All references','TAA','TAG','TGA')))
# Each cell/stratum contains at most one observation per species.
cortest<-function(d) {
 fit<-if(nrow(d)>=4 && n_distinct(d$gc_percent)>1 && n_distinct(d$probability)>1)
   suppressWarnings(cor.test(d$gc_percent,d$probability,method='spearman',exact=FALSE)) else NULL
 tibble(n_species=nrow(d),rho=if(is.null(fit)) NA_real_ else unname(fit$estimate),p=if(is.null(fit)) NA_real_ else fit$p.value)
}
stat<-prob %>% group_by(denominator,reference_group,predicted_class) %>% group_modify(~cortest(.x)) %>% ungroup() %>%
 mutate(q=p.adjust(p,'BH'),label=sprintf('n=%d; rho=%.2f; q=%.3g',n_species,rho,q))
within<-prob %>% group_by(denominator,reference_group,predicted_class,division) %>% group_modify(~cortest(.x)) %>% ungroup() %>%
 mutate(q=p.adjust(p,'BH'),label=sprintf('n=%d; rho=%.2f; q=%.3g',n_species,rho,q))
write_tsv(stat,file.path(output,'GC_correlations_all_tests.tsv'))
write_tsv(within,file.path(output,'GC_correlations_within_division_tests.tsv'))
for(den in c('canonical_only','all_sites')) {
 d<-filter(prob,denominator==den);st<-filter(stat,denominator==den)
 yl<-if(den=='canonical_only') 'P(predicted codon | predicted canonical, reference)' else 'P(predicted codon | reference); all predictions in denominator'
 export(ggplot(d,aes(gc_percent,probability))+geom_point(aes(color=division),alpha=.7)+
   geom_smooth(method='lm',formula=y~x,se=TRUE,color='black',linewidth=.45)+
   geom_text(data=st,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.02,vjust=1.2,size=2.8)+
   facet_grid(reference_group~predicted_class)+coord_cartesian(ylim=c(0,1))+
   labs(x='Species genome GC (%)',y=yl,subtitle='Rows: reference codon; columns: predicted codon. Spearman tests; OLS trend is descriptive.'),
   paste0('09_GC_predicted_',den),d,12,12,st)
 for(ref in levels(prob$reference_group)) {
  dd<-filter(d,reference_group==ref);ww<-filter(within,denominator==den,reference_group==ref)
  export(ggplot(dd,aes(gc_percent,probability))+geom_point(alpha=.7)+geom_smooth(method='lm',formula=y~x,linewidth=.4)+
   geom_text(data=ww,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.02,vjust=1.2,size=2.5)+
   facet_grid(division~predicted_class)+coord_cartesian(ylim=c(0,1))+
   labs(x='Species genome GC (%)',y=yl,title=paste('Reference:',ref)),
   paste0('10_GC_within_division_',den,'_',gsub(' ','_',ref)),dd,12,max(5,2.5*n_distinct(dd$division)),ww)
 }
}
# Equal-reference standardization distinguishes reference usage from preference.
standard<-prob %>% filter(reference_group!='All references') %>% group_by(short_name,division,denominator,predicted_class) %>%
 filter(n_distinct(reference_group)==3) %>%
 summarise(gc_percent=first(gc_percent),probability=mean(probability),.groups='drop')
ststd<-standard %>% group_by(denominator,predicted_class) %>% group_modify(~cortest(.x)) %>% ungroup() %>% mutate(q=p.adjust(p,'BH'),label=sprintf('rho=%.2f; q=%.3g',rho,q))
export(ggplot(standard,aes(gc_percent,probability))+geom_point(aes(color=division),alpha=.7)+geom_smooth(method='lm',formula=y~x,color='black')+
 geom_text(data=ststd,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.02,vjust=1.2,size=3)+
 facet_grid(denominator~predicted_class)+coord_cartesian(ylim=c(0,1))+
 labs(x='Species genome GC (%)',y='Mean conditional probability across TAA/TAG/TGA references',subtitle='Each reference codon weighted 1/3 within species'),
 '11_GC_equal_reference_weight',standard,12,7,ststd)
# Reference usage versus GC provides a direct composition comparison.
u<-usage %>% left_join(gcmeta,by='short_name') %>% filter(is.finite(gc_percent)) %>% mutate(probability=proportion)
ust<-u %>% group_by(reference_motif) %>% group_modify(~cortest(.x)) %>% ungroup() %>% mutate(q=p.adjust(p,'BH'))
export(ggplot(u,aes(gc_percent,proportion))+geom_point(aes(color=division),alpha=.7)+geom_smooth(method='lm',formula=y~x,color='black')+
 facet_wrap(~reference_motif,nrow=1)+labs(x='Species genome GC (%)',y='Reference codon frequency'), '12_GC_reference_usage',u,12,4,ust)
writeLines(capture.output(sessionInfo()),file.path(output,'R_sessionInfo.txt'))
writeLines(c(paste('Model label:',model_label),paste('Input:',normalizePath(input)),'Seed: 20261008','Shuffle replicates: 2000'),file.path(output,'run_info.txt'))
message('Finished TTS analysis: ',output)
