#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
a<-commandArgs(TRUE);if(length(a)!=2)stop('Provide metrics TSV and output directory')
d<-read.delim(a[1],na.strings='NA',check.names=FALSE);out<-a[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
cols<-c(Hornworts='#F8766D',Liverworts='#00BA38',Mosses='#619CFF')
d$division<-factor(d$division,levels=names(cols));d$key<-paste(d$section,d$metric,sep='/')
rows<-list()
for(k in unique(d$key))for(group in c('All species',names(cols))){
 x<-d[d$key==k & (group=='All species' | as.character(d$division)==group),];z<-x[!is.na(x$status)&x$status=='available',];f<-z$f1[is.finite(z$f1)]
 q<-if(length(f))quantile(f,c(.25,.5,.75)) else rep(NA_real_,3)
 rows[[paste(k,group)]]<-data.frame(metric=k,group=group,n_species=nrow(x),n_available=nrow(z),n_defined_F1=length(f),median_precision=if(any(is.finite(z$precision)))median(z$precision,na.rm=TRUE) else NA,median_recall=if(any(is.finite(z$recall)))median(z$recall,na.rm=TRUE) else NA,F1_Q1=q[1],F1_median=q[2],F1_Q3=q[3],F1_min=if(length(f))min(f) else NA,F1_max=if(length(f))max(f) else NA,n_F1_ge_0.5=sum(f>=.5),n_F1_ge_0.8=sum(f>=.8),row.names=NULL)
}
s<-do.call(rbind,rows);write.table(s,file.path(out,'metric_summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
keys<-c('base_level/cds','CDS/cds_ex','CDS/tx','coding_structure/intron','coding_structure/ichain','sites/Coding.Junc.Donor')
titles<-c('CDS bases','CDS segments','Complete CDS chain','Coding introns','Coding intron chains','Coding donor sites')
x<-d[d$key %in% keys & d$status=='available' & is.finite(d$f1),];x$Metric<-factor(x$key,levels=keys,labels=titles)
labels<-setNames(paste0(names(cols),'\n(n=',as.integer(table(factor(unique(d[c('short_name','division')])$division,levels=names(cols)))),')'),names(cols))
p<-ggplot(x,aes(division,f1,fill=division))+geom_boxplot(width=.55,outlier.shape=NA,alpha=.3)+geom_point(aes(colour=division),position=position_jitter(width=.16,height=0,seed=23),size=1.3,alpha=.65)+facet_wrap(~Metric,ncol=3)+scale_fill_manual(values=cols)+scale_colour_manual(values=cols)+scale_x_discrete(labels=labels)+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2))+labs(x=NULL,y='F1',title='Prediction agreement across 123 unseen species',subtitle='BryoCAD-Small / 4bryo; one point per species; all six metrics available for all 123 species',caption='Species are absent from head training (user-confirmed); base-model pretraining exposure is not established.\nReference agreement on selected whole-contig subsets (~50 Mb target), not whole-genome coverage.\nBoxes: median and IQR; whiskers: within 1.5 IQR. Four hornwort species; no significance tests.')+theme_bw(base_size=11)+theme(legend.position='none',panel.grid.minor=element_blank(),plot.caption=element_text(hjust=0,size=8.5))
for(ext in c('pdf','svg'))ggsave(file.path(out,paste0('summary_F1_boxplot.',ext)),p,width=12,height=8,device=if(ext=='svg')grDevices::svg else grDevices::pdf)
write.table(x,file.path(out,'summary_F1_boxplot.source.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
