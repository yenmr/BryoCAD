#!/usr/bin/env Rscript
# Rscript plot_comparisons.R comparison.metrics.tsv plots_by_species
suppressPackageStartupMessages(library(ggplot2))
a <- commandArgs(trailingOnly=TRUE)
if(length(a)!=2) stop('Provide metrics TSV and output directory')
d <- read.delim(a[1],check.names=FALSE,na.strings='NA')
out <- a[2]
for(sub in c('metrics','sections')) dir.create(file.path(out,sub),recursive=TRUE,showWarnings=FALSE)
d$Model <- factor(d$Model,levels=c('BryoCAD-Small','PlantCAD2-Small','Helixer'))
d$Head[d$Head=='-'] <- 'Baseline'
d$Head <- factor(d$Head,levels=c('2bryo','4bryo','5angio','Baseline'))
d$Species <- factor(d$Species,levels=c('Ath','Mpo','Ppa'))
d$Measure <- paste(d$section,d$metric,sep=' / ')
keys <- unique(d$Measure)
d$Measure <- factor(d$Measure,levels=keys)
cols <- c('BryoCAD-Small'='#94BDEC','PlantCAD2-Small'='#F5B18D',Helixer='#444444')
shapes <- c('2bryo'=16,'4bryo'=17,'5angio'=15,Baseline=18)
# Analytic iso-F1 curves: P = F*R/(2*R-F). No threshold sweep.
curves <- do.call(rbind,lapply(c(.2,.4,.6,.8,.9,.95),function(f) {
 r <- seq(f/(2-f),1,length.out=250)
 data.frame(recall=r,precision=pmin(1,f*r/(2*r-f)),F=f)
}))
fl <- data.frame(recall=.98,F=c(.2,.4,.6,.8,.9,.95))
fl$precision <- with(fl,F*recall/(2*recall-F));fl$label <- paste0('F1=',fl$F)
base_theme <- theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position='bottom',legend.box='vertical',strip.background=element_rect(fill='#EEEEEE'),plot.caption=element_text(hjust=0,size=9))
caption <- paste0('Structure policy: ',paste(unique(d$Policy),collapse=', '),'. Unavailable scores omitted; missing combinations are not zero.')
pr_plot <- function(x,title) {
 valid <- x[x$status=='available' & !is.na(x$precision) & !is.na(x$recall),]
 ggplot(valid,aes(recall,precision))+
 geom_path(data=curves,aes(recall,precision,group=F),inherit.aes=FALSE,colour='grey80',linewidth=.3)+
 geom_text(data=fl,aes(recall,precision,label=label),inherit.aes=FALSE,size=2.2,hjust=1,colour='grey55',vjust=-.4)+
 geom_point(aes(colour=Model,shape=Head),size=3.1,alpha=.95)+
 facet_grid(Measure~Species,drop=TRUE)+scale_colour_manual(values=cols,drop=FALSE)+scale_shape_manual(values=shapes,drop=FALSE)+
 scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1))+scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+coord_equal()+
 labs(title=title,x='Recall',y='Precision',caption=paste('Grey curves: iso-F1. One point per model/head, not a threshold PR curve.',caption,sep='\n'))+base_theme
}
bar_plot <- function(x,title) {
 ggplot(x[x$status=='available',],aes(Head,f1,fill=Model))+
 geom_col(position=position_dodge2(width=.8,preserve="single"),width=.7,na.rm=TRUE)+
 facet_grid(Measure~Species,drop=TRUE)+scale_fill_manual(values=cols,drop=FALSE)+
 scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+
 labs(title=title,x='Fine-tuning head',y='F1',caption=caption)+base_theme+
 theme(axis.text.x=element_text(angle=35,hjust=1))
}
save_pair <- function(p,name,h) {
 for(ext in c('pdf','svg')) ggsave(file.path(out,paste0(name,'.',ext)),p,width=12,height=h,limitsize=FALSE,device=if(ext=='svg') grDevices::svg else grDevices::pdf)
}
index <- data.frame(section=character(),metric=character(),PR=character(),F1_bar=character())
for(key in keys) {
 x <- droplevels(d[d$Measure==key,]);stem <- paste(x$section[1],x$metric[1],sep='__')
 save_pair(pr_plot(x,paste(key,'- precision / recall')),file.path('metrics',paste0(stem,'_PR')),5.6)
 save_pair(bar_plot(x,paste(key,'- F1')),file.path('metrics',paste0(stem,'_F1_bar')),5.4)
 index <- rbind(index,data.frame(section=x$section[1],metric=x$metric[1],PR=paste0('metrics/',stem,'_PR.pdf'),F1_bar=paste0('metrics/',stem,'_F1_bar.pdf')))
}
for(s in unique(d$section)) {
 x <- droplevels(d[d$section==s,]);n <- nlevels(x$Measure)
 save_pair(pr_plot(x,paste(s,'- precision / recall')),file.path('sections',paste0(s,'_PR')),3.45*n+1.8)
 save_pair(bar_plot(x,paste(s,'- F1')),file.path('sections',paste0(s,'_F1_bar')),2.5*n+1.8)
}
# Complete overview: one metric per page, all 20 metrics, same species columns.
for(kind in c('PR','F1_bar')) {
 pdf(file.path(out,paste0('overview_all_metrics_',kind,'.pdf')),width=12,height=5.6,onefile=TRUE)
 for(key in keys) {
  x <- droplevels(d[d$Measure==key,])
  print(if(kind=='PR') pr_plot(x,paste(key,'- precision / recall')) else bar_plot(x,paste(key,'- F1')))
 }
 dev.off()
}
write.table(index,file.path(out,'figure_index.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
message('Saved ',length(keys),' metric pairs, ',length(unique(d$section)),' section pairs, and two complete multipage PDFs; PDF/SVG only.')
