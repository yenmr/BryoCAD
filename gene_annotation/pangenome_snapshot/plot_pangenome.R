#!/usr/bin/env Rscript
# Rscript plot_pangenome.R pangenome.metrics.tsv plots
suppressPackageStartupMessages(library(ggplot2))
a<-commandArgs(TRUE);if(length(a)!=2)stop('Provide metrics TSV and output directory')
d<-read.delim(a[1],na.strings='NA',check.names=FALSE);out<-a[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
# plot_recovery.R uses ggplot's default discrete hue palette, alphabetic divisions.
cols<-c(Hornworts='#F8766D',Liverworts='#00BA38',Mosses='#619CFF')
if(any(!d$division %in% names(cols)))stop('Unknown division')
d$division<-factor(d$division,levels=names(cols))
meta<-unique(d[c('short_name','species','division','order')]);meta<-meta[order(meta$division,meta$order,meta$species),]
d$short_name<-factor(d$short_name,levels=rev(meta$short_name))
labels<-setNames(paste(meta$species,paste0(' [',meta$short_name,']')),meta$short_name)
curves<-do.call(rbind,lapply(c(.2,.4,.6,.8,.9,.95),function(f){r<-seq(f/(2-f),1,length.out=300);data.frame(recall=r,precision=pmin(1,f*r/(2*r-f)),F=f)}))
fl<-data.frame(recall=.98,F=c(.2,.4,.6,.8,.9,.95));fl$precision<-with(fl,F*recall/(2*recall-F));fl$label<-paste0('F1=',fl$F)
base<-theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position='bottom',plot.caption=element_text(hjust=0,size=9))
save_pair<-function(p,name,w,h){for(ext in c('pdf','svg'))ggsave(file.path(out,paste0(name,'.',ext)),p,width=w,height=h,limitsize=FALSE,device=if(ext=='svg')grDevices::svg else grDevices::pdf)}
summary<-list();plots<-list()
for(key in unique(paste(d$section,d$metric,sep='/'))){
 x<-d[paste(d$section,d$metric,sep='/')==key,];stem<-gsub('/','__',key,fixed=TRUE)
 valid<-x[x$status=='available'&is.finite(x$precision)&is.finite(x$recall),]
 counts<-table(factor(valid$division,levels=names(cols)))
 legend_labels<-setNames(paste0(names(cols),' (n=',as.integer(counts),')'),names(cols))
 sub<-sprintf('BryoCAD-Small / 4bryo; %d/%d species with defined precision and recall',nrow(valid),nrow(x))
 p<-ggplot(valid,aes(recall,precision))+geom_path(data=curves,aes(recall,precision,group=F),inherit.aes=FALSE,colour='grey80',linewidth=.35)+
 geom_text(data=fl,aes(recall,precision,label=label),inherit.aes=FALSE,size=2.6,hjust=1,colour='grey55',vjust=-.4)+
 geom_point(aes(colour=division),size=2.5,alpha=.75)+scale_colour_manual(values=cols,limits=names(cols),labels=legend_labels,drop=FALSE)+
 scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.2))+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2))+coord_equal()+
 labs(title=key,subtitle=sub,x='Recall',y='Precision',colour=NULL,caption='One point per species; grey curves are iso-F1, not threshold sweeps.\nSelected whole-contig subset (target ~50 Mb); structure policy: exon. NA is not zero.')+base
 if(!nrow(valid))p<-p+annotate('text',x=.5,y=.5,label='No evaluable species for this metric',size=4)
 save_pair(p,paste0(stem,'_PR'),7.5,7.5);plots[[key]]<-p
 b<-ggplot(x,aes(f1,short_name,fill=division))+geom_col(width=.8,na.rm=TRUE)+
 geom_text(data=x[!is.finite(x$f1),],aes(x=.02,label='NA'),hjust=0,size=2.1,colour='grey45')+
 scale_fill_manual(values=cols,limits=names(cols),drop=FALSE)+scale_y_discrete(labels=labels,drop=FALSE)+scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.2))+
 labs(title=key,subtitle='Species grouped by division, order, and name',x='F1',y=NULL,fill=NULL,caption='No gene filtering. Unavailable scores are shown as NA; they are not zero. Structure policy: exon.')+base+
 theme(axis.text.y=element_text(size=5.5),panel.grid.major.y=element_blank())
 save_pair(b,paste0(stem,'_F1_bar'),11,max(7,nrow(x)*.14+2))
 for(div in names(cols)){
  z<-x[x$division==div,];ok<-z$status=='available'&is.finite(z$f1)
  summary[[paste(key,div)]]<-data.frame(section=x$section[1],metric=x$metric[1],division=div,n_species=nrow(z),n_F1=sum(ok),n_PR=sum(valid$division==div),median_F1=if(any(ok))median(z$f1[ok]) else NA,mean_F1=if(any(ok))mean(z$f1[ok]) else NA)
 }
}
pdf(file.path(out,'all_metrics_PR.pdf'),width=7.5,height=7.5,onefile=TRUE);for(p in plots)print(p);dev.off()
write.table(do.call(rbind,summary),file.path(out,'division_summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
message('Saved all 20 metrics: individual PR and F1 barplots, PDF/SVG only.')
