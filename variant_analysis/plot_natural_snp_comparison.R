suppressPackageStartupMessages({library(data.table); library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE)
input_dir <- if(length(args)>=1) args[1] else '/Users/yenmr/project/BryoCAD/SNP'
out <- if(length(args)>=2) args[2] else 'output/pdf'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
colors <- c('PlantCAD2-Small'='#F5B18D', 'BryoCAD-Small'='#94BDEC')
files <- c('PlantCAD2-Small'='Marpol_135_accessions.PlantCAD2_small.scored.context.tsv.gz',
           'BryoCAD-Small'='Marpol_135_accessions.BryoCAD_small10.scored.context.tsv.gz')
cols <- c('chrom','pos','ref','alt','maf','zero_shot_score','fig4d_group')
read_model <- function(f) fread(cmd=paste('gzip -dc', shQuote(file.path(input_dir,f))), select=cols)
a <- read_model(files[1]); b <- read_model(files[2])
stopifnot(nrow(a)==nrow(b))
for(k in setdiff(cols,'zero_shot_score')) if(!identical(a[[k]],b[[k]])) stop('Inputs differ in ',k,'; align variants before comparison')
valid <- is.finite(a$zero_shot_score) & is.finite(b$zero_shot_score) & is.finite(a$maf) & a$fig4d_group %in% c('Missense','Synonymous','Noncoding')
a <- a[valid]; b <- b[valid]
message('Matched SNPs: ',nrow(a),'; excluded: ',sum(!valid))
a[,model:='PlantCAD2-Small']; b[,model:='BryoCAD-Small']
d <- rbindlist(list(a,b)); rm(a,b); gc()
d[,fig4d_group:=factor(fig4d_group,levels=c('Missense','Synonymous','Noncoding'))]
d[,model:=factor(model,levels=names(colors))]
d[,deleterious_score:=-zero_shot_score]
# Equivalent to dplyr::percent_rank: ties receive their minimum rank.
d[,rank_pct:=if(.N>1) (frank(deleterious_score,ties.method='min')-1)/(.N-1) else 0,by=.(model,fig4d_group)]
levels_top <- c('Top 50%','Top 10%','Top 1%','Top 0.5%')
d[,top_group:=fcase(rank_pct>=.995,'Top 0.5%',rank_pct>=.99,'Top 1%',rank_pct>=.9,'Top 10%',rank_pct>=.5,'Top 50%',default=NA_character_)]
s <- d[!is.na(top_group),.(n=.N,mean_maf=mean(maf),median_maf=median(maf),mean_deleterious_score=mean(deleterious_score)),by=.(model,fig4d_group,top_group)]
s[,top_group:=factor(top_group,levels=levels_top)]
cs <- rbindlist(lapply(seq_along(levels_top), function(i) d[rank_pct>=c(.5,.9,.99,.995)[i],.(n=.N,mean_maf=mean(maf),median_maf=median(maf),mean_deleterious_score=mean(deleterious_score)),by=.(model,fig4d_group)][,top_group:=levels_top[i]]))
cs[,top_group:=factor(top_group,levels=levels_top)]
fwrite(s,file.path(out,'model_comparison_fig4D_original_bins_summary_135.tsv'),sep='\t')
fwrite(cs,file.path(out,'model_comparison_fig4D_cumulative_summary_135.tsv'),sep='\t')
fwrite(d[,.(n=.N),by=.(model,fig4d_group)],file.path(out,'model_comparison_snp_counts_135.tsv'),sep='\t')
base_theme <- theme_classic(base_size=12)+theme(legend.position='top',legend.title=element_blank(),strip.background=element_rect(fill='#F3F3F3',color=NA),strip.text=element_text(face='bold'),plot.title=element_text(face='bold'),plot.caption=element_text(hjust=0,size=9,color='#555555'),axis.text.x=element_text(angle=30,hjust=1))
make_d <- function(dt,subtitle,caption) ggplot(dt,aes(top_group,mean_maf,color=model,group=model,shape=model))+geom_line(linewidth=.85)+geom_point(size=2.8)+scale_color_manual(values=colors)+scale_shape_manual(values=c(16,17))+facet_wrap(~fig4d_group,nrow=1)+labs(x='SNPs ranked within each category by each model',y='Mean MAF',title='MAF of model-prioritized SNPs | 135 accessions',subtitle=subtitle,caption=caption)+base_theme
save_plot <- function(p,stem,w,h) {
 ggsave(file.path(out,paste0(stem,'.pdf')),p,width=w,height=h,device=pdf)
 ggsave(file.path(out,paste0(stem,'.png')),p,width=w,height=h,dpi=220,bg='white')
}
save_plot(make_d(s,'Original script: mutually exclusive rank bins','Bins, left to right: 50-90%, 90-99%, 99-99.5%, and 99.5-100% of deleterious-score ranks.\nDeleterious score = -zero-shot score; identical SNPs and MAF values used for both models.'),'model_comparison_fig4D_original_bins_135',10.5,4.7)
save_plot(make_d(cs,'Cumulative selection: top 50%, 10%, 1%, and 0.5%','Each point includes all SNPs above the corresponding deleterious-score percentile.\nDeleterious score = -zero-shot score; identical SNPs and MAF values used for both models.'),'model_comparison_fig4D_cumulative_135',10.5,4.7)
# Aggregate common histogram breaks to avoid rendering millions of observations.
lims <- d[,.(lo=quantile(zero_shot_score,.001),hi=quantile(zero_shot_score,.999)),by=model]
breaks <- seq(min(d$zero_shot_score),max(d$zero_shot_score),length.out=301)
hist_df <- d[, {h<-hist(zero_shot_score,breaks=breaks,plot=FALSE); .(score=h$mids,count=h$counts)},by=model]
pC <- ggplot(hist_df,aes(score,count,color=model,fill=model))+geom_area(position='identity',alpha=.3)+geom_line(linewidth=.45)+geom_vline(data=lims,aes(xintercept=lo,color=model),linetype=2,linewidth=.6)+scale_color_manual(values=colors)+scale_fill_manual(values=colors)+coord_cartesian(xlim=c(min(lims$lo),max(lims$hi)))+labs(x='Zero-shot score',y='Number of SNPs per bin',title='Zero-shot score distributions | 135 accessions',subtitle='PlantCAD2-Small and continued-training BryoCAD-Small',caption='Identical SNP set; 300 shared bins. Dashed lines: model-specific 0.1st percentiles.\nDisplayed range spans the models\' 0.1st-99.9th percentiles; raw score scales may differ.')+base_theme
save_plot(pC,'model_comparison_fig4C_score_distribution_135',8,5)
fwrite(lims,file.path(out,'model_comparison_score_quantiles_135.tsv'),sep='\t')
print(s); print(cs)
writeLines(c(paste('Matched SNPs per model:',sum(valid)),paste('Excluded rows per model:',sum(!valid)),'Validated: chromosome, position, reference, alternative, MAF and category identical row-by-row.','Original-bin figure reproduces the supplied script rank logic; cumulative figure uses nested selections.','Colors from motif_recovery_with_small.svg: PlantCAD2-Small #F5B18D; BryoCAD-Small #94BDEC.'),file.path(out,'comparison_notes.txt'))
