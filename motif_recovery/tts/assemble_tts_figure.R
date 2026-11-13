#!/usr/bin/env Rscript
# Rscript assemble_tts_figure.R analyze_output figure_output [model_label] [canonical_only|all_sites]
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(readr);library(grid)})
a<-commandArgs(TRUE)
if(length(a)<2 || length(a)>4) stop('Usage: Rscript assemble_tts_figure.R analyze_output figure_output [model_label] [canonical_only|all_sites]')
input<-a[1];output<-a[2];model<-if(length(a)>=3) a[3] else ''
den<-if(length(a)>=4) a[4] else 'canonical_only'
if(!den %in% c('canonical_only','all_sites')) stop('Invalid denominator')
dir.create(output,recursive=TRUE,showWarnings=FALSE)
rd<-function(n) read_tsv(file.path(input,n),show_col_types=FALSE)
A<-rd('07b_reference_codon_distribution.tsv')
B<-rd('07c_reference_codon_recovery.tsv')
C<-rd('08c_shuffle_exact.tsv')
D<-rd(paste0('09_GC_predicted_',den,'.tsv')) %>% filter(reference_group=='TAA',predicted_class=='TAA')
DS<-rd(paste0('09_GC_predicted_',den,'_statistics.tsv')) %>% filter(reference_group=='TAA',predicted_class=='TAA')
if(!nrow(D) || nrow(DS)!=1) stop('Expected one TAA/TAA statistics row and nonempty source data')
lev<-c('Hornworts','Liverworts','Mosses');lev<-c(intersect(lev,unique(A$division)),sort(setdiff(unique(A$division),lev)))
colors<-setNames(grDevices::hcl.colors(length(lev),'Dark 3'),lev)
# Fixed bryophyte colors; additional angiosperm groups receive distinct colors.
fixed<-c(Hornworts='#CC6677',Liverworts='#228833',Mosses='#4477AA')
colors[intersect(names(fixed),names(colors))]<-fixed[intersect(names(fixed),names(colors))]
for(nm in c('A','B','C','D')) {
 d<-get(nm);d$division<-factor(d$division,levels=lev)
 if('reference_motif' %in% names(d)) d$reference_motif<-factor(d$reference_motif,levels=c('TAA','TAG','TGA'))
 assign(nm,d)
}
B$metric<-factor(B$metric,levels=c('Exact recovery','Canonical recovery'))
theme_set(theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
 strip.background=element_rect(fill='grey96'),plot.title=element_text(face='bold',size=11),
 plot.tag=element_text(face='bold',size=16),legend.position='bottom',plot.margin=margin(8,12,8,8)))
colscale<-function() scale_color_manual(values=colors,drop=FALSE)
pA<-ggplot(A,aes(reference_motif,proportion))+geom_boxplot(outlier.shape=NA,width=.55,fill='grey96')+
 geom_point(aes(color=division),position=position_jitter(width=.12,height=0,seed=7),alpha=.6,size=1.2)+
 facet_wrap(~division,nrow=1)+colscale()+coord_cartesian(ylim=c(0,1))+
 labs(tag='A',title='Reference stop-codon usage',x=NULL,y='Reference proportion')+theme(legend.position='none')
pB<-ggplot(B,aes(reference_motif,recovery_rate))+geom_boxplot(outlier.shape=NA,width=.55,fill='grey96')+
 geom_point(aes(color=division),position=position_jitter(width=.12,height=0,seed=7),alpha=.6,size=1.2)+
 facet_grid(metric~division)+colscale()+coord_cartesian(ylim=c(0,1))+
 labs(tag='B',title='Recovery by reference stop codon',x='Reference stop codon',y='Recovery rate')+theme(legend.position='none')
r<-range(c(C$expected_exact,C$observed_exact),na.rm=TRUE);pad<-max(diff(r)*.05,.001);r<-r+c(-pad,pad)
pC<-ggplot(C,aes(expected_exact,observed_exact,color=division))+geom_abline(slope=1,intercept=0,linetype=2,color='grey40')+
 geom_point(alpha=.8,size=1.5)+colscale()+coord_equal(xlim=r,ylim=r)+
 labs(tag='C',title='Exact recovery versus shuffled expectation',x='Expected exact recovery',y='Observed exact recovery',color=NULL)
ylab<-if(den=='canonical_only') 'P(predicted TAA | canonical prediction, ref TAA)' else 'P(predicted TAA | ref TAA)'
pD<-ggplot(D,aes(gc_percent,probability))+geom_point(aes(color=division),alpha=.8,size=1.5)+
 geom_smooth(method='lm',formula=y~x,color='black',linewidth=.5)+colscale()+coord_cartesian(ylim=c(0,1))+
 labs(tag='D',title='Genome GC and TAA prediction',subtitle=sprintf('Spearman rho = %.2f; BH q = %.3g; n = %d',DS$rho,DS$q,DS$n_species),
 x='Genome GC (%)',y=ylab,color=NULL)
# Native grid assembly avoids adding a patchwork/cowplot dependency.
width<-max(11,3.2*length(lev));height<-12.5
pdf(file.path(output,'TTS_manuscript_figure.pdf'),width=width,height=height,useDingbats=FALSE)
grid.newpage();pushViewport(viewport(layout=grid.layout(3,2,heights=unit(c(2.6,4.6,5.3),'null'))))
print(pA,vp=viewport(layout.pos.row=1,layout.pos.col=1:2),newpage=FALSE)
print(pB,vp=viewport(layout.pos.row=2,layout.pos.col=1:2),newpage=FALSE)
print(pC,vp=viewport(layout.pos.row=3,layout.pos.col=1),newpage=FALSE)
print(pD,vp=viewport(layout.pos.row=3,layout.pos.col=2),newpage=FALSE)
popViewport();dev.off()
for(nm in c('A','B','C','D')) write_tsv(get(nm),file.path(output,paste0('panel_',nm,'_source.tsv')))
write_tsv(DS,file.path(output,'panel_D_statistics.tsv'))
if(file.exists(file.path(input,'08c_shuffle_global.tsv'))) write_tsv(rd('08c_shuffle_global.tsv'),file.path(output,'panel_C_global_statistics.tsv'))
# Long-form plot-source table with distinct panel definitions; no mixed denominators.
combined<-bind_rows(
 A %>% transmute(panel='A',short_name,division,reference_motif=as.character(reference_motif),metric='reference_usage',x=as.character(reference_motif),y=proportion),
 B %>% transmute(panel='B',short_name,division,reference_motif=as.character(reference_motif),metric=as.character(metric),x=as.character(reference_motif),y=recovery_rate),
 C %>% transmute(panel='C',short_name,division,reference_motif=NA_character_,metric='observed_vs_expected_exact',x=as.character(expected_exact),y=observed_exact),
 D %>% transmute(panel='D',short_name,division,reference_motif='TAA',metric=paste0('predicted_TAA_',den),x=as.character(gc_percent),y=probability))
write_tsv(combined,file.path(output,'TTS_manuscript_source_long.tsv'))
writeLines(c(paste('Model:',model),
 'A: Species reference stop-codon proportions. Boxes show median and IQR; whiskers extend to 1.5 IQR; dots are species.',
 'B: Exact and canonical recovery conditional on reference codon, with all predictions in each denominator. Dots are species.',
 'C: Observed species exact recovery versus fixed-margin within-species shuffled expectation; dashed identity line. Both axes share limits.',
 paste('D: Reference=TAA, prediction=TAA; denominator =',den,'. Species genome GC; OLS line/95% mean CI; Spearman two-sided q retained from original GC analysis family.'),
 'D is a selected exploratory subpanel. It does not establish a general GC mechanism across all reference codons. Full GC panels should remain supplementary.',
 'Canonical-only means a stop codon was predicted, not necessarily the correct stop codon. Genome GC is not local-context GC.'),file.path(output,'figure_notes.txt'))
message('Wrote ',file.path(output,'TTS_manuscript_figure.pdf'))
