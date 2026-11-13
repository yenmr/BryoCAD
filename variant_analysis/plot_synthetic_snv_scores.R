library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
annotation_file <- if (length(args) >= 1) args[[1]] else "core_gene_snvs.annotation.tsv"
scored_vcf <- if (length(args) >= 2) args[[2]] else "core_gene_snvs.scored.vcf"
anno <- fread(annotation_file)

vcf <- fread(
    cmd=paste("grep -v", shQuote("^#"), shQuote(scored_vcf)),
    sep="\t",
    header=FALSE
)

colnames(vcf)[1:8] <- c(
    "chrom","pos","id","ref","alt",
    "qual","filter","info"
)

score <- vcf %>%
    transmute(
        chrom,
        pos=as.integer(pos),
        ref,
        alt,
        zero_shot_score=
            as.numeric(
                str_match(
                    info,
                    "plantCAD_zero_shot=([^;]+)"
                )[,2]
            )
    )

df <- left_join(
    anno,
    score,
    by=c("chrom","pos","ref","alt")
)

df$context <- factor(
  df$context,
  levels = c(
    "stop_gained_CDS",
    "start_lost_CDS",
    "splice_donor",
    "splice_acceptor",
    "stop_lost_CDS",
    "missense_CDS",
    "synonymous_CDS",
    "UTR",
    "promoter_1kb",
    "intron_or_gene_body"
  )
)

p <- ggplot(
    df,
    aes(context, zero_shot_score)
)+

geom_violin(
    trim=FALSE,
    scale="width",
    fill="#90CAF9",
    color="black",
    linewidth=0.3
)+

geom_boxplot(
    width=0.08,
    outlier.shape=NA,
    fill="white"
)+

geom_hline(
    yintercept=0,
    linetype=2,
    color="red"
)+

coord_cartesian(
    ylim=c(-12,4)
)+

theme_classic(base_size=14)+

theme(
    axis.text.x=
        element_text(
            angle=45,
            hjust=1
        ),
    panel.border=
        element_rect(
            colour="black",
            fill=NA
        )
)+

labs(
    x=NULL,
    y="BryoCAD zero-shot score"
)

ggsave(
    "BryoCAD_Fig4B_like.pdf",
    p,
    width=8,
    height=5
)

summary <- df %>%
group_by(context) %>%
summarise(

n=n(),

median=median(zero_shot_score),

mean=mean(zero_shot_score),

sd=sd(zero_shot_score)

)

write.table(
    summary,
    "BryoCAD_Fig4B_summary.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)


