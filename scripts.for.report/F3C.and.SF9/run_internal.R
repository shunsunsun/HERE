#' ---
#' title: "F3C and SF9 summary"
#' ---

library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("foreach")
library("magick")
library("scales")
source("./scripts/common/ggpubr.A4.R")

F3B.for.plot.dt <- fread(snakemake@input[["F3B_for_plot_dt_txt_gz_filename"]])

if (FALSE){
    F3B.for.plot.dt <- fread("report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F3B.for.plot.dt.txt.gz")
}


combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt <- fread(snakemake@input[["combined_gexpr_FPKM_pc_only_melt_with_phenotype_normal_sample_only_median_annotated_dt_txt_gz_filename"]])
if (FALSE) {
    combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt <- fread("result/S42_1__annotate_embryonic_genes/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/201221-fifth-phenotype-collection/combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt.txt.gz")
}


combined.dt <- merge(
    x=F3B.for.plot.dt[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell'), list(edited.in.at.least.one.stage=TRUE), list(stage, Gene_ID, Gene_Name, Annotation.class.pasted, Annotation.class.pasted.description, gene.occurrence.pct)],
    y=combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt[, list(gene.id, gene.name, is.maternal, is.8.cell, cluster)],
    by.x=c("Gene_ID", "Gene_Name"),
    by.y=c("gene.id", "gene.name"),
    all.x=TRUE, all.y=TRUE) %>%
    {.[is.na(edited.in.at.least.one.stage), edited.in.at.least.one.stage:=FALSE]}


INTERNAL.compute.phyper.enrichment.of.maternal.decay.in.edited <- function(data, type){
    temp.p <- data %>% {phyper(
                     q=.[list(type, "maternal.decay"), count],
                     m=.[list(c(type, "unedited"), "maternal.decay"), sum(count)],
                     n=.[list(c(type, "unedited"), "maternal.others"), sum(count)],
                     k=.[list(type), sum(count)],
                     lower.tail=FALSE
             )}
    return(temp.p)
}

INTERNAL.compute.phyper.enrichment.of.edited.in.maternal.decay <- function(data, type){
    temp.p <- data %>% {phyper(
                     q=.[list(type, "maternal.decay"), count],
                     m=.[list(c(type), c("maternal.decay", "maternal.others")), sum(count)],
                     n=.[list(c("unedited"), c("maternal.others", "maternal.others")), sum(count)],
                     k=.[list(c(type, "unedited"), "maternal.decay"), sum(count)],
                     lower.tail=FALSE
             )}
    return(temp.p)
}


recurrently.edited.genes.simple.result.list <- combined.dt %>%
    {.[gene.occurrence.pct >= 0.8 | is.na(gene.occurrence.pct), list(Gene_ID, Gene_Name, cluster, edited.in.at.least.one.stage, Annotation.class.pasted.description)]} %>%
    unique %>%
    {.[, list(count=.N), list(cluster, edited.in.at.least.one.stage, Annotation.class.pasted.description)]} %>%
    {.[, count.log10:=log10(count)]} %>%
    {.[cluster %in% c("", NA), cluster:="others"]} %>%
    {.[is.na(Annotation.class.pasted.description)==TRUE, Annotation.class.pasted.description:="unedited"]} %>%
    {data.table(.[cluster %in% c("maternal.decay", "maternal.others")])} %>%
    {.[, cluster.description:=c(
             "maternal.decay"="decay at 8-cell",
             "maternal.others"="others"
         )[cluster]]} %>%
    setkey(Annotation.class.pasted.description, cluster) %>%
    {list(
         data=.,
         stat.test=data.table(
             `.y.`="count",
             `group1`=c("exonic", "intronic", "mixed"),
             `group2`="unedited",
             y.position=.[, max(count.log10)] * c(1.15, 1.3, 1.45)
         )
     )}
F3C.left.panel.ggplot <- recurrently.edited.genes.simple.result.list %>%
    {
        .[["stat.test"]][, p.raw:=lapply(group1, function(temp.type){INTERNAL.compute.phyper.enrichment.of.maternal.decay.in.edited(.[["data"]], temp.type)})][, p.adj:=p.adjust(p.raw, method="BH")][, p.adj.to.plot:=scientific(p.adj)]
        .
    } %>%
    {
        ggbarplot(data=.[["data"]], x="Annotation.class.pasted.description", y="count.log10", fill="cluster.description", position=position_dodge()) +
            stat_pvalue_manual(data=.[["stat.test"]], label="p.adj.to.plot", size=3) +
            labs(x="", y="log10(# genes)", fill="") +
            theme_pubr(base_size=10) +
            theme(axis.text.x=element_text(angle=45, hjust=1)) +
            scale_y_continuous(breaks=waiver(), n.breaks=7, limits=c(0, max(.[["data"]][, count.log10]) * 1.5))
    }
F3C.right.panel.ggplot <- combined.dt %>%
    {.[gene.occurrence.pct >= 0.8 | is.na(gene.occurrence.pct), list(Gene_ID, Gene_Name, cluster, edited.in.at.least.one.stage, stage, Annotation.class.pasted.description)]} %>%
    unique %>%
    {.[, .N, list(cluster, edited.in.at.least.one.stage, stage, Annotation.class.pasted.description)]} %>%
    {dcast(., Annotation.class.pasted.description+edited.in.at.least.one.stage+stage~cluster, value.var="N", fill=0)} %>%
    {.[, pct.decay:=`maternal.decay`/(`maternal.decay` + `maternal.others`)]} %>%
    {.[, pct.zga:=`zygotic.genomic.activation`/(`zygotic.genomic.activation` + `maternal.decay` + `maternal.others` + `V1`)]} %>%
    {melt(., id.vars=c("stage", "Annotation.class.pasted.description"), measure.vars=c("pct.decay", "pct.zga"), variable.name="pct.type", value.name="pct.value")} %>%
    {.[, pct.value.to.plot:=pct.value * 100]} %>%
    {.[, pct.type.description:=c("pct.decay"="Pct( maternal decay | recurrently edited & maternal )", "pct.others"="Pct( ZGA | recurrently edited )")[pct.type]]} %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])][, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))]} %>%
    {
        ggline(.[is.na(Annotation.class.pasted.description)==FALSE][pct.type=='pct.decay'], x="stage.description.ordered", y="pct.value.to.plot", color="#F8766D") +
            geom_hline(yintercept=.[is.na(Annotation.class.pasted.description)==TRUE & pct.type=='pct.decay', pct.value.to.plot], color="#F8766D", linetype="dashed") +
            ##geom_hline(yintercept=.[is.na(Annotation.class.pasted.description)==TRUE & pct.type=='pct.zga', pct.value], color="#4DBBD5", linetype="dashed") +
            lims(y=c(0, 100)) + 
            facet_grid(~Annotation.class.pasted.description) +
            labs(x="", y="%", color="") +
            theme_pubr(base_size=10) +
            theme(axis.text.x=element_text(angle=45, hjust=1)) +
            guides(fill=guide_legend(nrow=1))}
F3C.ggplot <- ggarrange(plotlist=list(F3C.left.panel.ggplot, F3C.right.panel.ggplot), ncol=2, widths=c(1,2))
F3C.ggplot

saveRDS(F3C.ggplot, snakemake@output[["main_ggplot_RDS_filename"]])


ggsave.A4(filename=snakemake@output[["main_png_filename"]], plot=F3C.ggplot, width.r=0.9, height.r=0.25)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F3C.png", plot=F3C.ggplot, width.r=0.9, height.r=0.25)
}


print(F3C.ggplot)

image_read(snakemake@output[["main_png_filename"]]) %>% print




SF9.ggplot <- recurrently.edited.genes.simple.result.list %>%
    {
        .[["stat.test"]][, p.raw:=lapply(group1, function(temp.type){INTERNAL.compute.phyper.enrichment.of.edited.in.maternal.decay(.[["data"]], temp.type)})][, p.adj:=p.adjust(p.raw, method="BH")][, p.adj.to.plot:=scientific(p.adj)]
        .
    } %>%
    {
        ggbarplot(data=.[["data"]], x="Annotation.class.pasted.description", y="count.log10", fill="cluster.description", position=position_dodge()) +
            stat_pvalue_manual(data=.[["stat.test"]], label="p.adj.to.plot", size=3) +
            labs(x="", y="log10(# genes)", fill="") +
            theme_pubr(base_size=10) +
            theme(axis.text.x=element_text(angle=45, hjust=1)) +
            scale_y_continuous(breaks=waiver(), n.breaks=7, limits=c(0, max(.[["data"]][, count.log10]) * 1.5))
    }
## SF9.ggplot

saveRDS(SF9.ggplot, snakemake@output[["supp1_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(SF9.ggplot, "report/201218-fifth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF9.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["supp1_png_filename"]], plot=SF9.ggplot, width.r=0.4, height.r=0.4)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF9.png", plot=SF9.ggplot, width.r=0.4, height.r=0.4)
}


print(SF9.ggplot)

image_read(snakemake@output[["supp1_png_filename"]]) %>% print


