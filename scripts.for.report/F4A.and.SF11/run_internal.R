#' ---
#' title: "F4A and SF11 summary"
#' ---

library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("foreach")
library("magick")
library("scales")
source("./scripts/common/ggpubr.A4.R")

F4B.for.plot.dt <- fread(snakemake@input[["F4B_for_plot_dt_txt_gz_filename"]])

if (FALSE){
    F4B.for.plot.dt <- fread("report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/F4B.for.plot.dt.txt.gz")
}




combined.gexpr.FPKM.pc.only.melt.with.phenotype.GSE133854.sample.only.median.annotated.dt <- fread(snakemake@input[["combined_gexpr_FPKM_pc_only_melt_with_phenotype_GSE133854_sample_only_median_annotated_dt_txt_gz_filename"]])
if (FALSE) {
    combined.gexpr.FPKM.pc.only.melt.with.phenotype.GSE133854.sample.only.median.annotated.dt <- fread("result/S42_2__annotate_uniparental_embryonic_genes/210215-sixth-dataset/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/201221-fifth-phenotype-collection/combined.gexpr.FPKM.pc.only.melt.with.phenotype.GSE133854.sample.only.median.annotated.dt.txt.gz")
}


combined.gexpr.FPKM.pc.only.melt.with.phenotype.GSE133854.sample.only.median.annotated.dt[, disease.description:=disease][disease %in% c("", NA), disease.description:="biparental"]

combined.dt <- merge(
    x=F4B.for.plot.dt[, list(edited.in.at.least.one.stage=TRUE), list(stage, Gene_ID, Gene_Name, Annotation.class.pasted, Annotation.class.pasted.description, gene.occurrence.pct, disease.description)],
    y=combined.gexpr.FPKM.pc.only.melt.with.phenotype.GSE133854.sample.only.median.annotated.dt[, list(gene.id, gene.name, is.maternal, is.8.cell, cluster, disease.description)],
    by.x=c("Gene_ID", "Gene_Name", "disease.description"),
    by.y=c("gene.id", "gene.name", "disease.description"),
    all.x=TRUE, all.y=TRUE) %>%
    {.[is.na(edited.in.at.least.one.stage), edited.in.at.least.one.stage:=FALSE]}



INTERNAL.compute.phyper.enrichment.of.maternal.decay.in.edited.per.disease <- function(data, type, disease){
    data %>% {print(list(
        type=type,
        disease=disease,
        q=.[list(type, "maternal.decay", disease), count],
        m=.[list(c(type, "unedited"), "maternal.decay", disease), sum(count)],
        n=.[list(c(type, "unedited"), "maternal.others", disease), sum(count)],
        k=.[list(type, c("maternal.decay", "maternal.others"), disease), sum(count)]
    ))}
    temp.p <- data %>% {phyper(
                     q=.[list(type, "maternal.decay", disease), count],
                     m=.[list(c(type, "unedited"), "maternal.decay", disease), sum(count)],
                     n=.[list(c(type, "unedited"), "maternal.others", disease), sum(count)],
                     k=.[list(type, c("maternal.decay", "maternal.others"), disease), sum(count)],
                     lower.tail=FALSE
             )}
    return(temp.p)
}

INTERNAL.compute.phyper.enrichment.of.edited.in.maternal.decay.per.disease <- function(data, type, disease){
    temp.p <- data %>% {phyper(
                     q=.[list(type, "maternal.decay", disease), count],
                     m=.[list(c(type), c("maternal.decay", "maternal.others"), disease), sum(count)],
                     n=.[list(c("unedited"), c("maternal.others", "maternal.others"), disease), sum(count)],
                     k=.[list(c(type, "unedited"), "maternal.decay", disease), sum(count)],
                     lower.tail=FALSE
             )}
    return(temp.p)
}


recurrently.edited.genes.simple.result.list <- combined.dt %>%
    {.[gene.occurrence.pct >= 0.8 | is.na(gene.occurrence.pct)]} %>%
    {.[, list(Gene_ID, Gene_Name, disease.description, cluster, edited.in.at.least.one.stage, Annotation.class.pasted.description)]} %>%
    unique %>% ## collapse edit sites
    {.[cluster %in% c("", NA), cluster:="others"]} %>%
    {.[, list(count=.N), list(cluster, edited.in.at.least.one.stage, Annotation.class.pasted.description, disease.description)]} %>%
    {.[, count.log10:=log10(count)]} %>%
    {.[is.na(Annotation.class.pasted.description)==TRUE, Annotation.class.pasted.description:="unedited"]} %>%
    {data.table(.[cluster %in% c("maternal.decay", "maternal.others")])} %>%
    {.[, cluster.description:=c(
             "maternal.decay"="decay at 8-cell",
             "maternal.others"="others"
         )[cluster]]} %>%
    setkey(Annotation.class.pasted.description, cluster, disease.description) %>%
    {list(
         data=.,
         stat.test=data.table(
             `.y.`="count",
             `group1`=rep(c("exonic", "intronic", "mixed"), 3),
             `group2`="unedited",
             disease.description=c(rep("androgenetic", 3), rep("biparental", 3), rep("parthenogenetic", 3)),
             y.position=.[, max(count.log10)] * c(1.15, 1.3, 1.45)
         )
     )}


F4A.left.panel.ggplot <- recurrently.edited.genes.simple.result.list %>%
    {
        .[["stat.test"]][, p.raw:=mapply(group1, disease.description, FUN=function(temp.type, temp.disease.description){INTERNAL.compute.phyper.enrichment.of.maternal.decay.in.edited.per.disease(.[["data"]], temp.type, temp.disease.description)})][, p.adj:=p.adjust(p.raw, method="BH")][, p.adj.to.plot:=scientific(p.adj)]
        .
    } %>%
    {
        ggbarplot(data=.[["data"]], x="Annotation.class.pasted.description", y="count.log10", fill="cluster.description", position=position_dodge()) +
            stat_pvalue_manual(data=.[["stat.test"]], label="p.adj.to.plot", size=3) +
            labs(x="", y="log10(# genes)", fill="") +
            theme_pubr(base_size=10) +
            theme(axis.text.x=element_text(angle=45, hjust=1)) +
            scale_y_continuous(breaks=waiver(), n.breaks=7, limits=c(0, max(.[["data"]][, count.log10]) * 1.5)) +
            facet_grid(~disease.description)
    }
F4A.left.panel.ggplot

F4A.right.panel.ggplot <- combined.dt %>%
    {.[gene.occurrence.pct >= 0.8 | is.na(gene.occurrence.pct), list(Gene_ID, Gene_Name, cluster, edited.in.at.least.one.stage, stage, Annotation.class.pasted.description, disease.description)]} %>%
    unique %>%
    {.[, .N, list(cluster, edited.in.at.least.one.stage, stage, Annotation.class.pasted.description, disease.description)]} %>%
    {dcast(., Annotation.class.pasted.description+edited.in.at.least.one.stage+stage+disease.description~cluster, value.var="N", fill=0)} %>%
    {.[, pct.decay:=`maternal.decay`/(`maternal.decay` + `maternal.others`)]} %>%
    {.[, pct.zga:=`zygotic.genomic.activation`/(`zygotic.genomic.activation` + `maternal.decay` + `maternal.others` + `V1`)]} %>%
    {melt(., id.vars=c("stage", "Annotation.class.pasted.description", "disease.description"), measure.vars=c("pct.decay", "pct.zga"), variable.name="pct.type", value.name="pct.value")} %>%
    {.[, pct.value.to.plot:=pct.value * 100]} %>%
    {.[, pct.type.description:=c("pct.decay"="Pct( maternal decay | recurrently edited & maternal )", "pct.others"="Pct( ZGA | recurrently edited )")[pct.type]]} %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])][, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))]} %>%
    {
        ggline(.[is.na(Annotation.class.pasted.description)==FALSE][pct.type=='pct.decay'], x="stage.description.ordered", y="pct.value.to.plot", color="disease.description") +
            geom_hline(data=.[is.na(Annotation.class.pasted.description)==TRUE & pct.type=='pct.decay'] %>% {list(copy(.)[, Annotation.class.pasted.description:='exonic'], copy(.)[, Annotation.class.pasted.description:='intronic'], copy(.)[, Annotation.class.pasted.description:='mixed'])} %>% rbindlist, aes(yintercept=pct.value.to.plot, color=disease.description), linetype="dashed") +
            ##geom_hline(yintercept=.[is.na(Annotation.class.pasted.description)==TRUE & pct.type=='pct.zga', pct.value], color="#4DBBD5", linetype="dashed") +
            lims(y=c(0, 100)) + 
            facet_grid(~Annotation.class.pasted.description) +
            labs(x="", y="%", color="") +
            theme_pubr(base_size=10) +
            theme(axis.text.x=element_text(angle=45, hjust=1)) +
            guides(fill=guide_legend(nrow=1))}
F4A.right.panel.ggplot

F4A.ggplot <- ggarrange(plotlist=list(F4A.left.panel.ggplot, F4A.right.panel.ggplot), ncol=2, widths=c(1.5,2))
F4A.ggplot

saveRDS(F4A.ggplot, snakemake@output[["main_ggplot_RDS_filename"]])


ggsave.A4(filename=snakemake@output[["main_png_filename"]], plot=F4A.ggplot, width.r=0.9, height.r=0.25)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/F4A.png", plot=F4A.ggplot, width.r=0.9, height.r=0.25)
}


print(F4A.ggplot)

image_read(snakemake@output[["main_png_filename"]]) %>% print




SF11.ggplot <- recurrently.edited.genes.simple.result.list %>%
    {
        .[["stat.test"]][, p.raw:=mapply(group1, disease.description, FUN=function(temp.type, temp.disease.description){INTERNAL.compute.phyper.enrichment.of.edited.in.maternal.decay.per.disease(.[["data"]], temp.type, temp.disease.description)})][, p.adj:=p.adjust(p.raw, method="BH")][, p.adj.to.plot:=scientific(p.adj)]
        .
    } %>%
    {
        ggbarplot(data=.[["data"]], x="Annotation.class.pasted.description", y="count.log10", fill="cluster.description", position=position_dodge()) +
            stat_pvalue_manual(data=.[["stat.test"]], label="p.adj.to.plot", size=3) +
            labs(x="", y="log10(# genes)", fill="") +
            theme_pubr(base_size=10) +
            theme(axis.text.x=element_text(angle=45, hjust=1)) +
            scale_y_continuous(breaks=waiver(), n.breaks=7, limits=c(0, max(.[["data"]][, count.log10]) * 1.5)) +
            facet_grid(~disease.description)
    }
SF11.ggplot

saveRDS(SF11.ggplot, snakemake@output[["supp1_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(SF11.ggplot, "report/201218-fifth-dataset/201221-fifth-phenotype-collection/GSE133854.all/SF11.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["supp1_png_filename"]], plot=SF11.ggplot, width.r=0.4, height.r=0.4)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/SF11.png", plot=SF11.ggplot, width.r=0.4, height.r=0.4)
}


print(SF11.ggplot)

image_read(snakemake@output[["supp1_png_filename"]]) %>% print
