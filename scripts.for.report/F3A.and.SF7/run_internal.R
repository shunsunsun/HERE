#' ---
#' title: "F3A and SF7 summary"
#' ---

library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("foreach")
library("magick")
source("./scripts/common/ggpubr.A4.R")


subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt <- fread(snakemake@input[["subset_recurrent_edits_only_with_snpEff_annotation_on_valid_genes_only_dt_txt_gz_filename"]])

if (FALSE){
    subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt <- fread("result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")
}


subset.name <- snakemake@wildcards[["SUBSET_NAME"]]
if (FALSE){
    subset.name <- "all.normal.samples"
}

for.plot.dt <- subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt %>%
    {.[is.na(AC)==FALSE]} %>%
    {.[, list(Gene_ID, Gene_Name, SAMPLE, gse, stage, Annotation.class)]} %>% unique %>%
    {.[, list(gene.count=.N), list(SAMPLE, gse, stage, Annotation.class)]} %>%
    {.[, Annotation.class.description:=c('exonic.or.splicing.related'='exonic', 'purely.intronic'='intronic')[Annotation.class]]} %>%
    {.[stage %in% c('oocyte.GV', 'oocyte.MI', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell', 'morula', 'blastocyst.early', 'blastocyst.middle', 'blastocyst.late', 'trophoblast', 'TE', 'ICM', 'CTB', 'STB', 'MTB', 'epiblast', 'hypoblast')][stage %in% c('oocyte.GV', 'oocyte.MI', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell'), stage.interval:="till 8-cell"][is.na(stage.interval)==TRUE, stage.interval:="morula and afterwards"][, stage.interval.to.plot:=factor(stage.interval, levels=c("till 8-cell", "morula and afterwards"))]} %>%
    {.[, gse.combined:=gse][gse %in% c("GSE133854", "GSE36552", "GSE136447", "GSE130289", "GSE73211", "GSE119324", "GSE126488", "GSE62772", "GSE64417", "GSE95477", "GSE125616", "GSE65481")== FALSE, gse.combined:="others"][, gse.combined.sample.size.per.stage:=.N, list(stage, gse.combined, Annotation.class)][, gse.combined.with.Annotation.class.description.and.sample.size.per.stage:=paste(sep="", Annotation.class.description, ":", gse.combined, "/", gse.combined.sample.size.per.stage)] } %>%
    {.[, gse.combined.with.Annotation.class.description.and.sample.size.per.stage.ordered:=factor(gse.combined.with.Annotation.class.description.and.sample.size.per.stage, levels=sort(unique(gse.combined.with.Annotation.class.description.and.sample.size.per.stage)))]} %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])]}


F3A.ggplot <- for.plot.dt %>%
    {ggboxplot(data=., x="stage.description.ordered", y="gene.count", fill="Annotation.class.description")  +
         labs(x="", y="# genes affected by\nrecurrent edits per sample", fill="") +
         theme(axis.text.x=element_text(angle=45, hjust=1))}
F3A.ggplot



saveRDS(F3A.ggplot, snakemake@output[["main_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(F3A.ggplot, "report/201218-fifth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F3A.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["main_png_filename"]], plot=F3A.ggplot, width.r=0.45, height.r=0.4)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F3A.png", plot=F3A.ggplot, width.r=0.45, height.r=0.4)
}


print(F3A.ggplot)

image_read(snakemake@output[["main_png_filename"]]) %>% print



SF7.ggplot <- foreach(temp.stage.interval.to.plot=c("till 8-cell", "morula and afterwards"), temp.filename.label=c("till.8.cell", "morula.and.afterwards")) %do% {
    for.plot.dt %>%
        {.[stage.interval.to.plot==temp.stage.interval.to.plot]} %>%
        {ggboxplot(data=., x="gse.combined.with.Annotation.class.description.and.sample.size.per.stage.ordered", y="gene.count", fill="Annotation.class.description") + labs(x="type: dataset / sample size", y="# genes affected by\nrecurrent edits per sample", fill="") + theme_pubr(base_size=10) + facet_grid(~stage.description.ordered, scales="free")  + ggtitle(label=temp.stage.interval.to.plot) + theme(axis.text.x=element_text(angle=90, vjust=0.5))}
} %>%
    {ggarrange(plotlist=., nrow=2)}
SF7.ggplot


saveRDS(SF7.ggplot, snakemake@output[["supp1_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(SF7.ggplot, "report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF7.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["supp1_png_filename"]], plot=SF7.ggplot, width.r=1.2, height.r=0.8)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF7.png", plot=SF7.ggplot, width.r=1.2, height.r=0.8)
}


print(SF7.ggplot)

image_read(snakemake@output[["supp1_png_filename"]]) %>% print
