#' ---
#' title: "F2A and SF3 summary"
#' ---

library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("magick")
library("foreach")
source("./scripts/common/ggpubr.A4.R")


subset.site.recurrence.comparison.dt <- fread(snakemake@input[["subset_site_recurrence_comparison_dt_txt_gz_filename"]])

if (FALSE){
    subset.site.recurrence.comparison.dt <- fread("result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.site.recurrence.comparison.dt.txt.gz")
}

subset.recurrent.edits.only.dt <- fread(snakemake@input[["subset_recurrent_edits_only_dt_txt_gz_filename"]])

if (FALSE){
    subset.recurrent.edits.only.dt <- fread("result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.dt.txt.gz")
}


subset.name <- snakemake@wildcards[["SUBSET_NAME"]]
if (FALSE){
    subset.name <- "all.normal.samples"
}

if (subset.name == "all.normal.samples"){
    subset.site.recurrence.comparison.dt[, stage:=sub(pattern="@.*", replacement="", x=group)]
}



subset.site.recurrence.comparison.max.site.occurrence.dt <- subset.site.recurrence.comparison.dt %>%
    {.[, list(site.occurrence.max=max(site.occurrence.for.this.group)), list(CHROM, POS)]} %>%
    {.[, site.occurrence.max.clipped:=as.character(site.occurrence.max)][site.occurrence.max>=5, site.occurrence.max.clipped:=">=5"]} %>%
    {.[, site.occurrence.max.clipped.ordered:=factor(site.occurrence.max.clipped, levels=c(1:4, ">=5"))]} %>%
    {.[, list(count=.N), list(site.occurrence.max.clipped.ordered)]} %>%
    {.[, count.normalized:=count/1e5]}


F2A.upper.panel.dt <- subset.site.recurrence.comparison.dt %>%
    {.[, site.occurrence.for.this.group.clipped:=as.character(site.occurrence.for.this.group)][site.occurrence.for.this.group>=5, site.occurrence.for.this.group.clipped:=">=5"]} %>%
    {.[, site.occurrence.for.this.group.clipped.ordered:=factor(site.occurrence.for.this.group.clipped, levels=rev(c(1:4, ">=5")))]} %>%
    {.[, list(total.count=sum(.N)), list(stage, site.occurrence.for.this.group.clipped.ordered)][, total.percentage.to.plot:=total.count/sum(total.count)*100, list(stage)]} %>%
    {.[, panel:="upper"]}

F2A.lower.panel.dt <- subset.recurrent.edits.only.dt %>%
    {.[, list(CHROM, POS, stage)]} %>%
    unique %>%
    {.[, list(log10.count.recurrent.edits=log10(.N)), list(stage)]} %>%
    {.[, panel:="lower"]}

F2A.dt <- list(F2A.upper.panel.dt, F2A.lower.panel.dt) %>% rbindlist(fill=TRUE) %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])][, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))]} %>%
    {.[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell', 'morula', 'blastocyst.late', 'trophoblast', 'ICM', 'hESC', 'CTB', 'STB', 'EVT', 'MTB', 'epiblast', 'hypoblast')]}




F2A.ggplot <- F2A.dt %>%
    {
        ggplot() + 
            geom_bar(data=.[panel=="upper"], mapping=aes(x=stage.description.ordered, y=total.percentage.to.plot, fill=site.occurrence.for.this.group.clipped.ordered), color="white", stat="identity") +
            geom_bar(data=.[panel=="lower"], mapping=aes(x=stage.description.ordered, y=log10.count.recurrent.edits), stat="identity") +
            geom_hline(data=data.table(panel="lower", yintercept=log10(2000)), mapping=aes(yintercept=yintercept), color="red", linetype="dashed") +
            geom_text(data=data.table(panel="lower", stage.description.ordered=16, count.log10=log10(2000)+0.6, text="count\n=2000"), aes(x=stage.description.ordered, y=count.log10, label=text), color="red") +
            labs(x="", y="", fill="recurrence") +
            facet_grid(factor(c("upper"="recurrence distribution", "lower"="lg(# recurrent edits)")[panel])~., scales="free_y", switch="y") +
            theme_pubr(base_size=10) +
            scale_fill_manual(values=c("#E64B35", "#F39B7F", "#4DBBD5",  "#00A087", "#3C5488")) +
            theme(text=element_text(size=10), axis.text.x=element_text(angle=45, hjust=1))};
F2A.ggplot



saveRDS(F2A.ggplot, snakemake@output[["ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(F2A.ggplot, "report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F2A.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["png_filename"]], plot=F2A.ggplot, width.r=0.45, height.r=0.45)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F2A.png", plot=F2A.ggplot, width.r=0.45, height.r=0.45)
}

image_read(snakemake@output[["png_filename"]]) %>% print






SF3.ggplot <- subset.recurrent.edits.only.dt %>%
    {.[is.na(AC)==FALSE, list(CHROM, POS, SAMPLE, gse, stage)]} %>%
    unique %>%
    {.[, list(count.log10=log10(.N)), list(SAMPLE, gse, stage)]} %>%
    {.[stage %in% c('oocyte.GV', 'oocyte.MI', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell', 'morula', 'blastocyst.early', 'blastocyst.middle', 'blastocyst.late', 'trophoblast', 'TE', 'ICM', 'CTB', 'STB', 'MTB', 'epiblast', 'hypoblast')]} %>%
    {.[stage %in% c('oocyte.GV', 'oocyte.MI', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell'), stage.interval:="till 8-cell"][is.na(stage.interval)==TRUE, stage.interval:="morula and afterwards"][, stage.interval.to.plot:=factor(stage.interval, levels=c("till 8-cell", "morula and afterwards"))]} %>%
    {.[, gse.combined:=gse][gse %in% c("GSE133854", "GSE36552", "GSE136447", "GSE130289", "GSE73211", "GSE119324", "GSE126488", "GSE62772", "GSE64417", "GSE95477", "GSE125616", "GSE65481")== FALSE, gse.combined:="others"][, gse.combined.sample.size.per.stage:=.N, list(stage, gse.combined)][, gse.combined.with.sample.size.per.stage:=paste(sep="", gse.combined, "/", gse.combined.sample.size.per.stage)]} %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])][, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))]} %>%
    {
        foreach(temp.stage.interval.to.plot=c("till 8-cell", "morula and afterwards")) %do% {
            ggboxplot(data=.[stage.interval.to.plot==temp.stage.interval.to.plot], x="gse.combined.with.sample.size.per.stage", y="count.log10") + geom_hline(yintercept=log10(2000), color="red", linetype="dashed") + facet_grid(~stage.description.ordered, scales="free", space="free") + labs(x="dataset / sample size", y="log10(count) of highly\nrecurrent editing sites") + theme_pubr(base_size=10, legend="none") + scale_y_continuous(limits=c(0, 4.5)) + ggtitle(label=temp.stage.interval.to.plot) + theme(axis.text.x=element_text(angle=45, hjust=1))
        }
    } %>%
    {ggarrange(plotlist=., nrow=2)}
SF3.ggplot


saveRDS(SF3.ggplot, snakemake@output[["supp1_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(SF3.ggplot, "report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF3.ggplot.RDS")
}



ggsave.A4(filename=snakemake@output[["supp1_png_filename"]], plot=SF3.ggplot, width.r=1.2, height.r=0.9)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF3.png", plot=SF3.ggplot, width.r=1.2, height.r=0.9)
}



print(SF3.ggplot); image_read(snakemake@output[["supp1_png_filename"]]) %>% print

