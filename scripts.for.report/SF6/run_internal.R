#' ---
#' title: "SF6 summary"
#' ---

library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("foreach")
library("magick")
source("./scripts/common/ggpubr.A4.R")


subset.recurrent.edits.only.dt <- fread(snakemake@input[["subset_recurrent_edits_only_dt_txt_gz_filename"]])

if (FALSE){
    subset.recurrent.edits.only.dt <- fread("result/A02_3__check_recurrence_profile_for_a_subset_of_samples/201218-fifth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.dt.txt.gz")
}


REDIPortal.2021.hg38.editome.dt <- fread(snakemake@input[["REDIPortal_2021_hg38_editome_txt_gz_filename"]])

if (FALSE){
    REDIPortal.2021.hg38.editome.dt <- fread("external/REDIPortal/hg38/TABLE1_hg38.txt.gz")
}

subset.name <- snakemake@wildcards[["SUBSET_NAME"]]
if (FALSE){
    subset.name <- "all.normal.samples"
}

subset.recurrent.edits.only.dt[, detected.in.REDIPortal.2021.hg38:=FALSE][paste(sep="", CHROM, "_", POS) %in% REDIPortal.2021.hg38.editome.dt[, paste(sep="", Region, "_", Position)] == TRUE, detected.in.REDIPortal.2021.hg38:=TRUE]

fwrite(subset.recurrent.edits.only.dt, paste(sep="", snakemake@params[["result_directory"]], "/subset.recurrent.edits.only.with.REDIPortal.2021.hg38.status.dt.txt.gz"))


## discard unsequenced samples
for.plot.dt <- subset.recurrent.edits.only.dt %>%
    {.[is.na(AC)==FALSE, list(count=.N), list(SAMPLE, gse, stage, detected.in.REDIPortal.2021.hg38)]} %>%
    {.[, percentage:=count/sum(count), list(SAMPLE, gse, stage)][, percentage.to.plot:=percentage*100]} %>%
    {dcast(data=., formula=SAMPLE + gse + stage ~ detected.in.REDIPortal.2021.hg38, value.var="percentage.to.plot", fill=0)[, percentage.to.plot:=`FALSE`]} %>%
    {.[stage %in% c('oocyte.GV', 'oocyte.MI', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell', 'morula', 'blastocyst.early', 'blastocyst.middle', 'blastocyst.late', 'trophoblast', 'TE', 'ICM', 'CTB', 'STB', 'MTB', 'epiblast', 'hypoblast')][stage %in% c('oocyte.GV', 'oocyte.MI', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell'), stage.interval:="till 8-cell"][is.na(stage.interval)==TRUE, stage.interval:="morula and afterwards"][, stage.interval.to.plot:=factor(stage.interval, levels=c("till 8-cell", "morula and afterwards"))]} %>%
    {.[, gse.combined:=gse][gse %in% c("GSE133854", "GSE36552", "GSE136447", "GSE130289", "GSE73211", "GSE119324", "GSE126488", "GSE62772", "GSE64417", "GSE95477", "GSE125616", "GSE65481")== FALSE, gse.combined:="others"][, gse.combined.sample.size.per.stage:=.N, list(stage, gse.combined)][, gse.combined.with.sample.size.per.stage:=paste(sep="", gse.combined, "/", gse.combined.sample.size.per.stage)]} %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])][, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))]}


SF6.ggplot <- for.plot.dt %>%
    {ggboxplot(data=., x="stage.description.ordered", y="percentage.to.plot", fill="#00BFC4") + geom_hline(data=data.table(is.highly.recurrent.description="highly\nrecurrent"), aes(yintercept=20), color="#2E3436", linetype="dashed") + labs(x="", y="% not detected in\nREDIPortal hg38", fill="detected in REDIportal hg38") + theme_pubr(base_size=10) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous(limits=c(0, 100), breaks=seq(0, 100, 20)) + facet_grid(~stage.interval.to.plot, scales="free", space="free")}


saveRDS(SF6.ggplot, snakemake@output[["main_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(SF6.ggplot, "report/201218-fifth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF6.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["main_png_filename"]], plot=SF6.ggplot, width.r=0.45, height.r=0.25)
if (FALSE){
    ggsave.A4(filename="report/201218-fifth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF6.png", plot=SF6.ggplot, width.r=0.45, height.r=0.25)
}


print(SF6.ggplot)

image_read(snakemake@output[["main_png_filename"]]) %>% print


