#' ---
#' title: "F2B and SF4 summary"
#' ---

library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("foreach")
library("magick")
library("ggalluvial")
source("./scripts/common/ggpubr.A4.R")


subset.site.recurrence.comparison.CJ.dt <- fread(snakemake@input[["subset_site_recurrence_comparison_CJ_dt_txt_gz_filename"]])

if (FALSE){
    subset.site.recurrence.comparison.CJ.dt <- fread("result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.site.recurrence.comparison.CJ.dt.txt.gz")
}



subset.name <- snakemake@wildcards[["SUBSET_NAME"]]
if (FALSE){
    subset.name <- "all.normal.samples"
}



if (subset.name == "all.normal.samples"){
    subset.site.recurrence.comparison.CJ.dt[, stage:=sub(pattern="@.*", replacement="", x=group)]
}


sites.for.sankey.to.plot.dt <- subset.site.recurrence.comparison.CJ.dt %>%
    {.[, recurrence.type:=cut(site.occurrence.for.this.group, breaks=c(0, 1, 10, Inf), include.lowest=TRUE, right=FALSE) %>% {
        c("[0,1)"="0", "[1,10)"="[1,9]", "[10,Inf]"=">=10")[.]
    }]} %>%
    {.[, recurrence.type.ordered:=factor(recurrence.type, levels=c("0", "[1,9]", ">=10"))]} %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])][, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))]}

## only consider for each stage set all edits recurrent in at least one stage; otherwise we will run out of memory
F2B.dt <- foreach(
    temp.stage.set=list(
        c('oocyte.GV', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell'),
        c('8-cell', 'morula'),
        c('morula', 'blastocyst.late', 'trophoblast'),
        c('morula', 'blastocyst.late', 'ICM')
    ),
    temp.part=c(
        'till 8-cell',
        '8-cell to morula',
        'morula to TE',
        'morula to ICM'
    )) %do% {
        sites.for.sankey.to.plot.dt %>%
            {.[stage %in% temp.stage.set]} %>%
            {.[CHROM.and.POS %in% .[recurrence.type==">=10", CHROM.and.POS]][, part:=temp.part]}
    } %>% rbindlist %>% {.[, part.ordered:=factor(part, levels=unique(part))]}


F2B.upper.panel.ggplot <- F2B.dt %>%
    {ggplot(data=.[part=='till 8-cell'], mapping=aes(x = stage.description.ordered, stratum = recurrence.type.ordered, alluvium = CHROM.and.POS, fill=recurrence.type.ordered)) +
         scale_fill_brewer(type = "qual", palette = "Set2") +
         geom_flow(stat = "flow", color = "darkgray", na.rm=FALSE) +
         geom_stratum() +
         facet_wrap(~part.ordered, scales="free", nrow=1) + labs(x="", y="# editing sites", fill="# samples detected") + theme_pubr(base_size=10) + theme(axis.text.x=element_text(angle=45, hjust=1))}


F2B.lower.panel.ggplot <- F2B.dt %>%
    {ggplot(data=.[part %in% c('8-cell to morula', 'morula to TE', 'morula to ICM')], mapping=aes(x = stage.description.ordered, stratum = recurrence.type.ordered, alluvium = CHROM.and.POS, fill=recurrence.type.ordered)) +
         scale_fill_brewer(type = "qual", palette = "Set2") +
         geom_flow(stat = "flow", color = "darkgray", na.rm=FALSE) +
         geom_stratum() +
         facet_wrap(~part.ordered, scales="free", nrow=1) + labs(x="", y="# editing sites", fill="# samples detected") + theme_pubr(base_size=10) + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none")}


F2B.ggplot <- ggarrange(plotlist=list(F2B.upper.panel.ggplot, F2B.lower.panel.ggplot), ncol=1, heights=c(0.6, 0.4))

saveRDS(F2B.ggplot, snakemake@output[["main_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(F2B.ggplot, "report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F2B.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["main_png_filename"]], plot=F2B.ggplot, width.r=0.45, height.r=0.45)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F2B.png", plot=F2B.ggplot, width.r=0.45, height.r=0.45)
}

print(F2B.ggplot)

image_read(snakemake@output[["main_png_filename"]]) %>% print




SF4.ggplot <- foreach(
    temp.stage.set=list(
        c('trophoblast', 'CTB'),
        c('trophoblast', 'EVT'),
        c('trophoblast', 'STB'),
        c('trophoblast', 'MTB'),
        c('ICM', 'epiblast'),
        c('ICM', 'hypoblast')
    ),
    temp.part=c(
        'TE to CTB',
        'TE to EVT',
        'TE to STB',
        'TE to MTB',
        'ICM to epiblast',
        'ICM to hypoblast'
    )) %do% {
        sites.for.sankey.to.plot.dt %>%
            {.[stage %in% temp.stage.set]} %>%
            {.[CHROM.and.POS %in% .[recurrence.type==">=10", CHROM.and.POS]][, part:=temp.part]}
    } %>% rbindlist %>% {.[, part.ordered:=factor(part, levels=unique(part))]} %>%
    {ggplot(., aes(x = stage.description.ordered, stratum = recurrence.type.ordered, alluvium = CHROM.and.POS, fill=recurrence.type.ordered)) + scale_fill_brewer(type = "qual", palette = "Set2") + geom_flow(stat = "flow", color = "darkgray", na.rm=FALSE) + geom_stratum() + facet_wrap(~part.ordered, scales="free", nrow=2) + labs(x="", y="count of editing sites", fill="# samples detected") + theme_pubr(base_size=10) + theme(axis.text.x=element_text(angle=45, hjust=1))}

saveRDS(SF4.ggplot, snakemake@output[["supp1_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(SF4.ggplot, "report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF4.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["supp1_png_filename"]], plot=SF4.ggplot, width.r=0.9, height.r=0.6)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF4.png", plot=SF4.ggplot, width.r=0.9, height.r=0.6)
}

print(SF4.ggplot)

image_read(snakemake@output[["supp1_png_filename"]]) %>% print
