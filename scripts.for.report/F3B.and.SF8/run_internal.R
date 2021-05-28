#' ---
#' title: "F3B and SF8 summary"
#' ---

library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("foreach")
library("magick")
library("ggalluvial")
source("./scripts/common/ggpubr.A4.R")


subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt <- fread(snakemake@input[["subset_recurrent_edits_only_with_snpEff_annotation_on_valid_genes_only_dt_txt_gz_filename"]])

if (FALSE){
    subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt <- fread("result/A02_3__check_recurrence_profile_for_a_subset_of_samples/201218-fifth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")
}


subset.name <- snakemake@wildcards[["SUBSET_NAME"]]
if (FALSE){
    subset.name <- "all.normal.samples"
}

for.plot.dt <- subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt %>%
    {.[is.na(AC)==FALSE]} %>%
    {.[, list(Annotation.class.pasted=paste(collapse=";", Annotation.class %>% sort %>% unique)), list(Gene_ID, Gene_Name, SAMPLE, gse, stage)]} %>%
    ## compute % occurrence of each gene x Annotation.class.pasted
    ## note here that for total sample number per stage, `unique` must be used instead of `.N`
    {.[, total.sample.count:=length(unique(SAMPLE)), list(stage)]} %>%
    {.[, list(gene.occurrence=.N), list(Gene_ID, Gene_Name, stage, total.sample.count, Annotation.class.pasted)]} %>%
    {.[, gene.occurrence.pct:=gene.occurrence/total.sample.count]} %>%
    ## finish computing
    {.[, Annotation.class.pasted.description:=c('exonic.or.splicing.related'='exonic', 'purely.intronic'='intronic', 'exonic.or.splicing.related;purely.intronic'='mixed')[Annotation.class.pasted]]} %>%
    {.[stage %in% c('oocyte.GV', 'oocyte.MI', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell', 'morula', 'blastocyst.early', 'blastocyst.middle', 'blastocyst.late', 'trophoblast', 'TE', 'ICM', 'CTB', 'STB', 'MTB', 'epiblast', 'hypoblast')][stage %in% c('oocyte.GV', 'oocyte.MI', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell'), stage.interval:="till 8-cell"][is.na(stage.interval)==TRUE, stage.interval:="morula and afterwards"][, stage.interval.to.plot:=factor(stage.interval, levels=c("till 8-cell", "morula and afterwards"))]} %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])]}


## test
"
for.plot.dt[Gene_ID=='ENSG00000167065.13' & stage=='oocyte.MII']
        stage            Gene_ID Gene_Name total.sample.count
1: oocyte.MII ENSG00000167065.13    DUSP18                 30
2: oocyte.MII ENSG00000167065.13    DUSP18                 30
                      Annotation.class.pasted gene.occurrence
1: exonic.or.splicing.related;purely.intronic              24
2:                 exonic.or.splicing.related               6
   gene.occurrence.pct Annotation.class.pasted.description stage.interval
1:                 0.8                                both    till 8-cell
2:                 0.2                         E or S only    till 8-cell
   stage.interval.to.plot stage.description          category   stage.AtoG
1:            till 8-cell      oocyte (MII) Gamete (5 stages) oocyte (MII)
2:            till 8-cell      oocyte (MII) Gamete (5 stages) oocyte (MII)
   category.AtoG stage.description.ordered
1:        Gamete              oocyte (MII)
2:        Gamete              oocyte (MII)

subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt[is.na(AC)==FALSE][Gene_Name=='DUSP18' & stage=='oocyte.MII'][, list(Annotation.class, SAMPLE)] %>% unique %>% {.[, list(Annotation.class.pasted=paste(collapse=';', Annotation.class)), SAMPLE]} %>% {.[, .N, Annotation.class.pasted]}

                      Annotation.class.pasted  N
1: purely.intronic;exonic.or.splicing.related 24
2:                 exonic.or.splicing.related  6

" %>% invisible

fwrite(for.plot.dt, snakemake@output[["for_plot_dt_txt_gz_filename"]])

F3B.ggplot <- for.plot.dt %>%
    ## pick those with gene occurrence pct >= 80% samples
    ## the resulting `dt` conforms to alluvial form (i.e., <=1 record for each gene x stage)
    {.[gene.occurrence.pct>=0.8]} %>%
    {ggplot(.[stage.interval.to.plot=='till 8-cell'], aes(x = stage.description.ordered, stratum = Annotation.class.pasted.description, alluvium = Gene_ID, fill=Annotation.class.pasted.description)) + scale_fill_brewer(type = "qual", palette = "Set2") + geom_flow(stat = "flow", color = "darkgray", na.rm=FALSE) + geom_stratum() + labs(x="", y="# genes recurrently affected\nby recurrent edits", fill="") + theme_pubr(base_size=10) + theme(axis.text.x=element_text(angle=45, hjust=1))}




saveRDS(F3B.ggplot, snakemake@output[["main_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(F3B.ggplot, "report/201218-fifth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F3B.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["main_png_filename"]], plot=F3B.ggplot, width.r=0.45, height.r=0.4)
if (FALSE){
    ggsave.A4(filename="report/201218-fifth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F3B.png", plot=F3B.ggplot, width.r=0.45, height.r=0.4)
}


print(F3B.ggplot)

image_read(snakemake@output[["main_png_filename"]]) %>% print


SF8.ggplot <- for.plot.dt %>%
    {.[stage.interval.to.plot=='till 8-cell']} %>%
    {.[, list(Annotation.class.pasted.description.mixed=paste(collapse=";", Annotation.class.pasted.description %>% sort)), list(stage.description.ordered, Gene_ID)]} %>%
    {.[, Annotation.class.pasted.description.mixed.ordered:=factor(Annotation.class.pasted.description.mixed, levels=c("exonic", "intronic", "mixed", "exonic;mixed", "intronic;mixed", "exonic;intronic;mixed"))]} %>%
    {ggplot(., aes(x = stage.description.ordered, stratum = Annotation.class.pasted.description.mixed.ordered, alluvium = Gene_ID, fill=Annotation.class.pasted.description.mixed.ordered)) + scale_fill_brewer(type = "qual", palette = "Set2") + geom_flow(stat = "flow", color = "darkgray", na.rm=FALSE) + geom_stratum() + labs(x="", y="# genes affected\n(either recurrently or stochastically)\nby recurrent edits", fill="") + theme_pubr(base_size=10) + theme(axis.text.x=element_text(angle=45, hjust=1))}




saveRDS(SF8.ggplot, snakemake@output[["supp1_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(SF8.ggplot, "report/201218-fifth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF8.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["supp1_png_filename"]], plot=SF8.ggplot, width.r=0.9, height.r=0.8)
if (FALSE){
    ggsave.A4(filename="report/201218-fifth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF8.png", plot=SF8.ggplot, width.r=0.9, height.r=0.8)
}


print(SF8.ggplot)

image_read(snakemake@output[["supp1_png_filename"]]) %>% print
