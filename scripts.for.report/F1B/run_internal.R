#' ---
#' title: "F1B summary"
#' ---

library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("foreach")
library("magick")
library("scales")
library("iterators")
source("./scripts/common/ggpubr.A4.R")

## first row, sample profile (sequenced only; could have no edits)

phenotype.output.at.gsm.level.dt <- fread(snakemake@input[["phenotype_output_at_gsm_level_dt_filename"]])
if (FALSE) { 
    phenotype.output.at.gsm.level.dt <- fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt")
}

sample.sequenced.dt <- rbindlist(lapply(readLines(paste(sep="", "external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/", snakemake@wildcards[['DATASET_COLLECTION_NAME']])), function(temp.DATASET_RNA_EDITING_NAME){return(unique(fread(paste(sep="", "external/DATASET_RNA_EDITING_NAME_DIRECTORY/", temp.DATASET_RNA_EDITING_NAME))))}))

if (FALSE) {
    sample.sequenced.dt <- rbindlist(lapply(readLines(paste(sep="", "external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/", "210215-sixth-dataset")), function(temp.DATASET_RNA_EDITING_NAME){return(unique(fread(paste(sep="", "external/DATASET_RNA_EDITING_NAME_DIRECTORY/", temp.DATASET_RNA_EDITING_NAME))))}))
}

phenotype.output.at.gsm.level.sequenced.only.dt <- phenotype.output.at.gsm.level.dt[gsm %in% sample.sequenced.dt[, SAMPLE_NAME]]

sample.profile.dt <- phenotype.output.at.gsm.level.sequenced.only.dt %>%
    {.[, list(metric="log2.sample.count", value=log2(.N)), list(is.normal, stage)]}


'
temp.sample.profile.ggplot <- ggbarplot(data=sample.profile.to.plot.dt, x="stage.description.ordered", y="log2.count", fill="is.normal.description", position=position_dodge(), palette="npg") +
    geom_hline(yintercept=log2(100), linetype="dashed") +
    labs(x="stage", y="log2(total number of samples) ", fill="") +
    facet_grid(~category.ordered, scales="free", space="free") +
    theme(text=element_text(size=10), axis.text.x=element_text(angle=90, vjust=0.5)); temp.sample.profile.ggplot
'
## second row, count of edits identified

merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt <- fread(snakemake@input[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_txt_gz_filename"]])

if (FALSE) {
    merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")
}

total.count.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt[, list(ID, stage, is.normal)] %>% unique %>% {.[, list(metric="log10.total.count.of.editing.sites", value=log10(.N)), list(stage, is.normal)]}

'
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])][, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))]} %>%
    {ggbarplot(., x="stage.description.ordered", y="total.count.of.editing.sites", fill="is.normal.description", position=position_dodge(), pallete="npg") +
         labs(x="stage", y="# editing sites", fill="") +
         facet_grid(~category.ordered, scales="free", space="free") +
         theme(text=element_text(size=10), axis.text.x=element_text(angle=90, vjust=0.5))}; temp.total.count.ggplot
'

## third row, A-to-G ratio

merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt <- fread(snakemake@input[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_with_event_summary_dt_txt_gz_filename"]])

if (FALSE) {
    merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt.txt.gz")
}

A.to.G.ratio.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt[, list(count=.N), list(SAMPLE, stage, is.normal, disease, treatment, maternal.age, developmental.day, cell.line, event.summary)] %>%
    {.[, total.count:=sum(count), list(SAMPLE)][, percentage:=count/total.count][, percentage.to.plot:=percentage*100]} %>%
    {.[event.summary %in% c('A>G', 'A>G;T>C'), list(metric="A.to.G.ratio", value=sum(percentage.to.plot)), list(SAMPLE, stage, is.normal)]}

'
    {ggboxplot(., x="stage.description.ordered", y="percentage.to.plot.combined", fill="is.normal.description", pallete="npg") + geom_hline(yintercept=80, color="red", linetype="dashed")  + labs(x="stage", y="% A>G of all variants", fill="") + facet_grid(~category.ordered, scales="free", space="free") + scale_y_continuous(breaks=seq(0, 100, 20)) + theme(text=element_text(size=10), axis.text.x=element_text(angle=90, vjust=0.5))}; temp.A.to.G.ratio.ggplot
'

## fourth row, Alu ratio


Alu.ratio.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt[, list(count=.N), list(SAMPLE, stage, is.normal, disease, treatment, maternal.age, developmental.day, cell.line, SUBSET)] %>%
    {.[, total.count:=sum(count), list(SAMPLE)][, percentage:=count/total.count][, SUBSET.description:=factor(c("Alu"="rep.\nAlu", "RepNOTAlu"="rep.\nothers", "nonRep"="not\nrep.")[SUBSET], levels=c("rep.\nAlu", "rep.\nothers", "not\nrep."))][, percentage.to.plot:=percentage*100]} %>%
    {.[SUBSET=='Alu', list(SAMPLE, stage, is.normal, metric="Alu.ratio", value=percentage.to.plot)]}




mixed.dt <- list(
    data.table(sample.profile.dt, SAMPLE=NA), data.table(total.count.dt, SAMPLE=NA), A.to.G.ratio.dt, Alu.ratio.dt
) %>% rbindlist(use.names=TRUE) %>%
    {.[, is.normal.description:=c("Abnormal samples", "Normal samples")[is.normal+1]]} %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description.no.newline, levels=temp.stage.dt[, stage.description.no.newline])][, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))]}

F1B.ggplot <- ggplot() +
    geom_bar(data=mixed.dt[metric %in% c("log2.sample.count", "log10.total.count.of.editing.sites")], mapping=aes(x=stage.description.ordered, y=value, fill=is.normal.description), stat="identity", position="dodge") +
    geom_boxplot(data=mixed.dt[metric %in% c("A.to.G.ratio", "Alu.ratio")], mapping=aes(x=stage.description.ordered, y=value, fill=is.normal.description)) +
   geom_blank(data=data.table(
                  metric=rep(c("log2.sample.count", "log10.total.count.of.editing.sites", "A.to.G.ratio", "Alu.ratio"), 2),
                  value=c(0, 0, 0, 0, 10, 5, 100, 100)
              ), aes(y=value)) + 
    geom_hline(
        data=data.table(
            yintercept.to.plot=c(log2(100), 80, 80),
            metric=c("log2.sample.count", "A.to.G.ratio", "Alu.ratio")
        ),
        mapping=aes(yintercept=yintercept.to.plot),
        color="red", linetype="dashed") +
    facet_grid(c("log2.sample.count"="log2(# samples)", "log10.total.count.of.editing.sites"="lg(# editing sites)", "A.to.G.ratio"="% A-to-G", "Alu.ratio"="% Alu")[metric] %>% factor(levels=c("log2(# samples)", "lg(# editing sites)", "% A-to-G", "% Alu"))~category.ordered, scales="free", space="free_x", switch="y") +
    labs(x="", y="", fill="") +
    theme_pubr(base_size=10) + theme(axis.text.x=element_text(angle=90, vjust=0.5)) 

    
saveRDS(F1B.ggplot, snakemake@output[["main_ggplot_RDS_filename"]])


ggsave.A4(filename=snakemake@output[["main_png_filename"]], plot=F1B.ggplot, width.r=0.55, height.r=0.7)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F1B.png", plot=F1B.ggplot, width.r=0.55, height.r=0.7)
}


print(F1B.ggplot)

image_read(snakemake@output[["main_png_filename"]]) %>% print

