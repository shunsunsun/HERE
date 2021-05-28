#' ---
#' title: "SF1A and SF1B summary"
#' ---

library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("foreach")
library("magick")
library("scales")
source("./scripts/common/ggpubr.A4.R")


merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt <- fread(snakemake@input[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_with_event_summary_dt_txt_gz_filename"]])

if (FALSE) {
    merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt <- fread("result/S71_5__filter_for_A_to_G_sites_for_control/210203-GSE144296.A375-RNA-with-DNA-37-37/210203-GSE144296.A375-RNA-with-DNA-37-37/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt.txt.gz")
}

subset.name <- snakemake@wildcards[["SUBSET_NAME"]]
if (FALSE){
    subset.name <- "all.samples"
}

if (subset.name != "all.samples"){
    stop("subset.name must be 'all.samples'")
}

event.summary.per.sample.with.phenotype.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt[, list(count=.N), list(SAMPLE, stage, is.normal, event.summary)]


SF1A.ggplot <- event.summary.per.sample.with.phenotype.dt %>%
    {.[grepl(";", event.summary)==FALSE][, total.count:=sum(count), list(SAMPLE)][, percentage:=count/total.count][, percentage.to.plot:=percentage*100]} %>%
    {.[event.summary %in% c('A>G', 'A>G;T>C'), list(percentage.to.plot.combined=sum(percentage.to.plot)), list(SAMPLE, stage, is.normal)]} %>%
    {ggboxplot(., x="stage", y="percentage.to.plot.combined", pallete="npg")  + labs(x="", y="% A>G of\nall definite\nvariants", fill="") + scale_y_continuous(breaks=seq(0, 100, 20), limits=c(0, 100))}
SF1A.ggplot





ggsave.A4(filename=snakemake@output[["SF1A_png_filename"]], plot=SF1A.ggplot, width.r=0.2, height.r=0.15)
if (FALSE){
    ggsave.A4(filename="report/210203-GSE144296.A375-RNA-with-DNA-37-37/210203-GSE144296.A375-RNA-with-DNA-37-37/all.samples/SF1A.png", plot=SF1A.ggplot, width.r=0.2, height.r=0.15)
}


print(SF1A.ggplot)

image_read(snakemake@output[["SF1A_png_filename"]]) %>% print




SF1B.ggplot <- event.summary.per.sample.with.phenotype.dt %>%
    {.[, total.count:=sum(count), list(SAMPLE)][, percentage:=count/total.count][, percentage.to.plot:=percentage*100]} %>%
    {.[event.summary %in% c('A>G', 'A>G;T>C'), list(percentage.to.plot.combined=sum(percentage.to.plot)), list(SAMPLE, stage, is.normal)]} %>%
    {ggboxplot(., x="stage", y="percentage.to.plot.combined", pallete="npg")  + labs(x="", y="% A>G of\nall variants", fill="") + scale_y_continuous(breaks=seq(0, 100, 20), limits=c(0, 100))}
SF1B.ggplot



ggsave.A4(filename=snakemake@output[["SF1B_png_filename"]], plot=SF1B.ggplot, width.r=0.2, height.r=0.15)
if (FALSE){
    ggsave.A4(filename="report/210203-GSE144296.A375-RNA-with-DNA-37-37/210203-GSE144296.A375-RNA-with-DNA-37-37/all.samples/SF1B.png", plot=SF1B.ggplot, width.r=0.2, height.r=0.15)
}


print(SF1B.ggplot)

image_read(snakemake@output[["SF1B_png_filename"]]) %>% print

