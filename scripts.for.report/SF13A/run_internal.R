#' ---
#' title: "SF13A summary"
#' ---

library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
source("./scripts/common/ggpubr.A4.R")


merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt <- fread(snakemake@input[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_merged_with_coverage_dt_txt_gz_filename"]])

not.used.variable <- '
merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt <- fread("result/S52_3__mark_unsequenced_editing_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt.txt.gz")
'


edit.and.sample.unsequenced.status.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt %>%
    {.[, list(count.of.sequenced.samples.for.this.sample.group=.N), list(stage, is.normal, total.sample.count.for.this.sample.group, CHROM, POS)]} %>%
    {.[, count.of.unsequenced.samples.for.this.sample.group:=total.sample.count.for.this.sample.group - count.of.sequenced.samples.for.this.sample.group]} %>%
    {.[, percentage.of.unsequenced.samples.for.this.sample.group:=count.of.unsequenced.samples.for.this.sample.group/total.sample.count.for.this.sample.group]} %>%
    {.[, count.of.unsequenced.samples.for.this.sample.group.discretized:={temp.a <- NA; if (length(unique(count.of.unsequenced.samples.for.this.sample.group)) > 1) {temp.a <- cut(count.of.unsequenced.samples.for.this.sample.group, breaks=c(seq(0, total.sample.count.for.this.sample.group, length.out=11), Inf), right=FALSE)} else temp.a <- paste(sep="", "{", count.of.unsequenced.samples.for.this.sample.group, "}"); temp.a}, list(stage, is.normal, total.sample.count.for.this.sample.group)]} %>%
    {.[, percentage.of.unsequenced.samples.for.this.sample.group.discretized:={temp.a <- NA; if (length(unique(percentage.of.unsequenced.samples.for.this.sample.group)) > 1) {temp.a <- cut(percentage.of.unsequenced.samples.for.this.sample.group, breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, Inf), right=FALSE)} else temp.a <- paste(sep="", "{", percentage.of.unsequenced.samples.for.this.sample.group, "}"); temp.a}, list(stage, is.normal, total.sample.count.for.this.sample.group)]} ## cannot use seq(0, 1, length.out=11) here; otherwise the resulting percentage regions are not consistent with those of the counts


fwrite(edit.and.sample.unsequenced.status.dt, snakemake@output[["edit_and_sample_unsequenced_status_dt_txt_gz_filename"]])
'
fwrite(edit.and.sample.unsequenced.status.dt, "report/A01_7__check_unsequenced_status_profile/210215-sixth-dataset/201221-fifth-phenotype-collection/edit.and.sample.unsequenced.status.dt.txt.gz")
' %>% as.null

edit.and.sample.unsequenced.status.edit.histogram.per.sample.group.dt <- edit.and.sample.unsequenced.status.dt[, list(count.of.editing.sites=.N), list(stage, is.normal, count.of.unsequenced.samples.for.this.sample.group.discretized, percentage.of.unsequenced.samples.for.this.sample.group.discretized, total.sample.count.for.this.sample.group)] %>% {.[, rbindlist(list(.SD, data.table(count.of.unsequenced.samples.for.this.sample.group.discretized=paste(sep="", "{", total.sample.count.for.this.sample.group, "}"), percentage.of.unsequenced.samples.for.this.sample.group.discretized="{1}", count.of.editing.sites={edit.and.sample.unsequenced.status.dt[, list(CHROM, POS)] %>% unique %>% nrow} - .SD[, sum(count.of.editing.sites)]))), list(stage, is.normal, total.sample.count.for.this.sample.group)]}


fwrite(edit.and.sample.unsequenced.status.edit.histogram.per.sample.group.dt, snakemake@output[["edit_and_sample_unsequenced_status_edit_histogram_per_sample_group_dt_txt_gz_filename"]])
'
fwrite(edit.and.sample.unsequenced.status.edit.histogram.per.sample.group.dt, "result/A01_7__check_unsequenced_status_profile/210215-sixth-dataset/201221-fifth-phenotype-collection/edit.and.sample.unsequenced.status.edit.histogram.per.sample.group.dt.txt.gz")
' %>% as.null

edit.and.sample.unsequenced.status.edit.histogram.per.sample.group.to.plot.dt <- edit.and.sample.unsequenced.status.edit.histogram.per.sample.group.dt  %>%
    {.[, `:=`(is.normal.description=c("Abnormal samples", "Normal samples")[is.normal+1])]} %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])][, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))]}


## simple percentage distribution plot

SF13A.ggplot <- ggbarplot(data=edit.and.sample.unsequenced.status.edit.histogram.per.sample.group.to.plot.dt[is.normal==TRUE & percentage.of.unsequenced.samples.for.this.sample.group.discretized=="[0.9,1)"], x="stage.description.ordered", y="count.of.editing.sites", palette="npg")  +
    labs(x="stage", y="count of editing sites detected in >=1\nbut unsequenced in >= 90% normal samples") +
    theme(axis.text.x=element_text(angle=45, hjust=1))
SF13A.ggplot



saveRDS(SF13A.ggplot, snakemake@output[["main_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(SF13A.ggplot, "report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF13A.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["main_png_filename"]], plot=SF13A.ggplot, width.r=0.9, height.r=0.45)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF13A.png", plot=SF13A.ggplot, width.r=0.9, height.r=0.45)
}


print(SF13A.ggplot)
