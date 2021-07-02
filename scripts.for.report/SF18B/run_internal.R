#' ---
#' title: "SF18B summary"
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

site.recurrence.comparison.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt %>%
    {.[, count.of.sequenced.samples.for.this.sample.group:=.N, list(stage, is.normal, total.sample.count.for.this.sample.group, CHROM, POS)][is.na(site.occurrence)==FALSE, list(CHROM, POS, stage, is.normal, total.sample.count.for.this.sample.group, site.occurrence, count.of.sequenced.samples.for.this.sample.group)]} %>%
    unique %>%
    {.[, `:=`(site.frequency.wrt.sequenced.samples=site.occurrence/count.of.sequenced.samples.for.this.sample.group, site.frequency.wrt.all.samples=site.occurrence/total.sample.count.for.this.sample.group)]}

fwrite(site.recurrence.comparison.dt, snakemake@output[["site_recurrence_comparison_dt_txt_gz_filename"]])




site.recurrence.diff.to.plot.dt <- site.recurrence.comparison.dt[site.occurrence >= 5 & site.frequency.wrt.all.samples<1][, `:=`(site.frequency.diff.to.plot=site.frequency.wrt.sequenced.samples*100 - site.frequency.wrt.all.samples*100)]


options(scipen = 999)

SF18B.ggplot <- gghistogram(data=site.recurrence.diff.to.plot.dt, x="site.frequency.diff.to.plot", palette="npg") + labs(x="increase in site recurrence\nif discarding unsequenced sites per sample(%)", y="count of such editing sites\n(log10 scale)") + scale_y_log10() + geom_hline(yintercept=1000, color="red", linetype="dashed")
SF18B.ggplot






saveRDS(SF18B.ggplot, snakemake@output[["main_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(SF18B.ggplot, "report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF18B.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["main_png_filename"]], plot=SF18B.ggplot, width.r=0.9, height.r=0.3)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF18B.png", plot=SF18B.ggplot, width.r=0.9, height.r=0.3)
}


print(SF18B.ggplot)
