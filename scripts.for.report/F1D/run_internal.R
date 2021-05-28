#' ---
#' title: "F1D summary"
#' ---

library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("foreach")
library("magick")
library("scales")
source("./scripts/common/ggpubr.A4.R")

merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt <- fread(snakemake@input[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_txt_gz_filename"]])

if (FALSE) {
    merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt <- fread("result/S51_5__filter_for_A_to_G_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")
}


for.plot.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt %>%
    {.[POS==37519170 & CHROM=='chr20']} %>%
    {.[is.normal==TRUE & is.na(AC)==FALSE]} %>%
    {.[, list(CHROM, POS, SAMPLE, gse, stage, AF, AC, AN)]} %>%
    {.[stage %in% c("oocyte.GV", "oocyte.MII", "zygote", "2-cell", "4-cell", "8-cell", "morula", "hESC", "trophoblast", "ICM")]} %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=c("gene", "genomic DNA", "unedited RNA", "edited RNA", "unedited AA", "edited AA", temp.stage.dt[, stage.description]))][, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))]}

F1D.ggplot <- for.plot.dt %>%
    {
        ggplot(., aes(x=paste(sep="", CHROM, ":", POS), y=SAMPLE, fill=cut(x=AF, breaks=c(0, 0.1, 0.2, 0.4, 1)))) +
            geom_tile() +
            labs(x="", y="", fill="editing level") +
            geom_blank(data=data.table(POS=seq(37519140, 37519200, 4), SAMPLE="GSM3928470", stage.description.ordered=factor("morula", levels=levels(.[, stage.description.ordered])), CHROM="chr20", AF=as.numeric(NA))) +
            facet_grid(stage.description.ordered~., scales="free", switch="y") +
            theme_pubr() +
            theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.background = element_rect(color="black", fill="white"), strip.text.y.left = element_text(angle = 0), plot.margin=margin(t=75,r=0,b=0,l=0), legend.position=c(0.3, 1.1)) +
            guides(fill=guide_legend(nrow=1, title.hjust=0.5)) +
            scale_x_discrete(breaks=c("chr20:37519170")) +
            geom_text(data=data.table(CHROM="chr20", POS=37519170, SAMPLE=0, text=c("BLCAP", "--GCAATACAT->", "<-CGUUAUGUA--", "<-CGUU UGUA--", "G", "<--C--Y--M---", "<--C-- --M---", "C"), color=c("black", "black", "black", "black", "red", "black", "black", "red"), stage.description.ordered=factor(c("gene", "genomic DNA", "unedited RNA", "edited RNA", "edited RNA", "unedited AA", "edited AA", "edited AA"), levels=levels(.[, stage.description.ordered])), AF=as.numeric(NA)), mapping=aes(x=paste(sep="", CHROM, ":", POS), y=SAMPLE, label=text, color=color), family='mono') + scale_color_manual(values=c("black", "red")) + guides(color="none")
    }
F1D.ggplot





saveRDS(F1D.ggplot, snakemake@output[["main_ggplot_RDS_filename"]])


ggsave.A4(filename=snakemake@output[["main_png_filename"]], plot=F1D.ggplot, width.r=0.4, height.r=0.55)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F1D.png", plot=F1D.ggplot, width.r=0.4, height.r=0.55)
}


print(F1D.ggplot)

image_read(snakemake@output[["main_png_filename"]]) %>% print

