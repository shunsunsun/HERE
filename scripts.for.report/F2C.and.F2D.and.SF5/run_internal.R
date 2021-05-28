#' ---
#' title: "F2C and F2D and SF5 summary"
#' ---

library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("foreach")
library("magick")
library("scales")
source("./scripts/common/ggpubr.A4.R")


subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt <- fread(snakemake@input[["subset_recurrent_edits_only_with_snpEff_annotation_on_valid_genes_only_dt_txt_gz_filename"]])

if (FALSE){
    subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt <- fread("result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")
}

final.count.dt <- fread(snakemake@input[["final_count_csv_filename"]])[, list(dataset, name, count.5.prime.UTR, count.3.prime.UTR, count.CDS, count.intron)]

if (FALSE){
    final.count.dt <- fread("result/A02_7__combine_summaries_of_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/210215-sixth-dataset/201221-fifth-phenotype-collection/final.count.csv")[, list(dataset, name, count.5.prime.UTR, count.3.prime.UTR, count.CDS, count.intron)]
}



for.plot.dt <- subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt %>%
    {.[is.na(AC)==FALSE]} %>%
    {.[, list(SAMPLE, gse, stage, Annotation.pasted, Annotation.class, CHROM, POS)]} %>% unique %>%
    ## annotate Annotation.subset
    {.[, Annotation.subset:="intron"]} %>%
    {.[grepl("splice", Annotation.pasted) == TRUE, Annotation.subset:="splice"]} %>%
    {.[grepl("3_prime_UTR", Annotation.pasted) == TRUE, Annotation.subset:="3_prime_UTR"]} %>%
    {.[grepl("5_prime_UTR", Annotation.pasted) == TRUE, Annotation.subset:="5_prime_UTR"]} %>%
    ## 5'UTR premature start codon is considered in the 5'UTR subset
    {.[grepl("(missense|synonymous|stop)", Annotation.pasted) == TRUE, Annotation.subset:="CDS"]} %>%
    {.[,
       list(
           .SD[, list(count=.N), list(metric=Annotation.class)],
           .SD[Annotation.subset != 'splice', list(count=.N), list(metric=Annotation.subset)]
       ) %>% rbindlist,
       list(SAMPLE, gse, stage)]} %>%
    {dcast(., gse+stage+SAMPLE~metric, value.var="count", fill=0)} %>%
    {.[, ES.to.I.ratio:=`exonic.or.splicing.related` / (`purely.intronic` + 1)][, log2.ES.to.I.ratio:=pmax(-3, log2(ES.to.I.ratio))]} %>%
    {.[, `3.prime.UTR.to.5.prime.UTR.or.CDS.ratio`:=`3_prime_UTR` / (`5_prime_UTR` + `CDS` + 1)][, log2.3.prime.UTR.to.5.prime.UTR.or.CDS.ratio:=pmax(-3, log2(`3.prime.UTR.to.5.prime.UTR.or.CDS.ratio`))]} %>%
    {melt(data=., id.vars=c("SAMPLE", "gse", "stage"), measure.vars=c("log2.ES.to.I.ratio", "log2.3.prime.UTR.to.5.prime.UTR.or.CDS.ratio"), variable.name="ratio.type", value.name="ratio.value")} %>%
    {.[, ratio.type.description:=c("log2.ES.to.I.ratio"="exonic v.s.\nintronic", "log2.3.prime.UTR.to.5.prime.UTR.or.CDS.ratio"="3'-UTR v.s.\n 5'-UTR or CDS")[ratio.type]]} %>%
    {.[stage %in% c('oocyte.GV', 'oocyte.MI', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell', 'morula', 'blastocyst.early', 'blastocyst.middle', 'blastocyst.late', 'trophoblast', 'TE', 'ICM', 'CTB', 'STB', 'MTB', 'epiblast', 'hypoblast')][stage %in% c('oocyte.GV', 'oocyte.MI', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell'), stage.interval:="till 8-cell"][is.na(stage.interval)==TRUE, stage.interval:="morula and afterwards"][, stage.interval.to.plot:=factor(stage.interval, levels=c("till 8-cell", "morula and afterwards"))]} %>%
    {.[, gse.combined:=gse][gse %in% c("GSE133854", "GSE36552", "GSE136447", "GSE130289", "GSE73211", "GSE119324", "GSE126488", "GSE62772", "GSE64417", "GSE95477", "GSE125616", "GSE65481")== FALSE, gse.combined:="others"][, gse.combined.sample.size.per.stage:=.N, list(stage, gse.combined, ratio.type, ratio.type.description)][, gse.combined.with.sample.size.per.stage:=paste(sep="", gse.combined, "/", gse.combined.sample.size.per.stage)]} %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])]}


F2C.and.F2D.upper.panel.ggplot <- for.plot.dt %>%
    {ggboxplot(data=., x="stage.description.ordered", y="ratio.value", fill="#00BFC4") + theme_pubr(base_size=10) +
         geom_hline(data.table(y.intercept=c(0, 2.5), ratio.type="log2.ES.to.I.ratio"), mapping=aes(yintercept=y.intercept), color="red", linetype="dashed") +
         geom_hline(data.table(y.intercept=c(2.5, 9), ratio.type="log2.3.prime.UTR.to.5.prime.UTR.or.CDS.ratio"), mapping=aes(yintercept=y.intercept), color="red", linetype="dashed") +
         theme(axis.text.x=element_text(angle=45, hjust=1)) +
         facet_grid(c("log2.ES.to.I.ratio"="exonic v.s.\nintronic", "log2.3.prime.UTR.to.5.prime.UTR.or.CDS.ratio"="3'-UTR v.s.\n 5'-UTR, CDS")[ratio.type]~stage.interval.to.plot, scales="free", space="free_x") +
         labs(
             x="",
             ## y="log2(# dominant edits \n/ (1 + # other edits))\n(lower-clipped by -3)"
             y="enrichment score"
         ) +
         theme(panel.border=element_rect(color="black", fill=NA))
    }




sample.stats.dt <- subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt %>%
    ## remove all unedited sites
    {.[is.na(AC)==FALSE]} %>%
    ## consider Alu sites only
    {.[SUBSET == "Alu"]} %>%
    ## collapse records to site x sample x Annotation.pasted
    {.[, list(CHROM, POS, Annotation.pasted, SAMPLE, group, stage)]} %>%
    unique %>%
    ## annotate Annotation.subset
    {.[, Annotation.subset:="intron"]} %>%
    {.[grepl("splice", Annotation.pasted) == TRUE, Annotation.subset:="splice"]} %>%
    {.[grepl("3_prime_UTR", Annotation.pasted) == TRUE, Annotation.subset:="3_prime_UTR"]} %>%
    {.[grepl("5_prime_UTR", Annotation.pasted) == TRUE, Annotation.subset:="5_prime_UTR"]} %>%
    ## 5'UTR premature start codon is considered in the 5'UTR subset
    {.[grepl("(missense|synonymous|stop)", Annotation.pasted) == TRUE, Annotation.subset:="CDS"]} %>%
    ## remove all sites with splicing annotations
    {.[Annotation.subset!="splice"]} %>%
    {.[,
       data.table(
           RE.count.exon=sum(Annotation.subset %in% c("3_prime_UTR", "5_prime_UTR", "CDS")),
           RE.count.intron=sum(Annotation.subset == "intron"),
           RE.count.3.prime.UTR=sum(Annotation.subset=="3_prime_UTR"),
           RE.count.5.prime.UTR=sum(Annotation.subset=="5_prime_UTR"),
           RE.count.CDS=sum(Annotation.subset=="CDS"),
           RE.count.5.prime.UTR.or.CDS=sum(Annotation.subset %in% c("5_prime_UTR", "CDS")),
           RE.count.total=.N
       ),
       list(SAMPLE, group, stage)]} %>%
    {merge(x=., y=final.count.dt[, list(name, count.5.prime.UTR, count.3.prime.UTR, count.CDS, count.intron, count.exon=count.5.prime.UTR + count.3.prime.UTR + count.CDS, count.5.prime.UTR.or.CDS=count.5.prime.UTR + count.CDS)], by.x="SAMPLE", by.y="name", all.x=TRUE, all.y=FALSE)} %>%
    {.[, phyper.exon.enriched.in.RE:=phyper(q=RE.count.exon, m=count.exon, n=count.intron, k=RE.count.total, lower.tail=FALSE)]} %>%
    {.[, phyper.3.prime.UTR.enriched.in.RE:=phyper(q=RE.count.3.prime.UTR, m=count.3.prime.UTR, n=count.5.prime.UTR.or.CDS, k=RE.count.exon, lower.tail=FALSE)]} %>%
    ## remove invalid cases
    {.[RE.count.exon==0, phyper.exon.enriched.in.RE:=NA]} %>%
    {.[RE.count.3.prime.UTR==0, phyper.3.prime.UTR.enriched.in.RE:=NA]}


fwrite(sample.stats.dt, snakemake@output[["sample_stats_dt_txt_gz_filename"]])

F2C.and.F2D.lower.panel.ggplot <- sample.stats.dt %>%
    {melt(data=., id.vars=c("SAMPLE", "group", "stage"), measure.vars=c("phyper.exon.enriched.in.RE", "phyper.3.prime.UTR.enriched.in.RE"), variable.name="pvalue.type", value.name="pvalue.value")} %>%
    ## first select stages and then adjust the p-values
    {.[stage %in% c('oocyte.GV', 'oocyte.MI', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell', 'morula', 'blastocyst.early', 'blastocyst.middle', 'blastocyst.late', 'trophoblast', 'TE', 'ICM', 'CTB', 'STB', 'MTB', 'epiblast', 'hypoblast')][stage %in% c('oocyte.GV', 'oocyte.MI', 'oocyte.MII', 'zygote', '2-cell', '4-cell', '8-cell'), stage.interval:="till 8-cell"][is.na(stage.interval)==TRUE, stage.interval:="morula and afterwards"][, stage.interval.to.plot:=factor(stage.interval, levels=c("till 8-cell", "morula and afterwards"))]} %>%
    {.[, padj:=p.adjust(p=pvalue.value, method="BH")]} %>%
    {.[, padj.log10.negative:= -1 * log10(padj)]} %>%
    {.[, padj.log10.negative.clipped:=padj.log10.negative][padj < 1e-5, padj.log10.negative.clipped:=-1 * log10(1e-5)]} %>%
    {.[, pvalue.type.description:=c("phyper.exon.enriched.in.RE"="exonic\nenrichment", "phyper.3.prime.UTR.enriched.in.RE"="3'-UTR\nenrichment")[pvalue.type]]} %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])][, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))]} %>%
    {
        ggboxplot(data=., x="stage.description.ordered", y="padj.log10.negative.clipped", fill="#00BFC4") +
            geom_hline(aes(yintercept=-1 * log10(0.01)), color="#2E3436", linetype="dashed") +
            labs(x="", y="-log10(p-value)\n(upper-clipped by 5)", fill="") +
            theme_pubr(base_size=10) +
            theme(axis.text.x=element_text(angle=45, hjust=1)) +
            facet_grid(pvalue.type.description~stage.interval.to.plot, scales="free", space="free") +
            theme(panel.border=element_rect(color="black", fill=NA))
    }

F2C.and.F2D.ggplot <- ggarrange(plotlist=list(F2C.and.F2D.upper.panel.ggplot, F2C.and.F2D.lower.panel.ggplot), ncol=1, align="v")

saveRDS(F2C.and.F2D.ggplot, snakemake@output[["main_ggplot_RDS_filename"]])


ggsave.A4(filename=snakemake@output[["main_png_filename"]], plot=F2C.and.F2D.ggplot, width.r=0.9, height.r=0.50)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F2C.and.F2D.png", plot=F2C.and.F2D.ggplot, width.r=0.9, height.r=0.50)
}


print(F2C.and.F2D.ggplot)

image_read(snakemake@output[["main_png_filename"]]) %>% print





SF5.ggplot <- foreach(temp.stage.interval.to.plot=c("till 8-cell", "morula and afterwards"), temp.filename.label=c("till.8.cell", "morula.and.afterwards")) %do% {
    for.plot.dt %>%
        {.[stage.interval.to.plot==temp.stage.interval.to.plot]} %>%
        {ggboxplot(data=., x="gse.combined.with.sample.size.per.stage", y="ratio.value", fill="#00BFC4") + theme_pubr(base_size=10) + geom_hline(yintercept=2.5, color="blue", linetype="dashed") + theme_pubr(base_size=10)  + labs(x="dataset / sample size", y="log2(# E or S edits / # I edits)\nacross all genes per sample") + lims(y=c(-4.0, 10.0)) + ggtitle(label=temp.stage.interval.to.plot) + 
         geom_hline(data.table(y.intercept=c(0, 2.5), ratio.type="log2.ES.to.I.ratio"), mapping=aes(yintercept=y.intercept), color="red", linetype="dashed") +
         geom_hline(data.table(y.intercept=c(2.5, 9), ratio.type="log2.3.prime.UTR.to.5.prime.UTR.or.CDS.ratio"), mapping=aes(yintercept=y.intercept), color="red", linetype="dashed") +
         theme(axis.text.x=element_text(angle=45, hjust=1)) +
             facet_grid(c("log2.ES.to.I.ratio"="exonic v.s.\nintronic", "log2.3.prime.UTR.to.5.prime.UTR.or.CDS.ratio"="3'-UTR v.s.\n 5'-UTR, CDS")[ratio.type]~stage.description.ordered, scales="free", space="free_x") + labs(x="", y="enrichment score") +
             theme(panel.border=element_rect(color="black", fill=NA))}
} %>%
    {ggarrange(plotlist=., nrow=2)}
SF5.ggplot


saveRDS(SF5.ggplot, snakemake@output[["supp1_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(SF5.ggplot, "report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF5.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["supp1_png_filename"]], plot=SF5.ggplot, width.r=0.9, height.r=0.8)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF5.png", plot=SF5.ggplot, width.r=1.2, height.r=0.8)
}


print(SF5.ggplot)

image_read(snakemake@output[["supp1_png_filename"]]) %>% print
