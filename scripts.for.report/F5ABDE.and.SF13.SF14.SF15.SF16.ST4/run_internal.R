#' ---
#' title: "F5ABDE and SF13 SF14 SF15 SF16 ST4 summary"
#' ---

library("data.table")
library("magrittr")
library("ggpubr")
library("readxl")
library("foreach")
library("magick")
library("clusterProfiler")
library("scales")
library("iterators")
library("org.Hs.eg.db")
source("./scripts/common/ggpubr.A4.R")


ALS.subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt <- fread(snakemake@input[["ALS_subset_recurrent_edits_only_with_snpEff_annotation_on_valid_genes_only_dt_txt_gz_filename"]])


if (FALSE){
    ALS.subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt <- fread("result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")
}


subset.site.recurrence.comparison.CJ.dt <- fread(snakemake@input[["subset_site_recurrence_comparison_CJ_dt_txt_gz_filename"]])

if (FALSE){
    subset.site.recurrence.comparison.CJ.dt <- fread("result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/subset.site.recurrence.comparison.CJ.dt.txt.gz")
}

subset.dt <- fread(snakemake@input[["subset_dt_txt_gz_filename"]], select=c("CHROM", "POS", "disease", "SAMPLE", "AF", "stage"))


if (FALSE){
    subset.dt <- fread("result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/subset.dt.txt.gz", select=c("CHROM", "POS", "disease", "SAMPLE", "AF", "stage"))
}


subset.site.recurrence.comparison.CJ.dt[, `:=`(
    stage=sub(pattern="(.*)@.*", replacement="\\1", x=group),
    disease=sub(pattern=".*@(.*)", replacement="\\1", x=group)
)]

subset.site.recurrence.comparison.CJ.dcast.dt <- dcast(subset.site.recurrence.comparison.CJ.dt, CHROM + POS + stage ~ disease, value.var="site.occurrence.for.this.group", fill=-1)



phenotype.dt <- fread(snakemake@input[["phenotype_output_at_gsm_level_dt_filename"]])
if (FALSE){
    phenotype.dt <- fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt")
}
phenotype.subset.dt <- phenotype.dt[gse=="GSE133854"] %>% setnames("gsm", "SAMPLE")


## F5A

annotation.mapping.vector <- c(
    "3_prime_UTR_variant"="3'-UTR variant",
    "intron_variant"="intron variant",
    "5_prime_UTR_variant"="5'-UTR variant",
    "missense_variant"="missense\nvariant",
    "synonymous_variant"="synonymous variant",
    "splice_region_variant"="splice region variant",
    "5_prime_UTR_premature_start_codon_gain_variant"="5'-UTR premature\nstart codon gain variant",
    "stop_lost"="stop\nlost",
    "splice_acceptor_variant"="splice acceptor variant")
variant.of.interest.dt <- merge(x=subset.site.recurrence.comparison.CJ.dcast.dt[, list(CHROM, POS, stage, androgenetic, biparental, parthenogenetic)], y=ALS.subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt[stage %in% c("zygote", "2-cell", "4-cell", "8-cell", "morula")][Annotation.corrected %in% c("missense_variant", "stop_lost", "5_prime_UTR_premature_start_codon_gain_variant", "missense_variant&splice_region_variant")][, list(CHROM, POS, stage, Gene_Name, Annotation.corrected)] %>% unique, by=c("CHROM", "POS", "stage"), all=FALSE)[
    (stage %in% c("zygote") & ((biparental >= 7*0.5 & (androgenetic <=1 | parthenogenetic <=1 )))), pattern:="lost in zygote"
    ][
        (stage %in% c("2-cell") & ((biparental >= 10*0.5 & (androgenetic <=2 | parthenogenetic <=2 )))), pattern:="lost in 2-cell"
    ][
        (stage %in% c("4-cell") & ((biparental >= 21*0.5 & (androgenetic <=2 | parthenogenetic <=2 )))), pattern:="lost in 4-cell"
    ][
        (stage %in% c("8-cell") & ((biparental >= 34*0.5 & (androgenetic <=2 | parthenogenetic <=2 )))), pattern:="lost in 8-cell"
    ][
        (stage %in% c("morula") & ((biparental >= 48*0.5 & (androgenetic <=2 | parthenogenetic <=2 )))), pattern:="lost in morula"
] %>% 
    {
        .[
            is.na(pattern)==FALSE
        ][
          , Annotation.corrected.ordered := factor(annotation.mapping.vector[Annotation.corrected])
        ]
    }
F5A.for.plot.dt  <- merge(x=variant.of.interest.dt[, list(CHROM, POS, Gene_Name, Annotation.corrected.ordered, pattern, stage)] %>% unique, y=subset.dt[, list(CHROM, POS, disease, SAMPLE, AF, stage)], by=c("CHROM", "POS", "stage"), all=FALSE)[disease=="", disease:="biparental"][, disease.ordered:=factor(c("biparental"="BI", "androgenetic"="AG", "parthenogenetic"="PG")[disease], levels=c("BI", "AG", "PG"))][, site.name:=paste(sep="", CHROM, ":", POS, "_", Gene_Name)]
F5A.plot.list <- foreach(temp.pattern=paste(sep="", "lost in ", c("zygote", "2-cell", "4-cell")), temp.stage=c("zygote", "2-cell", "4-cell")) %do% {
    temp.F5A.for.plot.subset.dt <- F5A.for.plot.dt[pattern==temp.pattern]
    temp.F5A.subset.ggplot <- ggplot(temp.F5A.for.plot.subset.dt, aes(x=SAMPLE, y=site.name, fill=AF)) + geom_tile() + facet_grid(Annotation.corrected.ordered~disease.ordered, scales="free", space="free") + theme_pubr(base_size=10) + theme(panel.background=element_rect(fill="grey70"), legend.position="top", axis.text.x=element_blank()) + scale_fill_gradientn(limits=c(0, 1), colours=c("white", "lightpink", "darkred"), values=c(0, 0.1,  1),  na.value="white") + labs(x=paste(sep="", "GSE133854 ", temp.stage, " samples"), y=temp.pattern, fill="editing level")
    if (temp.pattern != "lost in zygote") {
        temp.F5A.subset.ggplot <- temp.F5A.subset.ggplot + theme(legend.position="none")
    }
    temp.F5A.subset.ggplot
}
F5A.ggplot <- ggarrange(plotlist=F5A.plot.list, ncol=1, align="v", heights=c(6, 5, 5))
ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/F5A.png", plot=F5A.ggplot, width.r=0.5, height.r=0.4)


## F5B

setkey(ALS.subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt, CHROM, POS)
F5B.for.plot.dt <- merge(
    x=ALS.subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt[F5A.for.plot.dt[, list(CHROM, POS)]][is.na(AF)==FALSE][, list(CHROM, POS, Gene_Name, Annotation.corrected, SAMPLE, gse, stage)] %>% unique %>% {.[, list(sample.count=.N), list(CHROM, POS, Gene_Name, Annotation.corrected, gse, stage)]},
    y=F5A.for.plot.dt[, list(CHROM, POS, Gene_Name, stage)] %>% unique,
    by=c("CHROM", "POS", "Gene_Name", "stage"),
    all=FALSE) %>%
    {
        merge(
            x=.,
            y=phenotype.dt[is.normal==TRUE][, list(total.sample.count=.N), list(gse, stage)],
            by=c("gse", "stage"),
            all.x=TRUE, all.y=FALSE
            )
    } %>% {
        .[, sample.pct:=sample.count / total.sample.count * 100]
        .[, site.name:=paste(sep="", CHROM, ":", POS, "_", Gene_Name)]
    } %>% {
        setkey(., site.name, stage, gse)
        temp.dt <- merge(x=.[CJ(site.name, stage, gse, unique=TRUE)], y=.[, list(site.name, stage)] %>% unique, by=c("site.name", "stage"), all=FALSE)
        setnafill(temp.dt, type="const", fill=0, cols="sample.pct")
        ## temp.dt[, site.name:=sub(pattern="(.*)@.*", replacement="\\1", x=site.name.and.stage)]
        ## temp.dt[, stage:=sub(pattern=".*@(.*)", replacement="\\1", x=site.name.and.stage)]
    } %>% {
        .[, Annotation.corrected:=Annotation.corrected %>% na.omit %>% unique, list(site.name)]
        .[, Annotation.corrected.ordered := factor(annotation.mapping.vector[Annotation.corrected])]
    }
F5B.ggplot <- ggplot(F5B.for.plot.dt, aes(x=site.name, y=sample.pct, fill=gse)) + geom_bar(stat='identity', position="dodge") + coord_flip() + facet_grid(factor(stage, levels=c("zygote", "2-cell", "4-cell")) + Annotation.corrected.ordered~., scales="free_y", space="free") + theme_pubr(base_size=10) + theme(plot.margin=margin(2, 0, 0, 0, "cm"), legend.position=c(0.2, 1.15)) + labs(x="", y="% normal samples detected", fill="") + guides(fill = guide_legend(nrow=2))
ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/F5B.png", plot=F5B.ggplot, width.r=0.46, height.r=0.4)


## F5D

F5D.for.plot.dt  <- subset.dt[CHROM=="chr19" & POS==4531743, list(CHROM, POS, disease, SAMPLE, AF, stage)] %>%
    {.[stage != 'oocyte.MII' & disease=="", disease:="biparental"][stage == 'oocyte.MII', `:=`(disease="NA", stage="oocyte (MII)")][, disease.ordered:=factor(c("NA"="NA", "biparental"="BI", "androgenetic"="AG", "parthenogenetic"="PG")[disease], levels=c("", "BI", "AG", "PG"))]}
F5D.plot.list <- foreach(temp.stage=c("oocyte (MII)", "zygote", "2-cell", "4-cell", "8-cell", "morula"), temp.sample.class=c("oocytes (MII)", "zygotes", "2-cells", "4-cells", "8-cells", "morulae")) %do% {
    temp.F5D.for.plot.subset.dt <- F5D.for.plot.dt[stage==temp.stage]
    temp.F5D.subset.ggplot <- ggplot(temp.F5D.for.plot.subset.dt, aes(x=SAMPLE, y="1", fill=AF)) + geom_tile() + facet_grid(~disease.ordered, scales="free", space="free") + theme_pubr(base_size=10) + theme(panel.background=element_rect(fill="grey70"), legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_fill_gradientn(limits=c(0, 1), colours=c("white", "lightpink", "darkred"), values=c(0, 0.1,  1),  na.value="white") + labs(x=paste(sep="", "GSE133854 ", temp.sample.class), y="", fill="editing level at\nchr19:4531743_\nPLIN5\n(sequenced\nsamples\nonly)")
    if (FALSE) {
        ##temp.F5D.subset.ggplot <- temp.F5D.subset.ggplot
    }
    temp.F5D.subset.ggplot
}
F5D.ggplot <- ggarrange(plotlist=F5D.plot.list, ncol=1, heights=c(2, 2, 2, 2, 2, 2), common.legend=TRUE, legend="left")
ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/F5D.png", plot=F5D.ggplot, width.r=0.4, height.r=0.4)



## F5E
setkey(ALS.subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt, CHROM, POS)
F5E.for.plot.dt <- ALS.subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt[data.table(CHROM="chr19", POS=4531743)][is.na(AF)==FALSE][, list(CHROM, POS, Gene_Name, Annotation.corrected, SAMPLE, gse, stage)] %>% unique %>% {.[, list(sample.count=.N), list(CHROM, POS, gse, stage)]} %>%
    {
        merge(
            x=.,
            y=phenotype.dt[is.normal==TRUE][, list(total.sample.count=.N), list(gse, stage)],
            by=c("gse", "stage"),
            all.x=TRUE, all.y=FALSE
            )
    } %>% {
        .[, sample.pct:=sample.count / total.sample.count * 100]
    } %>% {
        setkey(., stage, gse)
        temp.dt <- merge(x=.[CJ(stage, gse, unique=TRUE)], y=.[, list(stage)] %>% unique, by=c("stage"), all=FALSE)
        setnafill(temp.dt, type="const", fill=0, cols="sample.pct")
    } %>% {
        .[stage=="oocyte.MII", stage:="oocyte (MII)"]
    }
F5E.ggplot <- ggplot(F5E.for.plot.dt, aes(x="1", y=sample.pct, fill=gse)) + geom_bar(stat='identity', position="dodge") + coord_flip() + facet_grid(factor(stage, levels=c("oocyte (MII)", "zygote", "2-cell", "4-cell", "8-cell", "morula")) ~., scales="free_y", space="free") + theme_pubr(base_size=10) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + labs(x="", y="% normal samples detected", fill="") + guides(fill = guide_legend(nrow=2)) + ggtitle(label="chr19:4531743_PLIN5", subtitle="(displaying only those stages where this edit is recurrent)")
ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/F5E.png", plot=F5E.ggplot, width.r=0.46, height.r=0.4)


## SF13, SF14, SF15, SF16

variant.of.interest.full.dt <- merge(x=subset.site.recurrence.comparison.CJ.dcast.dt[, list(CHROM, POS, stage, androgenetic, biparental, parthenogenetic)], y=ALS.subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt[stage %in% c("zygote", "2-cell", "4-cell", "8-cell", "morula")][, list(CHROM, POS, stage, Gene_Name, Annotation.corrected)] %>% unique, by=c("CHROM", "POS", "stage"), all=FALSE)[
    (stage %in% c("zygote") & ((biparental >= 7*0.5 & (androgenetic <=1 | parthenogenetic <=1 )))), pattern:="lost in zygote"
    ][
        (stage %in% c("2-cell") & ((biparental >= 10*0.5 & (androgenetic <=2 | parthenogenetic <=2 )))), pattern:="lost in 2-cell"
    ][
        (stage %in% c("4-cell") & ((biparental >= 21*0.5 & (androgenetic <=2 | parthenogenetic <=2 )))), pattern:="lost in 4-cell"
    ][
        (stage %in% c("8-cell") & ((biparental >= 34*0.5 & (androgenetic <=2 | parthenogenetic <=2 )))), pattern:="lost in 8-cell"
    ][
        (stage %in% c("morula") & ((biparental >= 48*0.5 & (androgenetic <=2 | parthenogenetic <=2 )))), pattern:="lost in morula"
] %>% 
    {
        .[
            is.na(pattern)==FALSE
        ][
          , Annotation.corrected.ordered := factor(annotation.mapping.vector[Annotation.corrected])
        ]
    }
SF13.SF14.SF15.SF16.for.plot.dt  <- merge(x=variant.of.interest.full.dt[, list(CHROM, POS, Gene_Name, Annotation.corrected.ordered, pattern, stage)] %>% unique, y=subset.dt[, list(CHROM, POS, disease, SAMPLE, AF, stage)], by=c("CHROM", "POS", "stage"), all=FALSE)[disease=="", disease:="biparental"][, disease.ordered:=factor(c("biparental"="BI", "androgenetic"="AG", "parthenogenetic"="PG")[disease], levels=c("BI", "AG", "PG"))][, site.name:=paste(sep="", CHROM, ":", POS, "_", Gene_Name)]
SF13.SF14.SF15.SF16.plot.list <- foreach(temp.pattern=paste(sep="", "lost in ", c("zygote", "2-cell", "4-cell", "8-cell", "morula")), temp.stage=c("zygote", "2-cell", "4-cell", "8-cell", "morula")) %do% {
    temp.SF13.SF14.SF15.SF16.for.plot.subset.dt <- SF13.SF14.SF15.SF16.for.plot.dt[pattern==temp.pattern]
    temp.SF13.SF14.SF15.SF16.subset.ggplot <- ggplot(temp.SF13.SF14.SF15.SF16.for.plot.subset.dt, aes(x=SAMPLE, y=site.name, fill=AF)) + geom_tile() + facet_grid(Annotation.corrected.ordered~disease.ordered, scales="free", space="free") + theme_pubr(base_size=10) + theme(panel.background=element_rect(fill="grey70"), legend.position="top", axis.text.x=element_blank()) + scale_fill_gradientn(limits=c(0, 1), colours=c("white", "lightpink", "darkred"), values=c(0, 0.1,  1),  na.value="white") + labs(x=paste(sep="", "GSE133854 ", temp.stage, " samples"), y=temp.pattern, fill="editing level")
    if (temp.pattern != "lost in zygote") {
        temp.SF13.SF14.SF15.SF16.subset.ggplot <- temp.SF13.SF14.SF15.SF16.subset.ggplot + theme(legend.position="none")
    }
    temp.SF13.SF14.SF15.SF16.subset.ggplot
}
ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/SF13.png", plot=SF13.SF14.SF15.SF16.plot.list[[1]], width.r=1, height.r=3)
ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/SF14.png", plot=SF13.SF14.SF15.SF16.plot.list[[2]], width.r=1, height.r=3)
ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/SF15.png", plot=SF13.SF14.SF15.SF16.plot.list[[3]], width.r=1, height.r=2.5)
ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/SF16.png", plot=ggarrange(plotlist=SF13.SF14.SF15.SF16.plot.list[4:5], ncol=1, align="v", heights=c(5, 2)), width.r=1, height.r=0.4)




## ST4

annotation.mapping.vector <- c(
    "3_prime_UTR_variant"="3'-UTR variant",
    "intron_variant"="intron variant",
    "5_prime_UTR_variant"="5'-UTR variant",
    "missense_variant"="missense\nvariant",
    "synonymous_variant"="synonymous variant",
    "splice_region_variant"="splice region variant",
    "5_prime_UTR_premature_start_codon_gain_variant"="5'-UTR premature\nstart codon gain variant",
    "stop_lost"="stop\nlost",
    "splice_acceptor_variant"="splice acceptor variant")
variant.of.interest.full.dt <- merge(x=subset.site.recurrence.comparison.CJ.dcast.dt[, list(CHROM, POS, stage, androgenetic, biparental, parthenogenetic)], y=ALS.subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt[stage %in% c("zygote", "2-cell", "4-cell", "8-cell", "morula")][, list(CHROM, POS, stage, Gene_Name, Annotation.corrected)] %>% unique, by=c("CHROM", "POS", "stage"), all=FALSE)[
    (stage %in% c("zygote") & ((biparental >= 7*0.5 & (androgenetic <=1 | parthenogenetic <=1 )))), pattern:="lost in zygote"
    ][
        (stage %in% c("2-cell") & ((biparental >= 10*0.5 & (androgenetic <=2 | parthenogenetic <=2 )))), pattern:="lost in 2-cell"
    ][
        (stage %in% c("4-cell") & ((biparental >= 21*0.5 & (androgenetic <=2 | parthenogenetic <=2 )))), pattern:="lost in 4-cell"
    ][
        (stage %in% c("8-cell") & ((biparental >= 34*0.5 & (androgenetic <=2 | parthenogenetic <=2 )))), pattern:="lost in 8-cell"
    ][
        (stage %in% c("morula") & ((biparental >= 48*0.5 & (androgenetic <=2 | parthenogenetic <=2 )))), pattern:="lost in morula"
] %>% 
    {
        .[
            is.na(pattern)==FALSE
        ][
          , Annotation.corrected.ordered := factor(annotation.mapping.vector[Annotation.corrected])
        ]
    }
ST4.raw.dt  <- merge(x=variant.of.interest.full.dt[, list(CHROM, POS, Gene_Name, Annotation.corrected, pattern, stage)] %>% unique, y=subset.dt[, list(CHROM, POS, disease, SAMPLE, AF, stage)], by=c("CHROM", "POS", "stage"), all=FALSE)[disease=="", disease:="biparental"][, list(CHROM, POS, Gene_Name, Annotation.corrected, pattern, SAMPLE, stage, disease, AF)][, sequenced:=TRUE]
ST4.dt <- ST4.raw.dt[, {
    temp.rest.SAMPLEs.vector <- setdiff(phenotype.subset.dt[stage == stage.cache & disease == disease.cache, SAMPLE], .SD[, SAMPLE])
    if (length(temp.rest.SAMPLEs.vector) == 0){
        .SD
    } else {
        list(.SD, data.table(SAMPLE=temp.rest.SAMPLEs.vector, AF=as.numeric(NA), sequenced=FALSE)) %>% rbindlist(use.names=TRUE)}
}, list(CHROM, POS, Gene_Name, Annotation.corrected, pattern, stage.cache=stage, disease.cache=disease)]
fwrite(ST4.dt, "report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/ST4.csv", na="NA")
