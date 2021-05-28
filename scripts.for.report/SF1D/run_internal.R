#' ---
#' title: "SF1D summary"
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


RNA.merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt <- fread(snakemake@input[["RNA_merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_txt_gz_filename"]])

if (FALSE) {
    RNA.merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt <- fread("result/S71_5__filter_for_A_to_G_sites_for_control/210203-GSE144296.A375-RNA-with-DNA-37-37/210203-GSE144296.A375-RNA-with-DNA-37-37/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz")
}


subset.name <- snakemake@wildcards[["SUBSET_NAME"]]
if (FALSE){
    subset.name <- "all.samples"
}

if (subset.name != "all.samples"){
    stop("subset.name must be 'all.samples'")
}


DNA.dataset.collection.name <- snakemake@wildcards[["DNA_DATASET_COLLECTION_NAME"]]
if (FALSE) {
    DNA.dataset.collection.name <- "210203-GSE144296.A375-DNA"
}

if (DNA.dataset.collection.name == "210203-GSE144296.A375-DNA"){

    GSE144296.paired.sequenced.info.dt <- fread("external/NCBI.SRA.MetaData/GSE144296.txt")[Cell_Line=="A375"][, cell_ID_occurrence:=.N, list(cell_ID)][cell_ID_occurrence==2] %>%
        dcast(formula=cell_ID ~ LibrarySource, value.var="Sample Name") %>%
        {.[TRANSCRIPTOMIC %in% RNA.merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt[, SAMPLE]]}
    
    all.DNA.alignments.dt <- foreach(temp.row.dt=iter(GSE144296.paired.sequenced.info.dt, by="row")) %do% {
        cat(date(), "reading DNA alignment for sample ", temp.row.dt[1, GENOMIC], "\n")
        temp.DNA.alignment.bcf.filename <- paste(sep="", "result/S15_1__get_sample_RNA_editing_sites_v3/paired-37-37/210203-GSE144296.A375-DNA-37-37/", temp.row.dt[1, GENOMIC], "/__merged__/DNTRSeq-DNA-trimming/hg38.fa/32/bwa-index-default/DNA/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.bcf")
        temp.DNA.alignment.dt <- fread(cmd=paste(sep="", "bcftools view ", temp.DNA.alignment.bcf.filename, " | grep -v '^##'"), select=1:7) %>% setnames("#CHROM", "CHROM")
        data.table(temp.DNA.alignment.dt, temp.row.dt)
    } %>% rbindlist

    RNA.DNA.comparison.dt <- merge(x=RNA.merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt, y=all.DNA.alignments.dt[, list(CHROM, POS, SAMPLE=TRANSCRIPTOMIC, is.DNA.detected=TRUE)], by=c("CHROM", "POS", "SAMPLE"), all.x=TRUE, all.y=FALSE) %>% {.[is.na(is.DNA.detected)==TRUE, is.DNA.detected:=FALSE]}

    fwrite(RNA.DNA.comparison.dt, snakemake@output[["RNA_DNA_comparison_dt_txt_gz_filename"]])

    all.raw.RNA.alignments.dt <- foreach(temp.row.dt=iter(GSE144296.paired.sequenced.info.dt, by="row")) %do% {
        cat(date(), "reading DNA alignment for sample ", temp.row.dt[1, GENOMIC], "\n")
        temp.RNA.alignment.bcf.filename <- paste(sep="", "result/S15_1__get_sample_RNA_editing_sites_v3/paired-37-37/210203-GSE144296.A375-RNA-37-37/", temp.row.dt[1, TRANSCRIPTOMIC], "/__merged__/DNTRSeq-RNA-trimming/hg38.fa/32/bwa-index-10.1038_nmeth.2330/32/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.bcf")
        temp.RNA.alignment.dt <- fread(cmd=paste(sep="", "bcftools view ", temp.RNA.alignment.bcf.filename, " | grep -v '^##'"), select=1:7) %>% setnames("#CHROM", "CHROM")
        data.table(temp.RNA.alignment.dt, temp.row.dt)
    } %>% rbindlist

    raw.RNA.DNA.comparison.dt <- merge(x=all.raw.RNA.alignments.dt[, list(CHROM, POS, SAMPLE=TRANSCRIPTOMIC)], y=all.DNA.alignments.dt[, list(CHROM, POS, SAMPLE=TRANSCRIPTOMIC, is.DNA.detected=TRUE)], by=c("CHROM", "POS", "SAMPLE"), all.x=TRUE, all.y=FALSE) %>% {.[is.na(is.DNA.detected)==TRUE, is.DNA.detected:=FALSE]}

    fwrite(raw.RNA.DNA.comparison.dt, snakemake@output[["raw_RNA_DNA_comparison_dt_txt_gz_filename"]])

    SF1D.ggplot <- list(
        data.table(raw.RNA.DNA.comparison.dt[, list(CHROM, POS, SAMPLE, is.DNA.detected)], RNA.edit.type="raw"),
        data.table(RNA.DNA.comparison.dt[, list(CHROM, POS, SAMPLE, is.DNA.detected)], RNA.edit.type="filtered")
    ) %>% rbindlist %>%
    {.[, list(count=.N), list(SAMPLE, is.DNA.detected, RNA.edit.type)]} %>%
    setkey(SAMPLE, is.DNA.detected, RNA.edit.type) %>%
    {.[CJ(SAMPLE, is.DNA.detected, RNA.edit.type, unique=TRUE)]} %>%
    {.[is.na(count)==TRUE, count:=0]} %>%
    {.[, count.log2.corrected:=log2(count+0.5)]} %>%
    {.[, is.DNA.detected.description:=c("might be real edits as not detected by DNA sequencing", "possible false positives as detected by DNA sequencing")[is.DNA.detected+1]]} %>%
    {
        ggbarplot(data=., x="SAMPLE", y="count.log2.corrected", fill="is.DNA.detected.description", position=position_dodge()) +
            coord_flip() +
            labs(x="SAMPLE", y="log2(# RNA edits detected + 0.5)\n(should be -1 if no edits were detected)", fill="") +
            facet_grid(~RNA.edit.type, scales="free_x") +
            guides(fill=guide_legend(nrow=2))
    }
    SF1D.ggplot
    
} else {
    stop(paste0("DNA.dataset.collection.name ", DNA.dataset.collection.name, " is not supported"))
}



saveRDS(SF1D.ggplot, snakemake@output[["main_ggplot_RDS_filename"]])


ggsave.A4(filename=snakemake@output[["main_png_filename"]], plot=SF1D.ggplot, width.r=0.9, height.r=0.8)
if (FALSE){
    ggsave.A4(filename="report/210203-GSE144296.A375-RNA-with-DNA-37-37/210203-GSE144296.A375-RNA-with-DNA-37-37/all.samples/SF1D.png", plot=SF1D.ggplot, width.r=0.9, height.r=0.8)
}


print(SF1D.ggplot)

image_read(snakemake@output[["main_png_filename"]]) %>% print




