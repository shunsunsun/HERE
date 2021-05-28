#' ---
#' title: "SF2 summary"
#' ---


library("ggpubr")
library("stringr")
library("magrittr")
library("foreach")
library("data.table")
library("readxl")
library("iterators")

library("gtable")
library("grid")

source("./scripts/common/ggpubr.A4.R")




phenotype.output.at.gsm.level.dt <- fread(snakemake@input[["phenotype_output_at_gsm_level_dt_filename"]])
not.used.variable <- '
phenotype.output.at.gsm.level.dt <- fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt")
'


temp.combinations.dt <- data.table(
    type.and.dataset.and.sample=c(
        "/paired-100-100/201101-GSE95477-full20-100-100/GSM2514775",
        "/paired-100-100/201101-GSE95477-full20-100-100/GSM2514772",
        "/paired-125-125/200919-GSE71318-full48-125-125/GSM1833283",
        "/paired-100-100/200902-GSE101571-full-100-100/GSM2706261",
        "/single-100/201104-GSE36552-full124-100/GSM922156",
        "/paired-125-125/200902-GSE101571-full-125-125/GSM2706243",   
        "/paired-150-150/200924-GSE133854-all296-150-150/GSM3928473"
    ),
    index.parameter=c(
        100-5,
        100-5,
        125-5,
        100-5,
        100-5,
        125-5,
        150-5
    )) %>%
    {.[, gsm:=sub(pattern=".*/(GSM[0-9]+)$", replacement="\\1", x=type.and.dataset.and.sample)]} %>%
    {merge(x=., y=phenotype.output.at.gsm.level.dt, by="gsm", all.x=TRUE, all.y=FALSE)} %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=c("gene", "genomic DNA", "unedited RNA", "edited RNA", "unedited AA", "edited AA", temp.stage.dt[, stage.description]))][, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))]}


## remove reads with < and/or > only in this window
temp.dt <- foreach(temp.row.dt=iter(temp.combinations.dt, by="row")) %do% { system(paste(sep="", "samtools tview -d Text -p chr20:37519160 result/S15_1__get_sample_RNA_editing_sites_v3/", temp.row.dt[1, type.and.dataset.and.sample], "/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/", temp.row.dt[1, index.parameter], "/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/recal.bam.subset/BLCAP.Y2C.u10.d10/recal.subset.bam | tail -n +4 | cut -c 1-21 | grep -v -P '^[<>]+$'  "), intern=TRUE) %>% toupper %>% str_split(pattern="") %>% {do.call(rbind, .)} %>% as.data.table %>% setnames(37519160:37519180 %>% as.character) %>% {.[`37519170`!=" "][order(`37519170`)]} %>% {data.table(index=1:nrow(.), .)} %>% melt(id.vars="index") %>% {.[, position:=as.numeric(as.character(variable))][, value:=factor(value, levels=c("A", "C", "G", "T"))]} %>% {data.table(temp.row.dt, .)} } %>% rbindlist %>%
    {.[, final.facet:=paste(sep="\n", stage.description.ordered, gse, gsm)]} %>%
    {.[, final.facet.ordered:=factor(final.facet, c("genomic\nDNA", "unedited\nRNA", "edited\nRNA", "unedited\nprotein", "edited\nprotein", final.facet[order(stage.description.ordered)] %>% unique))]}

SF2.ggplot <- temp.dt %>% {
    temp.total.ggplot <- ggplot() + geom_tile(data=.[value==" "], mapping=aes(x=position, y=index), fill="white") + geom_tile(data=.[value!=" "], mapping=aes(x=position, y=index), fill="grey50") +   geom_tile(data=.[position=="37519170"], mapping=aes(x=position, y=index, fill=value)) + scale_fill_manual(values=c("pink", "orange", "blue", "brown"), drop=FALSE) + facet_grid(final.facet.ordered~., scales="free") + theme_pubr() +
        geom_segment(
            data=data.table(final.facet.ordered=factor("genomic\nDNA", levels=levels(.[, final.facet.ordered]))),
            mapping=aes(x=37519159, xend=37519181, y=1, yend=1), arrow=arrow(), color="grey80") + 
        geom_text(
            data=data.table(position=37519160:37519180, index=1, text="CTGGAGGCAATACATGATCTC" %>% str_split(pattern="") %>% unlist, final.facet.ordered=factor("genomic\nDNA", levels=levels(.[, final.facet.ordered]))),
            mapping=aes(x=position, y=index, label=text), color="black") + 
        geom_segment(
            data=data.table(final.facet.ordered=factor(c("unedited\nRNA", "edited\nRNA", "unedited\nprotein", "edited\nprotein"), levels=levels(.[, final.facet.ordered]))),
            mapping=aes(x=37519181, xend=37519159, y=1, yend=1), arrow=arrow(), color="grey80") + 
        geom_text(
            data=data.table(position=37519160:37519180, index=1, text="GACCUCCGUUAUGUACUAGAG" %>% str_split(pattern="") %>% unlist, final.facet.ordered=factor("unedited\nRNA", levels=levels(.[, final.facet.ordered]))),
            mapping=aes(x=position, y=index, label=text), color="black") + 
        geom_text(
            data=data.table(position=37519160:37519180, index=1, text="GACCUCCGUUGUGUACUAGAG" %>% str_split(pattern="") %>% unlist, final.facet.ordered=factor("edited\nRNA", levels=levels(.[, final.facet.ordered])), color=c(rep("black", 10), "red", rep("black", 10))),
            mapping=aes(x=position, y=index, label=text, color=color)) +
        geom_text(
            data=data.table(position=37519160:37519180, index=1, text=" Q  L  C  Y  M xxxxxx" %>% str_split(pattern="") %>% unlist, final.facet.ordered=factor("unedited\nprotein", levels=levels(.[, final.facet.ordered]))),
            mapping=aes(x=position, y=index, label=text), color="black") + 
        geom_text(
            data=data.table(position=37519160:37519180, index=1, text=" Q  L  C  C  M xxxxxx" %>% str_split(pattern="") %>% unlist, final.facet.ordered=factor("edited\nprotein", levels=levels(.[, final.facet.ordered])), color=c(rep("black", 10), "red", rep("black", 10))),
            mapping=aes(x=position, y=index, label=text, color=color)) +
        scale_color_manual(values=c("black", "red")) +
        theme(panel.background = element_rect(color="black", fill="white"), strip.text.y.right = element_text(angle = 0)) + guides(color=FALSE) +
        labs(x="position at chr20", y="# reads\t\t\t\t\t\t\t\t\t", fill="read base at chr20:37519170") +
        scale_x_continuous(n.breaks=3)
    ggplotGrob(temp.total.ggplot) %>% gtable_filter("axis-l-[12345]{1}$", invert=TRUE) %>% ggarrange
}

ggsave.A4(filename=snakemake@output[["main_png_filename"]], plot=SF2.ggplot, width.r=0.6, height.r=0.9)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF2.png", plot=SF2.ggplot, width.r=0.6, height.r=0.9)
}
