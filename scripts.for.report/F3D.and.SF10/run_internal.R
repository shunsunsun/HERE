#' ---
#' title: "F3D and SF10 summary"
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
source("./scripts/common/ggpubr.A4.R")



F3B.for.plot.dt <- fread(snakemake@input[["F3B_for_plot_dt_txt_gz_filename"]])

if (FALSE){
    F3B.for.plot.dt <- fread("report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F3B.for.plot.dt.txt.gz")
}

temp.GO.directory <- paste(sep="", snakemake@params[["result_directory"]], "/F3D.and.SF10.temp/GO")
if (FALSE) {
    temp.GO.directory <- "report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F3D.and.SF10.temp/GO"
}

dir.create(paste(sep="", temp.GO.directory), recursive=TRUE)

GO.result.dt <- F3B.for.plot.dt %>%
    {.[gene.occurrence.pct>=0.8]} %>%
    {.[stage %in% c('oocyte.GV', 'oocyte.MII', 'zygote', "2-cell", '4-cell', '8-cell')]} %>%
    setkey(stage, Annotation.class.pasted) %>%
    {.[, Gene_ID.no.version:=sub(pattern="\\.[0-9]*$", replacement="", x=Gene_ID)]} %>%
    {foreach(
         temp.keys=.[, list(stage, Annotation.class.pasted, Annotation.class.pasted.description)] %>%
             unique %>%
             as.data.frame(stringsAsFactors=FALSE) %>%
             merge(y=data.frame(temp.ontology=c("BP", "MF", "CC"), stringsAsFactors=FALSE)) %>%
             data.table %>%
             iter(by='row')) %do%
         {
             cat(date(), " Processing ", temp.keys[1, stage], " x ", temp.keys[1, Annotation.class.pasted], " x ", temp.keys[1, temp.ontology], "\n")
             temp.path <- paste(sep="", temp.GO.directory, "/", temp.keys[1, stage], "_", temp.keys[1, Annotation.class.pasted], "_", temp.keys[1, temp.ontology], "_enrichGO.result.RDS")
             ## print(temp.path)
             if (file.info(temp.path)$size %in% c(NA, 0)==TRUE){
                 temp.enrichGO.result <- enrichGO(.[temp.keys[, list(stage, Annotation.class.pasted)], Gene_ID.no.version], OrgDb='org.Hs.eg.db', keyType="ENSEMBL", ont=temp.keys[1, temp.ontology], pvalueCutoff=0.01)
                 print(temp.enrichGO.result)
                 saveRDS(temp.enrichGO.result, temp.path)
             }
             readRDS(temp.path)@result %>% data.table(temp.keys) %>% {.[pvalue<0.1]} ## increases statistical power
         }
    }  %>% rbindlist %>%
    {.[, p.adjusted.BH:=p.adjust(p=pvalue, method="BH")]} %>%
    {.[p.adjusted.BH<0.1]} %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])][, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))]}

fwrite(GO.result.dt, snakemake@output[["GO_result_dt_txt_gz_filename"]])


F3D.ggplot <- GO.result.dt %>%
    {.[, p.adjusted.BH.cut:=cut(x=p.adjusted.BH, breaks=c(0, 0.01, 0.05, 0.1, 1))]} %>%
    {.[p.adjusted.BH<0.05]} %>%
    {.[, Description.occurrence:=.N, Description][Description.occurrence>=3]} %>%
    {.[, Description.renamed:= Description]} %>%
    {ggplot(., aes(x="a", y = Description.renamed, fill=p.adjusted.BH.cut))  + geom_tile()  + scale_fill_manual(values=seq_gradient_pal("red", "grey90")(c(0, 2/3, 1))) + facet_grid(temp.ontology~Annotation.class.pasted.description + stage.description.ordered, scale="free_y", space="free") + labs(x="", y="", fill="BH-adjusted p-value") + theme_pubr(base_size=10) + theme(axis.text.x=element_blank(), axis.text.y=element_text(lineheight=0.6)) + guides(fill=guide_legend(nrow=1))}
F3D.ggplot

saveRDS(F3D.ggplot, snakemake@output[["main_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(F3D.ggplot, "report/201218-fifth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F3D.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["main_png_filename"]], plot=F3D.ggplot, width.r=0.9, height.r=0.32)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F3D.png", plot=F3D.ggplot, width.r=0.9, height.r=0.32)
}


print(F3D.ggplot)

image_read(snakemake@output[["main_png_filename"]]) %>% print


SF10.ggplot <- GO.result.dt %>%
    {.[, p.adjusted.BH.cut:=cut(x=p.adjusted.BH, breaks=c(0, 0.001, 0.01, 0.05, 0.1, 1))]} %>%
    {.[p.adjusted.BH<0.05]} %>%
    {ggplot(., aes(x=stage.description.ordered, y = Description, fill=p.adjusted.BH.cut))  + geom_tile()  + scale_fill_manual(values=seq_gradient_pal("red", "grey90")(c(0, 2/3, 1))) + facet_grid(temp.ontology~Annotation.class.pasted.description, scale="free_y", space="free") + labs(x="", y="", fill="BH-adjusted p-value") + theme_pubr(base_size=10) + theme(axis.text.x=element_text(angle=45, hjust=1)) + guides(fill=guide_legend(nrow=1))}
SF10.ggplot



saveRDS(SF10.ggplot, snakemake@output[["supp1_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(SF10.ggplot, "report/201218-fifth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF10.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["supp1_png_filename"]], plot=SF10.ggplot, width.r=1.3, height.r=1.3)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF10.png", plot=SF10.ggplot, width.r=1.3, height.r=1.3)
}


print(SF10.ggplot)

image_read(snakemake@output[["supp1_png_filename"]]) %>% print
