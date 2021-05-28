#' ---
#' title: "F4C summary"
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

GO.result.dt <- fread(snakemake@input[["F4B_GO_result_dt_txt_gz_filename"]])
if (FALSE){
    GO.result.dt <- fread("report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/F4A.and.F4AS.GO.result.dt.txt.gz")
}

## no >=3 stage terms for biparental-only terms
F4C.ggplot <- GO.result.dt %>%
    {.[, p.adjusted.BH.cut:=cut(x=p.adjusted.BH, breaks=c(0, 0.001, 0.01, 0.05, 0.1, 1))]} %>%
    {.[p.adjusted.BH<0.1]} %>%
    {.[, .SD %>% copy %>% {.[(Description %in% .SD[disease.description %in% c("parthenogenetic", "androgenetic"), Description]) == FALSE, disease.specific.Description:=TRUE]}, list(stage)][disease.specific.Description==TRUE]} %>%
    {.[p.adjusted.BH<0.05]} %>%
    {.[, Description.renamed:=
             Description %>%
             {sub(pattern="mitochondrion localization, ", replacement="mitochondrion localization,\n", x=.)}
       ]} %>%
    {.[, Description.renamed.ordered:=factor(Description.renamed, .[order(Annotation.class.pasted.description), Description.renamed %>% unique])]} %>%
    ## reset stage order
    {.[, `:=`(stage.description=NULL, stage.description.ordered=NULL)]; temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])]} %>%
    {ggplot(.[temp.ontology=='BP'], aes(x=disease.description, y = Description.renamed.ordered, fill=disease.description))  + geom_tile()  + facet_grid(temp.ontology~Annotation.class.pasted.description+stage.description.ordered, scale="free_y", space="free") + labs(x="", y="", fill="") + theme_pubr(base_size=10) + theme(axis.text.x=element_blank(), axis.text.y=element_text(lineheight=0.5)) + scale_fill_manual(values=c("#00BA38"))}
F4C.ggplot

saveRDS(F4C.ggplot, snakemake@output[["main_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(F4C.ggplot, "report/201218-fifth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F4C.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["main_png_filename"]], plot=F4C.ggplot, width.r=0.9, height.r=0.35)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/F4C.png", plot=F4C.ggplot, width.r=0.9, height.r=0.35)
}


print(F4C.ggplot)

image_read(snakemake@output[["main_png_filename"]]) %>% print

