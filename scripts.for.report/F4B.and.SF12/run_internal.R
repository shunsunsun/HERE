#' ---
#' title: "F4B and SF12 summary"
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

subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt <- fread(snakemake@input[["subset_recurrent_edits_only_with_snpEff_annotation_on_valid_genes_only_dt_txt_gz_filename"]])

if (FALSE){
    subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt <- fread("result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz")
}

for.plot.dt <- subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt %>%
    {.[is.na(AC)==FALSE]} %>%
    {.[, list(Annotation.class.pasted=paste(collapse=";", Annotation.class %>% sort %>% unique)), list(Gene_ID, Gene_Name, SAMPLE, gse, group, stage)]} %>%
    {.[, disease.description:=sub(pattern=".*@", replacement="", x=group)]} %>%
    ## compute % occurrence of each gene x Annotation.class.pasted
    ## note here that for total sample number per stage, `unique` must be used instead of `.N`
    {.[, total.sample.count:=length(unique(SAMPLE)), list(group, stage, disease.description)]} %>%
    {.[, list(gene.occurrence=.N), list(Gene_ID, Gene_Name, group, stage, disease.description, total.sample.count, Annotation.class.pasted)]} %>%
    {.[, gene.occurrence.pct:=gene.occurrence/total.sample.count]} %>%
    ## finish computing
    {.[, Annotation.class.pasted.description:=c('exonic.or.splicing.related'='exonic', 'purely.intronic'='intronic', 'exonic.or.splicing.related;purely.intronic'='mixed')[Annotation.class.pasted]]} %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])]}

fwrite(for.plot.dt, snakemake@output[["F4B_for_plot_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(for.plot.dt, "report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/F4B.for.plot.dt.txt.gz")
}


temp.GO.directory <- paste(sep="", snakemake@params[["result_directory"]], "/F4B.and.SF12.temp/GO")
if (FALSE) {
    temp.GO.directory <- "report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/F4B.and.SF12.temp/GO"
}

dir.create(paste(sep="", temp.GO.directory), recursive=TRUE)

GO.result.dt <- for.plot.dt %>%
    {.[gene.occurrence.pct>=0.8]} %>%
    setkey(stage, disease.description, Annotation.class.pasted) %>%
    {.[, Gene_ID.no.version:=sub(pattern="\\.[0-9]*$", replacement="", x=Gene_ID)]} %>%
    {foreach(
         temp.keys=.[, list(stage, disease.description, Annotation.class.pasted, Annotation.class.pasted.description)] %>%
             unique %>%
             as.data.frame(stringsAsFactors=FALSE) %>%
             merge(y=data.frame(temp.ontology=c("BP", "MF", "CC"), stringsAsFactors=FALSE)) %>%
             data.table %>%
             iter(by='row')) %do%
         {
             cat(date(), " Processing ", temp.keys[1, stage], " x ", temp.keys[1, disease.description], " x ", temp.keys[1, Annotation.class.pasted], " x ", temp.keys[1, temp.ontology], "\n")
             temp.path <- paste(sep="", temp.GO.directory, "/", temp.keys[1, stage], "_", temp.keys[1, disease.description], "_", temp.keys[1, Annotation.class.pasted], "_", temp.keys[1, temp.ontology], "_enrichGO.result.RDS")
             ## print(temp.path)
             if (file.info(temp.path)$size %in% c(NA, 0)==TRUE){
                 temp.enrichGO.result <- enrichGO(.[temp.keys[, list(stage, disease.description, Annotation.class.pasted)], Gene_ID.no.version], OrgDb='org.Hs.eg.db', keyType="ENSEMBL", ont=temp.keys[1, temp.ontology], pvalueCutoff=0.01)
                 print(temp.enrichGO.result)
                 saveRDS(temp.enrichGO.result, temp.path)
             }
             temp.enrichGO.result <- readRDS(temp.path)
             if (is.null(temp.enrichGO.result) == TRUE){
                 NULL
             } else {
                 temp.enrichGO.result@result %>% data.table(temp.keys) %>% {.[pvalue<0.1]} ## increases statistical power
             }
         }
    }  %>% rbindlist %>%
    {.[, p.adjusted.BH:=p.adjust(p=pvalue, method="BH")]} %>%
    {.[p.adjusted.BH<0.1]} %>%
    {temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])][, category.ordered:=factor(category, levels=unique(temp.stage.dt[, category]))]}

fwrite(GO.result.dt, snakemake@output[["F4B_GO_result_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(GO.result.dt, "report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/F4B.GO.result.dt.txt.gz")
}

F4B.ggplot <- GO.result.dt %>%
    {.[, p.adjusted.BH.cut:=cut(x=p.adjusted.BH, breaks=c(0, 0.001, 0.01, 0.05, 0.1, 1))]} %>%
    {.[p.adjusted.BH<0.1]} %>%
    {.[, .SD %>% copy %>% {.[(Description %in% .SD[disease.description=="biparental", Description]) == FALSE, disease.specific.Description:=TRUE]}, list(stage)][disease.specific.Description==TRUE]} %>%
    {.[, Description.occurrence:=.N, list(Description, disease.description, Annotation.class.pasted.description)][Description.occurrence>=3]} %>%
    {.[, Description.renamed:=
             Description %>%
             {sub(pattern="^cytoplasmic pattern recognition receptor ", replacement="cytoplasmic pattern recognition receptor\n", x=.)} %>%
             {sub(pattern="viral-induced cytoplasmic pattern recognition receptor", replacement="viral-incuded cytoplasmic\npattern recognition receptor", x=.)} %>%
             {sub(pattern="oxidoreductase activity, ", replacement="oxidoreductase activity,\n", x=.)} %>%
             {sub(pattern="stress-activated protein kinase", replacement="stress-activated protein\nkinase", x=.)}
       ]} %>%
    {.[, Description.renamed.ordered:=factor(Description.renamed, .[order(Annotation.class.pasted.description), Description.renamed %>% unique])]} %>%
    ## reset stage order
    {.[, `:=`(stage.description=NULL, stage.description.ordered=NULL)]; temp.stage.dt <- read_excel("./manuscript/table_for_processing.xlsx", sheet="stage.properties") %>% data.table; merge(x=., y=temp.stage.dt, by.x="stage", by.y="stage", all.x=TRUE, all.y=FALSE)[, stage.description.ordered:=factor(stage.description, levels=temp.stage.dt[, stage.description])]} %>%
    {ggplot(.[temp.ontology=='BP'], aes(x=disease.description, y = Description.renamed.ordered, fill=disease.description))  + geom_tile()  + scale_fill_manual(values=c("#F8766D", "#00BFC4")) + facet_grid(temp.ontology~Annotation.class.pasted.description+stage.description.ordered, scale="free_y", space="free") + labs(x="", y="", fill="") + theme_pubr(base_size=10) + theme(axis.text.x=element_blank(), axis.text.y=element_text(lineheight=0.5)) + guides(fill=guide_legend(nrow=1)) + scale_fill_manual(values=c("#F8766D", "#619CFF"))} ## "#00BA38" for green
F4B.ggplot

saveRDS(F4B.ggplot, snakemake@output[["main_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(F4B.ggplot, "report/201218-fifth-dataset/201221-fifth-phenotype-collection/all.normal.samples/F4B.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["main_png_filename"]], plot=F4B.ggplot, width.r=0.9, height.r=0.45)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/F4B.png", plot=F4B.ggplot, width.r=0.9, height.r=0.45)
}


print(F4B.ggplot)

image_read(snakemake@output[["main_png_filename"]]) %>% print


human.pc.gene.info.dt <- fread(cmd=paste(sep="", ' grep -v "#" ', 'external/reference.gene.annotation/GENCODE.annotation/32/gencode.annotation.gtf', ' | grep -P "\tgene\t" | grep -P \'gene_type "protein_coding"\' | cut -f 9 | sed -E -e \'s@gene_id "([^;]+)";.*gene_name "([^;]+)";.*@\\1\\t\\2@\' '), header=FALSE, col.names=c("gene.id", "gene.name") ) %>% {.[, gene.id.no.version:=sub(pattern="\\.[0-9]+", replacement="", x=gene.id)]}

GO.result.filtered.dt <- GO.result.dt %>%
    {.[, p.adjusted.BH.cut:=cut(x=p.adjusted.BH, breaks=c(0, 0.001, 0.01, 0.05, 0.1, 1))]} %>%
    {.[p.adjusted.BH<0.1]} %>%
    {.[, .SD %>% copy %>% {.[(Description %in% .SD[disease.description=="biparental", Description]) == FALSE, disease.specific.Description:=TRUE]}, list(stage)][disease.specific.Description==TRUE]} %>%
    {.[, Description.occurrence:=.N, list(Description, disease.description, Annotation.class.pasted.description)][Description.occurrence>=3]} %>%
    {.[, data.table(gene.id=strsplit(x=geneID, split="/")[[1]]), list(disease.description, Annotation.class.pasted.description, stage.description.ordered, Description)]} %>%
    {merge(x=., y=human.pc.gene.info.dt[, list(gene.id.no.version, gene.name)], by.x="gene.id", by.y="gene.id.no.version", all.x=TRUE, all.y=FALSE)}


SF12.ggplot <- GO.result.filtered.dt %>%
    {.[Description %in% .[gene.name=='PTK2B', Description]]} %>%
    {.[, Description.type:="others"][grepl("calcium", Description), Description.type:="calcium-related"][grepl("(MAPK|JUN|JNK)", Description), Description.type:="JUN/MAPK-related"]} %>%
    {ggplot(., aes(x=gene.name, y = Description, color=Description.type))  + geom_point()  + facet_grid(stage.description.ordered~disease.description, scale="free", space="free") + labs(x="Genes associated with PTK2B-annotated GO terms\n(Labelled with points only when the term is enriched)", y="Enriched PTK2B-annotated GO terms", color="Type of GO terms") + theme_pubr(base_size=10) + theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_text(lineheight=0.5)) + guides(fill=guide_legend(nrow=1)) + scale_color_manual(values=c("red", "blue", "black"))}


    



saveRDS(SF12.ggplot, snakemake@output[["supp1_ggplot_RDS_filename"]])

if (FALSE){
    saveRDS(SF12.ggplot, "report/201218-fifth-dataset/201221-fifth-phenotype-collection/all.normal.samples/SF12.ggplot.RDS")
}

ggsave.A4(filename=snakemake@output[["supp1_png_filename"]], plot=SF12.ggplot, width.r=1.2, height.r=0.8)
if (FALSE){
    ggsave.A4(filename="report/210215-sixth-dataset/201221-fifth-phenotype-collection/GSE133854.all/SF12.png", plot=SF12.ggplot, width.r=1.2, height.r=0.8)
}


print(SF12.ggplot)

image_read(snakemake@output[["supp1_png_filename"]]) %>% print
