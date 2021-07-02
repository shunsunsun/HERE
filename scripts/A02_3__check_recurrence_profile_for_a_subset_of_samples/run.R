library("data.table")
library("magrittr")
source("./scripts/common/logger.R")

merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt <- fread(snakemake@input[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_merged_with_coverage_dt_txt_gz_filename"]])

if (FALSE){
    merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt <- fread("result/S52_3__mark_unsequenced_editing_sites/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt.txt.gz")
}

#### 20210622 test
## merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt
## GSE133854 uniparental disomy
## GSE65481 (zygote viability good/bad prediction)
## ICM v.s. hESC (customized) v.s. hESC(long-term cell line; can be extracted from all.normal.samples)
## post-implantation, developmental day (focus on individual recurrent edits with abnormally high signals; can be extracted from all.normal.samples)
## GSE95477 (oocyte age)
## GSE101571 amanitin
## GSE100118 (blastocyst OCT4 CRISPR knockdown)


## phenotype.dt <- fread("result/S21_1__merge_phenotype_tables/201221-fifth-phenotype-collection/phenotype.output.at.gsm.level.dt.txt")
#### END of 20210622 test

subset.name <- snakemake@wildcards[['SUBSET_NAME']]

if (FALSE){
    subset.name <- "all.normal.samples"
    subset.name <- "GSE133854.all"
    subset.name <- "GSE100118.all.and.all.normal.blastocysts"
    subset.name <- "GSE95477.all"
}

report.expr(merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt %>% dim)


## start subsetting

subset.dt <- NULL

if (subset.name == "all.normal.samples"){
    subset.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt[is.normal==TRUE]
    subset.dt[, group:=paste(sep="", stage, "@", is.normal)]
} else if (subset.name == "GSE133854.all"){
    subset.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt[gse=="GSE133854"]
    subset.dt[, group:=paste(sep="", stage, "@", disease)]
    subset.dt[disease %in% c("", NA) == TRUE, group:=paste(sep="", stage, "@", "biparental")]
} else if (subset.name == "GSE65481.all"){
    subset.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt[gse=="GSE65481"]
    subset.dt[, group:=paste(sep="", stage, "@", disease)]
} else if (subset.name == "GSE95477.all"){
    subset.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt[gse=="GSE95477"]
    subset.dt[, group:=paste(sep="", stage, "@", maternal.age>=35)]
} else if (subset.name == "GSE100118.all.and.all.normal.blastocysts"){
    subset.dt <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt[gse=="GSE100118" | (is.normal==TRUE & grepl("blastocyst", stage))]
    subset.dt[, group:="none"]
    subset.dt[grepl("CRISPR", treatment), group:="case"]
    subset.dt[grepl("Cas9", treatment), group:="control"]
} ##else if (subset.name == "normal.samples.from.the.largest.two.datasets.from.two.labs"){
    ## temp. <- merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt[is.normal==TRUE & gse %in% c("")]
## }

if (is.null(subset.dt) == TRUE){
    stop(paste0("Unsupported subset.name ", subset.name, "\n"))
}


## generate datasets needed by downstream analyses

## 1. subset.dt
fwrite(subset.dt, snakemake@output[["subset_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(subset.dt, paste(sep="", "result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/subset.dt.txt.gz"))
}

## 2. recurrence profile (for F2A)

subset.site.recurrence.comparison.dt <- subset.dt %>%
    {.[is.na(AC)==FALSE,
       list(site.occurrence.for.this.group=.N),
       list(CHROM, POS, group)]}

fwrite(subset.site.recurrence.comparison.dt, snakemake@output[["subset_site_recurrence_comparison_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(subset.site.recurrence.comparison.dt, paste(sep="", "result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/subset.site.recurrence.comparison.dt.txt.gz"))
}

## 3. recurrent edits (site + group only; for F2B)

if (subset.name == 'all.normal.samples'){
    subset.site.recurrence.comparison.recurrent.edits.only.dt <- subset.site.recurrence.comparison.dt[site.occurrence.for.this.group >= 10]
} else if (subset.name == 'GSE133854.all'){
    subset.site.recurrence.comparison.recurrent.edits.only.dt <- subset.site.recurrence.comparison.dt %>%
        {.[(grepl('zygote', group)==TRUE & site.occurrence.for.this.group >= 3) |
           (grepl('2-cell', group)==TRUE & site.occurrence.for.this.group >= 5) |
           (group == "4-cell@biparental" & site.occurrence.for.this.group >= 10) |
           (group %in% c("4-cell@parthenogenetic", "4-cell@androgenetic") == TRUE & site.occurrence.for.this.group >= 5) |
           (grepl('(8-cell|morula)', group) == TRUE & site.occurrence.for.this.group >= 10), ]}
} else if (subset.name == 'GSE100118.all.and.all.normal.blastocysts'){
    subset.site.recurrence.comparison.recurrent.edits.only.dt <- subset.site.recurrence.comparison.dt[site.occurrence.for.this.group >= 10]
} else if (subset.name == 'GSE65481.all'){
    subset.site.recurrence.comparison.recurrent.edits.only.dt <- subset.site.recurrence.comparison.dt[site.occurrence.for.this.group >= 10]
} else if (subset.name == 'GSE95477.all'){
    subset.site.recurrence.comparison.recurrent.edits.only.dt <- subset.site.recurrence.comparison.dt[site.occurrence.for.this.group >= 4]
}


fwrite(subset.site.recurrence.comparison.recurrent.edits.only.dt, snakemake@output[["subset_site_recurrence_comparison_recurrent_edits_only_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(subset.site.recurrence.comparison.recurrent.edits.only.dt, paste(sep="", "result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/subset.site.recurrence.comparison.recurrent.edits.only.dt.txt.gz"))
}

## 4. subset.dt with recurrent edits only (for supp. fig. of F2B; F2D and supp. fig. of F2D (preprocessing) -- per-dataset distribution)
## here the recurrency is group-specific (i.e., if an editing site is recurrent in group A but not in group B, then the final subset.dt will only have this site kept in group A samples but not group B samples)

subset.recurrent.edits.only.dt <- merge(x=subset.dt, y=subset.site.recurrence.comparison.recurrent.edits.only.dt, by=c("CHROM", "POS", "group"), all.x=FALSE)

fwrite(subset.recurrent.edits.only.dt, snakemake@output[["subset_recurrent_edits_only_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(subset.recurrent.edits.only.dt, paste(sep="", "result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/subset.recurrent.edits.only.dt.txt.gz"))
}

## 5. subset.dt with recurrence comparison info (for F2C and its supp.)
## must use CHROM.and.POS rather than CHROM, POS; otherwise all combinations of CHROM and POS will appear

subset.site.recurrence.comparison.CJ.dt <- subset.site.recurrence.comparison.dt %>%
    {.[, list(CHROM.and.POS=paste(sep="", CHROM, "_", POS), group, site.occurrence.for.this.group)]} %>%
    setkey(CHROM.and.POS, group) %>%
    {.[CJ(CHROM.and.POS, group, unique=TRUE)]} %>%
    setnafill(type="const", fill=0, cols="site.occurrence.for.this.group") %>%
    {.[, `:=`(
         CHROM=sub(pattern="_.*", replacement="", x=CHROM.and.POS),
         POS=sub(pattern=".*_", replacement="", x=CHROM.and.POS) %>% as.integer
     )]}

fwrite(subset.site.recurrence.comparison.CJ.dt, snakemake@output[["subset_site_recurrence_comparison_CJ_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(subset.site.recurrence.comparison.CJ.dt, paste(sep="", "result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/subset.site.recurrence.comparison.CJ.dt.txt.gz"))
}

## 6. subset.dt with recurrent edits only plus gene annotations (for F2E and supp. fig. of F2E)

### 6.1. annotations for subsetted & valid sites (falls onto a gene, only A>G edits, only on protein-coding transcripts)

#### 6.1.1. all subsetted sites
merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt <- fread(snakemake@input[["merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_variant_only_snpEff_annotation_dt_txt_gz_filename"]])

if (FALSE) {
    merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt <- fread("result/S51_6__get_snpEff_annotation_subset_of_filtered_result/210215-sixth-dataset/201221-fifth-phenotype-collection/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt.txt.gz")
}

snpEff.annotation.for.subset.recurrent.edits.dt <- merge(
    x=subset.recurrent.edits.only.dt[, list(CHROM, POS)] %>% unique,
    y=merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt[, list(CHROM, POS, Annotation, Annotation_Impact, Gene_Name, Gene_ID, Feature_Type, Feature_ID, Transcript_BioType, Rank, HGVS.c, HGVS.p, `cDNA.pos / cDNA.length`, `CDS.pos / CDS.length`, `AA.pos / AA.length`, Distance, `ERRORS / WARNINGS / INFO`, event)],
    by=c("CHROM", "POS"), all.x=FALSE, all.y=FALSE)

fwrite(snpEff.annotation.for.subset.recurrent.edits.dt, snakemake@output[["snpEff_annotation_for_subset_recurrent_edits_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(snpEff.annotation.for.subset.recurrent.edits.dt, paste(sep="", "result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/snpEff.annotation.for.subset.recurrent.edits.dt.txt.gz"))
}

## 6.1.2. valid subsetted sites, transcript level
snpEff.annotation.for.subset.recurrent.edits.on.valid.transcripts.dt <- snpEff.annotation.for.subset.recurrent.edits.dt %>%
    {.[Annotation %in% c("downstream_gene_variant", "intergenic_region", "upstream_gene_variant") == FALSE]} %>%
    {.[event == 'A>G']} %>%
    {.[Transcript_BioType=='protein_coding']}

fwrite(snpEff.annotation.for.subset.recurrent.edits.on.valid.transcripts.dt, snakemake@output[["snpEff_annotation_for_subset_recurrent_edits_on_valid_transcripts_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(snpEff.annotation.for.subset.recurrent.edits.on.valid.transcripts.dt, paste(sep="", "result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/snpEff.annotation.for.subset.recurrent.edits.on.valid.transcripts.dt.txt.gz"))
}


## 6.1.3. valid subsetted sites, gene level (with Annotation collapsed)
snpEff.annotation.for.subset.recurrent.edits.on.valid.genes.dt <- snpEff.annotation.for.subset.recurrent.edits.on.valid.transcripts.dt %>%
    {.[, list(CHROM, POS, Annotation, Gene_Name, Gene_ID)]} %>%
    unique %>%
    {
        cat(date(), "generating Annotation.pasted ...\n");
        .[, Annotation.pasted:=paste(collapse=";", Annotation %>% sort %>% unique), list(CHROM, POS, Gene_Name, Gene_ID)]
    } %>%
    {.[, Annotation.corrected:="intron_variant"]} %>%
    {
        cat(date(), "looping over terms...\n");
        for (term in c("synonymous_variant", "stop_retained_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", 'missense_variant', "splice_region_variant", "splice_acceptor_variant", "splice_donor_variant", "5_prime_UTR_premature_start_codon_gain_variant", "stop_lost", "start_lost")){
            cat(date(), "  updating term: ", term, "\n");
            .[grepl(term, Annotation.pasted), Annotation.corrected:=term]
        }
        .
    } %>%
    {
        .[, Annotation.class:="exonic.or.splicing.related"]
        .[Annotation.corrected=='intron_variant', Annotation.class:="purely.intronic"]
    }

## each CHROM+POS+Gene_Name+Gene_ID has only one possible Annotation class
if (FALSE) {
    snpEff.annotation.for.subset.recurrent.edits.on.valid.genes.dt[, list(CHROM, POS, Gene_Name, Gene_ID, Annotation.class)] %>% unique %>% {.[, .N, list(CHROM, POS, Gene_Name, Gene_ID)][N!=1]}
'
    Empty data.table (0 rows and 5 cols): CHROM,POS,Gene_Name,Gene_ID,N
'
}

fwrite(snpEff.annotation.for.subset.recurrent.edits.on.valid.genes.dt, snakemake@output[["snpEff_annotation_for_subset_recurrent_edits_on_valid_genes_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(snpEff.annotation.for.subset.recurrent.edits.on.valid.genes.dt, paste(sep="", "result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/snpEff.annotation.for.subset.recurrent.edits.on.valid.genes.dt.txt.gz"))
}



## 6.1.4. valid subsetted sites, gene level (with Annotation collapsed) and merged with per-sample stats
#### Note: the plotting should discard site x Gene x SAMPLE combination when the site is unsequenced in that sample
subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt <- merge(
    x=subset.recurrent.edits.only.dt,
    y=snpEff.annotation.for.subset.recurrent.edits.on.valid.genes.dt,
    by=c("CHROM", "POS"), all=FALSE)


if (FALSE) {
    subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt[is.na(AC)==FALSE, .N, list(SAMPLE, stage, Annotation.class)] %>% {dcast(., SAMPLE+stage~Annotation.class, value.var="N")[, log2.ES.I.ratio:=log2(exonic.or.splicing.related/purely.intronic)]} %>% {.[, quantile(log2.ES.I.ratio, na.rm=TRUE) %>% t %>% data.table, list(stage)]}
'
              stage         0%         25%        50%      75%     100%
 1:          zygote  1.9343238  2.03493379  2.0747678 2.094367 2.276245
 2:          2-cell  1.5482677  1.70147599  1.7349735 1.787952 1.916952
 3:          4-cell  1.1053530  1.72582504  1.8061988 1.898120 2.194549
 4:          8-cell  0.7136958  1.32325329  1.4712873 1.669692 2.095924
 5:          morula -0.4150375  0.73696559  1.1595053 1.691071 3.502500
 6:            hESC  0.0000000  0.85533621  1.3683874 2.155278 3.584963
 7:      zygote.2PN  1.8505422  1.85540222  1.8611019 1.873287 1.930472
 8:      oocyte.MII  1.5785058  1.73568786  1.8053333 1.854991 2.213481
 9: blastocyst.late -1.0000000 -1.00000000 -1.0000000 0.000000 0.000000
10:             ICM -2.5849625  0.06875176  0.5849625 1.057739 3.700440
11:       oocyte.GV  1.6102055  1.73406013  1.7483122 1.771719 1.814769
12:     trophoblast -3.5849625 -1.00000000  0.0000000 1.000000 3.523562
13:             STB -1.0000000  0.93586966  1.3066613 1.758992 4.000000
14:             CTB -0.3625701  1.45943162  1.8443491 2.442943 5.700440
15:             MTB  1.5849625  3.11485790  3.5849625 3.671570 3.807355
16:        epiblast -3.0000000  1.28950662  1.8365013 2.751250 5.169925
17:       hypoblast -1.0000000  1.24144607  1.3219281 1.688722 3.169925
18:             EVT  0.7776076  1.93359749  2.1888353 2.844070 4.066089
'
}


fwrite(subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt, snakemake@output[["subset_recurrent_edits_only_with_snpEff_annotation_on_valid_genes_only_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt, paste(sep="", "result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz"))
}



## 7. subset.dt plus gene annotations (for sample-level comparison plot)

#### 7.1. all subsetted sites

snpEff.annotation.for.subset.dt <- merge(
    x=subset.dt,
    y=merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt[, list(CHROM, POS, Annotation, Annotation_Impact, Gene_Name, Gene_ID, Feature_Type, Feature_ID, Transcript_BioType, Rank, HGVS.c, HGVS.p, `cDNA.pos / cDNA.length`, `CDS.pos / CDS.length`, `AA.pos / AA.length`, Distance, `ERRORS / WARNINGS / INFO`, event)],
    by=c("CHROM", "POS"), all.x=FALSE, all.y=FALSE, allow.cartesian=TRUE)

fwrite(snpEff.annotation.for.subset.dt, snakemake@output[["snpEff_annotation_for_subset_dt_txt_gz_filename"]])
if (FALSE){
    fwrite(snpEff.annotation.for.subset.dt, paste(sep="", "result/A02_3__check_recurrence_profile_for_a_subset_of_samples/210215-sixth-dataset/201221-fifth-phenotype-collection/", subset.name, "/snpEff.annotation.for.subset.dt.txt.gz"))
}

