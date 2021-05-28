import pandas

for temp_key in ['threads_generate_vcf', 'threads_annotate_with_snpEff', 'threads_generate_bed']:
    if temp_key not in config.keys():
        config[temp_key] = 20

for temp_key in ['snpEff_Xmx']:
    if temp_key not in config.keys():
        config[temp_key] = "-Xmx150G"



def __get_dataset_name_collection(dataset_collection_name, prefix_for_dataset_collection_name_directory="external/DATASET_COLLECTION_NAME_DIRECTORY/"):
    dataset_name_collection = []
    with open(prefix_for_dataset_collection_name_directory + dataset_collection_name, "r") as temp_fileobject:
        dataset_name_collection = [temp_filename.strip() for temp_filename in temp_fileobject.readlines() if temp_filename != '\n']
    return dataset_name_collection


def __get_parameters_collection_for_RNA_editing_calling(dataset_collection_name, prefix_for_dataset_name_directory="external/DATASET_RNA_EDITING_NAME_DIRECTORY/", prefix_for_dataset_collection_name_directory="external/DATASET_RNA_EDITING_COLLECTION_NAME_DIRECTORY/"):
    dataset_name_collection = __get_dataset_name_collection(dataset_collection_name, prefix_for_dataset_collection_name_directory=prefix_for_dataset_collection_name_directory)
    final_collection = []
    for temp_dataset_name in dataset_name_collection:
        temp_dataset_DataFrame = pandas.read_csv( prefix_for_dataset_name_directory + temp_dataset_name, sep=",")
        temp_row_collection = [ temp_row for temp_rowindex, temp_row in temp_dataset_DataFrame.iterrows()]
        final_collection = final_collection + temp_row_collection
    return final_collection





rule A02_3__check_recurrence_profile_for_a_subset_of_samples:
    input:
        flag_S52_3="result/S52_3__mark_unsequenced_editing_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_merged_with_coverage_dt_txt_gz_filename="result/S52_3__mark_unsequenced_editing_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt.txt.gz",
        flag_S51_6="result/S51_6__get_snpEff_annotation_subset_of_filtered_result/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_variant_only_snpEff_annotation_dt_txt_gz_filename="result/S51_6__get_snpEff_annotation_subset_of_filtered_result/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.variant.only.snpEff.annotation.dt.txt.gz"
    output:
        flag_A02_3=touch("result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/finished"),
        subset_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.dt.txt.gz",
        subset_site_recurrence_comparison_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.site.recurrence.comparison.dt.txt.gz",
        subset_site_recurrence_comparison_recurrent_edits_only_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.site.recurrence.comparison.recurrent.edits.only.dt.txt.gz",
        subset_recurrent_edits_only_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.recurrent.edits.only.dt.txt.gz",
        subset_site_recurrence_comparison_CJ_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.site.recurrence.comparison.CJ.dt.txt.gz",
        snpEff_annotation_for_subset_recurrent_edits_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/snpEff.annotation.for.subset.recurrent.edits.dt.txt.gz",
        snpEff_annotation_for_subset_recurrent_edits_on_valid_transcripts_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/snpEff.annotation.for.subset.recurrent.edits.on.valid.transcripts.dt.txt.gz",
        snpEff_annotation_for_subset_recurrent_edits_on_valid_genes_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/snpEff.annotation.for.subset.recurrent.edits.on.valid.genes.dt.txt.gz",
        subset_recurrent_edits_only_with_snpEff_annotation_on_valid_genes_only_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz"
    script:
        "./scripts/A02_3__check_recurrence_profile_for_a_subset_of_samples/run.R"


rule sA02_5__generate_all_transcribable_Alu_A_vcf____step01__generate_Alu_only_vcf:
    input:
        contigs_fasta_filename="external/contigs/{CONTIGS_FASTA_FILENAME}",
        repeatmasker_repFamily_Alu_bed_filename="external/UCSC.Table.Browser.repeatmasker/repFamily.Alu/repeatmasker.bed"
    output:
        flag_sA02_5_step01=touch("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/finished.step01__generate_Alu_only_vcf")
    params:
        result_directory="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/"
    threads:
        config['threads_generate_vcf']
    script:
        "./scripts/sA02_5__generate_all_transcribable_Alu_A_vcf/step01__generate_Alu_only_vcf/run.R"


rule sA02_5__generate_all_transcribable_Alu_A_vcf____step02__annotate_a_single_Alu_only_vcf_with_snpEff:
    input:
        flag_sA02_5_step01="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/finished.step01__generate_Alu_only_vcf"
    output:
        flag_sA02_5_step02=touch("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.step02__annotate_a_single_Alu_only_vcf_with_snpEff.{CHROMOSOME}"),
        single_contig_Alu_only_snpEff_annotated_vcf_gz_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.vcf.gz"
    params:
        single_contig_Alu_only_vcf_gz_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{CHROMOSOME}.Alu.only.vcf.gz",
        snpEff_config_filename="result/x16_2__annotate_merged_sites_from_a_dataset_collection/step01__annotate_with_snpEff/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/snpEff.config",
        snpEff_database_name="{CONTIGS_FASTA_FILENAME}.GENCODE.{GENCODE_VERSION}",
        snpEff_Xmx=config['snpEff_Xmx']
    threads:
        config['threads_annotate_with_snpEff']
    shell:
        """
        snpEff ann {params.snpEff_Xmx} -verbose -no-intergenic -no-upstream -no-downstream -config {params.snpEff_config_filename} {params.snpEff_database_name} {params.single_contig_Alu_only_vcf_gz_filename} | bcftools view -o {output.single_contig_Alu_only_snpEff_annotated_vcf_gz_filename} -Oz
        sleep 1
        tabix {output.single_contig_Alu_only_snpEff_annotated_vcf_gz_filename}
        """


rule sA02_5__generate_all_transcribable_Alu_A_vcf____step03__extract_annotated_records_from_a_single_Alu_only_vcf:
    input:
        flag_sA02_5_step02="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.step02__annotate_a_single_Alu_only_vcf_with_snpEff.{CHROMOSOME}",
        single_contig_Alu_only_snpEff_annotated_vcf_gz_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.vcf.gz"
    output:
        flag_sA02_5_step03=touch("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.step03__extract_annotated_records_from_a_single_Alu_only_vcf.{CHROMOSOME}"),
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_vcf_gz_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.vcf.gz"
    threads:
        config['threads_generate_vcf']
    shell:
        """
        bcftools view --threads {threads} --output-file {output.single_contig_Alu_only_snpEff_annotated_annotated_records_only_vcf_gz_filename} -Oz --include 'INFO/ANN ~ "protein_coding.*A>G"' {input.single_contig_Alu_only_snpEff_annotated_vcf_gz_filename}
        sleep 1
        tabix {output.single_contig_Alu_only_snpEff_annotated_annotated_records_only_vcf_gz_filename}
        """


rule sA02_5__generate_all_transcribable_Alu_A_vcf____step04__split_a_single_annotated_records_of_Alu_only_vcf_into_subsets:
    input:
        flag_sA02_5_step03="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.step03__extract_annotated_records_from_a_single_Alu_only_vcf.{CHROMOSOME}",
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_vcf_gz_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.vcf.gz"
    output:
        flag_sA02_5_step04=touch("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.step04__split_a_single_annotated_records_of_Alu_only_vcf_into_subsets.{CHROMOSOME}"),
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_intron_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.intron.bed",
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_3_prime_UTR_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.3.prime.UTR.bed",
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_5_prime_UTR_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.5.prime.UTR.bed",
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_CDS_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.CDS.bed"
    threads:
        config['threads_generate_bed']
    script:
        "./scripts/sA02_5__generate_all_transcribable_Alu_A_vcf/step04__split_a_single_annotated_records_of_Alu_only_vcf_into_subsets/run.R"


temp_chr_collection = ['chr' + str(temp_item) for temp_item in list(range(1, 23)) + ['X', 'Y']]

rule sA02_5__generate_all_transcribable_Alu_A_vcf____step05__combine_all_subsets:
    input:
        flag_sA02_5_step04_collection=expand("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.step04__split_a_single_annotated_records_of_Alu_only_vcf_into_subsets.{CHROMOSOME}", CHROMOSOME=temp_chr_collection, allow_missing=True),
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_intron_bed_filename_collection=expand("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.intron.bed", CHROMOSOME=temp_chr_collection, allow_missing=True),
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_3_prime_UTR_bed_filename_collection=expand("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.3.prime.UTR.bed", CHROMOSOME=temp_chr_collection, allow_missing=True),
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_5_prime_UTR_bed_filename_collection=expand("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.5.prime.UTR.bed", CHROMOSOME=temp_chr_collection, allow_missing=True),
        single_contig_Alu_only_snpEff_annotated_annotated_records_only_CDS_bed_filename_collection=expand("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/{CHROMOSOME}.contigs.Alu.only.snpEff.annotated.annotated.records.only.CDS.bed", CHROMOSOME=temp_chr_collection, allow_missing=True)
    output:
        flag_sA02_5_step04=touch("result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/finished.step05__combine_all_subsets"),
        all_contigs_Alu_only_snpEff_annotated_annotated_records_only_intron_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.intron.combined.bed",
        all_contigs_Alu_only_snpEff_annotated_annotated_records_only_3_prime_UTR_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.3.prime.UTR.combined.bed",
        all_contigs_Alu_only_snpEff_annotated_annotated_records_only_5_prime_UTR_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.5.prime.UTR.combined.bed",
        all_contigs_Alu_only_snpEff_annotated_annotated_records_only_CDS_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/{CONTIGS_FASTA_FILENAME}/{GENCODE_VERSION}/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.CDS.combined.bed"
    threads:
        1
    shell:
        """
        cat {input.single_contig_Alu_only_snpEff_annotated_annotated_records_only_intron_bed_filename_collection} > {output.all_contigs_Alu_only_snpEff_annotated_annotated_records_only_intron_combined_bed_filename}
        cat {input.single_contig_Alu_only_snpEff_annotated_annotated_records_only_3_prime_UTR_bed_filename_collection} > {output.all_contigs_Alu_only_snpEff_annotated_annotated_records_only_3_prime_UTR_combined_bed_filename}
        cat {input.single_contig_Alu_only_snpEff_annotated_annotated_records_only_5_prime_UTR_bed_filename_collection} > {output.all_contigs_Alu_only_snpEff_annotated_annotated_records_only_5_prime_UTR_combined_bed_filename}
        cat {input.single_contig_Alu_only_snpEff_annotated_annotated_records_only_CDS_bed_filename_collection} > {output.all_contigs_Alu_only_snpEff_annotated_annotated_records_only_CDS_combined_bed_filename}
        """



rule A02_5__check_transcribable_Alu_A_coverage_of_merged_bam:
    input:
        flag_sA02_5_step04="result/sA02_5__generate_all_transcribable_Alu_A_vcf/hg38.fa/32/finished.step05__combine_all_subsets",
        ## all_contigs_Alu_only_snpEff_annotated_annotated_records_only_intron_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/hg38.fa/32/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.intron.combined.bed",
        all_contigs_Alu_only_snpEff_annotated_annotated_records_only_3_prime_UTR_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/hg38.fa/32/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.3.prime.UTR.combined.bed",
        all_contigs_Alu_only_snpEff_annotated_annotated_records_only_5_prime_UTR_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/hg38.fa/32/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.5.prime.UTR.combined.bed",
        all_contigs_Alu_only_snpEff_annotated_annotated_records_only_CDS_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/hg38.fa/32/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.CDS.combined.bed",
        ## flag_S15_1__step02__part04="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/{INDEXER_PARAMETERS}/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/finished.step02__call_variants____part04__correct_coordinates",
        alignment_sorted_withRG_dedup_converted_bq_sorted_without_splicing_junction_SN_bam_filename="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/{INDEXER_PARAMETERS}/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.bam"
    output:
        flag_A02_5=touch("result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/finished"),
        ## alignment_merged_depth_on_transcribable_Alu_A_intron_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.intron.bed.txt.gz",
        alignment_merged_depth_on_transcribable_Alu_A_5_prime_UTR_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.5.prime.UTR.bed.txt.gz",
        alignment_merged_depth_on_transcribable_Alu_A_3_prime_UTR_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.3.prime.UTR.bed.txt.gz",
        alignment_merged_depth_on_transcribable_Alu_A_CDS_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.CDS.bed.txt.gz"
    threads:
        1
    shell:
        """
        samtools depth -d 1 -b {input.all_contigs_Alu_only_snpEff_annotated_annotated_records_only_5_prime_UTR_combined_bed_filename} {input.alignment_sorted_withRG_dedup_converted_bq_sorted_without_splicing_junction_SN_bam_filename} | gzip > {output.alignment_merged_depth_on_transcribable_Alu_A_5_prime_UTR_bed_txt_gz_filename}
        samtools depth -d 1 -b {input.all_contigs_Alu_only_snpEff_annotated_annotated_records_only_3_prime_UTR_combined_bed_filename} {input.alignment_sorted_withRG_dedup_converted_bq_sorted_without_splicing_junction_SN_bam_filename} | gzip > {output.alignment_merged_depth_on_transcribable_Alu_A_3_prime_UTR_bed_txt_gz_filename}
        samtools depth -d 1 -b {input.all_contigs_Alu_only_snpEff_annotated_annotated_records_only_CDS_combined_bed_filename} {input.alignment_sorted_withRG_dedup_converted_bq_sorted_without_splicing_junction_SN_bam_filename} | gzip > {output.alignment_merged_depth_on_transcribable_Alu_A_CDS_bed_txt_gz_filename}
        """



rule BA02_5__check_transcribable_Alu_A_coverage_of_merged_bam:
    input:
        lambda wildcards: ["result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME  + "/" + str(temp_row.INDEXER_PARAMETERS) +  "/finished" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    output:
        flag_BA02_5=touch("result/BA02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished")




rule A02_5__check_transcribable_Alu_A_coverage_of_merged_bam____patch01__intron:
    input:
        flag_sA02_5_step04="result/sA02_5__generate_all_transcribable_Alu_A_vcf/hg38.fa/32/finished.step05__combine_all_subsets",
        all_contigs_Alu_only_snpEff_annotated_annotated_records_only_intron_combined_bed_filename="result/sA02_5__generate_all_transcribable_Alu_A_vcf/hg38.fa/32/all.contigss.Alu.only.snpEff.annotated.annotated.records.only.intron.combined.bed",
        ## flag_S15_1__step02__part04="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/{INDEXER_PARAMETERS}/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/finished.step02__call_variants____part04__correct_coordinates",
        alignment_sorted_withRG_dedup_converted_bq_sorted_without_splicing_junction_SN_bam_filename="result/S15_1__get_sample_RNA_editing_sites_v3/{TYPE}/{DATASET}/{SAMPLE}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/{INDEXER_PARAMETERS}/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/alignment.sorted.withRG.dedup.converted.bq.sorted.without.splicing.junction.SN.bam"
    output:
        flag_A02_5_patch01=touch("result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/finished.patch01__intron"),
        alignment_merged_depth_on_transcribable_Alu_A_intron_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.intron.bed.txt.gz"
    threads:
        1
    shell:
        """
        samtools depth -d 1 -b {input.all_contigs_Alu_only_snpEff_annotated_annotated_records_only_intron_combined_bed_filename} {input.alignment_sorted_withRG_dedup_converted_bq_sorted_without_splicing_junction_SN_bam_filename} | gzip > {output.alignment_merged_depth_on_transcribable_Alu_A_intron_bed_txt_gz_filename}
        """

rule BA02_5__check_transcribable_Alu_A_coverage_of_merged_bam____patch01__intron:
    input:
        lambda wildcards: ["result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME  + "/" + str(temp_row.INDEXER_PARAMETERS) +  "/finished.patch01__intron" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    output:
        flag_B52_1=touch("result/BA02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished.patch01__intron")


rule A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam:
    input:
        flag_A02_5="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/finished",
        ## alignment_merged_depth_on_transcribable_Alu_A_intron_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.intron.bed.txt.gz",
        alignment_merged_depth_on_transcribable_Alu_A_5_prime_UTR_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.5.prime.UTR.bed.txt.gz",
        alignment_merged_depth_on_transcribable_Alu_A_3_prime_UTR_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.3.prime.UTR.bed.txt.gz",
        alignment_merged_depth_on_transcribable_Alu_A_CDS_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.CDS.bed.txt.gz"
    output:
        flag_A02_6=touch("result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/finished"),
        count_5_prime_UTR_txt_filename="result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/count.5.prime.UTR.txt",
        count_3_prime_UTR_txt_filename="result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/count.3.prime.UTR.txt",
        count_CDS_txt_filename="result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/count.CDS.txt"
    threads:
        1
    shell:
        """
        zcat {input.alignment_merged_depth_on_transcribable_Alu_A_5_prime_UTR_bed_txt_gz_filename} | awk '{{if ($3 > 0) print}}' | wc -l > {output.count_5_prime_UTR_txt_filename}
        zcat {input.alignment_merged_depth_on_transcribable_Alu_A_3_prime_UTR_bed_txt_gz_filename} | awk '{{if ($3 > 0) print}}' | wc -l > {output.count_3_prime_UTR_txt_filename}
        zcat {input.alignment_merged_depth_on_transcribable_Alu_A_CDS_bed_txt_gz_filename} | awk '{{if ($3 > 0) print}}' | wc -l > {output.count_CDS_txt_filename}
        """

rule BA02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam:
    input:
        flags_collection=lambda wildcards: ["result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME  + "/" + str(temp_row.INDEXER_PARAMETERS) +  "/finished" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    output:
        flag_BA02_5=touch("result/BA02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished")



rule A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam____patch01__intron:
    input:
        flag_A02_5="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/finished",
        alignment_merged_depth_on_transcribable_Alu_A_intron_bed_txt_gz_filename="result/A02_5__check_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/alignment.merged.depth.on.transcribable.Alu.A.intron.bed.txt.gz"
    output:
        flag_A02_6__patch01=touch("result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/finished.patch01__intron"),
        count_intron_txt_filename="result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{TYPE}/{DATASET}/{SAMPLE}/{INDEXER_PARAMETERS}/count.intron.txt"
    threads:
        1
    shell:
        """
        zcat {input.alignment_merged_depth_on_transcribable_Alu_A_intron_bed_txt_gz_filename} | awk '{{if ($3 > 0) print}}' | wc -l > {output.count_intron_txt_filename}
        """

rule BA02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam____patch01__intron:
    input:
        flags_collection=lambda wildcards: ["result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME  + "/" + str(temp_row.INDEXER_PARAMETERS) +  "/finished.patch01__intron" for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    output:
        flag_BA02_6__patch01=touch("result/BA02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished.patch01__intron")


rule A02_7__combine_summaries_of_distribution_of_transcribable_Alu_A_coverage_of_merged_bam:
    input:
        flag_BA02_6="result/BA02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        flag_BA02_6__patch01="result/BA02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished.patch01__intron"
    output:
        flag_A02_7=touch("result/A02_7__combine_summaries_of_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished"),
        final_count_csv_filename="result/A02_7__combine_summaries_of_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/final.count.csv"
    params:
        lambda wildcards: [ [temp_row.DATASET_NAME, temp_row.SAMPLE_NAME] + ["result/A02_6__summarize_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/" + wildcards.DATASET_COLLECTION_NAME + "/" + wildcards.DATASET_PHENOTYPE_COLLECTION_NAME + "/" + temp_row.TYPE + "/" + temp_row.DATASET_NAME + '/' + temp_row.SAMPLE_NAME  + "/" + str(temp_row.INDEXER_PARAMETERS) +  "/count." + Annotation_subset + ".txt" for Annotation_subset in ["5.prime.UTR", "3.prime.UTR", "CDS", "intron"]] for temp_row in __get_parameters_collection_for_RNA_editing_calling(wildcards.DATASET_COLLECTION_NAME)]
    run:
        dataset_and_sample_and_counts_collection=params[0]
        all_records_list = []
        used_samples = []
        for dataset, sample, count_5_prime_UTR_filename, count_3_prime_UTR_filename, count_CDS_filename, count_intron_filename in dataset_and_sample_and_counts_collection:
            if sample in used_samples:
                print("The sample " + sample + " has been used")
                continue
            used_samples.append(sample)
            temp_counts_list = []
            for temp_filename in [count_5_prime_UTR_filename, count_3_prime_UTR_filename, count_CDS_filename, count_intron_filename]:
                temp_count = "NA"
                with open(temp_filename, "r") as temp_f:
                    temp_count=temp_f.readline().strip()
                temp_counts_list.append(temp_count)
            all_records_list.append([dataset, sample] + temp_counts_list)
        final_pd = pandas.DataFrame(all_records_list)
        final_pd.columns = ["dataset", "name", "count.5.prime.UTR", "count.3.prime.UTR", "count.CDS", "count.intron"]
        final_pd.to_csv(output.final_count_csv_filename)

