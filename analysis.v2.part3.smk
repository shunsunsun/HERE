if 'large_data_threads' not in config.keys():
    config['large_data_threads'] = 10

rule SF1A_and_SF1B:
    input:
        flag_S71_5="result/S71_5__filter_for_A_to_G_sites_for_control/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_with_event_summary_dt_txt_gz_filename="result/S71_5__filter_for_A_to_G_sites_for_control/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt.txt.gz"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF1A.and.SF1B.finished"),
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF1A.and.SF1B.summary.html",
        SF1A_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF1A.png",
        SF1B_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF1B.png"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    threads:
        1
    script:
        "./scripts.for.report/SF1A.and.SF1B/run.R"


rule SF1D:
    input:
        DNA_flag_B15_1__step02__part06="result/B15_1__get_sample_RNA_editing_sites_v3/{DNA_DATASET_COLLECTION_NAME}/__merged__/DNTRSeq-DNA-trimming/hg38.fa/32/bwa-index-default/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/finished.step02__call_variants____part06__really_call_variants",
        RNA_flag_S71_5="result/S71_5__filter_for_A_to_G_sites_for_control/{RNA_DATASET_COLLECTION_NAME}/{RNA_DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        RNA_merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_txt_gz_filename="result/S71_5__filter_for_A_to_G_sites_for_control/{RNA_DATASET_COLLECTION_NAME}/{RNA_DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz"
    output:
        flag=touch("report/{RNA_DATASET_COLLECTION_NAME}/{RNA_DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{DNA_DATASET_COLLECTION_NAME}/SF1D.finished"),
        raw_RNA_DNA_comparison_dt_txt_gz_filename="report/{RNA_DATASET_COLLECTION_NAME}/{RNA_DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{DNA_DATASET_COLLECTION_NAME}/raw.RNA.DNA.comparison.dt.txt.gz",
        RNA_DNA_comparison_dt_txt_gz_filename="report/{RNA_DATASET_COLLECTION_NAME}/{RNA_DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{DNA_DATASET_COLLECTION_NAME}/RNA.DNA.comparison.dt.txt.gz",
        summary_html_filename="report/{RNA_DATASET_COLLECTION_NAME}/{RNA_DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{DNA_DATASET_COLLECTION_NAME}/SF1D.summary.html",
        main_ggplot_RDS_filename="report/{RNA_DATASET_COLLECTION_NAME}/{RNA_DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{DNA_DATASET_COLLECTION_NAME}/SF1D.ggplot.RDS",
        main_png_filename="report/{RNA_DATASET_COLLECTION_NAME}/{RNA_DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{DNA_DATASET_COLLECTION_NAME}/SF1D.png"
    params:
        result_directory="report/{RNA_DATASET_COLLECTION_NAME}/{RNA_DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/{DNA_DATASET_COLLECTION_NAME}/"
    threads:
        1
    script:
        "./scripts.for.report/SF1D/run.R"


rule F1B:
    input:
        flag_S21_1="result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        phenotype_output_at_gsm_level_dt_filename="result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/phenotype.output.at.gsm.level.dt.txt",
        flag_S51_5="result/S51_5__filter_for_A_to_G_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_txt_gz_filename="result/S51_5__filter_for_A_to_G_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_with_event_summary_dt_txt_gz_filename="result/S51_5__filter_for_A_to_G_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.with.event.summary.dt.txt.gz"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F1B.finished"),
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F1B.summary.html",
        main_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F1B.ggplot.RDS",
        main_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F1B.png"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    threads:
        1
    script:
        "./scripts.for.report/F1B/run.R"



rule F1D:
    input:
        flag_S51_5="result/S51_5__filter_for_A_to_G_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_dt_txt_gz_filename="result/S51_5__filter_for_A_to_G_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.dt.txt.gz"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F1D.finished"),
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F1D.summary.html",
        main_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F1D.ggplot.RDS",
        main_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F1D.png"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    threads:
        1
    script:
        "./scripts.for.report/F1D/run.R"


rule SF2:
    input:
        flag_B15_1__step03="result/B15_1__get_sample_RNA_editing_sites_v3/{DATASET_COLLECTION_NAME}/__merged__/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/bwa-index-10.1038_nmeth.2330/bwa-aln-samsepe/none/GATK-3.6.0/none/151/common_all/recal.bam.subset/BLCAP.Y2C.u10.d10/finished.step03__extract_bam_subset",
        flag_S21_1="result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        phenotype_output_at_gsm_level_dt_filename="result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/phenotype.output.at.gsm.level.dt.txt"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF2.finished"),
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF2.summary.html",
        main_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF2.png"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    threads:
        1
    script:
        "./scripts.for.report/SF2/run.R"


        
        
rule F2A_and_SF3:
    input:
        subset_site_recurrence_comparison_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.site.recurrence.comparison.dt.txt.gz",
        subset_recurrent_edits_only_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.recurrent.edits.only.dt.txt.gz"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F2A.and.SF3.finished"),
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F2A.and.SF3.summary.html",
        ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F2A.ggplot.RDS",
        png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F2A.png",
        supp1_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF3.ggplot.RDS",
        supp1_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF3.png"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    script:
        "./scripts.for.report/F2A.and.SF3/run.R"


rule F2B_and_SF4:
    input:
        subset_site_recurrence_comparison_CJ_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.site.recurrence.comparison.CJ.dt.txt.gz"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F2B.and.SF4.finished"),
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F2B.and.SF4.summary.html",
        main_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F2B.ggplot.RDS",
        main_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F2B.png",
        supp1_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF4.ggplot.RDS",
        supp1_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF4.png"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    script:
        "./scripts.for.report/F2B.and.SF4/run.R"


rule F2C_and_F2D_and_SF5:
    input:
        subset_recurrent_edits_only_with_snpEff_annotation_on_valid_genes_only_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz",
        flag_A02_7="result/A02_7__combine_summaries_of_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        final_count_csv_filename="result/A02_7__combine_summaries_of_distribution_of_transcribable_Alu_A_coverage_of_merged_bam/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/final.count.csv"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F2C.and.F2D.and.SF5.finished"),
        sample_stats_dt_txt_gz_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F2C.and.F2D.and.SF5.sample.stats.dt.txt.gz",
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F2C.and.F2D.and.SF5.summary.html",
        main_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F2C.and.F2D.ggplot.RDS",
        main_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F2C.and.F2D.png",
        supp1_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF5.ggplot.RDS",
        supp1_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF5.png"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    threads:
        1
    script:
        "./scripts.for.report/F2C.and.F2D.and.SF5/run.R"


rule SF6:
    input:
        subset_recurrent_edits_only_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.recurrent.edits.only.dt.txt.gz",
        REDIPortal_2021_hg38_editome_txt_gz_filename="external/REDIPortal/hg38/TABLE1_hg38.txt.gz"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF6.finished"),
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF6.summary.html",
        main_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF6.ggplot.RDS",
        main_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF6.png"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    script:
        "./scripts.for.report/SF6/run.R"



rule F3A_and_SF7:
    input:
        subset_recurrent_edits_only_with_snpEff_annotation_on_valid_genes_only_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3A.and.SF7.finished"),
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3A.and.SF7.summary.html",
        main_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3A.ggplot.RDS",
        main_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3A.png",
        supp1_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF7.ggplot.RDS",
        supp1_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF7.png"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    script:
        "./scripts.for.report/F3A.and.SF7/run.R"


rule F3B_and_SF8:
    input:
        subset_recurrent_edits_only_with_snpEff_annotation_on_valid_genes_only_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3B.and.SF8.finished"),
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3B.and.SF8.summary.html",
        for_plot_dt_txt_gz_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3B.for.plot.dt.txt.gz",
        main_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3B.ggplot.RDS",
        main_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3B.png",
        supp1_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF8.ggplot.RDS",
        supp1_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF8.png"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    threads:
        1
    script:
        "./scripts.for.report/F3B.and.SF8/run.R"


rule F3C_and_SF9:
    input:
        F3B_and_SF8_flag="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3B.and.SF8.finished",
        F3B_for_plot_dt_txt_gz_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3B.for.plot.dt.txt.gz",
        flag_S42_1="result/S42_1__annotate_embryonic_genes/{DATASET_COLLECTION_NAME}/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        combined_gexpr_FPKM_pc_only_melt_with_phenotype_normal_sample_only_median_annotated_dt_txt_gz_filename="result/S42_1__annotate_embryonic_genes/{DATASET_COLLECTION_NAME}/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/{DATASET_PHENOTYPE_COLLECTION_NAME}/combined.gexpr.FPKM.pc.only.melt.with.phenotype.normal.sample.only.median.annotated.dt.txt.gz"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3C.and.SF9.finished"),
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3C.and.SF9.summary.html",
        main_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3C.ggplot.RDS",
        main_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3C.png",
        supp1_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF9.ggplot.RDS",
        supp1_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF9.png"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    threads:
        1
    script:
        "./scripts.for.report/F3C.and.SF9/run.R"


rule F3D_and_SF10:
    input:
        F3B_and_SF8_flag="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3B.and.SF8.finished",
        F3B_for_plot_dt_txt_gz_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3B.for.plot.dt.txt.gz"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3D.and.SF10.finished"),
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3D.and.SF10.summary.html",
        GO_result_dt_txt_gz_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3D.GO.result.dt.txt.gz",
        main_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3D.ggplot.RDS",
        main_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F3D.png",
        supp1_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF10.ggplot.RDS",
        supp1_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF10.png"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    threads:
        1
    script:
        "./scripts.for.report/F3D.and.SF10/run.R"


rule F4B_and_SF12:
    input:
        subset_recurrent_edits_only_with_snpEff_annotation_on_valid_genes_only_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz",
        subset_site_recurrence_comparison_CJ_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.site.recurrence.comparison.CJ.dt.txt.gz",
        subset_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.dt.txt.gz"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F4B.and.SF12.finished"),
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F4B.summary.html",
        F4B_for_plot_dt_txt_gz_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F4B.for.plot.dt.txt.gz",
        F4B_GO_result_dt_txt_gz_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F4B.GO.result.dt.txt.gz",
        main_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F4B.ggplot.RDS",
        main_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F4B.png",
        supp1_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF12.ggplot.RDS",
        supp1_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF12.png"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    script:
        "./scripts.for.report/F4B.and.SF12/run.R"


rule F4A_and_SF11:
    input:
        F4B_for_plot_dt_txt_gz_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F4B.for.plot.dt.txt.gz",
        flag_S42_2="result/S42_2__annotate_uniparental_embryonic_genes/{DATASET_COLLECTION_NAME}/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/201221-fifth-phenotype-collection/finished",
        combined_gexpr_FPKM_pc_only_melt_with_phenotype_GSE133854_sample_only_median_annotated_dt_txt_gz_filename="result/S42_2__annotate_uniparental_embryonic_genes/{DATASET_COLLECTION_NAME}/auto-detect-and-cut-adapter-by-trim-galore-and-select-reads-with-base-quality-no-smaller-than-25-by-fastp/hg38.fa/32/STAR-expression/__sample_dependent__/STAR-expression/default/stringtie/none/201221-fifth-phenotype-collection/combined.gexpr.FPKM.pc.only.melt.with.phenotype.GSE133854.sample.only.median.annotated.dt.txt.gz"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F4A.and.SF11.finished"),
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F4A.and.SF11.summary.html",
        main_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F4A.ggplot.RDS",
        main_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F4A.png",
        supp1_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF11.ggplot.RDS",
        supp1_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF11.png"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    threads:
        1
    script:
        "./scripts.for.report/F4A.and.SF11/run.R"



rule F4C:
    input:
        F4B_GO_result_dt_txt_gz_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F4B.GO.result.dt.txt.gz"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F4C.finished"),
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F4C.summary.html",
        main_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F4C.ggplot.RDS",
        main_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F4C.png"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    threads:
        1
    script:
        "./scripts.for.report/F4C/run.R"


rule F5ABDE_and_SF13_SF14_SF15_SF16_ST4:
    input:
        ALS_subset_recurrent_edits_only_with_snpEff_annotation_on_valid_genes_only_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/all.normal.samples/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz",
        subset_recurrent_edits_only_with_snpEff_annotation_on_valid_genes_only_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.recurrent.edits.only.with.snpEff.annotation.on.valid.genes.only.dt.txt.gz",
        subset_site_recurrence_comparison_CJ_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.site.recurrence.comparison.CJ.dt.txt.gz",
        subset_dt_txt_gz_filename="result/A02_3__check_recurrence_profile_for_a_subset_of_samples/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/subset.dt.txt.gz",
        phenotype_output_at_gsm_level_dt_filename="result/S21_1__merge_phenotype_tables/{DATASET_PHENOTYPE_COLLECTION_NAME}/phenotype.output.at.gsm.level.dt.txt"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F5ABDE.and.SF13.SF14.SF15.SF16.ST4.finished"),
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F5ABDE.and.SF13.SF14.SF15.SF16.ST4.summary.html",
        F5A_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F5A.png",
        F5B_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F5B.png",
        F5D_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F5D.png",
        F5E_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/F5E.png",
        SF13_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF13.png",
        SF14_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF14.png",
        SF15_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF15.png",
        SF16_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF16.png",
        ST4_csv_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/ST4.csv"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    threads:
        1
    script:
        "./scripts.for.report/F5ABDE.and.SF13.SF14.SF15.SF16.ST4/run.R"



rule SF18A:
    input:
        flag_S52_3="result/S52_3__mark_unsequenced_editing_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_merged_with_coverage_dt_txt_gz_filename="result/S52_3__mark_unsequenced_editing_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt.txt.gz"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF18A.finished"),
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF18A.summary.html",
        edit_and_sample_unsequenced_status_dt_txt_gz_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/edit.and.sample.unsequenced.status.dt.txt.gz",
        edit_and_sample_unsequenced_status_edit_histogram_per_sample_group_dt_txt_gz_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/edit.and.sample.unsequenced.status.edit.histogram.per.sample.group.dt.txt.gz",
        main_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF18A.ggplot.RDS",
        main_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF18A.png"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    threads:
        config["large_data_threads"]
    script:
        "./scripts.for.report/SF18A/run.R"


rule SF18B:
    input:
        flag_S52_3="result/S52_3__mark_unsequenced_editing_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/finished",
        merged_long_disjoint_with_population_without_potential_polymorphism_with_enough_read_support_with_phenotype_sequenced_samples_only_with_enough_sample_support_A_to_G_only_merged_with_coverage_dt_txt_gz_filename="result/S52_3__mark_unsequenced_editing_sites/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/merged.long.disjoint.with.population.without.potential.polymorphism.with.enough.read.support.with.phenotype.sequenced.samples.only.with.enough.sample.support.A.to.G.only.merged.with.coverage.dt.txt.gz"
    output:
        flag=touch("report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF18B.finished"),
        summary_html_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF18B.summary.html",
        site_recurrence_comparison_dt_txt_gz_filename="result/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/site.recurrence.comparison.dt.txt.gz",
        main_ggplot_RDS_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF18B.ggplot.RDS",
        main_png_filename="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/SF18B.png"
    params:
        result_directory="report/{DATASET_COLLECTION_NAME}/{DATASET_PHENOTYPE_COLLECTION_NAME}/{SUBSET_NAME}/"
    threads:
        config["large_data_threads"]
    script:
        "./scripts.for.report/SF18B/run.R"
