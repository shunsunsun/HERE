library("rmarkdown")

render(input="scripts.for.report/F2B.and.SF4/run_internal.R", output_format="html_document", output_file="F2B.and.SF4.summary.html", output_dir=snakemake@params[["result_directory"]], knit_root_dir=getwd())
