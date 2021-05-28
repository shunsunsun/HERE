library("rmarkdown")

render(input="scripts.for.report/F2A.and.SF3/run_internal.R", output_format="html_document", output_file="F2A.and.SF3.summary.html", output_dir=snakemake@params[["result_directory"]], knit_root_dir=getwd())
