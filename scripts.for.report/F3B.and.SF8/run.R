library("rmarkdown")

render(input="scripts.for.report/F3B.and.SF8/run_internal.R", output_format="html_document", output_file="F3B.and.SF8.summary.html", output_dir=snakemake@params[["result_directory"]], knit_root_dir=getwd())
