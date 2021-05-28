library("rmarkdown")

render(input="scripts.for.report/F3C.and.SF9/run_internal.R", output_format="html_document", output_file="F3C.and.SF9.summary.html", output_dir=snakemake@params[["result_directory"]], knit_root_dir=getwd())
