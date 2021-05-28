library("rmarkdown")

render(input="scripts.for.report/F4B.and.SF12/run_internal.R", output_format="html_document", output_file="F4B.summary.html", output_dir=snakemake@params[["result_directory"]], knit_root_dir=getwd())
