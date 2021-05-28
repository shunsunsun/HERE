library("rmarkdown")

render(input="scripts.for.report/SF13B/run_internal.R", output_format="html_document", output_file="SF13B.summary.html", output_dir=snakemake@params[["result_directory"]], knit_root_dir=getwd())
