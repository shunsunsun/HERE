library("rmarkdown")

render(input="scripts.for.report/SF18A/run_internal.R", output_format="html_document", output_file="SF18A.summary.html", output_dir=snakemake@params[["result_directory"]], knit_root_dir=getwd())
