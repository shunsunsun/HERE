library("rmarkdown")

render(input="scripts.for.report/SF2/run_internal.R", output_format="html_document", output_file="SF2.summary.html", output_dir=snakemake@params[["result_directory"]], knit_root_dir=getwd())
