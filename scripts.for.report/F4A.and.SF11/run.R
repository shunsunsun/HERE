library("rmarkdown")

render(input="scripts.for.report/F4A.and.SF11/run_internal.R", output_format="html_document", output_file="F4A.and.SF11.summary.html", output_dir=snakemake@params[["result_directory"]], knit_root_dir=getwd())
