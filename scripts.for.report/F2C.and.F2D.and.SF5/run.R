library("rmarkdown")

render(input="scripts.for.report/F2C.and.F2D.and.SF5/run_internal.R", output_format="html_document", output_file="F2C.and.F2D.and.SF5.summary.html", output_dir=snakemake@params[["result_directory"]], knit_root_dir=getwd())
