library("rmarkdown")

render(input="scripts.for.report/F5ABDE.and.SF13.SF14.SF15.SF16.ST4/run_internal.R", output_format="html_document", output_file="F5ABDE.and.SF13.SF14.SF15.SF16.ST4.summary.html", output_dir=snakemake@params[["result_directory"]], knit_root_dir=getwd())
