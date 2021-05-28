#!/usr/bin/env Rscript

# USAGE: render.R notebook.Rmd 
args = commandArgs(trailingOnly=TRUE)

# work around rmarkdown bug: https://github.com/rstudio/rmarkdown/issues/1508
rmdfile = args[1]
rmdfile_orig = paste0(rmdfile, ".orig")
system(sprintf("mv '%s' '%s'", rmdfile, rmdfile_orig)) 
system(sprintf("cp -L '%s' '%s'", rmdfile_orig, rmdfile))

nxfvars = list(nxfvars = yaml::read_yaml('.params.yml'))
rmarkdown::render(rmdfile, params = nxfvars)
