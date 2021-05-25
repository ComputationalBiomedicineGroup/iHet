#!/usr/bin/bash

# USAGE: render.R notebook.Rmd report.html
args = commandArgs(trailingOnly=TRUE)
nxfvars = list(nxfvars = yaml::read_yaml('.params.yml'))
rmarkdown::render(args[1], params = nxfvars, output_file=args[2])
