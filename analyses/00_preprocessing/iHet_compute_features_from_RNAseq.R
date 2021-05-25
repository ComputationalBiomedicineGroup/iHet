rm(list=ls())

source("iHet_functions.R")

cancertype <- "NSCLC"

## NSCLC data

load("../RData/NSCLC_expr_data_sel.RData")

NSCLCobj.sel.easier <- getEASIERinfo(NSCLCobj.sel, 
                                     cancertype = cancertype)

save(NSCLCobj.sel.easier, 
     file="../RData/NSCLC_easier_sel.RData")


## GTEx data

load("../RData/GTEx_expr_data.RData")

GTExobj.easier <- getEASIERinfo(GTExobj, 
                                 cancertype = cancertype)

save(GTExobj.easier, 
     file="../RData/GTEx_easier.RData")



