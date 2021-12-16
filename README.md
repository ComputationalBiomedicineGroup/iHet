# Immune Heterogeneity (iHet) in NSCLC

This repository holds the code for all analyses related to

> t.b.a

All analyses are integrated into a nextflow pipeline and all dependencies are packaged as singularity containers.
The pipeline consists of the following subworkflows

 * **MOFA**. Runs unsupervised multi-omics factorial analysis (MOFA) of bulk RNA-seq data from Jia, Sharma, TracerX and pan-cancer datasets from TCGA
 * **single-cell**. TODO

## Usage

To reproduce the analyses, run

```
#TODO
nextflow run main.nf
```

## Input data

### Bulk RNA-seq
 * Gene expression data for each dataset
 * Tumor mutational burden data for each dataset

### scRNA-seq 
 * TODO


## Output documentation

The analysis pipeline generates the following directory structure:

```
10_mofa
  11_easier
  12_prepare_mofa_data
  13_run_mofa
  14_mofa_analysis
20_annotate_cell_types
  ... # TODO
```

In this section, we describe the directories and their contents in more detail.

### 11_easier

Results of the [easier](https://github.com/olapuentesantana/easier) package to obtain
cell-type fraction, pathway- and transcription factor estimates. Contains one
`.rds` file per dataset group, which contains a list of lists of data frames.
There is one data frame for each dataset and each modality. 
Contained modalities: 

```
"count"    "tpm"      "response" "pathway"  "tf"       "cellfrac" "immresp"  
```

### 12_prepare_mofa_data

Result of an Rmarkdown notebook preparing the input data into a format compatible
with MOFA. Splits up the data into all individual dataset and creates bootstrap datasets
for each dataset. A rendered version of the notebook is in the main directory, all
generated files are in the `artifacts` directory:  

* `data_all_tidy.rds`: All modalities and datasets from 11_easier merged into a single, "tidy" data frame. Some features are renamed, some are removed and some are merged.  
* `mofa_*.rds`: The tidy data from above split up by dataset, with some additional scaling. 
* `mofa_boot_*.rds`: Same structure as `mofa_*.rds`, but with datasets resampled by bootstrapping. Each file represents a different boostrap. 

### 13_run_mofa

Results of MOFA after applying it to the datasets generated in the previous step. 
Contains `hdf5` datasets for each dataset holding the mofa models. 

### 14_mofa_analysis

Results of an Rmarkdown notebook with the analysis of the MOFA results. A rendered version of the 
notebook is in the main directory, all generated filesl are in the `artifacts` directory: 

* `median_factors.rds`: Median factors across all boostraps for each dataset.
  Contains a list of dataframes (one for each dataset). 
* `median_weights.rds`: Median weights (for each feature) across all bootstraps for each dataset. Contains a 
  list of data frames (one for each dataset)
* `*_factor_correlations{,_pvalues}.tsv`: Pearson correlations between F1 and F1-F3 between the
  different datasets and the associated p-values. 
* `plots/`: Various plots


