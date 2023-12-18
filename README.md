# Immune Heterogeneity (iHet) in NSCLC

This repository holds the code for all analyses related to

> t.b.a

All analyses are integrated into a nextflow pipeline and all dependencies are packaged as singularity containers.
The pipeline consists of the following subworkflows

 * **MOFA**. Run unsupervised multi-omics factorial analysis (MOFA) of bulk RNA-seq data from Jia, Sharma, TRACERx and 
   pan-cancer datasets from TCGA
 * **single-cell**. Create a custom subset of the [Lung Cancer Atlas](https://doi.org/10.1016/j.ccell.2022.10.008), 
   re-annotate myeloid subtypes and perform transcription factor and pathway analysis. 

## Launching the workflow

### Prerequisites
* [Nextflow](https://www.nextflow.io/index.html#GetStarted), version 22.04.5
* [Singularity/Apptainer](https://apptainer.org/), version 3.7 or higher (tested with 3.7.0-1.el7)
* A machine with at least 200GB of memory

### Obtaining data

Before launching the workflow, you need to obtain input data and singularity containers from zenodo.
First of all, clone this repository:

```bash
git clone https://github.com/ComputationalBiomedicineGroup/iHet.git
cd iHet
```

Then, within the repository, download the data archives and extract then to the corresponding directories:

```bash
# singularity containers
curl TODO

# input data
curl TODO
```

Additionally, you can obtain the pre-computed results without running the workflow using
```bash
curl TODO
```

Note that some results depend on the [TRACERx data](https://pubmed.ncbi.nlm.nih.gov/31591602/) (EGAS00001003458, EGAD00001003206)
which is not publicly available. The workflow is configured, by default, to run without these data.

Briefly, the input data contains
 * Gene expression data for each dataset
 * Tumor mutational burden data for each dataset
 * Single-cell data ([LuCA version 2022.05.10](https://doi.org/10.5281/zenodo.6411868))

### Run nextflow

```bash
# newer versions of nextflow are incompatible with the workflow. By setting this variable
# the correct version will be used automatically.
export NXF_VER=22.04.5

nextflow run main.nf --outdir data/results
```

## Structure of this repository

* `analyses`: Place for e.g. jupyter/rmarkdown notebooks, gropued by their respective (sub-)workflows.
* `containers`: place for singularity image files. Not part of the git repo and gets created by the download command.
* `data`: place for input data and results in different subfolders. Gets populated by the download commands and by running the workflows.
* `lib`: custom libraries and helper functions
* `modules`: nextflow DSL2.0 modules
* `subworkflows`: nextflow subworkflows
* `tables`: contains static content that should be under version control (e.g. manually created tables)


## Output documentation

The analysis pipeline generates the following directory structure:

```
10_mofa
  11_easier
  12_prepare_mofa_data
  13_run_mofa
  14_mofa_analysis
  15_iHet_predictions
20_annotate_cell_types
  21_subset_atlas
  22_annotate_myeloid
  23_tf_pw
  29_overview_plots
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

### 15_iHet_predictions

TODO

### 21_subset_atlas 

Subset the full LuCA atlas to only contain primary tumor samples from either LUAD or LUSC. Creates two h5ad files:
 * `adata_m.h5ad`: Subset of myeloid cell-types
 * `adata_nsclc.h5ad`: The full custom subset of the atlas. 

### 22_annotate_myeloid

Re-annotate myeloid cell-types at a better resolution than LuCA. Generates updated h5ad files:
 * `adata_nsclc_reannotated.h5ad`
 * `adata_myeloid_reannotated.h5ad`

### 23_tf_pw

Execute Dorothea and Progeny on the single-cell data. Generates heatmaps and summary tables. 

### 29_overview_plots

Generate single-cell related figures for the manuscript. 

## Contact

For reproducibility issues or any other requests regarding single-cell data analysis, please use the [issue tracker](https://github.com/ComputationalBiomedicineGroup/iHet/issues). For anything else, you can reach out to the corresponding author(s) as indicated in the manuscript.

