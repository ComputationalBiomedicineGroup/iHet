# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python [conda env:conda-2021-nsclc-heterogeneity-scanpy2]
#     language: python
#     name: conda-env-conda-2021-nsclc-heterogeneity-scanpy2-py
# ---

# %%
# %load_ext autoreload
# %autoreload 2 

# %%
import scanpy as sc
from nxfvars import nxfvars
import pandas as pd
from threadpoolctl import threadpool_limits
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from matplotlib.backends.backend_pdf import PdfPages

sc.settings.set_figure_params(figsize=(4, 4))
import scanpy_helpers as sh

# %%
path_adata_m = nxfvars.get(
    "adata_m", "../../data/results/20_single-cell/annotate_myeloid/adata_myeloid_reannotated.h5ad"
)
path_adata_nsclc = nxfvars.get(
    "adata_nsclc", "../../data/results/20_single-cell/annotate_myeloid/adata_nsclc_reannotated.h5ad"
)
cpus = nxfvars.get("cpus", 16)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/ihet")

# %%
threadpool_limits(cpus)

# %%
adata_m = sc.read_h5ad(path_adata_m)
adata_nsclc = sc.read_h5ad(path_adata_nsclc)

# %%
adatas = {"m": {}, "nsclc": {}}

# %%
adatas["m"]["progeny"] = sh.compare_groups.compute_scores.run_progeny(adata_m)
adatas["m"]["dorothea"] = sh.compare_groups.compute_scores.run_dorothea(adata_m)
adatas["nsclc"]["progeny"] = sh.compare_groups.compute_scores.run_progeny(adata_nsclc)
adatas["nsclc"]["dorothea"] = sh.compare_groups.compute_scores.run_dorothea(adata_nsclc)

# %%
sc.pl.umap(adatas["nsclc"]["dorothea"], cmap="coolwarm", color=["E2F1", "E2F2", "E2F3", "E2F4"], vmin=-2, vmax=2)

# %%
tfs_of_interest = [
    "SPI1", "NFKB1", "STAT1", "MYCN", "E2F3", "E2F2", "TFDP1", "ZNF263", "MYC", "E2F4"    
]
pws_of_interest = [
    "TNFa", "NFkB", "Trail", "TGFb", "WNT", "Androgen"
]

# %%
pb_progeny = sh.pseudobulk.pseudobulk(
    adatas["nsclc"]["progeny"], groupby=["patient", "cell_type"], aggr_fun=np.mean
)

# %% [markdown]
# Note: another round of z-scoreing after pseudobulking hardly makes any difference. 

# %%
fig  = sc.pl.matrixplot(pb_progeny, var_names=pws_of_interest, groupby="cell_type", cmap="coolwarm", swap_axes=True, vmin=-2.5, vmax=2.5, return_fig=True)
fig.savefig(f"{artifact_dir}/heatmap_progeny.svg")

# %%
dorothea = sh.pseudobulk.pseudobulk(
    adatas["nsclc"]["dorothea"], groupby=["patient", "cell_type"], aggr_fun=np.mean
)

# %%
fig = sc.pl.matrixplot(pb_dorothea, var_names=tfs_of_interest, groupby="cell_type", cmap="coolwarm", swap_axes=True, vmin=-2.5, vmax=2.5, return_fig=True)
fig.savefig(f"{artifact_dir}/heatmap_dorothea.svg")
