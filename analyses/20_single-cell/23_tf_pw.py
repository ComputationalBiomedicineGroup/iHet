# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
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

sc.settings.set_figure_params(figsize=(4, 4))
import scanpy_helpers as sh

# %%
path_adata_m = nxfvars.get(
    "adata_m", "../../data/results/20_single-cell/subset_atlas/adata_m.h5ad"
)
path_adata_nsclc = nxfvars.get("adata_nsclc", "../../data/results/20_single-cell/subset_atlas/adata_nsclc.h5ad")
cpus = nxfvars.get("cpus", 16)

# %%
threadpool_limits(cpus)

# %%
adata_m = sc.read_h5ad(path_adata_m)
adata_nsclc = sc.read_h5ad(path_adata_nsclc)

# %%
adatas = {"m": {}, "nsclc": {}}

# %%
adatas["m"]["progeny"] = sh.compare_groups.compute_scores.run_progeny(adata_m)

# %%
sh.compare_groups.compute_scores.run_dorothea(adata_m)
