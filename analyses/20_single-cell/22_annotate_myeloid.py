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
sc.settings.set_figure_params(figsize=(4, 4))

# %%
path_adata_m = nxfvars.get("adata_m", "../../data/results/20_single-cell/subset_atlas/adata_m.h5ad")

# %%
adata_m = sc.read_h5ad(path_adata_m)

# %% [markdown]
# # Re-annotate myeloid cells
#
# Create a more fine-graind annotaiton of the myeloid cells than in the atlas

# %%
sc.pl.umap(adata_m, color=["leiden", "cell_type"])

# %%
