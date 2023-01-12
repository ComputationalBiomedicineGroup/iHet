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

# %% [markdown]
# # Overview plots and statistics

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
from matplotlib.backends.backend_pdf import PdfPages

sc.settings.set_figure_params(figsize=(4, 4))
import scanpy_helpers as sh

# %%
path_adata_m = nxfvars.get(
    "adata_m",
    "../../data/results/20_single-cell/annotate_myeloid/adata_myeloid_reannotated.h5ad",
)
path_adata_nsclc = nxfvars.get(
    "adata_nsclc",
    "../../data/results/20_single-cell/annotate_myeloid/adata_nsclc_reannotated.h5ad",
)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/ihet/single-cell")

# %%
# !mkdir -p {artifact_dir}

# %%
adata_m = sc.read_h5ad(path_adata_m)
adata_nsclc = sc.read_h5ad(path_adata_nsclc)

# %% [markdown]
# ## Dataset stats

# %%
adata_nsclc.shape

# %%
adata_nsclc.obs["patient"].nunique() 

# %%
adata_nsclc.obs["study"].nunique()

# %%
adata_nsclc.obs["dataset"].nunique()

# %%
adata_nsclc.obs["cell_type_coarse"].nunique()

# %%
adata_nsclc.obs["cell_type_coarse"].unique().tolist()

# %%
adata_nsclc.obs["cell_type"].nunique()

# %% [markdown]
# ## UMAP plots

# %%
with plt.rc_context({"figure.figsize": (7, 7)}):
    sc.pl.umap(
        adata_nsclc,
        color="cell_type_tumor",
        legend_loc="on data",
        legend_fontoutline=2,
        size=500000 / adata_nsclc.shape[0],
    )

# %%
# %%capture
for scope, tmp_adata in {"m": adata_m, "nsclc": adata_nsclc}.items():
    size = 500000 / tmp_adata.shape[0]
    with plt.rc_context({"figure.figsize": (7, 7), "figure.dpi": 300}):
        fig = sc.pl.umap(
            tmp_adata,
            color="cell_type_coarse",
            legend_loc="on data",
            legend_fontoutline=2,
            size=size,
            show=False,
            frameon=False,
            return_fig=True,
        )
        fig.savefig(
            f"{artifact_dir}/umap_cell_type_coarse_{scope}.svg",
            bbox_inches="tight",
            dpi=1200,
        )

        fig = sc.pl.umap(
            tmp_adata,
            color="cell_type",
            legend_loc="on data",
            legend_fontoutline=0.5,
            size=size,
            legend_fontsize=7,
            show=False,
            frameon=False,
            return_fig=True,
        )
        fig.savefig(
            f"{artifact_dir}/umap_cell_type_fine_{scope}.svg",
            bbox_inches="tight",
            dpi=1200,
        )

        for color in ["platform", "condition", "dataset"]:
            fig = sc.pl.umap(
                tmp_adata,
                color=color,
                size=size,
                show=False,
                return_fig=True,
                frameon=False,
            )
            fig.savefig(
                f"{artifact_dir}/umap_{color}_{scope}.svg", bbox_inches="tight", dpi=600
            )

# %%
