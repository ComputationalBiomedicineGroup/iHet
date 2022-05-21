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
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/ihet")

# %%
adata_m = sc.read_h5ad(path_adata_m)
adata_nsclc = sc.read_h5ad(path_adata_nsclc)

# %%
# %%capture
for scope, tmp_adata in {"m": adata_m, "nsclc": adata_nsclc}.items():
    size = 500000 / tmp_adata.shape[0]
    with PdfPages(f"{artifact_dir}/overview_{scope}_umap.pdf") as pdf:
        with plt.rc_context({"figure.figsize": (10, 10), "figure.dpi": 300}):
            sc.pl.umap(
                tmp_adata,
                color="cell_type_coarse",
                legend_loc="on data",
                legend_fontoutline=2,
                size=size,
                show=False,
            )
            pdf.savefig(bbox_inches="tight")

            sc.pl.umap(
                tmp_adata,
                color="cell_type",
                legend_loc="on data",
                legend_fontoutline=1,
                size=size,
                legend_fontsize=4,
                show=False,
            )
            pdf.savefig(bbox_inches="tight")

            for color in ["platform", "condition", "dataset"]:
                sc.pl.umap(tmp_adata, color=color, size=size, show=False)
                pdf.savefig(bbox_inches="tight")

# %%
