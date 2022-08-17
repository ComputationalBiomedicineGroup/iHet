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
# %%capture
for scope, tmp_adatas in adatas.items():
    for tool, tmp_ad in tmp_adatas.items():
        with PdfPages(f"{artifact_dir}/{tool}_{scope}_umap.pdf") as pdf:
            for var in tmp_ad.var_names:
                with plt.rc_context({"figure.figsize": (6, 6), "figure.dpi": 300}):
                    sc.pl.umap(
                        tmp_ad,
                        color=var,
                        cmap="coolwarm",
                        vmax=3,
                        vmin=-3,
                        show=False,
                        size=500000 / tmp_ad.shape[0],
                    )
                    pdf.savefig(bbox_inches="tight")

# %%
# %%capture
for scope, tmp_adatas in adatas.items():
    for tool, tmp_ad in tmp_adatas.items():
        tmp_pb = sh.pseudobulk.pseudobulk(
            tmp_ad, groupby=["patient", "cell_type"], aggr_fun=np.mean
        )
        tmp_pb_z = tmp_pb.copy()
        tmp_pb_z.X = scipy.stats.zscore(tmp_pb_z.X, axis=0)
        fig = sc.pl.matrixplot(
            tmp_pb_z,
            groupby="cell_type",
            var_names=tmp_pb.var_names,
            cmap="coolwarm",
            swap_axes=True,
            dendrogram=True,
            return_fig=True,
            show=False,
            vmin=-2.5,
            vmax=2.5
        )
        fig.savefig(
            f"{artifact_dir}/{tool}_{scope}_pseudobulk_heatmap_clustered_zscore.pdf", bbox_inches="tight"
        )
        fig = sc.pl.matrixplot(
            tmp_pb,
            groupby="cell_type",
            var_names=tmp_pb.var_names,
            cmap="coolwarm",
            swap_axes=True,
            dendrogram=True,
            return_fig=True,
            show=False,
        )
        fig.savefig(
            f"{artifact_dir}/{tool}_{scope}_pseudobulk_heatmap_clustered.pdf", bbox_inches="tight"
        )
        fig = sc.pl.matrixplot(
            tmp_pb,
            groupby="cell_type",
            var_names=tmp_pb.var_names,
            cmap="coolwarm",
            swap_axes=True,
            dendrogram=False,
            return_fig=True,
            show=False,
        )
        fig.savefig(f"{artifact_dir}/{tool}_{scope}_pseudobulk_heatmap.pdf", bbox_inches="tight")

# %%
