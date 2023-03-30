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
import seaborn as sns

sc.settings.set_figure_params(figsize=(4, 4))
import scanpy_helpers as sh

# %%
path_adata_m = nxfvars.get(
    "adata_m",
    "../../data/results/20_single_cell/22_annotate_myeloid/artifacts/adata_myeloid_reannotated.h5ad",
)
path_adata_nsclc = nxfvars.get(
    "adata_nsclc",
    "../../data/results/20_single_cell/22_annotate_myeloid/artifacts/adata_nsclc_reannotated.h5ad",
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
with plt.rc_context({"figure.figsize": (6, 6)}):
    sc.pl.umap(
        adatas["nsclc"]["dorothea"],
        cmap="coolwarm",
        color=["E2F1", "E2F2", "E2F3", "E2F4"],
        vmin=-2,
        vmax=2,
        sort_order=False,
    )

# %%
tfs_of_interest = [
    "E2F2",
    "E2F3",
    "E2F4",
    "ETS1",
    "MYC",
    "MYCN",
    "NFKB1",
    "RELA",
    "RFX5",
    "SPI1",
    "STAT1",
    "STAT2",
    "TFDP1",
    "ZNF263",
]
pws_of_interest = ["Androgen", "JAK-STAT", "NFkB", "TGFb", "TNFa", "Trail", "WNT"]

# %%
pb_progeny = sh.pseudobulk.pseudobulk(
    adatas["nsclc"]["progeny"], groupby=["patient", "cell_type"], aggr_fun=np.mean
)
pb_progeny.obs["cell_type"] = pb_progeny.obs["cell_type"].astype(
    adatas["nsclc"]["progeny"].obs["cell_type"].dtype
)

# %% [markdown]
# Note: another round of z-scoreing after pseudobulking hardly makes any difference. 

# %%
fig = sc.pl.matrixplot(
    pb_progeny,
    var_names=pws_of_interest,
    groupby="cell_type",
    cmap="coolwarm",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/heatmap_progeny.svg")

# %%
pb_dorothea = sh.pseudobulk.pseudobulk(
    adatas["nsclc"]["dorothea"], groupby=["patient", "cell_type"], aggr_fun=np.mean
)
pb_dorothea.obs["cell_type"] = pb_dorothea.obs["cell_type"].astype(
    adatas["nsclc"]["dorothea"].obs["cell_type"].dtype
)

# %%
fig = sc.pl.matrixplot(
    pb_dorothea,
    var_names=tfs_of_interest,
    groupby="cell_type",
    cmap="coolwarm",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/heatmap_dorothea.svg")

# %% [markdown]
# ## Including myeloid subclusters

# %%
pb_progeny = sh.pseudobulk.pseudobulk(
    adatas["nsclc"]["progeny"], groupby=["patient", "cell_type_macro"], aggr_fun=np.mean
)
pb_progeny.obs["cell_type_macro"] = pb_progeny.obs["cell_type_macro"].astype(
    adatas["nsclc"]["progeny"].obs["cell_type_macro"].dtype
)

# %%
pb_dorothea = sh.pseudobulk.pseudobulk(
    adatas["nsclc"]["dorothea"],
    groupby=["patient", "cell_type_macro"],
    aggr_fun=np.mean,
)
pb_dorothea.obs["cell_type_macro"] = pb_dorothea.obs["cell_type_macro"].astype(
    adatas["nsclc"]["dorothea"].obs["cell_type_macro"].dtype
)

# %%
fig = sc.pl.matrixplot(
    pb_progeny[
        pb_progeny.obs["cell_type_macro"].isin(adata_m.obs["cell_type_macro"]), :
    ],
    var_names=pws_of_interest,
    groupby="cell_type_macro",
    cmap="coolwarm",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/heatmap_progeny_tam.svg")

# %%
pb_progeny.layers["zscore"] = scipy.stats.zscore(pb_progeny.X, axis=1)

# %%
fig = sc.pl.matrixplot(
    pb_progeny[
        pb_progeny.obs["cell_type_macro"].isin(adata_m.obs["cell_type_macro"]), :
    ],
    var_names=pws_of_interest,
    groupby="cell_type_macro",
    cmap="coolwarm",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
    return_fig=True,
    layer="zscore"
)
fig.savefig(f"{artifact_dir}/heatmap_progeny_tam_zscore.svg")

# %%
fig = sc.pl.matrixplot(
    pb_dorothea[pb_dorothea.obs["cell_type_macro"].isin(adata_m.obs["cell_type_macro"]), :],
    var_names=tfs_of_interest,
    groupby="cell_type_macro",
    cmap="coolwarm",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/heatmap_dorothea_tam.svg")

# %%
pb_dorothea.layers["zscore"] = scipy.stats.zscore(pb_dorothea.X, axis=1)

# %%
fig = sc.pl.matrixplot(
    pb_dorothea[pb_dorothea.obs["cell_type_macro"].isin(adata_m.obs["cell_type_macro"]), :],
    var_names=tfs_of_interest,
    groupby="cell_type_macro",
    cmap="coolwarm",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
    return_fig=True,
    layer="zscore"
)
fig.savefig(f"{artifact_dir}/heatmap_dorothea_tam_zscore.svg")

# %%
