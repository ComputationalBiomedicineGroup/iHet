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
pb_progeny = sh.pseudobulk.pseudobulk(adatas["nsclc"]["progeny"], groupby=["patient", "cell_type"], aggr_fun=np.mean)
pb_progeny.obs["cell_type"] = pb_progeny.obs["cell_type"].astype(adatas["nsclc"]["progeny"].obs["cell_type"].dtype)

# %% [markdown]
# Note: another round of z-scoreing after pseudobulking hardly makes any difference.

# %%
fig = sc.pl.matrixplot(
    pb_progeny,
    var_names=pws_of_interest,
    groupby="cell_type",
    cmap="coolwarm",
    swap_axes=False,
    vmin=-2.5,
    vmax=2.5,
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/heatmap_progeny.svg")


# %%
def _export_df(pb):
    df = pb.to_df()
    df["cell_type"] = pb.obs["cell_type"]
    return df.groupby("cell_type").agg("mean").T


_export_df(pb_progeny).to_csv(f"{artifact_dir}/table_progeny.csv")

# %%
pb_dorothea = sh.pseudobulk.pseudobulk(adatas["nsclc"]["dorothea"], groupby=["patient", "cell_type"], aggr_fun=np.mean)
pb_dorothea.obs["cell_type"] = pb_dorothea.obs["cell_type"].astype(adatas["nsclc"]["dorothea"].obs["cell_type"].dtype)

# %%
fig = sc.pl.matrixplot(
    pb_dorothea,
    var_names=tfs_of_interest,
    groupby="cell_type",
    cmap="coolwarm",
    swap_axes=False,
    vmin=-2.5,
    vmax=2.5,
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/heatmap_dorothea.svg")

_export_df(pb_dorothea).to_csv(f"{artifact_dir}/table_dorothea.csv")

# %% [markdown]
# ## Including myeloid subclusters

# %%
pb_progeny = sh.pseudobulk.pseudobulk(adatas["nsclc"]["progeny"], groupby=["patient", "cell_type"], aggr_fun=np.mean)
pb_progeny.obs["cell_type"] = pb_progeny.obs["cell_type"].astype(adatas["nsclc"]["progeny"].obs["cell_type"].dtype)

# %%
pb_dorothea = sh.pseudobulk.pseudobulk(
    adatas["nsclc"]["dorothea"],
    groupby=["patient", "cell_type"],
    aggr_fun=np.mean,
)
pb_dorothea.obs["cell_type"] = pb_dorothea.obs["cell_type"].astype(adatas["nsclc"]["dorothea"].obs["cell_type"].dtype)

# %%
fig = sc.pl.matrixplot(
    pb_progeny[pb_progeny.obs["cell_type"].isin(adata_m.obs["cell_type"]), :],
    var_names=pws_of_interest,
    groupby="cell_type",
    cmap="coolwarm",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/heatmap_progeny_myeloid.svg")

# %%
pb_progeny.layers["zscore"] = scipy.stats.zscore(pb_progeny.X, axis=1)

# %%
fig = sc.pl.matrixplot(
    pb_progeny[pb_progeny.obs["cell_type"].isin(adata_m.obs["cell_type"]), :],
    var_names=pws_of_interest,
    groupby="cell_type",
    cmap="coolwarm",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
    return_fig=True,
    layer="zscore",
)
fig.savefig(f"{artifact_dir}/heatmap_progeny_myeloid_zscore.svg")

# %%
fig = sc.pl.matrixplot(
    pb_dorothea[pb_dorothea.obs["cell_type"].isin(adata_m.obs["cell_type"]), :],
    var_names=tfs_of_interest,
    groupby="cell_type",
    cmap="coolwarm",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/heatmap_dorothea_myeloid.svg")

# %%
pb_dorothea.layers["zscore"] = scipy.stats.zscore(pb_dorothea.X, axis=1)

# %%
fig = sc.pl.matrixplot(
    pb_dorothea[pb_dorothea.obs["cell_type"].isin(adata_m.obs["cell_type"]), :],
    var_names=tfs_of_interest,
    groupby="cell_type",
    cmap="coolwarm",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
    return_fig=True,
    layer="zscore",
)
fig.savefig(f"{artifact_dir}/heatmap_dorothea_myeloid_zscore.svg")

# %% [markdown]
# ## UMAP plots of selected features

# %%
fig, axs = plt.subplots(2, 4, figsize=(18, 8))
axs = iter(axs.flatten())
for pw, ax in zip(["TNFa", "NFkB", "TGFb"], axs):
    sc.pl.umap(
        adatas["nsclc"]["progeny"],
        color=pw,
        ax=ax,
        show=False,
        cmap="coolwarm",
        vmin=-2.5,
        vmax=2.5,
        title=f"PW: {pw}",
        frameon=False,
    )
for tf, ax in zip(["SPI1", "MYCN", "RFX5", "E2F4", "TFDP1"], axs):
    sc.pl.umap(
        adatas["nsclc"]["dorothea"],
        color=tf,
        ax=ax,
        show=False,
        cmap="coolwarm",
        vmin=-2.5,
        vmax=2.5,
        title=f"TF: {tf}",
        frameon=False,
    )
fig.savefig(
    f"{artifact_dir}/umap_cell_selected_features.svg",
    bbox_inches="tight",
    dpi=600,
)
plt.show()

# %%
