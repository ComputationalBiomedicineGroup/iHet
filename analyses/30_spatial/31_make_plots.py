# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.4
# ---

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import scanpy as sc
from nxfvars import nxfvars
import pandas as pd
from threadpoolctl import threadpool_limits
import matplotlib.pyplot as plt
import decoupler as dc
import numpy as np

sc.settings.set_figure_params(figsize=(4, 4))

# %%
path_visium = nxfvars.get(
    "visium_slide",
    "/home/sturm/projects/2021/nsclc_heterogeneity/data/01_processed/spatial/lungCancerVisium/outs",
)
path_cell2location_estimates = nxfvars.get(
    "cell2location",
    "/home/sturm/projects/2021/nsclc_heterogeneity/data/01_processed/spatial/cell2location_estimates.csv",
)
cpus = nxfvars.get("cpus", 16)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/ihet")

# %%
threadpool_limits(cpus)

# %%
adata = sc.read_visium(path_visium)
cell2location_estimates = pd.read_csv(path_cell2location_estimates, index_col=0).drop(
    columns=["in_tissue", "array_row", "array_col"]
)
adata = adata[cell2location_estimates.index,]
adata.obs = adata.obs.join(cell2location_estimates)
adata.var_names_make_unique()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# %%
progeny = dc.get_progeny(organism="human", top=500)
dc.run_wmean(adata, progeny, use_raw=False, times=0, min_n=5, verbose=True)

# %%
# Aggregate cell_types
adata.obs["Myeloid cells"] = adata.obs.loc[
    :,
    [
        "c2l_Monocyte.classical",
        "c2l_cDC2",
        "c2l_Macrophage.alveolar",
        "c2l_cDC2.CD1A.",
        "c2l_TAM",
        "c2l_DC.mature",
        "c2l_Monocyte.non.classical",
        "c2l_cDC1",
    ],
].sum(axis=1)

adata.obs["Fibroblasts"] = adata.obs.loc[
    :,
    [
        "c2l_Fibroblast.adventitial",
        "c2l_Fibroblast.alveolar",
        "c2l_Fibroblast.peribronchial",
    ],
].sum(axis=1)


# %%
# Compute log ratios (same as spacedeconv plot_comparison)
adata.obs["Fibroblasts/Monocytes"] = np.log(adata.obs["c2l_Fibroblast.peribronchial"] + 1) - np.log(
    adata.obs["c2l_Monocyte.classical"] + 1
)
adata.obs["TGFb/TNFa"] = np.log(adata.obsm["wmean_estimate"]["TGFb"] + 1) - np.log(
    adata.obsm["wmean_estimate"]["TNFa"] + 1
)
adata.obs["TGFb/NFkB"] = np.log(adata.obsm["wmean_estimate"]["TGFb"] + 1) - np.log(
    adata.obsm["wmean_estimate"]["NFkB"] + 1
)

# %%
fig = sc.pl.spatial(adata, color="Myeloid cells", cmap="inferno", alpha_img=0, size=1.5, frameon=False, return_fig=True)
fig.savefig(f"{artifact_dir}/myeloid.svg", bbox_inches="tight", dpi=600)
fig = sc.pl.spatial(adata, color="Fibroblasts", cmap="inferno", alpha_img=0, size=1.5, frameon=False, return_fig=True)
fig.savefig(f"{artifact_dir}/fibroblasts.svg", bbox_inches="tight", dpi=600)
fig = sc.pl.spatial(
    adata,
    color="c2l_T.cell.CD8",
    cmap="inferno",
    alpha_img=0,
    size=1.5,
    frameon=False,
    title="CD8+ T cells",
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/t_cd8.svg", bbox_inches="tight", dpi=600)
fig = sc.pl.spatial(
    adata,
    color="c2l_Monocyte.classical",
    cmap="inferno",
    alpha_img=0,
    size=1.5,
    frameon=False,
    title="Monocytes classical",
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/monocytes_classical.svg", bbox_inches="tight", dpi=600)
fig = sc.pl.spatial(
    adata,
    color="c2l_Fibroblast.peribronchial",
    cmap="inferno",
    alpha_img=0,
    size=1.5,
    frameon=False,
    title="Fibroblasts peribronchial",
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/fibroblasts_peribronchial.svg", bbox_inches="tight", dpi=600)
fig = sc.pl.spatial(
    adata, color="c2l_TAM", cmap="inferno", alpha_img=0, size=1.5, frameon=False, title="TAMs", return_fig=True
)
fig.savefig(f"{artifact_dir}/TAMs.svg", bbox_inches="tight", dpi=600)

# %%
fig = sc.pl.spatial(
    adata,
    color="Fibroblasts/Monocytes",
    cmap="PRGn_r",
    alpha_img=0,
    size=1.5,
    frameon=False,
    title="CD8+ T cells",
    vmin=-0.5,
    vmax=0.5,
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/ratio_fibroblasts_monocytes.svg", bbox_inches="tight", dpi=600)

# %%
fig = sc.pl.spatial(
    adata,
    color="TGFb/TNFa",
    cmap="PRGn_r",
    alpha_img=0,
    size=1.5,
    frameon=False,
    title="CD8+ T cells",
    vmin=-0.2,
    vmax=0.2,
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/ratio_tbfb_tnfa.svg", bbox_inches="tight", dpi=600)


# %%
fig = sc.pl.spatial(
    adata,
    color="TGFb/NFkB",
    cmap="PRGn_r",
    alpha_img=0,
    size=1.5,
    frameon=False,
    title="CD8+ T cells",
    vmin=-0.2,
    vmax=0.2,
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/ratio_tbfb_nfkb.svg", bbox_inches="tight", dpi=600)

# %%
