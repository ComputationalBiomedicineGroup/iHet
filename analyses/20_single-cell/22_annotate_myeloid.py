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
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from threadpoolctl import threadpool_limits

sc.settings.set_figure_params(figsize=(4, 4))
import scanpy_helpers as sh

# %%
path_adata_m = nxfvars.get(
    "adata_m",
    "../../data/results/20_single_cell/21_subset_atlas/artifacts/adata_m.h5ad",
)
path_adata_nsclc = nxfvars.get(
    "adata_nsclc",
    "../../data/results/20_single_cell/21_subset_atlas/artifacts/adata_nsclc.h5ad",
)
path_hlca_markers = nxfvars.get(
    "hlca_markers",
    "../../tables/gene_annotations/hlca_cell_type_signatures.csv",
)
cpus = nxfvars.get("cpus", 16)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/ihet/single-cell/")

# %%
threadpool_limits(cpus)

# %%
adata_m = sc.read_h5ad(path_adata_m)
adata_nsclc = sc.read_h5ad(path_adata_nsclc)

# %%
ah = sh.annotation.AnnotationHelper()

# %%
# based on Human Lung Cell Atlas
ah2 = sh.annotation.AnnotationHelper(markers=pd.read_csv(path_hlca_markers))

# %% [markdown]
# # Re-annotate myeloid cells
#
# Create a more fine-graind annotaiton of the myeloid cells than in the atlas

# %%
sc.pl.umap(adata_m, color=["leiden", "cell_type"])

# %%
sc.pl.umap(adata_m, color=["dataset", "study"], wspace=1.5)

# %%
ah.plot_umap(adata_m, filter_cell_type=["DC", "macro", "mono", "div"])

# %%
sc.pl.umap(adata_m, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ah.plot_dotplot(adata_m)

# %%
ah2.score_cell_types(adata_m)

# %%
ah2.plot_dotplot_scores(adata_m)

# %%
ah.annotate_cell_types(
    adata_m,
    {
        "Macrophage alveolar": [0],
        "DC mature": [10],
        "dividing": [6],
        "Monocyte classical": [2],
        "Monocyte non-classical": [5],
        "cDC1": [9],
        "cDC2": [4],
        "7": [7],
        "8": [8],
        "Macrophage": [1, 3, 11],
    },
)

# %% [markdown]
# ## Dividing
# -> subdivide into alveolar/ other

# %%
adata_div = adata_m[adata_m.obs["cell_type"] == "dividing", :].copy()

# %%
ah.reprocess_adata_subset_scvi(adata_div, use_rep="X_scANVI", leiden_res=0.5)

# %%
sc.pl.umap(adata_div, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ah.plot_umap(adata_div, filter_cell_type=["DC", "macro", "mono", "div"])

# %%
ah.annotate_cell_types(
    adata_div,
    {
        "cDC dividing": [0, 9],
        "Macrophage alveolar dividing": [5, 1],
        "Macrophage dividing": [4, 2, 3, 8, 6, 7],
    },
)

# %%
ah.integrate_back(adata_m, adata_div)

# %% [markdown]
# ## cDC2
# -> subdivide into two subclusters

# %%
adata_cdc2 = adata_m[adata_m.obs["cell_type"] == "cDC2", :]

# %%
ah.reprocess_adata_subset_scvi(adata_cdc2, use_rep="X_scANVI", leiden_res=0.5)

# %%
sc.pl.umap(adata_cdc2, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
sc.pl.umap(adata_cdc2, color=["CD1C", "CD1A", "FCER1A"])

# %%
sc.pl.umap(adata_cdc2, color=["C1QA", "C1QB", "C1QC"])

# %%
sc.pl.dotplot(adata_cdc2, groupby="leiden", var_names=["CD1C", "CD1A", "FCER1A"])

# %%
ah.annotate_cell_types(
    adata_cdc2,
    {
        "cDC2": [0, 1, 2, 3, 5, 7],
        "cDC2 CD1A+": [4, 6],
    },
)

# %%
ah.integrate_back(adata_m, adata_cdc2)

# %% [markdown]
# ## Macrophage
# Unsupervised clusters

# %%
adata_macro = adata_m[adata_m.obs["cell_type"] == "Macrophage", :].copy()

# %%
ah.reprocess_adata_subset_scvi(adata_macro, use_rep="X_scANVI", leiden_res=0.3)

# %%
sc.pl.umap(adata_macro, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
sc.pl.umap(adata_macro, color=["platform"])

# %%
ah.annotate_cell_types(
    adata_macro,
    {"TAM-4": [3], "TAM-3": [1], "TAM-1": [0], "TAM-2": [2]},
)

# %%
sc.tl.rank_genes_groups(adata_macro, groupby="cell_type", method="wilcoxon")

# %%
sc.tl.filter_rank_genes_groups(
    adata_macro,
    min_in_group_fraction=0.2,
    max_out_group_fraction=0.2,
)

# %%
fig = sc.pl.rank_genes_groups_dotplot(
    adata_macro, n_genes=10, show=False, return_fig=True, min_logfoldchange=1.5
)
fig.savefig(f"{artifact_dir}/dotplot_tam.svg", dpi=600, bbox_inches="tight")

# %%
fig = sc.pl.umap(
    adata_macro,
    color="cell_type",
    legend_loc="on data",
    legend_fontoutline=2,
    show=False,
    return_fig=True,
    frameon=False,
)
fig.savefig(f"{artifact_dir}/umap_tam.svg", dpi=600, bbox_inches="tight")

# %%
ah.integrate_back(adata_m, adata_macro)

# %% [markdown]
# ## 7/8 what are those? 
# -> maybe use HLCA markers...

# %%
sc.pl.umap(
    adata_m, color=["n_genes_by_counts", "total_counts", "platform"], vmax=[6000, 20000]
)

# %% [markdown]
# --> looks like empty droplets

# %%
adata_m.obs["cell_type"] = adata_m.obs["cell_type"].astype(str)
adata_m.obs.loc[
    adata_m.obs["cell_type"].isin(["7", "8"]), "cell_type"
] = "potential empty droplets"

# %%
sc.pl.umap(adata_m, color="cell_type")

# %%
adata_m.obs["cell_type_macro"] = adata_m.obs["cell_type"]
adata_m.obs["cell_type"] = adata_m.obs["cell_type"].str.replace(
    "TAM(-\d+)", "TAM", regex=True
)

# %%
ah.integrate_back(adata_nsclc, adata_m)

# %%
adata_nsclc.obs["cell_type_macro"] = adata_nsclc.obs["cell_type"]
ah.integrate_back(adata_nsclc, adata_m, variable="cell_type_macro")

# %% [markdown]
# ## Exclude empty droplets

# %%
adata_m = adata_m[adata_m.obs["cell_type"] != "potential empty droplets", :].copy()

# %%
sc.pl.umap(adata_m, color="cell_type_macro")

# %%
adata_nsclc = adata_nsclc[
    adata_nsclc.obs["cell_type"] != "potential empty droplets", :
].copy()

# %%
sc.pl.umap(adata_nsclc, color="cell_type")

# %% [markdown]
# ## Fix cell-type names

# %%
fix_cell_type_names = lambda x: {
    "Neutrophils": "Neutrophil",
    "transitional club/AT2": "Transitional club/AT2",
    "Tumor cells": "Tumor cell",
    "stromal dividing": "Stromal dividing",
}.get(x, x)
adata_nsclc.obs["cell_type_macro"] = adata_nsclc.obs["cell_type_macro"].map(
    fix_cell_type_names
)
adata_nsclc.obs["cell_type"] = adata_nsclc.obs["cell_type"].map(fix_cell_type_names)
adata_m.obs["cell_type"] = adata_m.obs["cell_type"].map(fix_cell_type_names)

# %%
make_coarse_cell_types = {
    "Plasma cell": "Plasma cell",
    "Plasma cell dividing": "Plasma cell",
    "B cell": "B cell",
    "B cell dividing": "B cell",
    "T cell CD4": "T cell",
    "T cell regulatory": "T cell",
    "T cell CD8": "T cell",
    "T cell dividing": "T cell",
    "NK cell": "NK cell",
    "Mast cell": "Mast cell",
    "pDC": "pDC",
    "DC mature": "cDC",
    "cDC1": "cDC",
    "cDC2": "cDC",
    "cDC2 CD1A+": "cDC",
    "cDC dividing": "cDC",
    "Monocyte classical": "Macrophage/Monocyte",
    "Monocyte non-classical": "Macrophage/Monocyte",
    "Macrophage alveolar": "Macrophage/Monocyte",
    "Macrophage alveolar dividing": "Macrophage/Monocyte",
    "Macrophage dividing": "Macrophage/Monocyte",
    "TAM": "Macrophage/Monocyte",
    "Neutrophil": "Neutrophil",
    "Alveolar cell type 1": "Epithelial cell",
    "Alveolar cell type 2": "Epithelial cell",
    "Ciliated": "Epithelial cell",
    "Club": "Epithelial cell",
    "Transitional club/AT2": "Epithelial cell",
    "ROS1+ healthy epithelial": "Epithelial cell",
    "Tumor cell": "Epithelial cell",
    "Endothelial cell arterial": "Stromal cell",
    "Endothelial cell capillary": "Stromal cell",
    "Endothelial cell lymphatic": "Stromal cell",
    "Endothelial cell venous": "Stromal cell",
    "Fibroblast adventitial": "Stromal cell",
    "Fibroblast alveolar": "Stromal cell",
    "Fibroblast peribronchial": "Stromal cell",
    "Mesothelial": "Stromal cell",
    "Pericyte": "Stromal cell",
    "Smooth muscle cell": "Stromal cell",
    "Stromal dividing": "Stromal cell",
}
adata_nsclc.obs["cell_type"] = adata_nsclc.obs["cell_type"].astype(
    pd.CategoricalDtype(categories=make_coarse_cell_types)
)
adata_nsclc.obs["cell_type_coarse"] = (
    adata_nsclc.obs["cell_type"].map(make_coarse_cell_types)
    # keep order by using dict.fromkeys
    .astype(
        pd.CategoricalDtype(categories=dict.fromkeys(make_coarse_cell_types.values()))
    )
)
adata_m.obs["cell_type_coarse"] = (
    adata_m.obs["cell_type"].map(make_coarse_cell_types)
    # keep order by using dict.fromkeys
    .astype(
        pd.CategoricalDtype(categories=dict.fromkeys(make_coarse_cell_types.values()))
    )
)

# %%
assert not np.sum(pd.isnull(adata_nsclc.obs["cell_type"]))
assert not np.sum(pd.isnull(adata_nsclc.obs["cell_type_coarse"]))
assert not np.sum(pd.isnull(adata_m.obs["cell_type"]))
assert not np.sum(pd.isnull(adata_m.obs["cell_type_coarse"]))

# %% [markdown]
# ## Reprocess datasets

# %%
sc.pp.neighbors(adata_m, use_rep="X_scANVI")

# %%
sc.tl.umap(adata_m)

# %%
sc.pl.umap(adata_m, color="cell_type")

# %%
sc.pp.neighbors(adata_nsclc, use_rep="X_scANVI")

# %%
sc.tl.umap(adata_nsclc, init_pos="X_umap")

# %%
sc.pl.umap(adata_nsclc, color="cell_type")

# %%
sc.pl.umap(adata_nsclc, color="cell_type_coarse")

# %% [markdown]
# # Save reprocessed anndatas

# %%
adata_nsclc.write_h5ad(f"{artifact_dir}/adata_nsclc_reannotated.h5ad")
adata_m.write_h5ad(f"{artifact_dir}/adata_myeloid_reannotated.h5ad")

# %% [markdown]
# ### save for cellxgene

# %%
adata_nsclc_cellxgene = adata_nsclc.raw.to_adata()

# %%
adata_m_cellxgene = adata_m.raw.to_adata()

# %%
adata_nsclc_cellxgene.write_h5ad(
    f"{artifact_dir}/adata_nsclc_reannotated_cellxgene.h5ad", compression="lzf"
)
adata_m_cellxgene.write_h5ad(
    f"{artifact_dir}/adata_myeloid_reannotated_cellxgene.h5ad", compression="lzf"
)

# %%
