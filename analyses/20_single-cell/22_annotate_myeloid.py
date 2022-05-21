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
ah = sh.annotation.AnnotationHelper()

# %%
# based on Human Lung Cell Atlas
ah2 = sh.annotation.AnnotationHelper(
    markers=pd.read_csv(
        nxfvars.get(
            "hlca_markers",
            "../../tables/gene_annotations/hlca_cell_type_signatures.csv",
        )
    )
)

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
        "monocyte conventional": [2],
        "monocyte non-conventional": [5],
        "cDC1": [9],
        "cDC2": [4],
        "7": [7],
        "8": [8],
        "Macrophage": [1, 3, 11],
    },
)

# %% [markdown]
# ```
# DC_Lagerhans --> CD1C+CD1A+ DCs
# Macro_M0 --> SLAMF9 macro
# Macro_M1_like --> CD163+ macro
# Macro_M2_like_alevolar --> MARCO+ macro
# myeloid_TF_low --> Other
#
# For the other cells, we might just refine a bit the final nomenclature.
#
# There are also a couple of refinements in the clusters that we could perform:
#
# Restrict the NK cell cluster, which now wrongly contains some CD8 + T cells (expressing CD3 and TCRalpha chain genes)
# The cDC2 cluster should be made by CD1C+ DC also expressing FCER1A.
# The  DC_Lagerhans cluster should correspond to CD1C+CD1A+ DCs. In brief, the CD1C+ blob should cover this and the cDC2 subset
# ```

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
# -> subdivide into two subclusters

# %%
adata_macro = adata_m[adata_m.obs["cell_type"] == "Macrophage", :].copy()

# %%
ah.reprocess_adata_subset_scvi(adata_macro, use_rep="X_scANVI", leiden_res=0.5)

# %%
sc.pl.umap(adata_macro, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
sc.pl.umap(adata_macro, color=["SLAMF9", "CD163", "MARCO", "platform"])

# %%
ah.plot_umap(adata_macro, filter_cell_type=["macro", "mono"])

# %%
ah.annotate_cell_types(
    adata_macro,
    {
        "Macrophage SLAMF9+": [2],
        "Macrophage MARCO+": [4, 3],
        "Macrophage CD163+": [0, 5, 1, 6, 7],
    },
)

# %%
ah.integrate_back(adata_m, adata_macro)

# %% [markdown]
# ## 7/8 what are those? 
# -> maybe use HLCA markers...

# %%
sc.pl.umap(adata_m, color=["n_genes_by_counts", "total_counts", "platform"], vmax=[6000, 20000])

# %% [markdown]
# --> looks like empty droplets

# %%
adata_m.obs["cell_type"] = adata_m.obs["cell_type"].astype(str)
adata_m.obs.loc[adata_m.obs["cell_type"].isin(["7", "8"]), "cell_type"] = "potential empty droplets"

# %%
sc.pl.umap(adata_m, color="cell_type")

# %%
ah.integrate_back(adata_nsclc, adata_m)

# %% [markdown]
# ## Exclude empty droplets

# %%
adata_m = adata_m[adata_m.obs["cell_type"] != "potential empty droplets", :].copy()

# %%
sc.pl.umap(adata_m, color="cell_type")

# %%
adata_nsclc = adata_nsclc[adata_nsclc.obs["cell_type"] != "potential empty droplets", :].copy()

# %%
sc.pl.umap(adata_nsclc, color="cell_type")

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
sc.tl.umap(adata_nsclc)

# %%