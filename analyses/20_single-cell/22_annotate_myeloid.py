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

sc.settings.set_figure_params(figsize=(4, 4))
import scanpy_helpers as sh

# %%
path_adata_m = nxfvars.get(
    "adata_m", "../../data/results/20_single_cell/21_subset_atlas/artifacts/adata_m.h5ad"
)
path_adata_nsclc = nxfvars.get(
    "adata_nsclc", "../../data/results/20_single_cell/21_subset_atlas/artifacts/adata_nsclc.h5ad"
)
path_hlca_markers = nxfvars.get(
    "hlca_markers",
    "../../tables/gene_annotations/hlca_cell_type_signatures.csv",
)
cpus = nxfvars.get("cpus", 16)
artifact_dir = nxfvars.get(
    "artifact_dir", "/home/sturm/Downloads/ihet/single-cell/"
)

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
ah.reprocess_adata_subset_scvi(adata_macro, use_rep="X_scANVI", leiden_res=0.3)

# %%
sc.pl.umap(adata_macro, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
sc.tl.rank_genes_groups(adata_macro, groupby="leiden", method="wilcoxon")

# %%
sc.tl.filter_rank_genes_groups(adata_macro, min_in_group_fraction=0.2, max_out_group_fraction=0.2, min_fold_change=1)

# %%
sc.pl.rank_genes_groups_dotplot(adata_macro, min_logfoldchange=1, n_genes=20)

# %%
pb_macro = sh.pseudobulk.pseudobulk(adata_macro, groupby=["patient", "leiden"])

# %%
sc.pp.normalize_total(pb_macro)
sc.pp.log1p(pb_macro)

# %%
sc.tl.rank_genes_groups(pb_macro, groupby="leiden", method="wilcoxon")

# %%
sc.pl.rank_genes_groups_matrixplot(
    pb_macro,
    standard_scale=None,
    min_logfoldchange=2,
    values_to_plot="logfoldchanges",
    cmap="bwr",
    n_genes=20,
)

# %%
sc.tl.filter_rank_genes_groups(pb_macro, min_fold_change=2, min_in_group_fraction=0, max_out_group_fraction=1)

# %%
pd.set_option("display.max_rows", 300)

# %%
{k: pd.Series(v).dropna().tolist()[:20] for k, v in pd.DataFrame(pb_macro.uns["rank_genes_groups_filtered"]["names"]).to_dict(orient='list').items()}

# %%
sc.pl.rank_genes_groups_matrixplot(
    pb_macro,
    standard_scale=None,
    min_logfoldchange=2,
    cmap="viridis",
    n_genes=20,
)

# %%
sc.pl.rank_genes_groups_matrixplot(
    pb_macro,
    standard_scale=None,
    min_logfoldchange=2,
    values_to_plot="log10_pvals_adj",
    cmap="viridis",
    n_genes=20,
    vmax=6
)

# %%
sc.pl.umap(adata_macro, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
sc.pl.umap(adata_macro, color=["MCEMP1", "S100A8", "S100A4", "FBP1", "MARCO"], ncols=5)

# %%
sc.pl.umap(adata_macro, color=["CTSB", "SLAMF9", "SPP1", "NMB"], ncols=5)

# %%
sc.pl.umap(adata_macro, color=["MS4A6A", "F13A1", "CD74"])

# %%
sc.pl.umap(adata_macro, color=["APOE", "APOC1", "PLD3", "ACP5", "CCL18"])

# %%
sc.pl.umap(adata_macro, color=["platform"])

# %%
ah.plot_umap(adata_macro, filter_cell_type=["macro", "mono"])

# %%
ah.annotate_cell_types(
    adata_macro,
    {
        "Macrophage SPP1-hi": [3],
        "Macrophage MARCO-hi": [1],
        "Macrophage CD74-hi": [0],
        "Macrophage CCL18-hi": [2]
    },
)

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
ah.integrate_back(adata_nsclc, adata_m)

# %% [markdown]
# ## Exclude empty droplets

# %%
adata_m = adata_m[adata_m.obs["cell_type"] != "potential empty droplets", :].copy()

# %%
sc.pl.umap(adata_m, color="cell_type")

# %%
adata_nsclc = adata_nsclc[
    adata_nsclc.obs["cell_type"] != "potential empty droplets", :
].copy()

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
sc.pl.umap(adata_nsclc, color="cell_type")

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
