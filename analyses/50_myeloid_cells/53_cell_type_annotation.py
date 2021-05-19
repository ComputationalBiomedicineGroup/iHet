# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: SSH apollo-15 2020-nsclc_heterogeneity-scanpy
#     language: ''
#     name: rik_ssh_apollo_15_2020nsclc_heterogeneityscanpy
# ---

# %%
# %load_ext autoreload
# %autoreload 2
import scanpy as sc
import scvi
import numpy as np
import pandas as pd
from scanpy_helpers.annotation import AnnotationHelper
from warnings import filterwarnings

sc.settings.verbosity = 3

filterwarnings("ignore", category=FutureWarning)
sc.set_figure_params(figsize=(5, 5))

# %%
scvi.__version__

# %%
cell_type_scores_dir = "../../data/50_myeloid_cells/10_cell_type_scoring/"
output_dir = "../../data/50_myeloid_cells/20_cell_type_annotation/"
deploy_dir = "../../deploy/Ab1MxvlS9c56BaA0Z94J/figures/"

# %%
adata = sc.read_h5ad("../../data/50_myeloid_cells/01_myeloid_subset/adata.h5ad")
adata_scvi = sc.read_h5ad(
    "../../data/50_myeloid_cells/01_myeloid_subset/adata_scvi.h5ad"
)

# %% [markdown]
# The two adatas are identical, excapt `adata_scvi` contains the counts in `X`, `adata` the log-normalized counts. 

# %%
assert adata.shape == adata_scvi.shape

# %% [markdown]
# # Cell-type annotation
# In this notebook, we perform manual cell-type 
# annotation based on unsupervised clustering, marker genes, and cell-type scores
# obtained from the deconvolution methods in the previous notebook. 

# %% [markdown]
# ### Run scVI
#
# We use scVI to generate a batch-corrected, latent representation of the myeloid single-cells. 
# Since we also want to use scVI for differential expression testing, I do not subset it to the
# most highly variable cells. Since `n_cell` $\approx$ `n_genes`, this is ok according to the authors. 
# I also visually compared the result to using the subset only, and it changed only slightly. 

# %%
scvi.data.setup_anndata(adata_scvi, batch_key="sample")
scvi.data.view_anndata_setup(adata_scvi)

# %% [markdown]
# This would be the code to re-train the model. Since scVI is not 100% reproducible (and in particular not reproducible 
# on another system), I stored and load the pre-computed model instead. 

# %%
# model = scvi.model.SCVI(adata_scvi, use_cuda=True)
# model.train()
# model.save(f"{output_dir}/scvi_model_all_genes")

# %%
# %%capture
# Load the pre-computed model.
model = scvi.model.SCVI.load(f"{output_dir}/scvi_model_all_genes", adata_scvi)

# %%
adata.obsm["X_scVI"] = model.get_latent_representation(adata_scvi)

# %%
adata.obsm["X_scVI"]

# %% [markdown]
# ### Compute UMAP based on scVI

# %%
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, random_state=0)

# %%
sc.pl.umap(adata, color=["dataset", "sample"], legend_loc="none")

# %%
sc.pl.umap(
    adata, color=["cell_type", "leiden"], legend_loc="on data", legend_fontoutline=2
)

# %%
adata.obs["n_genes"] = np.sum(adata.X != 0, axis=1)

# %% [markdown]
# ## Annotate cell types

# %% [markdown]
# ### plot cell-type markers
#
# Cell-type markers have been curated from
#  * [Schupp et al](https://www.biorxiv.org/content/10.1101/2020.10.21.347914v1)
#  * [Lambrechts et al.](https://doi.org/10.1038/s41591-018-0096-5)
#  * [Madissoon et al. ](https://doi.org/10.1186/s13059-019-1906-x)
#  * [Alshetaiwi et al.](https://immunology.sciencemag.org/content/5/44/eaay6017)
#
# and are summarized in [this google sheet](https://docs.google.com/spreadsheets/d/1beW-9oeM31P50NFvNLVvsdXlfh_tOjsmIrnwV2ZlxDU/edit#gid=966184166). 

# %%
ah = AnnotationHelper()

# %%
ah.plot_umap(adata, filter_cell_type=["Macro", "Mono", "DC", "MDSC", "Div"])

# %% [markdown]
# ## Load results form SingleR, xCell and CIBERSORT
#
# Also tested quanTIseq, but does not perform well for this use case. 

# %% [markdown]
# ### SingleR scores

# %%
singler_res = pd.read_csv(f"{cell_type_scores_dir}/singler_res.csv")

# %%
adata.obs["singler_macrophages"] = singler_res["Macrophages"].values
adata.obs["singler_mono"] = singler_res["Monocytes"].values
adata.obs["singler_m1"] = singler_res["Macrophages M1"].values
adata.obs["singler_m2"] = singler_res["Macrophages M2"].values
adata.obs["singler_labels"] = singler_res["labels"].values

# %%
sc.pl.umap(
    adata,
    color=[
        "singler_macrophages",
        "singler_mono",
        "singler_m1",
        "singler_m2",
        "singler_labels",
    ],
)

# %% [markdown]
# ### xCell

# %%
xcell_res = pd.read_csv(f"{cell_type_scores_dir}/xcell_res.csv").set_index("cell_type")

# %%
for cell_type in xcell_res.index.values:
    adata.obs[f"xcell {cell_type}"] = xcell_res.loc[cell_type, :].values

# %%
sc.pl.umap(
    adata,
    color=[f"xcell {cell_type}" for cell_type in xcell_res.index.values] + ["leiden"],
)

# %% [markdown]
# ### CIBERSORT
#
# CIBERSORT failed for a few (~150) cells, but this appears random (possibly all signature genes = 0) 
# and not systematic for a certain subcluster. 

# %%
cibersort_res = pd.read_csv(f"{cell_type_scores_dir}/cibersort_res.csv").set_index(
    "cell_type"
)

# %%
cibersort_cell_types = [
    "Monocyte",
    "Macrophage M0",
    "Macrophage M1",
    "Macrophage M2",
    "Myeloid dendritic cell resting",
    "Myeloid dendritic cell activated",
]
for cell_type in cibersort_cell_types:
    adata.obs[f"cibersort {cell_type}"] = cibersort_res.loc[cell_type, :]

# %%
sc.pl.umap(
    adata, color=[f"cibersort {cell_type}" for cell_type in cibersort_cell_types]
)

# %% [markdown]
# ### Plot Transcription factors

# %%
# pathways_of_interest = [
#     f"PW:{x}" for x in ["JAK-STAT", "VEGF", "PI3K", "NFkB", "Trail", "TNFa"]
# ]

# %%
# tfs_of_interest = [
#     f"TF:{tf}" for tf in ["SPI1", "NFKB1", "STAT1", "MYC", "E2F4", "E2F2", "ZNF263"]
# ]

# %%
sc.pl.umap(
    adata,
    color=[
        "leiden",
        "cell_type",
        "TF:MYC",
        "TF:SPI1",
        "TF:STAT1",
        "TF:NFKB1",
        "TF:E2F4",
    ],
    ncols=3,
    wspace=0.6,
    legend_fontoutline=1,
    cmap="coolwarm",
)

# %% [markdown]
# ### Annotate cell types
#
# I annotate cell-types by manually assigning each leiden cluster a cell-type labels
# based on the plots of marker genes and deconvolution scores above. 
# Sometimes the border is ambiguous which warrants subclustering to refine it. 

# %%
sc.pl.umap(adata, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "mDC mature": [10],
    "cDC1": [8],
    "C9": [9],
    "C4": [4],
    "cDC2 CD1A+": [7],
    "ambiguous1": [0, 5, 3],
    "Macro CD163+": [1],
    "Macro MARCO+": [2],
    "Other": [6, 9],
}

# %%
ah.annotate_cell_types(adata, ct_map)

# %% [markdown]
# ### Subclustering of C4

# %%
adata_4 = adata[adata.obs["cell_type"] == "C4", :]

# %%
ah.reprocess_adata_subset_scvi(adata_4, leiden_res=0.3)

# %%
sc.pl.umap(adata_4, color=["leiden", "TF:E2F4"])

# %%
ah.plot_umap(adata_4, filter_cell_type=["Dividing", "cDC"])

# %%
ah.annotate_cell_types(
    adata_4,
    {
        "Myeloid dividing": [1],
        "cDC2": [0, 2],
    },
)

# %%
ah.integrate_back(adata, adata_4)

# %% [markdown]
# ### Subclustering of 'ambigous1'

# %%
adata_3 = adata[adata.obs["cell_type"] == "ambiguous1", :]

# %%
ah.reprocess_adata_subset_scvi(adata_3, leiden_res=0.5)

# %%
sc.pl.umap(
    adata_3,
    color=[
        "leiden",
        "cibersort Macrophage M0",
        "cibersort Monocyte",
        "xcell Monocyte",
        "APOE",
        "VCAN",
    ],
)

# %%
ah.plot_umap(adata_3, filter_cell_type=["MDSC", "Mono", "cDC2", "Macro"])

# %%
ct_map = {
    "cDC2": [0, 1],
    "Monocyte": [2, 4, 6, 5, 7],
    "Macro SLAMF9+": [3],
}

# %%
ah.annotate_cell_types(adata_3, ct_map)

# %% [markdown]
# ### Final cell-type annotation

# %%
ah.integrate_back(adata, adata_3)

# %%
sc.set_figure_params(figsize=(12, 10))
ax =sc.pl.umap(adata, color="cell_type", size=30, show=False)
fig = ax.get_figure()
fig.tight_layout()
fig.savefig(f"{deploy_dir}/umap_cell_type_myeloid.pdf")

# %%
sc.set_figure_params(figsize=(6,6))
sc.pl.umap(adata, color="dataset")

# %% [markdown]
# ## Export annotated anndata object

# %%
adata.obsm["dorothea"].join(adata.obsm["progeny"]).join(adata.obs.loc[:, ["patient", "condition", "dataset", "sex", "cell_type"]]).to_csv(f"{output_dir}/obs_dorothea_progeny.csv")

# %%
adata.write_h5ad("../../data/50_myeloid_cells/20_cell_type_annotation/adata_annotated.h5ad")
