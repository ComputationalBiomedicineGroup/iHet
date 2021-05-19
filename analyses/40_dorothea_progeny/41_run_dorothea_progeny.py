# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Python [conda env:conda-2021-nsclc_heterogeneity-scanpy]
#     language: python
#     name: conda-env-conda-2021-nsclc_heterogeneity-scanpy-py
# ---

# %%
import scanpy as sc
import logging

import anndata2ri
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Colormap, TwoSlopeNorm
import progeny
import dorothea

sc.set_figure_params(figsize=(4, 4))

# %%
adata = sc.read_h5ad(
    "../../data/20_annotate_all_cells/primary-tumor-annotated-integrated.h5ad"
)

# %%
sc.pl.umap(adata, color=["cell_type"])

# %%
sc.pl.umap(adata, color=["SPI1"])

# %% [markdown]
# ### Run progeny

# %%
model = progeny.getModel(organism="Human", top=1000)

# %%
progeny.run(adata, model, center=False, scale=True)

# %%
adata.obsm["progeny"].columns = [f"PW:{pw}" for pw in adata.obsm["progeny"].columns]

# %%
pathways_of_interest = [
    f"PW:{x}" for x in ["JAK-STAT", "VEGF", "PI3K", "NFkB", "Trail", "TNFa"]
]

# %%
adata.obs = adata.obs.join(
    adata.obsm["progeny"].loc[:, pathways_of_interest]
)

# %%
adata

# %% [markdown]
# ### Run dorothea

# %%
regulons = dorothea.load_regulons(
    [
        "A",
        "B",
    ],  # Which levels of confidence to use (A most confident, E least confident)
    organism="Human",  # If working with mouse, set to Mouse
)

# %%
dorothea.run(
    adata,
    regulons,
    center=True,  # Center gene expression by mean per cell
    num_perm=100,  # Simulate m random activities
    norm=True,  # Normalize by number of edges to correct for large regulons
    scale=True,  # Scale values per feature so that values can be compared across cells
    min_size=5,  # TF with less than 5 targets will be ignored
)

# %%
adata.obsm["dorothea"].columns = [f"TF:{tf}" for tf in adata.obsm["dorothea"].columns]

# %%
tfs_of_interest = [
    f"TF:{tf}" for tf in ["SPI1", "NFKB1", "STAT1", "MYC", "E2F4", "E2F2", "ZNF263"]
]

# %%
adata.obs = adata.obs.join(adata.obsm["dorothea"].loc[:, tfs_of_interest])

# %%
adata

# %% [markdown]
# ## Export adata

# %%
adata.write_h5ad(
    "../../data/40_dorothea_progeny/40_primary_tumor_dorothea_progeny.h5ad",
    compression="lzf",
)

# %% [markdown]
# ## Visualize
# ### Transcription factors

# %%
# fake adata object for stacked violin
dorothea_df = adata.obsm["dorothea"]
tmp_adata = sc.AnnData(X=dorothea_df.values.copy(), obs=adata.obs.loc[:, ["cell_type"]])
tmp_adata.var_names = dorothea_df.columns

# %%
sc.pp.scale(tmp_adata)

# %%
norm = TwoSlopeNorm(0, vmin=-5, vmax=5)
sc.pl.umap(adata, color=tfs_of_interest, ncols=3, cmap="coolwarm", norm=norm)

# %%
sc.pl.matrixplot(
    tmp_adata,
    groupby="cell_type",
    var_names=tfs_of_interest,
    swap_axes=True,
    dendrogram=True,
    cmap="bwr",
)

# %%
sc.pl.stacked_violin(
    tmp_adata,
    groupby="cell_type",
    var_names=tfs_of_interest,
    swap_axes=True,
    dendrogram=True,
    cmap="bwr",
)

# %% [markdown]
# ### Pathways

# %%
# fake adata object for stacked violin
progeny_df = adata.obsm["progeny"]
tmp_adata = sc.AnnData(
    X=progeny_df.values.copy(), obs=adata.obs.loc[:, ["cell_type"]], obsm=adata.obsm
)
tmp_adata.var_names = progeny_df.columns

# %%
sc.pp.scale(tmp_adata)

# %%
pathways_of_interest

# %%
progeny_df.loc[:, "PW:JAK-STAT"]

# %%
norm = TwoSlopeNorm(0, vmin=np.min(tmp_adata.X), vmax=np.max(tmp_adata.X))
sc.pl.umap(tmp_adata, color=pathways_of_interest, ncols=3, norm=norm, cmap="coolwarm")

# %%
sc.pl.matrixplot(
    tmp_adata,
    groupby="cell_type",
    var_names=pathways_of_interest,
    swap_axes=True,
    dendrogram=True,
    cmap="bwr",
)

# %%
sc.pl.stacked_violin(
    adata,
    groupby="cell_type",
    var_names=pathways_of_interest,
    swap_axes=True,
    dendrogram=True,
    cmap="bwr",
)

# %% [markdown]
# ### Genes

# %%
genes_of_interest = regulons.index[regulons["SPI1"] != 0]

# %%
genes_of_interest = list(genes_of_interest) + ["TREM2", "SPI1"]

# %%
genes_of_interest = list(set(genes_of_interest) & set(adata.var_names.values))

# %%
sc.pl.umap(adata, color=sorted(genes_of_interest), cmap="inferno")

# %%
sc.pl.matrixplot(
    adata,
    groupby="cell_type",
    var_names=sorted(genes_of_interest),
    swap_axes=True,
    dendrogram=True,
    cmap="bwr",
    standard_scale="var",
)
