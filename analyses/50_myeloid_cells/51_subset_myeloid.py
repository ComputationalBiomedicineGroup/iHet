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
import numpy as np
import pandas as pd
from warnings import filterwarnings
sc.settings.verbosity = 3
filterwarnings('ignore', category=FutureWarning)
sc.set_figure_params(figsize=(5, 5))

# %%
output_dir = "../../data/50_myeloid_cells/01_myeloid_subset/"

# %% [markdown]
# # Extract myeloid cells from integrated lung cancer dataset

# %%
adata_raw_counts = sc.read_h5ad(
    "../../data/10_processed/sc_rna_seq/merged_nsclc_heterogeneity.h5ad"
)

# %%
adata = sc.read_h5ad(
    "../../data/40_dorothea_progeny/40_primary_tumor_dorothea_progeny.h5ad"
)

# %%
adata.obs.to_csv(f"{output_dir}/obs.csv")

# %%
sc.pl.umap(adata, color="cell_type")

# %%
tmp_adata = adata[
    adata.obs["cell_type"].isin(
        [
            "Macrophage",
            "Monocyte conventional",
            "Monocyte non-conventional",
            "cDC1",
            "cDC2",
            "DC mature",
            "Macrophage alevolar",
        ]
    ),
    :,
].copy()
adata = tmp_adata

# %%
sc.pl.umap(adata, color="cell_type")

# %% [markdown]
# ## Export raw counts for analysis in R

# %%
pd.DataFrame(adata.X.T.toarray(), columns=adata.obs_names, index=adata.var_names).to_csv(f"{output_dir}/X.csv")

# %% [markdown]
# ## Export anndata with raw counts for analysis with scVI and DE methods

# %%
adata_scvi = adata.copy()
adata_scvi.X = adata_raw_counts[adata_scvi.obs_names, :].X.astype("int")

# %%
adata.write_h5ad(f"{output_dir}/adata.h5ad")
adata_scvi.write_h5ad(f"{output_dir}/adata_scvi.h5ad")
