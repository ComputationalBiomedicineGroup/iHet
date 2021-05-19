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
from scanpy_helpers.annotation import AnnotationHelper

sc.set_figure_params(figsize=(4,4))

# %%
adata_myeloid = sc.read_h5ad("../../data/50_myeloid_cells/20_cell_type_annotation/adata_annotated.h5ad")
adata = sc.read_h5ad("../../data/20_annotate_all_cells/primary-tumor-annotated-integrated.h5ad")

# %%
ah = AnnotationHelper()

# %%
adata.obs["cell_type"] = adata.obs["cell_type_coarse"]

# %%
sc.pl.umap(adata_myeloid, color="cell_type")

# %%
sc.pl.umap(adata, color="cell_type")

# %%
ah.integrate_back(adata, adata_myeloid)

# %%
adata = adata[adata.obs["cell_type"] != "other", :]

# %%
adata.write_h5ad("../../data/60_cellphonedb/01_prepare_input_data/adata.h5ad")

# %%
adata.obs.index.name = "Cell"
adata.obs.loc[:, "cell_type"].to_csv("../../data/60_cellphonedb/01_prepare_input_data/adata.obs.csv")
