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
from threadpoolctl import threadpool_limits

# %%
sc.settings.set_figure_params(figsize=(4, 4))

# %%
adata_path = nxfvars.get(
    "atlas_adata", "../../data/01_processed/sc_rna_seq/full_atlas_merged.h5ad"
)
artifact_dir = nxfvars.get("artifact_dir", "../../results/20_single-cell/subset_atlas")
threadpool_limits(nxfvars.get("cpus", 16))

# %%
adata = sc.read_h5ad(adata_path)

# %%
adata.shape

# %%
adata.obs["condition"].value_counts()

# %%
adata_cancer = adata[
    adata.obs["condition"].isin(["LUAD", "LUSC", "NSCLC_NOS"]), :
].copy()

# %%
adata_cancer

# %%
sc.pl.umap(adata_cancer, color="cell_type_coarse")

# %%
adata_m = adata_cancer[
    adata_cancer.obs["cell_type_coarse"].isin(["Macrophage/Monocyte", "cDC"]), :
].copy()

# %%
sc.pp.neighbors(adata_m, use_rep="X_scANVI")

# %%
sc.tl.umap(adata_m)

# %%
sc.tl.leiden(adata_m, resolution=0.5)

# %%
adata_cancer.write_h5ad(f"{artifact_dir}/adata_nsclc.h5ad")

# %%
adata_m.write_h5ad(f"{artifact_dir}/adata_m.h5ad")

# %%
