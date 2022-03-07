# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Python [conda env:.conda-2020-organoids-scanpy]
#     language: python
#     name: conda-env-.conda-2020-organoids-scanpy-py
# ---

# %%
# %load_ext autoreload
# %autoreload 2
import threadpoolctl
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
from scanpy_helpers.annotation import AnnotationHelper
from warnings import filterwarnings

sc.set_figure_params(figsize=(5, 5))
filterwarnings("ignore", category=FutureWarning)

# %%
adata = sc.read_h5ad("/data/datasets/SadeFeldman_2018/h5ad_tpm/SadeFeldman_2018.h5ad")

# %%
sns.displot(adata.X.data)

# %%
adata = adata[~adata.obs_names.str.contains("T_enriched"), :]

# %%
adata.obs["sample"] = [
    f"{patient}_{timepoint}"
    for patient, timepoint in zip(adata.obs["patient_id"], adata.obs["timepoint"])
]

# %%
adata.obs

# %%
sc.pp.highly_variable_genes(
    adata, flavor="cell_ranger", batch_key="sample", n_top_genes=3000
)

# %%
sc.tl.pca(adata)

# %%
sc.pp.neighbors(adata)

# %%
sc.tl.umap(adata)

# %%
sc.tl.leiden(adata)

# %%
sc.pl.umap(adata, color=["patient_id", "sample", "therapy", "timepoint", "response"])

# %%
ah = AnnotationHelper()

# %%
ah.plot_umap(adata)

# %%
sc.pl.umap(adata, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "B cell": [4],
    "Plasma cell": [9],
    "pDC": [10],
    "myeloid": [7, 8],
    "T dividing": [6],
    "T CD8": [0, 1, 5],
    "T reg": [2],
    "T CD4": [3],
}

# %%
ah.annotate_cell_types(adata, ct_map)

# %%
adata.obs

# %%
adata.write_h5ad('../../data/50_myeloid_cells/sade-feldman-melanoma.h5ad')

# %%
tmp_adata = adata[
    (adata.obs["therapy"] == "anti-PD1") & (adata.obs["timepoint"] == "Pre"), :
].copy()

# %%
cells_per_patient = (
    tmp_adata.obs.groupby("patient_id").size().reset_index(name="cell_count")
)

# %%
myeloid_cells_per_patient = (
    tmp_adata.obs.loc[tmp_adata.obs["cell_type"] == "myeloid", :]
    .groupby("patient_id")
    .size()
    .reset_index(name="myeloid_cell_count")
)

# %%
t_cells_per_patient = (
    tmp_adata.obs.loc[tmp_adata.obs["cell_type"].str.startswith("T "), :]
    .groupby("patient_id")
    .size()
    .reset_index(name="t_cell_count")
)

# %%
b_cells_per_patient = (
    tmp_adata.obs.loc[tmp_adata.obs["cell_type"].str.startswith("B "), :]
    .groupby("patient_id")
    .size()
    .reset_index(name="b_cell_count")
)

# %%
patient_table = (
    tmp_adata.obs.loc[:, ["patient_id", "response"]]
    .drop_duplicates()
    .reset_index(drop=True)
    .set_index("patient_id")
    .join(cells_per_patient.set_index("patient_id"))
    .join(myeloid_cells_per_patient.set_index("patient_id"))
    .join(t_cells_per_patient.set_index("patient_id"))
    .join(b_cells_per_patient.set_index("patient_id"))
)

# %%
patient_table["myeloid_ratio"] = (
    patient_table["myeloid_cell_count"] / patient_table["cell_count"]
)
patient_table["t_ratio"] = patient_table["t_cell_count"] / patient_table["cell_count"]
patient_table["b_ratio"] = patient_table["b_cell_count"] / patient_table["cell_count"]

# %%
sns.scatterplot(data=patient_table, x="myeloid_ratio", y="t_ratio")

# %%
sns.boxplot(data=patient_table, y="myeloid_ratio", x="response")
ax = sns.swarmplot(data=patient_table, y="myeloid_ratio", x="response", color="black")
ax.set_ylabel("ratio myeloid/total cells")

# %%
sns.boxplot(data=patient_table, y="t_ratio", x="response")
ax = sns.swarmplot(data=patient_table, y="t_ratio", x="response", color="black")
ax.set_ylabel("ratio T/total cells")

# %%
sns.boxplot(data=patient_table, y="b_ratio", x="response")
ax = sns.swarmplot(data=patient_table, y="b_ratio", x="response", color="black")
ax.set_ylabel("ratio B/total cells")
