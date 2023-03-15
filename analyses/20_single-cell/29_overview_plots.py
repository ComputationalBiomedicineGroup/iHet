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

# %% [markdown]
# # Overview plots and statistics

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import scanpy as sc
from nxfvars import nxfvars
import pandas as pd
from threadpoolctl import threadpool_limits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import altair as alt

sc.settings.set_figure_params(figsize=(4, 4))
import scanpy_helpers as sh

# %%
path_adata_m = nxfvars.get(
    "adata_m",
    "../../data/results/20_single_cell/22_annotate_myeloid/artifacts/adata_myeloid_reannotated.h5ad",
)
path_adata_nsclc = nxfvars.get(
    "adata_nsclc",
    "../../data/results/20_single_cell/22_annotate_myeloid/artifacts/adata_nsclc_reannotated.h5ad",
)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/ihet/single-cell")

# %%
# !mkdir -p {artifact_dir}

# %%
adata_m = sc.read_h5ad(path_adata_m)
adata_nsclc = sc.read_h5ad(path_adata_nsclc)

# %% [markdown]
# ## Dataset stats

# %%
adata_nsclc.obs["study"].unique().tolist()

# %%
adata_nsclc.shape

# %%
adata_nsclc.obs["patient"].nunique()

# %%
adata_nsclc.obs["study"].nunique()

# %%
adata_nsclc.obs["dataset"].nunique()

# %%
adata_nsclc.obs["cell_type_coarse"].nunique()

# %%
adata_nsclc.obs["cell_type_coarse"].unique().tolist()

# %%
adata_nsclc.obs["cell_type"].nunique()

# %% [markdown]
# ## Cell-type stats

# %%
n_myel = np.sum(
    adata_nsclc.obs["cell_type_coarse"].isin(["Macrophage/Monocyte", "cDC"])
)
n_myel, n_myel / adata_nsclc.shape[0]

# %%
n_t = np.sum(adata_nsclc.obs["cell_type_coarse"].isin(["T cell", "NK cell"]))
n_t, n_t / adata_nsclc.shape[0]

# %%
adata_nsclc.obs.groupby("cell_type_coarse").size().reset_index(name="n").assign(
    frac=lambda x: x["n"] / np.sum(x["n"])
).sort_values("frac", ascending=False)

# %% [markdown]
# ## UMAP plots

# %%
with plt.rc_context({"figure.figsize": (7, 7)}):
    sc.pl.umap(
        adata_nsclc,
        color="cell_type_tumor",
        legend_loc="on data",
        legend_fontoutline=2,
        size=500000 / adata_nsclc.shape[0],
    )

# %%
# %%capture
for scope, tmp_adata in {"m": adata_m, "nsclc": adata_nsclc}.items():
    size = 500000 / tmp_adata.shape[0]
    with plt.rc_context({"figure.figsize": (7, 7), "figure.dpi": 300}):
        fig = sc.pl.umap(
            tmp_adata,
            color="cell_type_coarse",
            legend_loc="on data",
            legend_fontoutline=2,
            size=size,
            show=False,
            frameon=False,
            return_fig=True,
        )
        fig.savefig(
            f"{artifact_dir}/umap_cell_type_coarse_{scope}.svg",
            bbox_inches="tight",
            dpi=1200,
        )

        fig = sc.pl.umap(
            tmp_adata,
            color="cell_type",
            legend_loc="on data",
            legend_fontoutline=0.5,
            size=size,
            legend_fontsize=7,
            show=False,
            frameon=False,
            return_fig=True,
        )
        fig.savefig(
            f"{artifact_dir}/umap_cell_type_fine_{scope}.svg",
            bbox_inches="tight",
            dpi=1200,
        )

        for color in ["platform", "condition", "dataset"]:
            fig = sc.pl.umap(
                tmp_adata,
                color=color,
                size=size,
                show=False,
                return_fig=True,
                frameon=False,
            )
            fig.savefig(
                f"{artifact_dir}/umap_{color}_{scope}.svg", bbox_inches="tight", dpi=600
            )

# %% [markdown]
# ## Dotplots

# %%
sc.pl.dotplot(
    adata_m,
    groupby="cell_type",
    var_names={
        "Macrophage SLAMF9/SPP1": ["SLAMF9", "SPP1", "MIF"],
        "Macrophage M1": ["CD163", "RNASE1", "CCL13"],
        "Macrophage M2": [
            "MARCO",
            "RBP4",
            "PPARG",
            "FABP4",
            "RBP4",
            "PCOLCE2",
            "CXCL9",
            "CXCL10",
            "CXCL11",
        ],
        "Macrophaghe": ["TREM2"],
        "DC mature": ["CCR7", "CCR7", "CD40", "RELB", "CD83", "CD274", "CD200"],
        "cDC1": ["CLEC9A", "XCR1", "IRF8", "BATF3"],
        "cDC2": ["CD1C", "FCER1A", "FCGR2B", "AXL"],
        "cDC2 CD1A+": ["CD1C", "CD1A", "CD207", "IL18"],
        "Monocytes": ["VCAN", "CD14", "FCGR3A", "LST1"],
    },
)

# %%
fig = sc.pl.dotplot(
    adata_m,
    groupby="cell_type_macro",
    var_names={
        "DC mature": ["CCR7", "CCL22", "CCR7", "CD40", "RELB", "LAMP3"],
        "Macrophages": ["APOE", "C1QB", "TREM2"],
        "Macrophage alveolar": ["FABP4"],
        "cDC1": ["CLEC9A", "XCR1", "IRF8", "CD274", "MS4A1"],
        "cDC2": ["CD1C", "FCER1A", "CLEC10A"],
        "cDC2 CD1A-hi": ["CD1A", "CD207"],
        "Monocytes": ["FCN1"],
        "Monocyte contentional": ["S100A12", "CD14", "VCAN"],
        "Monocyte non-conventional": ["FCGR3A", "LST1"],
        "dividing": ["CDK1", "MKI67"],
    },
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/myeloid_dotplot.svg")

# %%
sc.pl.umap(
    adata_nsclc,
    color=["CD14", "CD68", "CD163", "CSF1R", "FCGR3A", "CD84", "MS4A4A", "MS4A2"],
)

# %%
marker_dict = {
    "Alveolar cell type 1": ["AGER", "CLDN18"],
    "Alveolar cell type 2": ["SFTPC", "SFTPB", "SFTPA1"],
    "B cell": ["CD19", "CD79A", "MS4A1"],
    "Ciliated": ["PIFO", "FOXJ1", "HYDIN", "CFAP299"],
    "Club": ["SCGB3A1", "SCGB3A2"],
    "Endothelial cell": ["VWF", "CDH5", "SELE"],
    "Endothelial cell lymphatic": ["CCL21"],
    "Fibroblast": ["PDGFRA", "FAP", "COL1A1"],
    "Fibroblast adventitial": ["MFAP5", "SCARA5"],
    "Fibroblast alveolar": ["ITGA8", "SCN7A"],
    "Epithelial cell": ["EPCAM"],
    "Mast cells": ["TPSB2"],
    "Mesothelial": ["MSLN", "CALB2"],
    "Monocytic lineage": ["CD14", "CD68"],
    "NK cell": ["KLRD1", "GNLY"],
    "Neutrophils": ["FCGR3B", "CSF3R"],
    "pDC": ["IL3RA", "CLEC4C"],
    "Pericyte": ["COX4I2", "PDGFRB"],
    "Plasma cell": ["SDC1", "MZB1"],
    "Smooth muscle cell": ["TAGLN", "MYH11"],
    "T cell": ["CD3E"],
    "T cell CD4": ["CD4"],
    "T cell CD8": ["CD8A"],
    "T cell regulatory": ["FOXP3", "IL2RA", "CTLA4"],
    "Dividing": ["MKI67", "CDK1"],
    "Tumor cells LUAD": ["KRT7", "CD24"],
    "Tumor cells LUSC": ["SOX2", "KRT17"],
}
fig = sc.pl.dotplot(
    adata_nsclc, groupby="cell_type", var_names=marker_dict, return_fig=True
)
fig.savefig(f"{artifact_dir}/nsclc_dotplot.svg")

# %% [markdown]
# ## Patients per myeloid cell-type

# %%
m_counts = (
    adata_m.obs.groupby(["cell_type_macro", "dataset", "study", "patient"], observed=True)
    .size()
    .reset_index(name="n_cells")
    .groupby(["cell_type_macro", "dataset"])
    .apply(lambda x: x.assign(n_cells_cell_type_dataset=lambda k: k["n_cells"].sum()))
    .groupby("patient")
    .apply(lambda x: x.assign(n_cells_patient=lambda k: k["n_cells"].sum()))
).query("n_cells_patient >= 10")

# %%
patient_cell_type_combs = (
    m_counts.loc[:, ["dataset", "study", "patient"]]
    .merge(m_counts.loc[:, ["cell_type_macro"]], how="cross")
    .drop_duplicates()
)

# %%
tmp_df = m_counts.merge(
    patient_cell_type_combs,
    on=["dataset", "study", "patient", "cell_type_macro"],
    how="outer",
)

# %%
tmp_df2 = (
    tmp_df.groupby(["study", "cell_type_macro"], observed=True)
    .agg(n_cells=("n_cells", np.sum))
    .reset_index()
)

# %%
tmp_df2_study = (
    m_counts.groupby("study", observed=True)
    .agg(n_patients_ge_10_cells=("patient", lambda x: x.nunique()))
    .reset_index()
    .assign(x="#patients with ≥ 10 myeloid cells")
)

# %%
heatmp = (
    alt.Chart(tmp_df2)
    .mark_rect()
    .encode(
        x="cell_type_macro",
        y=alt.Y("study", axis=None),
        color=alt.Color("n_cells", scale=alt.Scale(scheme="inferno", reverse=True)),
    )
)
txt = (
    alt.Chart(tmp_df2)
    .mark_text()
    .encode(
        x="cell_type_macro",
        y=alt.Y("study", axis=None),
        text="n_cells",
        color=alt.condition(
            alt.datum.n_cells < 10000, alt.value("black"), alt.value("white")
        ),
    )
)
studies = (
    alt.Chart(tmp_df2_study)
    .mark_rect()
    .encode(
        y=alt.Y("study"),
        x=alt.X("x", title=None),
        color=alt.Color("study", scale=sh.colors.altair_scale("study"), legend=None),
    )
)
studies_txt = (
    alt.Chart(tmp_df2_study)
    .mark_text()
    .encode(y="study", x=alt.X("x", title=None), text="n_patients_ge_10_cells")
)

# %%
ch = ((studies + studies_txt) | (heatmp + txt).properties(width=500)).configure_concat(
    spacing=0
)
ch.display()

# %%
