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
import infercnvpy as cnv

sc.set_figure_params(figsize=(4, 4))
import numpy as np

# %%
adata = sc.read_h5ad(
    "../../data/40_dorothea_progeny/40_primary_tumor_dorothea_progeny.h5ad"
)

# %%
adata_m = sc.read_h5ad(
    "../../data/50_myeloid_cells/20_cell_type_annotation/adata_annotated.h5ad"
)

# %%
pathways_of_interest = [
    f"PW:{x}" for x in ["JAK-STAT", "VEGF", "PI3K", "NFkB", "Trail", "TNFa"]
]
tfs_of_interest = [
    f"TF:{tf}" for tf in ["SPI1", "NFKB1", "STAT1", "MYC", "E2F4", "E2F2", "ZNF263"]
]

# %%
sc.pl.dotplot(
    adata,
    groupby="cell_type_coarse",
    var_names={
        "B cell": ["MS4A1"],
        "Ciliated": ["PIFO"],
        "Endothelial cell": ["CDH5", "VWF"],
        "Endothelial cell lymphatic": ["CCL21"],
        "Epithelial cell (benign)": ["AGER", "SFTPA1", "SFTPC", "CLDN18"],
        "Mast cell": ["TPSB2"],
        "Myeloid cell": ["CD14"],
        "NK cell": ["CD160", "KLRB1"],
        "NKT cell": ["KLRG1"],
        "Plasma cell": ["MZB1"],
        "Stromal cell": ["TAGLN", "COL1A1"],
        "T cell": ["CD3G"],
        "T cell CD4": ["CD4"],
        "T cell CD8": ["CD8A"],
        "T cell dividing": ["CDK1"],
        "T reg": ["FOXP3", "IL2RA"],
        "Tumor cell": ["EPCAM"],
        "pDC": ["IL3RA", "CLEC4C"],
    },
)

# %%
sc.pl.dotplot(
    adata_m,
    groupby="cell_type",
    var_names={
        "Macro CD163+": ["CD163", "CCL13"],
        "Macro MARCO+": ["MARCO", "FABP4", "PCOLCE2"],
        "Macro SLAMF9+": ["SLAMF9"],
        "Monocyte": ["FCN1", "VCAN"],
        "Myeloid dividing": ["CDK1"],
        "cDC1": ["CLEC9A"],
        "cDC2": ["CD1C", "CLEC10A"],
        "cDC2 CD1A+": ["CD1A", "CD207"],
        "mDC mature": ["CCR7"],
    },
)

# %%
adata_m

# %% [markdown]
# ## TF plots (all cell_types)

# %%
adata_progeny = sc.AnnData(
    X=adata.obsm["progeny"], obs=adata.obs.loc[:, ["cell_type_coarse"]]
)

# %%
sc.pl.matrixplot(
    adata_progeny,
    var_names=adata_progeny.var_names,
    groupby="cell_type_coarse",
    cmap="bwr",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
)

# %%
sc.pl.dotplot(
    adata_progeny,
    var_names=adata_progeny.var_names,
    groupby="cell_type_coarse",
    cmap="bwr",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
)

# %%
adata_dorothea = sc.AnnData(
    X=adata.obsm["dorothea"], obs=adata.obs.loc[:, ["cell_type_coarse"]]
)

# %%
sc.pl.matrixplot(
    adata_dorothea,
    var_names=tfs_of_interest,
    groupby="cell_type_coarse",
    cmap="bwr",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
)

# %%
sc.pl.dotplot(
    adata_dorothea,
    var_names=tfs_of_interest,
    groupby="cell_type_coarse",
    cmap="bwr",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
)

# %%
sc.pl.matrixplot(
    adata_dorothea,
    var_names=adata_dorothea.var_names,
    groupby="cell_type_coarse",
    cmap="bwr",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
)

# %%
sc.pl.dotplot(
    adata_dorothea,
    var_names=adata_dorothea.var_names,
    groupby="cell_type_coarse",
    cmap="bwr",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
)

# %% [markdown]
# ## TF plots

# %%
adata_progeny = sc.AnnData(
    X=adata_m.obsm["progeny"], obs=adata_m.obs.loc[:, ["cell_type"]]
)

# %%
sc.pl.matrixplot(
    adata_progeny,
    var_names=adata_progeny.var_names,
    groupby="cell_type",
    cmap="bwr",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
)

# %%
sc.pl.dotplot(
    adata_progeny,
    var_names=adata_progeny.var_names,
    groupby="cell_type",
    cmap="bwr",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
)

# %%
adata_dorothea = sc.AnnData(
    X=adata_m.obsm["dorothea"], obs=adata_m.obs.loc[:, ["cell_type"]]
)

# %%
sc.pl.dotplot(
    adata_dorothea,
    var_names=tfs_of_interest,
    groupby="cell_type",
    cmap="bwr",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
)

# %%
sc.pl.matrixplot(
    adata_dorothea,
    var_names=adata_dorothea.var_names,
    groupby="cell_type",
    cmap="bwr",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
)

# %%
sc.pl.dotplot(
    adata_dorothea,
    var_names=adata_dorothea.var_names,
    groupby="cell_type",
    cmap="bwr",
    swap_axes=True,
    vmin=-2.5,
    vmax=2.5,
)

# %% [markdown]
# ## Experiments with paga

# %%
adata_dc = adata_m[adata_m.obs["cell_type"].str.contains("DC"), :].copy()

# %%
sc.pp.neighbors(adata_dc, use_rep="X_scVI")

# %%
sc.tl.diffmap(adata_dc)

# %%
sc.pp.neighbors(adata_dc, use_rep="X_diffmap")

# %%
sc.tl.leiden(adata_dc, resolution=0.5)

# %%
sc.tl.umap(adata_dc)

# %%
sc.pl.umap(adata_dc, color=["leiden"])

# %%
sc.tl.paga(adata_dc, groups="leiden")

# %%
adata_dc.uns["iroot"] = np.flatnonzero(adata_dc.obs["cell_type"] == "cDC2")[0]
sc.tl.dpt(adata_dc)

# %%
sc.pl.paga(adata_dc, color=["cell_type", "dpt_pseudotime"])

# %%
sc.tl.umap(adata_dc, init_pos="paga")

# %%
sc.pl.umap(
    adata_dc, color=["cell_type", "CD1C", "CD1A", "dpt_pseudotime"], wspace=0.4, ncols=2
)

# %%
sc.tl.draw_graph(adata_dc, init_pos="paga")

# %%
sc.pl.draw_graph(
    adata_dc, color=["cell_type", "CD1C", "CD1A", "dpt_pseudotime"], wspace=0.4, ncols=2
)

# %%
sc.pl.paga(adata_dc, color="dpt_pseudotime")

# %%
