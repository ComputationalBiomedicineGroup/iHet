# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.4
# ---

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import scanpy as sc
from nxfvars import nxfvars
import pandas as pd
from threadpoolctl import threadpool_limits
import matplotlib.pyplot as plt
import decoupler as dc
from IPython.display import Markdown, display

sc.settings.set_figure_params(figsize=(4, 4))

# %%
path_adata_m = nxfvars.get(
    "adata_m",
    "/scratch/nsclc_heterogeneity/results/20_single_cell/22_annotate_myeloid/artifacts/adata_myeloid_reannotated.h5ad",
)
path_adata_nsclc = nxfvars.get(
    "adata_nsclc",
    "/scratch/nsclc_heterogeneity/results/20_single_cell/22_annotate_myeloid/artifacts/adata_nsclc_reannotated.h5ad",
)
cpus = nxfvars.get("cpus", 16)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/ihet")

# %%
threadpool_limits(cpus)

# %%
adata_m = sc.read_h5ad(path_adata_m)
adata_nsclc = sc.read_h5ad(path_adata_nsclc)

# %%
# exclude these datasets because they contain patients that were already measured in another dataset (same study)
# Having repeated measures would be problematic because the statistics we use do not account for this.
excludelist = [
    "Leader_Merad_2021_10x_5p_v1_beads",
    "Leader_Merad_2021_10x_5p_v1_CD2",
    "Leader_Merad_2021_10x_3p_v2_beads_cite",
]
adata_m = adata_m[~adata_m.obs["dataset"].isin(excludelist)]
adata_nsclc = adata_nsclc[~adata_nsclc.obs["dataset"].isin(excludelist)]

sc_datasets = {"NSCLC full dataset": adata_nsclc, "NSCLC myeloid subset": adata_m}


# %%
references = {
    "progeny": dc.get_progeny(organism="human", top=500),
    "collectri": dc.get_collectri(organism="human", split_complexes=False),
}

# %%
tfs_of_interest = [
    "E2F4",
    "HOPX",
    "IRF1",
    "KMT2E",
    "MYC",
    "NFKB1",
    "REL",
    "RELA",
    "RFXAP",
    "SPI1",
    "STAT1",
    "ZNF143",
]
pws_of_interest = [
    "EGFR",
    "Hypoxia",
    "JAK-STAT",
    "NFkB",
    "PI3K",
    "TGFb",
    "TNFa",
    "Trail",
    "VEGF",
    "WNT",
    "p53",
]

# %% [markdown]
# ## Single-cell approach

# %%
for name, tmp_ad in sc_datasets.items():
    # .raw contains log-norm values
    for method, net in references.items():
        dc.run_wmean(tmp_ad, net, use_raw=True, times=0, min_n=4, verbose=True)
        tmp_ad.obsm[method] = tmp_ad.obsm["wmean_estimate"]
        print(f"completed {method}")


# %%
ad_pw = dc.get_acts(adata_nsclc, "progeny")
ad_tf = dc.get_acts(adata_nsclc, "collectri")
sc.pp.scale(ad_pw)
sc.pp.scale(ad_tf)
fig, axs = plt.subplots(2, 4, figsize=(18, 8))
axs = iter(axs.flatten())
for pw, ax in zip(["TNFa", "NFkB", "TGFb"], axs):
    sc.pl.umap(
        ad_pw,
        color=pw,
        ax=ax,
        show=False,
        cmap="coolwarm",
        vmin=-2.5,
        vmax=2.5,
        title=f"PW: {pw}",
        frameon=False,
    )
for tf, ax in zip(["SPI1", "MYCN", "RFXAP", "E2F4", "ZNF143"], axs):
    sc.pl.umap(
        ad_tf,
        color=tf,
        ax=ax,
        show=False,
        cmap="coolwarm",
        vmin=-2.5,
        vmax=2.5,
        title=f"TF: {tf}",
        frameon=False,
    )
fig.savefig(
    f"{artifact_dir}/umap_cell_selected_features.svg",
    bbox_inches="tight",
    dpi=600,
)
plt.show()

# %% [markdown]
# ## Pseudobulk approach


# %%
pb_res = {}
for name, tmp_ad in sc_datasets.items():
    tmp_pb_res = {}
    pb_res[name] = tmp_pb_res

    display(Markdown(f"## {name}"))

    # Make pseudobulk
    pb = dc.get_pseudobulk(
        tmp_ad, sample_col="patient", groups_col="cell_type", mode="sum", min_cells=10, min_counts=1000
    )
    # keep order
    pb.obs["cell_type"] = pd.Categorical(
        pb.obs["cell_type"],
        categories=[x for x in adata_nsclc.obs["cell_type"].cat.categories if x in pb.obs["cell_type"].values],
    )

    sc.pp.normalize_total(pb, target_sum=1e6)
    sc.pp.log1p(pb)

    for method, net in references.items():
        dc.run_wmean(pb, net, use_raw=False, times=0, min_n=4)
        pb.obsm[method] = pb.obsm["wmean_estimate"]
        tmp_ad_decoupler = dc.get_acts(pb, method)
        tmp_ad_decoupler.layers["scaled"] = sc.pp.scale(tmp_ad_decoupler.X.copy())
        # tmp_ad.layers["batch_corrected"] = tmp_ad.X.copy()
        # sc.pp.regress_out(tmp_ad, keys="dataset", layer="batch_corrected")
        # sc.pp.scale(tmp_ad, layer="batch_corrected")
        tmp_pb_res[method] = tmp_ad_decoupler

    display(Markdown("### Progeny"))

    fig = sc.pl.matrixplot(
        tmp_pb_res["progeny"],
        var_names=pws_of_interest,
        groupby="cell_type",
        cmap="coolwarm",
        vmin=-2.5,
        vmax=2.5,
        layer="scaled",
        swap_axes="myeloid" in name,
        title="PW scaled",
        return_fig=True,
    )
    fig.savefig(f"{artifact_dir}/matrixplot_{name}_progeny.svg", bbox_inches="tight")

    display(Markdown("### collecTRI"))

    fig = sc.pl.matrixplot(
        tmp_pb_res["collectri"],
        var_names=tfs_of_interest,
        groupby="cell_type",
        cmap="coolwarm",
        vmin=-2.5,
        vmax=2.5,
        layer="scaled",
        swap_axes="myeloid" in name,
        title="TF scaled",
        return_fig=True,
    )
    fig.savefig(f"{artifact_dir}/matrixplot_{name}_collectri.svg", bbox_inches="tight")

# %%
# # %% [markdown]
# # ---
# #
# # ## Comparison of approaches
# # ### Metacell approach

# # %%
# adata_metacells = sc.read_h5ad("/scratch/nsclc_heterogeneity/_ihet_metacells/iHet_metcells.h5ad")
# fix_cell_type_names = lambda x: {
#     "Neutrophils": "Neutrophil",
#     "transitional club/AT2": "Transitional club/AT2",
#     "Tumor cells": "Tumor cell",
#     "stromal dividing": "Stromal cycling",
#     "B cell dividing": "B cell cycling",
#     "Plasma cell dividing": "Plasma cell cycling",
#     "T cell dividing": "T cell cycling",
#     "cDC dividing": "cDC cycling",
#     "Macrophage dividing": "Macrophage cycling",
#     "Macrophage alveolar dividing": "Macrophage alveolar cycling",
# }.get(x, x)
# adata_metacells.obs["major_celltype"] = adata_metacells.obs["major_celltype"].map(fix_cell_type_names)
# adata_metacells.obs["major_celltype"] = pd.Categorical(
#     adata_metacells.obs["major_celltype"],
#     categories=[
#         x for x in adata_nsclc.obs["cell_type"].cat.categories if x in adata_metacells.obs["major_celltype"].values
#     ],
# )

# # %%
# adata_metacells.X = adata_metacells.layers["total_umis"].copy()
# sc.pp.normalize_total(adata_metacells, target_sum=1e6)
# sc.pp.log1p(adata_metacells)
# metacell_res = {}

# # %%
# for method, net in references.items():
#     dc.run_wmean(adata_metacells, net, use_raw=False, times=0, min_n=5)
#     adata_metacells.obsm[method] = adata_metacells.obsm["wmean_estimate"]
#     tmp_ad = dc.get_acts(adata_metacells, method)
#     tmp_ad.layers["scaled"] = sc.pp.scale(tmp_ad.X.copy())
#     metacell_res[method] = tmp_ad


# # %% [markdown]
# # ### Progeny

# # %%
# sc.pl.matrixplot(
#     metacell_res["progeny"],
#     var_names=pws_of_interest,
#     groupby="major_celltype",
#     cmap="coolwarm",
#     vmin=-2.5,
#     vmax=2.5,
#     layer="scaled",
#     swap_axes=True,
#     title="PW scaled",
# )

# # %% [markdown]
# # ### CollecTRI

# # %%
# sc.pl.matrixplot(
#     metacell_res["collectri"],
#     var_names=tfs_of_interest,
#     groupby="major_celltype",
#     cmap="coolwarm",
#     vmin=-2.5,
#     vmax=2.5,
#     layer="scaled",
#     swap_axes=True,
#     title="TF scaled",
# )

# # %% [markdown]
# # ## Comparison of approaches
# # ### single cell -> decoupler -> pseudobulk

# # %%
# # Single cell -> decoupler -> pseudobulk
# singlecell_res = {}
# for name, tmp_ad in {
#     # "NSCLC full dataset": adata_nsclc,
#     "NSCLC myeloid subset": adata_m,
# }.items():
#     for method in references:
#         # Make pseudobulk
#         pb = dc.get_pseudobulk(
#             dc.get_acts(tmp_ad, method),
#             sample_col="patient",
#             groups_col="cell_type",
#             mode="mean",
#             min_cells=10,
#             # decoupler pseudobulk is mainly intended for counts - need to disable checks and set all filtering options to 0
#             skip_checks=True,
#             min_counts=0,
#             min_prop=0,
#             min_smpls=0,
#         ).copy()
#         # keep order
#         pb.obs["cell_type"] = pd.Categorical(
#             pb.obs["cell_type"],
#             categories=[x for x in adata_nsclc.obs["cell_type"].cat.categories if x in pb.obs["cell_type"].values],
#         )
#         pb.layers["scaled"] = sc.pp.scale(pb.X.copy())
#         singlecell_res[method] = pb

#     display(Markdown(f"### Progeny"))

#     sc.pl.matrixplot(
#         singlecell_res["progeny"],
#         var_names=pws_of_interest,
#         groupby="cell_type",
#         cmap="coolwarm",
#         vmin=-2.5,
#         vmax=2.5,
#         layer="scaled",
#         swap_axes="myeloid" in name,
#         title="PW scaled",
#     )

#     display(Markdown("### collecTRI"))

#     sc.pl.matrixplot(
#         singlecell_res["collectri"],
#         var_names=tfs_of_interest,
#         groupby="cell_type",
#         cmap="coolwarm",
#         vmin=-2.5,
#         vmax=2.5,
#         layer="scaled",
#         swap_axes="myeloid" in name,
#         title="TF scaled",
#     )

# # %% [markdown]
# # ### correlations

# # %%
# import altair as alt  # noqa

# alt.data_transformers.enable("vegafusion")

# for tool in references:

#     ad_decoupler_pseudobulk = sc.get.aggregate(singlecell_res[tool], by="cell_type", func="mean", layer="scaled")
#     ad_pseudobulk_decoupler = sc.get.aggregate(
#         pb_res["NSCLC myeloid subset"][tool], by="cell_type", func="mean", layer="scaled"
#     )
#     ad_metacell_decoupler = sc.get.aggregate(metacell_res[tool], by="major_celltype", func="mean", layer="scaled")

#     agg_df = pd.concat(
#         [
#             ad_decoupler_pseudobulk.to_df(layer="mean").assign(approach="SC -> decoupler -> pseudobulk"),
#             ad_pseudobulk_decoupler.to_df(layer="mean").assign(approach="SC -> pseudobulk -> decoupler"),
#             ad_metacell_decoupler.to_df(layer="mean").assign(approach="SC -> metacell -> decoupler"),
#         ]
#     )
#     agg_df.index.name = "cell_type"
#     plot_df = (
#         agg_df.reset_index()
#         .melt(id_vars=["cell_type", "approach"])
#         .pivot(columns="approach", index=["cell_type", "variable"], values="value")
#         .reset_index()
#     )

#     # correlations
#     from scipy.stats import pearsonr

#     cor_singlecell_vs_pseudobulk = pearsonr(
#         plot_df["SC -> decoupler -> pseudobulk"], plot_df["SC -> pseudobulk -> decoupler"]
#     )
#     cor_singlecell_vs_metacells = pearsonr(
#         plot_df["SC -> decoupler -> pseudobulk"], plot_df["SC -> metacell -> decoupler"]
#     )
#     cor_pseudobulk_vs_metacells = pearsonr(
#         plot_df["SC -> pseudobulk -> decoupler"], plot_df["SC -> metacell -> decoupler"]
#     )

#     # single cell vs. pseudobulk
#     c1 = alt.Chart(plot_df).encode(
#         x="SC -> decoupler -> pseudobulk", y="SC -> pseudobulk -> decoupler", color="cell_type"
#     ).mark_point() + alt.Chart({"values": [{}]}).mark_text(align="left", baseline="top").encode(
#         x=alt.value(5),  # pixels from left
#         y=alt.value(5),  # pixels from top
#         text=alt.value(f"r={cor_singlecell_vs_pseudobulk[0]:.3f}"),
#     )

#     # single cell vs. metacell
#     c2 = alt.Chart(plot_df).encode(
#         x="SC -> decoupler -> pseudobulk", y="SC -> metacell -> decoupler", color="cell_type"
#     ).mark_point() + alt.Chart({"values": [{}]}).mark_text(align="left", baseline="top").encode(
#         x=alt.value(5),  # pixels from left
#         y=alt.value(5),  # pixels from top
#         text=alt.value(f"r={cor_singlecell_vs_metacells[0]:.3f}"),
#     )

#     # pseudobulk vs metacell
#     c3 = alt.Chart(plot_df).encode(
#         x="SC -> pseudobulk -> decoupler", y="SC -> metacell -> decoupler", color="cell_type"
#     ).mark_point() + alt.Chart({"values": [{}]}).mark_text(align="left", baseline="top").encode(
#         x=alt.value(5),  # pixels from left
#         y=alt.value(5),  # pixels from top
#         text=alt.value(f"r={cor_pseudobulk_vs_metacells[0]:.3f}"),
#     )

#     (c1 | c2 | c3).properties(title=tool).show()

# # %%
