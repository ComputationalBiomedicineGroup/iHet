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
# %load_ext autoreload
# %autoreload 2
import scanpy as sc
import scvi
import numpy as np
import pandas as pd
import scanpy_helpers.de as de
from scanpy_helpers.annotation import AnnotationHelper
from warnings import filterwarnings
from pygenesig.gini import GiniSignatureGenerator
import itertools
import matplotlib.pyplot as plt

sc.settings.verbosity = 2

filterwarnings("ignore", category=FutureWarning)
sc.set_figure_params(figsize=(5, 5))
import re

# %%
scvi.__version__

# %%
output_dir = "../../data/50_myeloid_cells/30_de_analysis"
deploy_dir = "../../deploy/Ab1MxvlS9c56BaA0Z94J/figures/"

# %%
adata = sc.read_h5ad(
    "../../data/50_myeloid_cells/20_cell_type_annotation/adata_annotated.h5ad"
)
adata_scvi = sc.read_h5ad(
    "../../data/50_myeloid_cells/01_myeloid_subset/adata_scvi.h5ad"
)

# %% [markdown]
# The two anndata objects are the same, except `adata` contains log-normalized counts and the cell-type annotations, whereas `adata_scvi` just contains the raw counts

# %%
assert adata.shape == adata_scvi.shape

# %%
adata.obs["cell_type"] = de._make_names(adata.obs["cell_type"])
adata._sanitize()

# %%
sc.pl.umap(adata, color="cell_type")

# %% [markdown]
# # DE analysis
#
# According to [this benchmark](https://pubmed.ncbi.nlm.nih.gov/29481549/), edgeR and MAST are among the best methods for 
# DE analysis of single cells. In any case, the authors recommend to include the number of detected genes into the model. 
# Moreover, the model needs to correct for batch effects arising from different patients and dataset:
#
#       `~ 1 + cell_type + sample + n_genes`. 
#
# I tried different approaches for DE analysis:
#
# * `edgeR` way too slow given the number of cells and categories. Aborted after 72h. 
# * `glmGamPoi` claims to be a faster version of edge R. Completed after 6 hours on a subset of 500 cells. Still way too slow. 
# * `scVI` provides a method for batch-corrected differential expression analysis based on the same autoencode used for integrating the data. It is very fast (few minutes) and their [preprint claims](https://www.biorxiv.org/content/10.1101/794289v1) superior performance over MAST and edgeR. 
# * `MAST`, while being slower than `scVI` is also based on a linear model. It can be massively parallelized and completed on the entire dataset in 2:30h on 44 cores. 
# * `Gini-Index` is not a DE-appraoch *per se*, but a useful approach to identify marker genes based on information gain
#
# In the following I will show the results of `scVI` and `MAST` and `Gini-Index`. 
#
#
# ## Results summary
#  * big overlap between scVI and gini. 
#  * MAST fold changes don't rank marker genes for specific clusters very high. Top hits tend to be markers for multiple clusters (e.g. Macrophages rather than M1/M2/M0). 

# %% [markdown]
# ## Run the DE methods
# ### scVI

# %%
adata_scvi.obs["cell_type"] = adata.obs["cell_type"]
scvi.data.setup_anndata(adata_scvi, batch_key="sample")
model = scvi.model.SCVI.load(
    f"../../data/50_myeloid_cells/20_cell_type_annotation/scvi_model_all_genes",
    adata_scvi,
)

# try to load from cache
try:
    scvi_res = pd.read_csv(f"{output_dir}/scvi_res.csv", index_col=0)
except FileNotFoundError:
    scvi_res = de.scvi(adata_scvi, model, groupby="cell_type")
    scvi_res.to_csv(f"{output_dir}/scvi_res.csv")

scvi_res.index.name = "gene_symbol"
scvi_res.reset_index(inplace=True)

# %% [markdown]
# ### MAST

# %%
# %%time
try:
    mast_res = pd.read_csv(f"{output_dir}/mast_res.csv", index_col=0)
except FileNotFoundError:
    mast_res = de.mast(
        adata,
        groupby="cell_type",
        groups="all",
        cofactors=["sample", "n_genes"],
        n_jobs=16,
        n_cores_per_job=4,
    )
    mast_res.to_csv(f"{output_dir}/mast_res.csv")

# %% [markdown]
# ### Gini-index

# %%
# %%time
# run signature generator without cutoff, can threshold later
sg = GiniSignatureGenerator(
    adata.X.toarray().astype(float).T,
    target=adata.obs["cell_type"].values,
    min_expr=0,
    min_gini=0,
    aggregate_fun=np.mean,
    max_rk=1,
    max_rel_rk=1,
)
gini_res = sg.get_rogini_format()

# %% [markdown]
# ## DE Results

# %% [markdown]
# ### scVI
#
# * scVI uses logBayesFactor as main output statistic. The log bayes factors are to be interpreted as described [here](https://docs.scvi-tools.org/en/stable/user_guide/notebooks/scVI_DE_worm.html#Interpreting-Bayes-factors). `proba_not_de` can be interpreted as a FDR-corrected p-value. 
# * a lot of the top hits that are very sparsely expressed in one cell-type and not at all in the other cell-types. While the signal is probably real, these genes are still bad candidate markers as they are hardly expressed. Below is an example of "TEX11" which is expressed in 8 cells in total - 7 of them are DC Lagerhans. 
# * Alternatively, I rank genes by their log-fold change require the minimum expression to be >= 0.5

# %% [markdown]
# #### TEX11 is a top hit, but a bad marker

# %%
adata.obs["_TEX11"] = (adata[:, "TEX11"].X.todense().A1 > 0).astype(str)

# %%
sc.pl.umap(
    adata,
    color="_TEX11",
    groups="True",
    size=[20 if x == "True" else 5 for x in adata.obs["_TEX11"]],
)

# %%
scvi_res

# %% [markdown]
# #### Results based on alternative scoring strategy
#  * rank by log-fold change, require minimal expression level

# %%
scvi_res["score"] = scvi_res["lfc_mean"] * (scvi_res["raw_normalized_mean1"] > 0.5)
# de_res["score"] = de_res["bayes_factor"] * (de_res["raw_normalized_mean2"] - de_res["raw_normalized_mean1"]) * -1

# %%
scvi_res["cell_type"] = [x.split("vs")[0].strip() for x in scvi_res["comparison"]]

# %%
de.de_res_to_anndata(
    adata,
    scvi_res,
    groupby="cell_type",
    gene_id_col="gene_symbol",
    pval_col="proba_not_de",
    pval_adj_col="proba_not_de",
    lfc_col="lfc_mean",
)

# %% [markdown]
# #### Top genes for each group
# score = log fold change

# %%
sc.pl.rank_genes_groups(adata, sharey=False)

# %%
ax = sc.pl.rank_genes_groups_dotplot(adata, show=False)
fig = ax["mainplot_ax"].get_figure()
fig.savefig(f"{deploy_dir}/dotplot_myeloid_markers_scvi.pdf", bbox_inches="tight")

# %% [markdown]
# On first glance, the markers make sense: 
#
# * CCR7 is the marker that was used for classifying mDC mature in the first place
# * CDK1 is a marker for dividing cells
# * Visual inspection of the top hits on UMAP are clearly restricted to a subset of cells (see "Plots" section at the end of the notebook) 

# %% [markdown]
# ### MAST
#
#  * Due to the immense statistical power arising from 15k cells, almost all differences between cells become statistically significant. In fact many comparisons hava a p-value of 0. 
#  * Ranking by p-value (only) is therefore not viable to find marker genes. Ranking by log-fold change instead
#  * Top MAST hits appear to be specific for several clusters, not a single one. 
#
#
# #### APOE is a top hit for Macrophage M0, but also M1 and M2

# %%
sc.pl.umap(adata, color=["APOE", "cell_type"])

# %%
# replace 0-pvalues with minimum nonzero pvalue.
mast_res["cell_type"] = mast_res["comparison"]
mast_res.loc[mast_res["coef"].isnull(), "coef"] = 0
mast_res.loc[mast_res["Pr(>Chisq)"] == 0, "Pr(>Chisq)"] = np.min(
    mast_res["Pr(>Chisq)"][mast_res["Pr(>Chisq)"] != 0]
)
mast_res["score"] = mast_res["coef"]

# %%
mast_res.sort_values("score")

# %%
de.de_res_to_anndata(
    adata,
    mast_res,
    groupby="cell_type",
    gene_id_col="primerid",
    pval_col="Pr(>Chisq)",
    lfc_col="coef",
    score_col="coef",
)

# %% [markdown]
# #### Top genes for each group
# score = log fold change

# %%
sc.pl.rank_genes_groups(adata, sharey=False)

# %%
sc.tl.dendrogram(adata, groupby="cell_type")
sc.pl.rank_genes_groups_dotplot(adata)

# %% [markdown]
# ### Gini-Index

# %%
gini_res = gini_res.sort_values(["CATEGORY", "GINI_IDX"], ascending=[True, False])

# %%
gini_res["GENE_SYMBOL"] = [adata.var_names[i] for i in gini_res["GENEID"]]

# %%
# merge MAST and gini results
gini_mast_df = gini_res.merge(
    mast_res,
    how="inner",
    left_on=["GENE_SYMBOL", "CATEGORY"],
    right_on=["primerid", "cell_type"],
)

# %% [markdown]
#  * VALUE = average gene expression in CATEGORY
#  * Pr(>Chisq) = p-value for that gene from the MAST analysis (CATEGORY vs rest) 
#  * coef = MAST log fold change

# %%
gini_mast_df

# %%
gini_res

# %%
gini_mast_df.to_csv(f"{output_dir}/gini_mast_res.csv")

# %%
gini_df_filtered = gini_mast_df.loc[
    (gini_mast_df["VALUE"] > 0.1)
    & (gini_mast_df["GINI_IDX"] > 0.6)
    & (gini_mast_df["Pr(>Chisq)"] < 0.01 / gini_mast_df.shape[0]),
    :,
]

# %%
gini_signatures = {
    cat: df.sort_values("VALUE", ascending=False)["GENE_SYMBOL"].values[:10]
    for cat, df in gini_df_filtered.groupby("CATEGORY")
}

# %%
var_group_labels = list(gini_signatures.keys())
c = -1
var_group_positions = [
    (c := c + 1, (c := c + len(group) - 1)) for group in gini_signatures.values()
]

# %%
ax = sc.pl.dotplot(
    adata,
    groupby="cell_type",
    var_names=list(itertools.chain.from_iterable(gini_signatures.values())),
    var_group_labels=var_group_labels,
    var_group_positions=var_group_positions,
    show=False
)
fig = ax["mainplot_ax"].get_figure()
fig.savefig(f"{deploy_dir}/dotplot_myeloid_markers_gini.pdf", bbox_inches="tight")

# %% [markdown]
# ### Create dataframe with all results

# %%
de_all_results = (
    gini_mast_df.merge(
        scvi_res.rename(columns={"score": "scvi_score"}).drop(columns=["comparison"]),
        left_on=["GENE_SYMBOL", "cell_type"],
        right_on=["gene_symbol", "cell_type"],
    )
    .drop(
        columns=[
            "primerid",
            "GENEID",
            "RANKING",
            "comparison",
            "CATEGORY",
            "gene_symbol",
            "proba_de",
            "scale1",
            "scale2",
            "lfc_median",
            "lfc_std",
            "lfc_min",
            "lfc_max",
            "is_de_fdr_0.05",
        ]
    )
    .rename(
        columns={
            "VALUE": "mean_expr_in_cell_type",
            "Pr(>Chisq)": "mast_pvalue",
            "coef": "mast_coef",
            "score": "mast_score",
            "proba_not_de": "scvi_adjusted_pvalue",
            "bayes_factor": "scvi_log_bayes_factor",
            "lfc_mean": "scvi_lfc_mean",
            "raw_mean1": "scvi_raw_mean1",
            "raw_mean2": "scvi_raw_mean2",
            "non_zeros_proportion1": "scvi_non_zeros_proportion1",
            "non_zeros_proportion2": "scvi_non_zeros_proportion2",
            "raw_normalized_mean1": "scvi_raw_normalized_mean1",
            "raw_normalized_mean2": "scvi_raw_normalized_mean2",
            "GENE_SYMBOL": "gene_symbol"
        }
    )
)

# %%
de_all_results.to_csv(f"{output_dir}/de_all_results.csv")
de_all_results.to_csv(f"{deploy_dir}/de_all_results.csv")

# %% [markdown]
# ### Compare gini against scVI

# %%
top10_scvi = set(
    list(
        scvi_res.sort_values("score", ascending=False)
        .groupby("cell_type")
        .apply(lambda x: x["gene_symbol"][:10])
        .reset_index()["gene_symbol"]
        .values
    )
)
top_gini = set(list(itertools.chain.from_iterable(gini_signatures.values())))
top_both = top10_scvi & top_gini
rest = set(de_all_results["gene_symbol"]) - top10_scvi - top_gini

# %% [markdown]
# Scatterplot of gini index against scVI Bayes Factor. Genes with a high gini index also have a high bayes factor. 

# %%
fig, ax = plt.subplots()
df_background = de_all_results.loc[de_all_results["gene_symbol"].isin(rest), :]
df_top10_scvi = de_all_results.loc[de_all_results["gene_symbol"].isin(top10_scvi - top_both), :]
df_top_gini = de_all_results.loc[de_all_results["gene_symbol"].isin(top_gini - top_both), :]
df_top_both = de_all_results.loc[de_all_results["gene_symbol"].isin(top_both), :]
background = ax.scatter(
    df_background["GINI_IDX"], df_background["scvi_log_bayes_factor"], c="lightgrey", s=3
)
background.set_label("other")
tmp_sc = ax.scatter(
    df_top_gini["GINI_IDX"], df_top_gini["scvi_log_bayes_factor"], c="#7570b3", s=5
)
tmp_sc.set_label("gini signature")
tmp_sc = ax.scatter(
    df_top10_scvi["GINI_IDX"], df_top10_scvi["scvi_log_bayes_factor"], c="#1b9e77", s=5
)
tmp_sc.set_label("top 10 scVI")
tmp_sc = ax.scatter(
    df_top_both["GINI_IDX"], df_top_both["scvi_log_bayes_factor"], c="#d95f02", s=5
)
tmp_sc.set_label("both")
red_line = ax.vlines(x=[0.6], ymin=0, ymax=5, color="red", lw=1)
red_line.set_label("min gini")
ax.legend(bbox_to_anchor=(1.05, 1))
ax.set_title("scVI vs gini index")
ax.set_xlabel("gini index")
ax.set_ylabel("scVI Bayes Factor")
ax.grid(False)

# %% [markdown]
# The same plot colored by the epression cutoff used for finding gini signature genes. 
# Genes witha high gini index tend to be expressed at a lower level. 

# %%
fig, ax = plt.subplots()
df_background = de_all_results.loc[de_all_results["mean_expr_in_cell_type"] < 0.1, :]
df_expr = de_all_results.loc[de_all_results["mean_expr_in_cell_type"] >= 0.1, :]
tmp_sc = ax.scatter(
    df_background["GINI_IDX"], df_background["scvi_log_bayes_factor"], c="lightgrey", s=3
)
tmp_sc.set_label("mean expr < 0.1")
tmp_sc = ax.scatter(df_expr["GINI_IDX"], df_expr["scvi_log_bayes_factor"], c="blue", s=5)
tmp_sc.set_label("mean expr >= 0.1")
red_line = ax.vlines(x=[0.6], ymin=0, ymax=5, color="red", lw=1)
red_line.set_label("min gini")
ax.legend(bbox_to_anchor=(1.05, 1))
ax.set_title("scVI vs gini index")
ax.set_xlabel("gini index")
ax.set_ylabel("scVI Bayes Factor")
ax.grid(False)

# %% [markdown]
# ## Plot marker genes
#

# %% [markdown]
# ### Cell-type markers
# Top markers tend to be in the top hits of both gini and scVI

# %%
print("DC Lagerhans")
sc.pl.umap(adata, color=["LTB", "FCER1A"], cmap="inferno")

print("MDSC")
sc.pl.umap(adata, color=["APOBEC3A", "LILRB2"], cmap="inferno")

print("Macro M0")
sc.pl.umap(adata, color=["SLAMF9", "CCL7", "SPP1"], cmap="inferno")

print("Macro M1")
sc.pl.umap(adata, color=["SELENOP", "CCL13"], cmap="inferno")

print("Macro M2 (alevolar)")
sc.pl.umap(adata, color=["FABP4", "PCOLCE2"], cmap="inferno")

print("Mono")
sc.pl.umap(adata, color=["CD300E", "VCAN", "FCN1"], cmap="inferno")

print("cDC1")
sc.pl.umap(adata, color=["CLEC9A", "DNASE1L3"], cmap="inferno")

print("cDC2")
sc.pl.umap(adata, color=["AREG", "CLEC10A", "CCL17"], cmap="inferno")

print("mDC_mature")
sc.pl.umap(adata, color=["BIRC3", "CCR7", "FSCN1"], cmap="inferno")

print("myeloid TF low")
sc.pl.umap(adata, color=["FCN3"], cmap="inferno")

print("myeloid dividing")
sc.pl.umap(adata, color=["PCLAF", "CDK1"], cmap="inferno")

# %% [markdown]
# ### Plot the cibersort/xCell scores by cell-type

# %%
for tmp_adata in [
    adata,
    adata[
        adata.obs["cell_type"].str.startswith("Macro"),
    ],
]:
    sc.pl.dotplot(
        tmp_adata,
        var_names=[
            "TF:MYC",
            "TF:SPI1",
            "TF:STAT1",
            "TF:NFKB1",
            "cibersort Macrophage M0",
            "cibersort Macrophage M1",
            "cibersort Macrophage M2",
            "cibersort Monocyte",
            "xcell Monocyte",
            "xcell Macrophage M1",
            "xcell Macrophage M2",
            "xcell Macrophage",
        ],
        groupby="cell_type",
        standard_scale="var",
        swap_axes=True,
    )

# %% [markdown]
# ## DE Analysis of CD1C+CD1A+ vs CD1C+ cluster

# %%
sc.pl.umap(adata, color="cell_type")

# %%
res = de.scvi(adata_scvi, model, groupby="cell_type", groups=["cDC2"], group2="cDC2_CD1A_")

# %%
res.to_csv(f"{output_dir}/cdc2_clusters.csv")

# %%
res.loc[(res["raw_normalized_mean1"] > 0.5) | (res["raw_normalized_mean1"] > 0.5), :].to_csv(f"{output_dir}/cdc2_clusters_filtered.csv")

# %%
