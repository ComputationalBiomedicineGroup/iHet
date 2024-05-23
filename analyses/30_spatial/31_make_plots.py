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

sc.settings.set_figure_params(figsize=(4, 4))

# %%
path_visium = nxfvars.get(
    "visium_slide",
    "/home/sturm/projects/2021/nsclc_heterogeneity/data/01_processed/spatial/lungCancerVisium/outs",
)
path_cell2location_estimates = nxfvars.get(
    "cell2location",
    "/home/sturm/projects/2021/nsclc_heterogeneity/data/01_processed/spatial/cell2location_estimates.csv",
)
cpus = nxfvars.get("cpus", 16)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/ihet")

# %%
threadpool_limits(cpus)

# %%
adata = sc.read_visium(path_visium)
cell2location_estimates = pd.read_csv(path_cell2location_estimates, index_col=0).drop(
    columns=["in_tissue", "array_row", "array_col"]
)
adata = adata[cell2location_estimates.index,]
adata.obs = adata.obs.join(cell2location_estimates)
# %%
