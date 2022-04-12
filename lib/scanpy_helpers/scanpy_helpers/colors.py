import altair as alt
import pandas as pd


def set_scale_anndata(adata, column, palette=None):
    if palette is None:
        palette = column

    adata._sanitize()

    tmp_cols = getattr(COLORS, palette)
    adata.uns[f"{column}_colors"] = [
        tmp_cols[cat] for cat in adata.obs[column].cat.categories
    ]


def altair_scale(variable, *, data=None, data_col=None, **kwargs):
    """
    Discrete color scale for altair based on our color definitions.

    Parameters:
    -----------
    variable
        name of the color scale
    data
        Data frame used for the chart. If specified, will only show values that actually occur in the data.
    data_col
        If specified, check this column in `data` instead of `variable`

    Returns
    -------
    Altair color scale
    """
    tmp_colors = getattr(COLORS, variable)
    if data is not None:
        data_col = variable if data_col is None else data_col
        tmp_colors = {k: tmp_colors[k] for k in sorted(data[data_col].unique())}

    return alt.Scale(
        domain=list(tmp_colors.keys()),
        range=list(tmp_colors.values()),
        **kwargs,
    )


def altair_scale_mpl(scheme, **kwargs):
    """
    Use a continuous color scheme from mpl with altair
    """
    from matplotlib import cm
    from matplotlib.colors import to_hex

    return alt.Scale(
        range=[to_hex(x) for x in cm.get_cmap(scheme, 1000)(range(1000))], **kwargs
    )


def plot_palette(variable):
    """Display a palette"""
    tmp_cols = getattr(COLORS, variable)
    return (
        alt.Chart(
            pd.DataFrame.from_dict(tmp_cols, orient="index", columns=["color"])
            .reset_index()
            .rename(columns={"index": variable})
        )
        .mark_rect(height=40, width=30)
        .encode(
            x=alt.X(variable),
            color=alt.Color(variable, scale=altair_scale(variable), legend=None),
        )
    )


def plot_all_palettes():
    return alt.vconcat(
        *[plot_palette(v) for v in COLORS.__dict__.keys() if not v.startswith("_")]
    ).resolve_scale(color="independent")


class COLORS:
    leiden = {
        str(k): v
        for k, v in enumerate(
            [
                "#1f77b4",
                "#ff7f0e",
                "#279e68",
                "#d62728",
                "#aa40fc",
                "#8c564b",
                "#e377c2",
                "#b5bd61",
            ]
        )
    }
    sex = {
        "male": "#80b1d3",
        "female": "#e41a1c",
        "unknown": "#dddddd",
    }
    tumor_stage = {
        "early": "#998ec3",
        "advanced": "#f1a340",
        "unknown": "#dddddd",
    }
    tumor_stage_verbose = {
        "early (I/II)": "#998ec3",
        "advanced (III/IV)": "#f1a340",
        "unknown": "#dddddd",
    }
    immune_infiltration = {
        "B": "#ff7f00",
        "T": "#1f78b4",
        "M": "#33a02c",
        "desert": "#999999",
        # "n/a": "#dddddd",
    }
    tumor_type = {
        "LUSC": "#beaed4",
        "LUAD": "#ffff99",
        "NSCLC NOS": "#999999",
        "LUAD EMT": "#386cb0",
        "LUAD NE": "#f0027f",
        "LUAD dedifferentiated": "#bf5b17",
    }
    condition = {
        "healthy_control": "#7fc97f",
        "non-cancer": "#7fc97f",
        "LUSC": "#beaed4",
        "COPD": "#fdc086",
        "LUAD": "#ffff99",
        "NSCLC NOS": "#999999",
        "NSCLC": "#999999",
        "NOS": "#999999",
        "other": "#DDDDDD",
        "LUAD EMT": "#386cb0",
        "LUAD NE": "#f0027f",
        "LUAD dedifferentiated": "#bf5b17",
    }
    origin = {
        "tumor_metastasis": "#a6761d",
        "metastases": "#a6761d",
        "other": "#DDDDDD",
        "normal": "#66a61e",
        "normal_adjacent": "#d95f02",
        "adj. normal": "#d95f02",
        "tumor_primary": "#e6ab02",
        "primary tumor": "#e6ab02",
        "nan": "#DDDDDD",
    }
    platform = {
        "10x": "#377eb8",
        "Smart-seq2": "#e41a1c",
        "DropSeq": "#984ea3",
        "BD-Rhapsody": "#4daf4a",
        "InDrop": "#ff7f00",
        "Singleron": "#ffff33",
    }
    cell_type_coarse = {
        "B cell": "#1f77b4",
        "Endothelial cell": "#ff7f0e",
        "Epithelial cell": "#279e68",
        "Macrophage/Monocyte": "#d62728",
        "Mast cell": "#aa40fc",
        "NK cell": "#8c564b",
        "Neutrophils": "#e377c2",
        "Plasma cell": "#b5bd61",
        "Stromal": "#17becf",
        "T cell": "#aec7e8",
        "cDC": "#ffbb78",
        "pDC": "#98df8a",
    }
    cell_type_major = {
        "Alveolar cell type 1": "#023fa5",
        "Alveolar cell type 2": "#7d87b9",
        "B cell": "#bec1d4",
        "Ciliated": "#d6bcc0",
        "Club": "#bb7784",
        "DC mature": "#8e063b",
        "Endothelial cell": "#4a6fe3",
        "Macrophage": "#8595e1",
        "Macrophage alveolar": "#b5bbe3",
        "Mast cell": "#e6afb9",
        "Monocyte": "#e07b91",
        "NK cell": "#d33f6a",
        "Neutrophils": "#11c638",
        "Plasma cell": "#8dd593",
        "Stromal": "#c6dec7",
        "T cell CD4": "#ead3c6",
        "T cell CD8": "#f0b98d",
        "T cell regulatory": "#ef9708",
        "Tumor cells": "#0fcfc0",
        "cDC1": "#9cded6",
        "cDC2": "#d5eae7",
        "other": "#f3e1eb",
        "pDC": "#f6c4e1",
        "transitional club/AT2": "#f79cd4",
    }
    cell_type_tumor = {
        "Hemoglobin+": "#1f77b4",
        "Hepatocytes": "#ff7f0e",
        "Tumor cells LUAD": "#279e68",
        "Tumor cells LUAD EMT": "#d62728",
        "Tumor cells LUAD MSLN": "#aa40fc",
        "Tumor cells LUAD NE": "#8c564b",
        "Tumor cells LUAD dedifferentiated": "#e377c2",
        "Tumor cells LUAD mitotic": "#b5bd61",
        "Tumor cells LUSC": "#17becf",
        "Tumor cells LUSC mitotic": "#aec7e8",
        "transitional club/AT2": "#ffbb78",
    }
    study = {
        "Adams_Kaminski_2020": "#023fa5",
        "Chen_Zhang_2020": "#7d87b9",
        "Goveia_Carmeliet_2020": "#bec1d4",
        "Guo_Zhang_2018": "#d6bcc0",
        "Habermann_Kropski_2020": "#bb7784",
        "He_Fan_2021": "#8e063b",
        "Kim_Lee_2020": "#4a6fe3",
        "Lambrechts_Thienpont_2018": "#8595e1",
        "Laughney_Massague_2020": "#e07b91",
        "Madissoon_Meyer_2020": "#d33f6a",
        "Leader_Merad_2021": "#11c638",
        "Maynard_Bivona_2020": "#8dd593",
        "Mayr_Schiller_2020": "#c6dec7",
        "Reyfman_Misharin_2018": "#ead3c6",
        "Travaglini_Krasnow_2020": "#f0b98d",
        "UKIM-V": "#ef9708",
        "Vieira_Teichmann_2019": "#0fcfc0",
        "Wu_Zhou_2021": "#9cded6",
        "Zilionis_Klein_2019": "#d5eae7",
    }

    scissor = {}
