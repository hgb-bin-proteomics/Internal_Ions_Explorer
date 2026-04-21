#!/usr/bin/env python3
#
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "pandas>=2.3.3",
#   "plotly[express]>=6.5.2",
# ]
# ///

# internal-ions-data-analysis
# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pandas as pd
import plotly.express as px


# Adding a column with combined frag types to the fragment centric dfs
def combined_frag_types(df: pd.DataFrame) -> pd.DataFrame:
    df["combined_frag_types"] = df.apply(
        lambda row: str(row["frag_type1"]).strip() + str(row["frag_type2"]).strip(),
        axis=1,
    )
    return df


# Removing fragments that are "nn"
def remove_nn(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["combined_frag_types"] != "nn"]  # pyright: ignore[reportReturnType]


# Binning
def frag_mz_binning(df: pd.DataFrame) -> pd.DataFrame:
    bins = [i for i in range(0, int(df["frag_mz"].max()) + 100, 100)]
    df["frag_mz_bins"] = pd.cut(df["frag_mz"], bins)
    return df.groupby("frag_mz_bins", as_index=False, observed=True)[
        "combined_frag_types"
    ].value_counts()  # pyright: ignore[reportReturnType]


# Proportional Distribution of m/z Values
def proportional_distribution_of_mz_values(fragment_centric_df: pd.DataFrame):
    df = fragment_centric_df.copy()
    df = combined_frag_types(df)
    df = remove_nn(df)
    vc = frag_mz_binning(df)
    vc["frag_mz_bins"] = vc["frag_mz_bins"].astype(str)
    fig = px.bar(
        vc,
        x="frag_mz_bins",
        y="count",
        color="combined_frag_types",
        barmode="group",
        labels={
            "frag_mz_bins": "Fragment M/Z Bins",
            "count": "Count",
            "combined_frag_types": "Fragment Types",
        },
        title="Proportional Distribution of m/z Values (Bin Size = 100)",
    )
    return fig


# Analyzing first and last residue of each ion
def add_first_and_last_aa(df: pd.DataFrame) -> pd.DataFrame:
    df["first_aa"] = df.apply(lambda row: str(row["frag_seq"]).strip()[0], axis=1)
    df["last_aa"] = df.apply(lambda row: str(row["frag_seq"]).strip()[-1], axis=1)
    return df


# Group and value counts by first amino acid
def group_by_first_aa(df: pd.DataFrame) -> pd.DataFrame:
    vc = df.groupby("first_aa", as_index=False, observed=True)[
        "combined_frag_types"
    ].value_counts()
    return vc  # pyright: ignore[reportReturnType]


# Group and value counts by last amino acid
def group_by_last_aa(df: pd.DataFrame) -> pd.DataFrame:
    vc = df.groupby("last_aa", as_index=False, observed=True)[
        "combined_frag_types"
    ].value_counts()
    return vc  # pyright: ignore[reportReturnType]


# First residue of internal and terminal ions
def first_residue_of_internal_and_terminal_ions(fragment_centric_df: pd.DataFrame):
    df = fragment_centric_df.copy()
    df = combined_frag_types(df)
    df = remove_nn(df)
    df = add_first_and_last_aa(df)
    vc = group_by_first_aa(df)
    fig = px.bar(
        vc,
        x="first_aa",
        y="count",
        color="combined_frag_types",
        barmode="group",
        labels={
            "first_aa": "Amino Acid",
            "count": "Count",
            "combined_frag_types": "Fragment Types",
        },
        title="First residue of internal and terminal ions",
    )
    return fig


# Last residue of internal and terminal ions
def last_residue_of_internal_and_terminal_ions(fragment_centric_df: pd.DataFrame):
    df = fragment_centric_df.copy()
    df = combined_frag_types(df)
    df = remove_nn(df)
    df = add_first_and_last_aa(df)
    vc = group_by_last_aa(df)
    fig = px.bar(
        vc,
        x="last_aa",
        y="count",
        color="combined_frag_types",
        barmode="group",
        labels={
            "last_aa": "Amino Acid",
            "count": "Count",
            "combined_frag_types": "Fragment Types",
        },
        title="Last residue of internal and terminal ions",
    )
    return fig


# Add top internal ion sequence lengths
def add_top_internal_seq_len(df: pd.DataFrame) -> pd.DataFrame:
    df["top1_internal_seq_len"] = df.apply(
        lambda row: len(str(row["top1_internal_seq"]).strip()), axis=1
    )
    df["top2_internal_seq_len"] = df.apply(
        lambda row: len(str(row["top2_internal_seq"]).strip()), axis=1
    )
    df["top3_internal_seq_len"] = df.apply(
        lambda row: len(str(row["top3_internal_seq"]).strip()), axis=1
    )
    return df


# Plot of the length of the top 1 internal ion
def density_plot_of_the_length_of_the_top_1_internal_ion(
    spectrum_centric_df: pd.DataFrame,
):
    df = spectrum_centric_df.copy()
    df = add_top_internal_seq_len(df)
    df = df[df["perc_internal"] > 0.0]
    df = df.dropna(subset="top1_internal_seq")  # pyright: ignore[reportCallIssue]
    vc = df.groupby("top1_internal_seq_len", as_index=False, observed=True)[
        "top1_internal_seq_len"
    ].value_counts()
    fig = px.bar(
        vc,
        x="top1_internal_seq_len",
        y="count",
        labels={
            "top1_internal_seq_len": "Sequence Length",
            "count": "Count",
        },
        title="Plot of the length of the top 1 internal ion",
    )
    return fig


# Plot of the length of the top 2 internal ion
def density_plot_of_the_length_of_the_top_2_internal_ion(
    spectrum_centric_df: pd.DataFrame,
):
    df = spectrum_centric_df.copy()
    df = add_top_internal_seq_len(df)
    df = df[df["perc_internal"] > 0.0]
    df = df.dropna(subset="top2_internal_seq")  # pyright: ignore[reportCallIssue]
    vc = df.groupby("top2_internal_seq_len", as_index=False, observed=True)[
        "top2_internal_seq_len"
    ].value_counts()
    fig = px.bar(
        vc,
        x="top2_internal_seq_len",
        y="count",
        labels={
            "top2_internal_seq_len": "Sequence Length",
            "count": "Count",
        },
        title="Plot of the length of the top 2 internal ion",
    )
    return fig


# Plot of the length of the top 3 internal ion
def density_plot_of_the_length_of_the_top_3_internal_ion(
    spectrum_centric_df: pd.DataFrame,
):
    df = spectrum_centric_df.copy()
    df = add_top_internal_seq_len(df)
    df = df[df["perc_internal"] > 0.0]
    df = df.dropna(subset="top3_internal_seq")  # pyright: ignore[reportCallIssue]  # pright: ignore[reportCallIssue]
    vc = df.groupby("top3_internal_seq_len", as_index=False, observed=True)[
        "top3_internal_seq_len"
    ].value_counts()
    fig = px.bar(
        vc,
        x="top3_internal_seq_len",
        y="count",
        labels={
            "top3_internal_seq_len": "Sequence Length",
            "count": "Count",
        },
        title="Plot of the length of the top 3 internal ion",
    )
    return fig


# Plot of the length of the top internal ions
def density_plot_of_the_length_of_the_top_internal_ions(
    spectrum_centric_df: pd.DataFrame,
):
    df = spectrum_centric_df.copy()
    df = df.melt(
        id_vars=list(
            set(df.columns.tolist())
            - {"top1_internal_seq", "top2_internal_seq", "top3_internal_seq"}
        ),
        value_vars=["top1_internal_seq", "top2_internal_seq", "top3_internal_seq"],
        var_name="top_n",
        value_name="internal_seq",
    )
    df["internal_seq_len"] = df.apply(
        lambda row: len(str(row["internal_seq"]).strip()), axis=1
    )
    df = df[df["perc_internal"] > 0.0]
    df = df.dropna(subset="internal_seq")  # pyright: ignore[reportCallIssue]
    vc = df.groupby("internal_seq_len", as_index=False, observed=True)[
        "top_n"
    ].value_counts()
    fig = px.bar(
        vc,
        x="internal_seq_len",
        y="count",
        color="top_n",
        barmode="group",
        labels={
            "internal_seq_len": "Sequence Length",
            "count": "Count",
            "top_n": "Top (X) Internal Ion",
        },
        title="Plot of the length of the top internal ions",
    )
    return fig
