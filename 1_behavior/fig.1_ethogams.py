# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 14:52:25 2025

@author: Mathieu_Thabault
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# Impot files and set directories
# csv_path = r"paste the path to the dataset here + comment out"
# output_folder = r"paste your path here + comment out" 


# Set colors for the plots
behaviour_colors = {
    "empty": "#FFFFFF",
    "supported_rearing": "#2E7D32",
    "unsupported_rearing": "#66BB6A",
    "digging": "#00897B",
    "grooming_long": "#D97904",
    "grooming_short": "#F2B705",
    "scratching": "#C2185B",
    "head_body_twitch": "#E64A19",
}


# Define function
def plot_group_ethogram(df, geno, age, behaviour_colors, mouse_order=None):
    subdf = df[(df["genotype"] == geno) & (df["age"] == age)]
    if subdf.empty:
        print(f"No data for group {geno}, {age}")
        return None
    if mouse_order is None:
        mice = (
            subdf[["video", "mouse"]]
            .drop_duplicates()
            .sort_values(["video", "mouse"])
        )
    else:
        mice = (
            subdf[["video", "mouse"]]
            .drop_duplicates()
            .set_index("mouse")
            .reindex(mouse_order)
            .dropna()
            .reset_index()
            .rename(columns={"index": "mouse"})
        )
    n_mice = len(mice)
    fig_width_in = 12 / 2.54
    fig_height_in = max(3, n_mice * 0.25)
    fig, ax = plt.subplots(figsize=(fig_width_in, fig_height_in))
    y_labels = []
    for i, row in enumerate(mice.itertuples(index=False)):
        vid = row.video
        m = row.mouse
        y_labels.append(f"{m} (v{vid})")
        bouts = (
            subdf[(subdf["mouse"] == m) & (subdf["video"] == vid)]
            .sort_values("start_time")
        )
        for _, bout in bouts.iterrows():
            beh = bout["behavior"]
            color = behaviour_colors.get(beh, "#000000")

            ax.barh(
                y=i,
                width=bout["end_time"] - bout["start_time"],
                left=bout["start_time"],
                height=0.8,
                color=color,
                edgecolor=color,
            )
    ax.set_yticks(range(n_mice))
    ax.set_yticklabels(y_labels, fontsize=7)
    ax.set_xlabel("Time (s)", fontsize=8)
    ax.tick_params(axis="x", labelsize=7, length=2)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.set_title(f"Ethogram – {geno.upper()} {age}", fontsize=10)
    plt.tight_layout()
    return fig


# Generate ethograms
if __name__ == "__main__":
    df = pd.read_csv(csv_path)
    df["age_num"] = df["age"].str.replace("w", "", regex=False).astype(int)
    wt_order = (
        df[df["genotype"] == "wt"]
        .sort_values(["video", "mouse"])
        ["mouse"]
        .unique()
    )
    ko_order = (
        df[df["genotype"] == "ko"]
        .sort_values(["video", "mouse"])
        ["mouse"]
        .unique()
    )
    fig_wt_10w = plot_group_ethogram(df, "wt", "10w", behaviour_colors, mouse_order=wt_order)
    fig_wt_20w = plot_group_ethogram(df, "wt", "20w", behaviour_colors, mouse_order=wt_order)
    fig_ko_10w = plot_group_ethogram(df, "ko", "10w", behaviour_colors, mouse_order=ko_order)
    fig_ko_20w = plot_group_ethogram(df, "ko", "20w", behaviour_colors, mouse_order=ko_order)

    # Create output directory if it does not already exists
    os.makedirs(output_folder, exist_ok=True)
    
    
    # Save ethograms as EPS for publicatio + PNG for presentation
    if fig_wt_10w:
        fig_wt_10w.savefig(f"{output_folder}/ethogram_wt_10w.eps", format="eps", bbox_inches="tight")
        fig_wt_10w.savefig(f"{output_folder}/ethogram_wt_10w.png", format="png", dpi=300, bbox_inches="tight")
    if fig_wt_20w:
        fig_wt_20w.savefig(f"{output_folder}/ethogram_wt_20w.eps", format="eps", bbox_inches="tight")
        fig_wt_20w.savefig(f"{output_folder}/ethogram_wt_20w.png", format="png", dpi=300, bbox_inches="tight")
    if fig_ko_10w:
        fig_ko_10w.savefig(f"{output_folder}/ethogram_ko_10w.eps", format="eps", bbox_inches="tight")
        fig_ko_10w.savefig(f"{output_folder}/ethogram_ko_10w.png", format="png", dpi=300, bbox_inches="tight")
    if fig_ko_20w:
        fig_ko_20w.savefig(f"{output_folder}/ethogram_ko_20w.eps", format="eps", bbox_inches="tight")
        fig_ko_20w.savefig(f"{output_folder}/ethogram_ko_20w.png", format="png", dpi=300, bbox_inches="tight")


    # For size harmonization of exported files, save legend with hex codes separately
    legend_handles = []
    for beh, col in behaviour_colors.items():
        if beh == "empty":
            continue
        label = f"{beh} ({col})"
        legend_handles.append(Patch(facecolor=col, edgecolor=col, label=label))
    fig_leg = plt.figure(figsize=(5, 3))
    ax_leg = fig_leg.add_subplot(111)
    ax_leg.axis("off")
    ax_leg.legend(
        handles=legend_handles,
        loc="center",
        fontsize=8,
        frameon=False,
        ncol=1,                # Set to 2 if you want a two-column legend
        columnspacing=1.2,
        handlelength=1.5,
        handletextpad=0.8
    )
    fig_leg.savefig(
        f"{output_folder}/ethogram_legend.png",
        format="png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close(fig_leg)
