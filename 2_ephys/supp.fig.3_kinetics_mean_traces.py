# -*- coding: utf-8 -*-
"""
Created on Mon May  5 10:44:18 2025

@author: Mathieu_Thabault
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Config parameters from Simplyfire detection parameters and kinetics fitting curves
seed = 4  # Exponent for rise sharpness
delay = 5  # Synaptic delay in ms
file_path = r"path"
save_dir = r"path"
os.makedirs(save_dir, exist_ok=True)

# Define colors per (genotype, age)
group_colors = {
    ("wt", "10w"): "lightgray",
    ("wt", "20w"): "lightgray",
    ("ko", "10w"): "#ccf6ff",
    ("ko", "20w"): "#6099fc"
}

mean_colors = {
    ("wt", "10w"): "black",
    ("wt", "20w"): "black",
    ("ko", "10w"): "#34d6fa",
    ("ko", "20w"): "#0213f7"
}

def compute_mean_params(df):
    stats = {
        "amplitude": (df["amp"].mean(), df["amp"].std()),
        "rise": (df["10_90_rise"].mean(), df["10_90_rise"].std()),
        "decay": (df["decay_const"].mean(), df["decay_const"].std()),
        "halfwidth": (df["halfwidth"].mean(), df["halfwidth"].std())
    }
    return stats

def plot_pscs(df, stats, seed, delay, title=None, genotype=None, age=None, save_basename=None):
    plt.figure(figsize=(4, 8))
    t = np.linspace(0, 100, 50000)

    # Determine colors based on genotype + age
    key = (genotype, age)
    color_individual = group_colors.get(key, "lightgray")
    color_mean = mean_colors.get(key, "black")

    # Plot individual PSCs
    for _, row in df.iterrows():
        A = row["amp"]
        tau_r = row["10_90_rise"]
        tau_d = row["decay_const"]

        t_shifted = t - delay
        t_shifted[t_shifted < 0] = 0

        psc = A * (np.exp(-t_shifted / tau_d) - np.exp(-t_shifted / tau_r)**seed)
        psc *= A / np.max(psc)

        plt.plot(t, psc, color=color_individual)

    # Mean IPSC
    amp_mean, _ = stats["amplitude"]
    rise_mean, _ = stats["rise"]
    decay_mean, _ = stats["decay"]

    psc_mean = amp_mean * (np.exp(-t_shifted / decay_mean) - np.exp(-t_shifted / rise_mean)**seed)
    psc_mean *= amp_mean / np.max(psc_mean)

    plt.plot(t, psc_mean, color=color_mean, linewidth=2, label="Mean PSC")

    plt.axis('off')

    # Scale bar
    scale_x = 20
    scale_y = 100
    x_ref = 75
    y_ref = 400

    plt.plot([x_ref, x_ref + scale_x], [y_ref, y_ref], color='black', linewidth=2)
    plt.plot([x_ref + scale_x, x_ref + scale_x], [y_ref, y_ref - scale_y], color='black', linewidth=2)

    plt.text(x_ref + scale_x / 2, y_ref + 3, f"{scale_x} ms", ha='center', va='bottom')
    plt.text(x_ref + scale_x + 2, y_ref - scale_y / 2, f"{scale_y} pA", ha='left', va='center')

    plt.xlim(0, 100)
    plt.ylim(-5, 600)

    if title:
        plt.text(50, 300, title, ha='center', va='center', fontsize=12)

    plt.tight_layout()

    # Save
    if save_basename:
        eps_path = os.path.join(save_dir, save_basename + ".eps")
        png_path = os.path.join(save_dir, save_basename + ".png")
        plt.savefig(eps_path, bbox_inches='tight')
        plt.savefig(png_path, dpi=300, bbox_inches='tight')
        print(f"Saved {eps_path} and {png_path}")

    plt.show()

# Run functions

df = pd.read_csv(file_path)
df = df[df["10_90_rise"] <= 6]

grouped = df.groupby(["genotype", "age"])

for (genotype, age), group_df in grouped:
    print(f"\n--- {genotype} - {age} ---")
    stats = compute_mean_params(group_df)

    for key, (mean, std) in stats.items():
        print(f"{key.capitalize()}: {mean:.2f} ± {std:.2f}")

    title = f"{genotype} - {age}"
    save_basename = f"representative_ipscs_kinetics_{genotype}_{age}"

    plot_pscs(group_df, stats, seed, delay, title=title, genotype=genotype, age=age, save_basename=save_basename)
