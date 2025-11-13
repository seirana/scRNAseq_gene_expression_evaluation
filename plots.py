#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: seirana
"""
import os
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
from scipy import sparse as sp
import matplotlib.pyplot as plt

ADATA_PATH = "/home/shashemi/Desktop/Data/Single Cell/sinlge cell datasets/scIBD_Colon.h5ad"  # change to your .h5ad path
adata = sc.read_h5ad(ADATA_PATH)

# Configuration
gene_name = "SPOCD1"
out_dir = "plots_spocd1_scatter"
os.makedirs(out_dir, exist_ok=True)

# Make sure gene exists
if gene_name not in adata.var_names:
    raise ValueError(f"Gene {gene_name} not found in adata.var_names")

gene_idx = adata.var_names.get_loc(gene_name)
cell_types = sorted(adata.obs["major_cluster"].astype(str).unique())
disease_col = "disease"

# Define a consistent disease order (if known)
desired_order = [
    "Healthy",
    "UC_non_inflamed",
    "UC_inflamed",
    "CD_inflamed",
    "Colitis_inflamed"
]

for ct in cell_types:
    ad_ct = adata[adata.obs["major_cluster"].astype(str) == ct]
    if ad_ct.n_obs == 0:
        continue

    # Extract expression vector for the gene
    col = ad_ct[:, gene_idx].X
    expr = col.toarray().ravel() if sp.issparse(col) else np.asarray(col).ravel()

    # Skip cell types where gene not expressed
    if np.all(expr == 0):
        print(f"[{ct}] SPOCD1 not expressed (all zeros) — skipping.")
        continue

    diseases = ad_ct.obs[disease_col].astype(str).values
    uniq_dis = np.unique(diseases)

    # Determine actual order (sorted intersection)
    order = [d for d in desired_order if d in uniq_dis] or sorted(uniq_dis)

    # Collect per-disease expression values
    grouped = {d: expr[diseases == d] for d in order}

    # --- Plot ---
    plt.figure(figsize=(7, 4))
    for i, d in enumerate(order):
        vals = grouped[d]
        x_jitter = np.random.normal(i + 1, 0.08, size=len(vals))  # small jitter
        plt.scatter(
            x_jitter,
            vals,
            alpha=0.5,
            s=8,
            label=d if i == 0 else "",
            color="royalblue"
        )

    plt.title(f"{gene_name} Expression per Disease\nCell type: {ct}")
    plt.xlabel("Disease")
    plt.ylabel("Expression level")
    plt.ylim(0, 5)
    plt.xticks(range(1, len(order) + 1), order, rotation=25, ha="right")
    plt.tight_layout()

    # Save plot
    out_path = os.path.join(out_dir, f"{gene_name}_{ct.replace(' ', '_')}.png")
    plt.savefig(out_path, dpi=150)
    plt.close()

    print(f"Saved scatter plot for {ct}: {out_path}")

print("All scatter plots generated successfully.")



#..............................................................................
#
#
#..............................................................................

# Configuration
gene_name = "SPOCD1"
out_dir = "plots_spocd1_box"
os.makedirs(out_dir, exist_ok=True)

# Ensure gene exists
if gene_name not in adata.var_names:
    raise ValueError(f"Gene {gene_name} not found in adata.var_names")

gene_idx = adata.var_names.get_loc(gene_name)
cell_types = sorted(adata.obs["major_cluster"].astype(str).unique())
disease_col = "disease"

# Define consistent disease order
desired_order = [
    "Healthy",
    "UC_non_inflamed",
    "UC_inflamed",
    "CD_inflamed",
    "Colitis_inflamed"
]

for ct in cell_types:
    ad_ct = adata[adata.obs["major_cluster"].astype(str) == ct]
    if ad_ct.n_obs == 0:
        continue

    # Extract expression for the gene
    col = ad_ct[:, gene_idx].X
    expr = col.toarray().ravel() if sp.issparse(col) else np.asarray(col).ravel()

    if np.all(expr == 0):
        print(f"[{ct}] SPOCD1 not expressed (all zeros) — skipping.")
        continue

    diseases = ad_ct.obs[disease_col].astype(str).values
    uniq_dis = np.unique(diseases)
    order = [d for d in desired_order if d in uniq_dis] or sorted(uniq_dis)

    # Collect stats
    data_per_dis = []
    means, mins, maxs, vars_ = [], [], [], []

    for d in order:
        vals = expr[diseases == d]
        data_per_dis.append(vals)
        means.append(np.mean(vals))
        mins.append(np.min(vals))
        maxs.append(np.max(vals))
        vars_.append(np.var(vals))

    # --- Plot ---
    plt.figure(figsize=(7, 4))
    # Boxplot (distribution)
    bp = plt.boxplot(data_per_dis, labels=order, showfliers=False, patch_artist=True)
    for box in bp["boxes"]:
        box.set(facecolor="lightgray", alpha=0.5)

    # Overlay mean (green dots)
    plt.scatter(range(1, len(order) + 1), means, color="green", s=50, label="Mean", zorder=3)

    # Overlay min (red triangles)
    plt.scatter(range(1, len(order) + 1), mins, color="red", marker="v", s=40, label="Min", zorder=3)

    # Overlay max (blue triangles)
    plt.scatter(range(1, len(order) + 1), maxs, color="blue", marker="^", s=40, label="Max", zorder=3)

    # Annotate variance
    for i, var in enumerate(vars_):
        plt.text(i + 1, 4.7, f"Var={var:.2f}", ha="center", va="top", fontsize=8, color="dimgray")

    plt.title(f"{gene_name} Expression Summary per Disease\nCell type: {ct}")
    plt.xlabel("Disease")
    plt.ylabel("Expression Level")
    plt.ylim(0, 10)
    plt.xticks(rotation=25, ha="right")
    plt.legend(frameon=False, loc="upper left")
    plt.tight_layout()

    # Save
    out_path = os.path.join(out_dir, f"{gene_name}_{ct.replace(' ', '_')}_boxplot.png")
    plt.savefig(out_path, dpi=150)
    plt.close()

    print(f"Saved boxplot with stats for {ct}: {out_path}")

print("All boxplots generated successfully.")
#------------------------------------------------------------------------------
#
#
#------------------------------------------------------------------------------
gene_name = "SPOCD1"
effect_sizes = []

cell_types = sorted(adata.obs["major_cluster"].astype(str).unique())
disease_col = "disease"

for ct in cell_types:
    ad_ct = adata[adata.obs["major_cluster"].astype(str) == ct]
    if ad_ct.n_obs == 0:
        continue

    col = ad_ct[:, adata.var_names.get_loc(gene_name)].X
    expr = col.toarray().ravel() if sp.issparse(col) else np.asarray(col).ravel()
    if np.all(expr == 0):
        continue

    diseases = ad_ct.obs[disease_col].astype(str).values
    uniq_dis = np.unique(diseases)

    means = [np.mean(expr[diseases == d]) for d in uniq_dis]
    eff_size = np.max(means) - np.min(means)
    effect_sizes.append({"cell_type": ct, "effect_size": eff_size})

# Convert to DataFrame
df = pd.DataFrame(effect_sizes).sort_values("effect_size", ascending=False)

# --- Plot ---
plt.figure(figsize=(8, 4))
bars = plt.bar(df["cell_type"], df["effect_size"], color="skyblue", edgecolor="black")

# Add thresholds for interpretation
plt.axhline(0.25, color="orange", linestyle="--", label="Medium threshold (0.25)")
plt.axhline(0.5, color="red", linestyle="--", label="High threshold (0.5)")

plt.title(f"Effect Size of {gene_name} Across Cell Types")
plt.ylabel("Effect Size (max(mean) - min(mean))")
plt.xticks(rotation=30, ha="right")
plt.ylim(0, 1)  # adjust depending on your range
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(f"Effect Size of {gene_name} Across Cell Types", dpi=150)
plt.show()


# Build a table of mean expression per disease per cell type
records = []
for ct in cell_types:
    ad_ct = adata[adata.obs["major_cluster"].astype(str) == ct]
    if ad_ct.n_obs == 0:
        continue

    col = ad_ct[:, adata.var_names.get_loc(gene_name)].X
    expr = col.toarray().ravel() if sp.issparse(col) else np.asarray(col).ravel()
    if np.all(expr == 0):
        continue

    diseases = ad_ct.obs[disease_col].astype(str).values
    uniq_dis = np.unique(diseases)
    for d in uniq_dis:
        records.append({
            "cell_type": ct,
            "disease": d,
            "mean_expression": np.mean(expr[diseases == d])
        })

heat_df = pd.DataFrame(records).pivot(index="cell_type", columns="disease", values="mean_expression")

# --- Heatmap ---
plt.figure(figsize=(8, 5))
sns.heatmap(heat_df, cmap="viridis", annot=True, fmt=".7f")
plt.title(f"Mean Expression of {gene_name} per Disease and Cell Type")
plt.xlabel("Disease")
plt.ylabel("Cell Type")
plt.tight_layout()
plt.savefig(f"Mean Expression of {gene_name} per Disease and Cell Type", dpi=150)
plt.show()
#------------------------------------------------------------------------------
#
#
#------------------------------------------------------------------------------

# -------------------------------------------------------------
# CONFIGURATION
# -------------------------------------------------------------
input_dir = "/home/shashemi/Desktop/Projects/antiTNF/output_SPOCD1/uc_inflamed_separated"
bonf_p = 0.05 / 161060  # Bonferroni corrected threshold
N_GENES = 16106         # <<< FIXED denominator
out_plot = "uc_inflamed_fraction_significant.png"

# -------------------------------------------------------------
# PROCESSING
# -------------------------------------------------------------
summary = []

files = [f for f in os.listdir(input_dir) if f.endswith("_UCinflamed.csv")]
if not files:
    raise FileNotFoundError(f"No *_UCinflamed.csv files found in {input_dir}")

for fname in sorted(files):
    path = os.path.join(input_dir, fname)
    try:
        df = pd.read_csv(path)
    except Exception as e:
        print(f"⚠️ Skipping {fname}: {e}")
        continue

    # Normalize column names
    cols = {c.lower(): c for c in df.columns}
    required = ["kw_p", "kw_fdr", "effect_size", "mw_p", "mw_fdr"]
    if not all(k in cols for k in required):
        print(f"⚠️ Missing required columns in: {fname}")
        continue

    # Apply filters
    mask = (
        (df[cols["kw_p"]] < bonf_p)
        & (df[cols["kw_fdr"]] < 0.05)
        & (df[cols["effect_size"]] > 0.5)
        & (df[cols["mw_p"]] < bonf_p)
        & (df[cols["mw_fdr"]] < 0.05)
    )

    n_sig = int(mask.sum())

    # NEW FRACTION: divide by 16106
    fraction = n_sig / N_GENES

    cell_type = fname.replace("_UCinflamed.csv", "")
    summary.append({
        "cell_type": cell_type,
        "significant": n_sig,
        "fraction": fraction
    })

# -------------------------------------------------------------
# PLOT RESULTS
# -------------------------------------------------------------
if not summary:
    raise RuntimeError("No valid data to plot.")

res_df = pd.DataFrame(summary).sort_values("fraction", ascending=False)

plt.figure(figsize=(10, 5))
ax = sns.barplot(data=res_df, x="cell_type", y="fraction", color="steelblue", edgecolor="black")
ax.set_ylim(0, 1)
ax.set_ylabel("Fraction of significant genes (out of 16106)")
ax.set_xlabel("Cell type")
ax.set_title("UC-inflamed specific genes per cell type\n"
             "(Strict multi-test filtering + effect size > 0.5)")

# Annotate bars with percentages
for i, v in enumerate(res_df["fraction"].values):
    ax.text(i, min(v + 0.02, 0.98), f"{v:.2%}", ha="center", va="bottom", fontsize=9)

plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(out_plot, dpi=150)
plt.show()

print(f"✅ Plot saved as: {out_plot}")
