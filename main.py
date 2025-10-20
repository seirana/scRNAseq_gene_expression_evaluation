#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: seirana hashemi
"""
import scanpy as sc
import pandas as pd
import numpy as np
import os

adata = sc.read_h5ad('./scIBD_Colon.h5ad')

# Ensure gene names are accessible
gene_names = adata.var_names

# Directory to save the output
output_dir = "cluster_gene_stats"
os.makedirs(output_dir, exist_ok=True)

# Loop over each disease type
for disease in adata.obs['disease'].unique():
    # Subset to that disease
    adata_disease = adata[adata.obs['disease'] == disease]
    
    # Create a DataFrame to collect results
    results = pd.DataFrame(index=gene_names)
    
    # Loop over each major cluster within this disease
    for cluster in adata_disease.obs['major_cluster'].unique(): #major_cluster
        subset = adata_disease[adata_disease.obs['major_cluster'] == cluster]
        X = subset.X
        
        # Convert to dense array if sparse
        if not isinstance(X, np.ndarray):
            X = X.toarray()
        
        # Compute statistics
        mean_vals = np.mean(X, axis=0)
        var_vals  = np.var(X, axis=0)
        min_vals  = np.min(X, axis=0)
        max_vals  = np.max(X, axis=0)
        
        # Add to DataFrame
        results[f"{cluster}_mean"] = mean_vals
        results[f"{cluster}_var"]  = var_vals
        results[f"{cluster}_min"]  = min_vals
        results[f"{cluster}_max"]  = max_vals
    
    # Save per disease
    out_path = os.path.join(output_dir, f"{disease}_gene_stats.csv")
    results.to_csv(out_path)
    print(f"Saved: {out_path}")

#..............................................................................
import pandas as pd
import os

input_dir = "cluster_gene_stats"
output_dir = "top10_genes_per_cluster"
os.makedirs(output_dir, exist_ok=True)

# Loop over each disease stats file
for file in os.listdir(input_dir):
    if not file.endswith("_gene_stats.csv"):
        continue
    
    disease = file.replace("_gene_stats.csv", "")
    df = pd.read_csv(os.path.join(input_dir, file), index_col=0)
    
    top_genes_summary = []
    
    # Identify cluster names from columns
    clusters = sorted(set([col.split('_')[0] for col in df.columns if col.endswith('_mean')]))
    
    for cluster in clusters:
        mean_col = f"{cluster}_mean"
        if mean_col not in df.columns:
            continue
        
        # Get top 10 genes with highest mean
        top_genes = df[mean_col].nlargest(10)
        top_df = pd.DataFrame({
            "cluster": cluster,
            "gene": top_genes.index,
            "mean_expression": top_genes.values
        })
        top_genes_summary.append(top_df)
    
    # Combine all clusters for this disease
    top_genes_summary = pd.concat(top_genes_summary, ignore_index=True)
    
    # Save summary
    out_path = os.path.join(output_dir, f"{disease}_top10_genes.csv")
    top_genes_summary.to_csv(out_path, index=False)
    print(f"Saved: {out_path}")
#..............................................................................
import pandas as pd
import os
import numpy as np

input_dir = "cluster_gene_stats"
output_dir = "disease_avg_expression"
os.makedirs(output_dir, exist_ok=True)

for file in os.listdir(input_dir):
    if not file.endswith("_gene_stats.csv"):
        continue
    
    disease = file.replace("_gene_stats.csv", "")
    df = pd.read_csv(os.path.join(input_dir, file), index_col=0)
    
    # Identify all mean columns for clusters
    mean_cols = [col for col in df.columns if col.endswith('_mean')]
    
    # Compute average of mean expression across clusters
    disease_avg = df[mean_cols].mean(axis=1)
    
    # Save
    out_df = pd.DataFrame({"gene": df.index, "avg_expression_across_clusters": disease_avg})
    out_path = os.path.join(output_dir, f"{disease}_avg_expression.csv")
    out_df.to_csv(out_path, index=False)
    print(f"Saved: {out_path}")
#..............................................................................
import pandas as pd
import os

input_dir = "cluster_gene_stats"
output_file = "disease_avg_expression_all.csv"

disease_dfs = []

for file in os.listdir(input_dir):
    if not file.endswith("_gene_stats.csv"):
        continue

    disease = file.replace("_gene_stats.csv", "")
    df = pd.read_csv(os.path.join(input_dir, file), index_col=0)

    # Get only mean columns
    mean_cols = [col for col in df.columns if col.endswith("_mean")]
    if not mean_cols:
        continue

    # Compute average of mean expression across clusters
    disease_avg = df[mean_cols].mean(axis=1)
    disease_df = pd.DataFrame({disease: disease_avg})
    disease_df.index.name = "gene"

    disease_dfs.append(disease_df)

# Merge all diseases on gene index
merged_df = pd.concat(disease_dfs, axis=1)

# Save one combined file
merged_df.to_csv(output_file)
print(f"Saved combined file: {output_file}")
#..............................................................................
import pandas as pd
import numpy as np
from scipy.stats import kruskal

# Example: adata.obs has 'disease', adata.X has expression
genes = adata.var_names
results = []

for gene_idx, gene in enumerate(genes):
    X = adata[:, gene_idx].X
    if not isinstance(X, np.ndarray):
        X = X.toarray().flatten()

    data = pd.DataFrame({
        "expression": X,
        "disease": adata.obs["disease"].values
    })

    # Kruskal–Wallis across all disease groups
    groups = [data.loc[data["disease"] == d, "expression"].values for d in data["disease"].unique()]
    stat, pval = kruskal(*groups)
    
    results.append({"gene": gene, "kruskal_p": pval})

results_df = pd.DataFrame(results).sort_values("kruskal_p")
results_df["FDR"] = np.minimum(1, results_df["kruskal_p"] * len(results_df) / (np.arange(1, len(results_df)+1)))
results_df.to_csv("gene_disease_significance.csv", index=False)
#..............................................................................
# from joblib import Parallel, delayed
# from scipy.stats import kruskal

# def kruskal_for_gene(gene_idx):
#     X = adata[:, gene_idx].X
#     if not isinstance(X, np.ndarray):
#         X = X.toarray().flatten()
#     df = pd.DataFrame({"expression": X, "disease": adata.obs["disease"].values})
#     groups = [df.loc[df["disease"] == d, "expression"].values for d in df["disease"].unique()]
#     stat, pval = kruskal(*groups)
#     return (adata.var_names[gene_idx], pval)

# results = Parallel(n_jobs=8)(delayed(kruskal_for_gene)(i) for i in range(adata.n_vars))
#..............................................................................
import os
import numpy as np
import pandas as pd
from scipy.stats import kruskal
from scipy import sparse as sp

# --------- config ----------
cluster_col = "major_cluster"
disease_col = "disease"
out_dir = "kruskal_per_cluster"
os.makedirs(out_dir, exist_ok=True)
# ---------------------------
import numpy as np

# Round adata.X to 3 decimal digits safely
if sp.issparse(adata.X):
    # convert to COO for in-place rounding
    X = adata.X.tocoo(copy=True)
    X.data = np.round(X.data, 3)
    adata.X = X.tocsr()  # keep efficient format
else:
    adata.X = np.round(adata.X, 3)

clusters = pd.Index(adata.obs[cluster_col].unique()).astype(str)

def extract_vector(mat, gene_idx):
    #Return 1D numpy array for column `gene_idx` from mat (dense or sparse).
    col = mat[:, gene_idx]
    if sp.issparse(col):
        return col.toarray().ravel()
    if isinstance(col, np.matrix):
        return np.asarray(col).ravel()
    return np.asarray(col).ravel()

for clus in clusters:
    ad_sub = adata[adata.obs[cluster_col].astype(str) == clus]
    if ad_sub.n_obs == 0:
        continue

    diseases = ad_sub.obs[disease_col].astype(str)
    uniq_diseases = diseases.unique().tolist()
    if len(uniq_diseases) < 2:
        continue

    masks = {d: (diseases.values == d) for d in uniq_diseases}
    results = []
    group_sizes = {d: int(np.sum(m)) for d, m in masks.items()}

    for gene_idx, gene in enumerate(adata.var_names):
        x = extract_vector(ad_sub.X, gene_idx)
        groups = [x[masks[d]] for d in uniq_diseases]
        groups = [g for g in groups if g.size > 0]
        if len(groups) < 2:
            continue

        flat = np.concatenate(groups)
        if np.allclose(flat, flat[0]):
            continue

        try:
            stat, pval = kruskal(*groups)
            results.append((gene, stat, pval))
        except Exception:
            continue

    if not results:
        continue

    df = pd.DataFrame(results, columns=["gene", "H_statistic", "p_value"]).sort_values("p_value", kind="mergesort")

    # Benjamini–Hochberg FDR
    m = len(df)
    ranks = np.arange(1, m + 1, dtype=float)
    bh = df["p_value"].values * m / ranks
    bh = np.minimum.accumulate(bh[::-1])[::-1]
    df["FDR"] = np.clip(bh, 0, 1)

    # Round all numerical columns to 3 digits
    df = df.round(3)

    # Save results
    out_path = os.path.join(out_dir, f"{clus}_gene_disease_significance.csv")
    df.to_csv(out_path, index=False)

    # Also save cell counts per disease
    pd.DataFrame.from_dict(group_sizes, orient="index", columns=["n_cells"])\
        .rename_axis("disease").reset_index()\
        .to_csv(os.path.join(out_dir, f"{clus}_group_sizes.csv"), index=False)

    print(f"Saved: {out_path}")