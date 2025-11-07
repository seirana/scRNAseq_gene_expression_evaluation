#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: seirana
"""
import scanpy as sc
import pandas as pd
import numpy as np
import os

adata = sc.read_h5ad('/work_beegfs/sukmb655/data_scDRS/sinlge cell datasets/scIBD_Colon.h5ad')

# Ensure gene names are accessible
gene_names = adata.var_names

# Directory to save the output
output_dir = "/work_beegfs/sukmb655/antiTNF/cluster_gene_stats"
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
    out_path = os.path.join(output_dir, f"{disease}_gene_stats_major.csv")
    results.to_csv(out_path)
    print(f"Saved: {out_path}")

#..............................................................................
import pandas as pd
import os

input_dir = "/work_beegfs/sukmb655/antiTNF/cluster_gene_stats"
output_dir = "/work_beegfs/sukmb655/antiTNF/top10_genes_per_cluster"
os.makedirs(output_dir, exist_ok=True)

# Loop over each disease stats file
for file in os.listdir(input_dir):
    if not file.endswith("_gene_stats_major.csv"):
        continue
    
    disease = file.replace("_gene_stats_major.csv", "")
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
    out_path = os.path.join(output_dir, f"{disease}_top10_genes_major.csv")
    top_genes_summary.to_csv(out_path, index=False)
    print(f"Saved: {out_path}")
#..............................................................................
import pandas as pd
import os
import numpy as np

input_dir = "/work_beegfs/sukmb655/antiTNF/cluster_gene_stats"
output_dir = "/work_beegfs/sukmb655/antiTNF/disease_avg_expression"
os.makedirs(output_dir, exist_ok=True)

for file in os.listdir(input_dir):
    if not file.endswith("_gene_stats_major.csv"):
        continue
    
    disease = file.replace("_gene_stats_major.csv", "")
    df = pd.read_csv(os.path.join(input_dir, file), index_col=0)
    
    # Identify all mean columns for clusters
    mean_cols = [col for col in df.columns if col.endswith('_mean')]
    
    # Compute average of mean expression across clusters
    disease_avg = df[mean_cols].mean(axis=1)
    
    # Save
    out_df = pd.DataFrame({"gene": df.index, "avg_expression_across_clusters": disease_avg})
    out_path = os.path.join(output_dir, f"{disease}_avg_expression_major.csv")
    out_df.to_csv(out_path, index=False)
    print(f"Saved: {out_path}")
#..............................................................................
import pandas as pd
import os

input_dir = "/work_beegfs/sukmb655/antiTNF/cluster_gene_stats"
output_file = "/work_beegfs/sukmb655/antiTNF/disease_avg_expression_all_major.csv"

disease_dfs = []

for file in os.listdir(input_dir):
    if not file.endswith("_gene_stats_major.csv"):
        continue

    disease = file.replace("_gene_stats_major.csv", "")
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
results_df.to_csv("gene_disease_significance_major.csv", index=False)
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
out_dir = "/work_beegfs/sukmb655/antiTNF/kruskal_per_cluster"
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
    df = df.round(10)

    # Save results
    out_path = os.path.join(out_dir, f"{clus}_gene_disease_significance_major.csv")
    df.to_csv(out_path, index=False)

    # Also save cell counts per disease
    pd.DataFrame.from_dict(group_sizes, orient="index", columns=["n_cells"])\
        .rename_axis("disease").reset_index()\
        .to_csv(os.path.join(out_dir, f"{clus}_group_sizes_major.csv"), index=False)

    print(f"Saved: {out_path}")
# -------------------------------------------------------------------------
# POST-HOC ANALYSIS AFTER KRUSKAL–WALLIS
# For each cluster:
#   - Identify genes with significant Kruskal–Wallis FDR (< 0.05)
#   - Run pairwise Mann–Whitney U tests between all disease groups
#   - Apply Benjamini–Hochberg correction for multiple testing
#   - Save detailed results per cluster
# -------------------------------------------------------------------------
from itertools import combinations
from scipy.stats import mannwhitneyu
import pandas as pd
import numpy as np
import os

cluster_col = "major_cluster"
disease_col = "disease"
out_dir = "/work_beegfs/sukmb655/antiTNF/kruskal_per_cluster"

for clus in pd.Index(adata.obs[cluster_col].unique()).astype(str):

    # Path to the Kruskal–Wallis results for this cluster
    kw_path = os.path.join(out_dir, f"{clus}_gene_disease_significance_major.csv")
    if not os.path.exists(kw_path):
        print(f"Skipping {clus}: no Kruskal–Wallis result found.")
        continue

    print(f"\nRunning post-hoc analysis for cluster: {clus}")

    # Read Kruskal–Wallis results
    df_kw = pd.read_csv(kw_path)

    # Filter to significant genes (you can change the FDR threshold if needed)
    sig_df = df_kw[df_kw["FDR"] < 0.05].copy()
    if sig_df.empty:
        print(f"No significant genes found in cluster {clus}. Skipping post-hoc tests.")
        continue

    # Subset AnnData to this cluster
    ad_sub = adata[adata.obs[cluster_col].astype(str) == clus]

    # Extract disease information
    diseases = ad_sub.obs[disease_col].astype(str)
    uniq_diseases = diseases.unique().tolist()
    if len(uniq_diseases) < 2:
        print(f"Only one disease in cluster {clus}, skipping.")
        continue

    # Build boolean masks for each disease (faster indexing)
    masks = {d: (diseases.values == d) for d in uniq_diseases}

    # Prepare disease pair combinations for Mann–Whitney
    disease_pairs = list(combinations(uniq_diseases, 2))
    pairwise_results = []

    # Loop over significant genes and perform pairwise Mann–Whitney U tests
    for _, row in sig_df.iterrows():
        gene = row["gene"]
        gene_idx = adata.var_names.get_loc(gene)

        # Extract expression vector for this gene
        col = ad_sub[:, gene_idx].X
        if hasattr(col, "toarray"):  # handle sparse
            x = col.toarray().ravel()
        else:
            x = np.asarray(col).ravel()

        # Loop over all disease pairs
        for d1, d2 in disease_pairs:
            g1 = x[masks[d1]]
            g2 = x[masks[d2]]

            # Skip pairs with very few cells (optional safeguard)
            if g1.size < 5 or g2.size < 5:
                continue

            # Mann–Whitney U test (two-sided)
            stat, pval = mannwhitneyu(g1, g2, alternative="two-sided")

            # Store results
            pairwise_results.append({
                "gene": gene,
                "disease_1": d1,
                "disease_2": d2,
                "mw_stat": stat,
                "mw_p": pval,
                "n1": g1.size,
                "n2": g2.size
            })

    # Adjust p-values (Benjamini–Hochberg FDR) and save
    if not pairwise_results:
        print(f"No valid pairwise comparisons for cluster {clus}.")
        continue

    pw_df = pd.DataFrame(pairwise_results).sort_values("mw_p")
    m = len(pw_df)
    ranks = np.arange(1, m + 1, dtype=float)
    bh = pw_df["mw_p"].values * m / ranks
    bh = np.minimum.accumulate(bh[::-1])[::-1]
    pw_df["mw_FDR"] = np.clip(bh, 0, 1)

    # Round for readability
    pw_df[["mw_stat", "mw_p", "mw_FDR"]] = pw_df[["mw_stat", "mw_p", "mw_FDR"]].round(10)

    # Save post-hoc results
    pw_path = os.path.join(out_dir, f"{clus}_pairwise_posthoc_major.csv")
    pw_df.to_csv(pw_path, index=False)

    print(f"Saved post-hoc results: {pw_path}")
#------------------------------------------------------------------------------
#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
from scipy import sparse as sp

# ----------------- CONFIG -----------------
disease_col = "disease"
out_dir = "/work_beegfs/sukmb655/antiTNF/kruskal_global"
os.makedirs(out_dir, exist_ok=True)
# ------------------------------------------


# 1) round adata.X to 3 decimals (sparse- and dense-safe)
if sp.issparse(adata.X):
    X = adata.X.tocoo(copy=True)
    X.data = np.round(X.data, 3)
    adata.X = X.tocsr()
else:
    adata.X = np.round(adata.X, 3)

# 2) get disease labels
diseases = adata.obs[disease_col].astype(str)
uniq_diseases = diseases.unique().tolist()

# need at least 2 diseases
if len(uniq_diseases) < 2:
    raise ValueError("Need at least two disease groups to run Kruskal–Wallis.")

# build boolean masks for each disease once (for speed)
masks = {d: (diseases.values == d) for d in uniq_diseases}

results = []
group_sizes = {d: int(np.sum(m)) for d, m in masks.items()}

# helper to extract a gene column as 1D array
def extract_vector(mat, gene_idx):
    col = mat[:, gene_idx]
    if sp.issparse(col):
        return col.toarray().ravel()
    if isinstance(col, np.matrix):
        return np.asarray(col).ravel()
    return np.asarray(col).ravel()

# 3) KRUSKAL–WALLIS FOR ALL GENES (global, no clusters)
for gene_idx, gene in enumerate(adata.var_names):
    x = extract_vector(adata.X, gene_idx)

    # make per-disease arrays
    groups = [x[masks[d]] for d in uniq_diseases]
    # drop empty groups (just in case)
    groups = [g for g in groups if g.size > 0]
    if len(groups) < 2:
        continue

    # skip genes that are completely constant across all diseases
    flat = np.concatenate(groups)
    if np.allclose(flat, flat[0]):
        continue

    try:
        stat, pval = kruskal(*groups)
        results.append((gene, stat, pval))
    except Exception:
        # numerical edge case, just skip
        continue

# turn into DataFrame
if not results:
    raise RuntimeError("No Kruskal–Wallis results were produced.")

df_kw = pd.DataFrame(results, columns=["gene", "H_statistic", "p_value"]).sort_values("p_value", kind="mergesort")

# 4) Benjamini–Hochberg FDR for KW
m = len(df_kw)
ranks = np.arange(1, m + 1, dtype=float)
bh = df_kw["p_value"].values * m / ranks
bh = np.minimum.accumulate(bh[::-1])[::-1]
df_kw["FDR"] = np.clip(bh, 0, 1)

# round numeric
df_kw = df_kw.round(10)

# save KW results
kw_path = os.path.join(out_dir, "gene_disease_significance_global_major.csv")
df_kw.to_csv(kw_path, index=False)
print(f"Saved Kruskal–Wallis results: {kw_path}")

# also save group sizes (so we remember how many cells per disease)
pd.DataFrame.from_dict(group_sizes, orient="index", columns=["n_cells"])\
    .rename_axis("disease").reset_index()\
    .to_csv(os.path.join(out_dir, "disease_group_sizes_major.csv"), index=False)


# -----------------------------------------------------------------
# 5) POST-HOC: pairwise Mann–Whitney U only for KW-significant genes
# -----------------------------------------------------------------
sig_df = df_kw[df_kw["FDR"] < 0.05].copy()
if sig_df.empty:
    print("No significant genes (KW FDR < 0.05). Skipping post-hoc tests.")
else:
    pairwise_results = []
    disease_pairs = [(d1, d2) for i, d1 in enumerate(uniq_diseases)
                               for d2 in uniq_diseases[i+1:]]

    for _, row in sig_df.iterrows():
        gene = row["gene"]
        gene_idx = adata.var_names.get_loc(gene)

        # expression vector for this gene (all cells)
        col = adata[:, gene_idx].X
        if hasattr(col, "toarray"):
            x = col.toarray().ravel()
        else:
            x = np.asarray(col).ravel()

        for d1, d2 in disease_pairs:
            g1 = x[masks[d1]]
            g2 = x[masks[d2]]

            # optional: skip very small groups
            if g1.size < 5 or g2.size < 5:
                continue

            stat, pval = mannwhitneyu(g1, g2, alternative="two-sided")

            pairwise_results.append({
                "gene": gene,
                "disease_1": d1,
                "disease_2": d2,
                "mw_stat": stat,
                "mw_p": pval,
                "n1": g1.size,
                "n2": g2.size
            })

    if pairwise_results:
        pw_df = pd.DataFrame(pairwise_results).sort_values("mw_p")

        # BH FDR across all pairwise tests
        m_pw = len(pw_df)
        ranks_pw = np.arange(1, m_pw + 1, dtype=float)
        bh_pw = pw_df["mw_p"].values * m_pw / ranks_pw
        bh_pw = np.minimum.accumulate(bh_pw[::-1])[::-1]
        pw_df["mw_FDR"] = np.clip(bh_pw, 0, 1)

        pw_df[["mw_stat", "mw_p", "mw_FDR"]] = pw_df[["mw_stat", "mw_p", "mw_FDR"]].round(10)

        pw_path = os.path.join(out_dir, "gene_disease_pairwise_posthoc_global_major.csv")
        pw_df.to_csv(pw_path, index=False)
        print(f"Saved post-hoc pairwise results: {pw_path}")
    else:
        print("No valid pairwise comparisons to save.")
#..............................................................................
#
# reprt the statisitac evaluation
#
#..............................................................................
import os
import pandas as pd

input_dir = "/work_beegfs/sukmb655/antiTNF/kruskal_per_cluster"
out_dir_A = os.path.join(input_dir, "filtered_mode_A_outlier")
out_dir_B = os.path.join(input_dir, "filtered_mode_B_outlier")
os.makedirs(out_dir_A, exist_ok=True)
os.makedirs(out_dir_B, exist_ok=True)

# helper: given a df of pairwise rows for ONE gene, find if one disease is sig vs all others
def find_single_outlier_disease(gene_df, sig_mask_col="sig"):
    """
    gene_df: rows of one gene, columns: disease_1, disease_2, and a boolean sig column
    returns: outlier_disease or None
    """
    # collect all diseases in this gene
    diseases = set(gene_df["disease_1"]).union(set(gene_df["disease_2"]))
    # for each disease, count sig comparisons involving it
    for d in diseases:
        # rows where this disease participates
        rows_d = gene_df[(gene_df["disease_1"] == d) | (gene_df["disease_2"] == d)]
        # rows among others (not including d)
        rows_others = gene_df[~((gene_df["disease_1"] == d) | (gene_df["disease_2"] == d))]

        # number of other diseases
        n_other = len(diseases) - 1

        # condition 1: this disease must be significant vs ALL others
        sig_count_d = rows_d[sig_mask_col].sum()

        if sig_count_d != n_other:
            continue

        # condition 2: among the other diseases, none of their pairwise comparisons should be significant
        if rows_others[sig_mask_col].any():
            continue

        # if both conditions satisfied, we found the outlier
        return d

    return None

# iterate over all cluster post-hoc files
for file in os.listdir(input_dir):
    if not file.endswith("_pairwise_posthoc_major.csv"):
        continue

    clus = file.replace("_pairwise_posthoc_major.csv", "")
    path = os.path.join(input_dir, file)
    df = pd.read_csv(path)

    # ------------------------------------------------------------------
    # MODE A: mw_p < 0.05 and mw_FDR < 0.05
    # ------------------------------------------------------------------
    df_A = df.copy()
    df_A["sig"] = (df_A["mw_p"] < 0.05) & (df_A["mw_FDR"] < 0.05)

    genes_A = []
    details_rows_A = []

    for gene, gdf in df_A.groupby("gene"):
        outlier = find_single_outlier_disease(gdf, sig_mask_col="sig")
        if outlier is not None:
            genes_A.append({"gene": gene, "disease_outlier": outlier})
            # keep only the significant rows for that gene (the 4 rows)
            details_rows_A.append(gdf[gdf["sig"]])

    if genes_A:
        # summary: gene + outlier disease
        summary_A = pd.DataFrame(genes_A)
        summary_path_A = os.path.join(out_dir_A, f"{clus}_summary_outlier_A.csv")
        summary_A.to_csv(summary_path_A, index=False)

        # details: all rows for all such genes (concatenated)
        details_A = pd.concat(details_rows_A, ignore_index=True)
        details_path_A = os.path.join(out_dir_A, f"{clus}_details_outlier_A.csv")
        details_A.to_csv(details_path_A, index=False)

        print(f"[Mode A] {clus}: {len(genes_A)} genes with single-disease outlier pattern")

    # ------------------------------------------------------------------
    # MODE B: mw_p < 0.05/161060 and mw_FDR < 0.05
    # ------------------------------------------------------------------
    strict_p = 0.05 / 161060
    df_B = df.copy()
    df_B["sig"] = (df_B["mw_p"] < strict_p) & (df_B["mw_FDR"] < 0.05)

    genes_B = []
    details_rows_B = []

    for gene, gdf in df_B.groupby("gene"):
        outlier = find_single_outlier_disease(gdf, sig_mask_col="sig")
        if outlier is not None:
            genes_B.append({"gene": gene, "disease_outlier": outlier})
            details_rows_B.append(gdf[gdf["sig"]])

    if genes_B:
        summary_B = pd.DataFrame(genes_B)
        summary_path_B = os.path.join(out_dir_B, f"{clus}_summary_outlier_B.csv")
        summary_B.to_csv(summary_path_B, index=False)

        details_B = pd.concat(details_rows_B, ignore_index=True)
        details_path_B = os.path.join(out_dir_B, f"{clus}_details_outlier_B.csv")
        details_B.to_csv(details_path_B, index=False)

        print(f"[Mode B] {clus}: {len(genes_B)} genes with single-disease outlier pattern")
