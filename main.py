#!/usr/bin/env python3
"""
@author: seiranai


Pipeline for per-cell-type disease comparisons on single-cell data.

Steps:
1. Count cells per (major_cluster, disease)
2. Per cell type: compute per-disease gene stats (min, max, mean, var)
3. Per cell type: Kruskal–Wallis per gene across diseases
4. Per cell type: effect size (max(mean) - min(mean)) and category
5. Per cell type: pairwise post-hoc Mann–Whitney tests + BH FDR
6. Per cell type: select genes where UC_inflamed is separated from all others
   under strict criteria, save per-cell-type file.

Assumptions:
- adata has .obs columns: 'major_cluster' and 'disease'
- adata.X can be dense or sparse
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse as sp
from scipy.stats import kruskal, mannwhitneyu
from itertools import combinations

# ---------------- CONFIG ----------------
ADATA_PATH = "data/input.h5ad"  # change to your .h5ad path
OUT_DIR = "outputs"
CELL_COUNTS_CSV = os.path.join(OUT_DIR, "cell_type_count_per_disease.csv")
STATS_DIR = os.path.join(OUT_DIR, "stats_per_celltype")
KW_DIR = os.path.join(OUT_DIR, "kw_per_celltype")
EFFECT_DIR = os.path.join(OUT_DIR, "effectsize_per_celltype")
POSTHOC_DIR = os.path.join(OUT_DIR, "posthoc_per_celltype")
UCSEP_DIR = os.path.join(OUT_DIR, "uc_inflamed_separated")

CLUSTER_COL = "major_cluster"
DISEASE_COL = "disease"
UC_DISEASE_NAME = "UC_inflamed"

# this is what you specified: 16,106 genes * 10 cell types
BONF_TOTAL_TESTS = 161060
STRICT_P = 0.05 / BONF_TOTAL_TESTS  # ≈ 3.1e-07
FDR_THRESH = 0.05
EFFECT_LOW = 0.25
# -----------------------------------------


def ensure_dirs():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(STATS_DIR, exist_ok=True)
    os.makedirs(KW_DIR, exist_ok=True)
    os.makedirs(EFFECT_DIR, exist_ok=True)
    os.makedirs(POSTHOC_DIR, exist_ok=True)
    os.makedirs(UCSEP_DIR, exist_ok=True)


def extract_vector(mat, gene_idx):
    """Return 1D numpy array for column gene_idx from dense or sparse matrix."""
    col = mat[:, gene_idx]
    if sp.issparse(col):
        return col.toarray().ravel()
    if isinstance(col, np.matrix):
        return np.asarray(col).ravel()
    return np.asarray(col).ravel()


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """Return BH-FDR for a 1D array of p-values."""
    m = len(pvals)
    order = np.argsort(pvals)
    ranks = np.arange(1, m + 1)
    fdr = pvals[order] * m / ranks
    # enforce monotonicity
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    fdr_out = np.empty_like(fdr)
    fdr_out[order] = np.clip(fdr, 0, 1)
    return fdr_out


def main():
    ensure_dirs()

    print(f"Loading AnnData from: {ADATA_PATH}")
    adata = sc.read_h5ad(ADATA_PATH)
    genes = adata.var_names.to_list()

    # -------------------------------------------------
    # Step 1: cell counts per (cell type, disease)
    # -------------------------------------------------
    print("Step 1: counting cells per (cell type, disease)...")
    counts = (
        adata.obs[[CLUSTER_COL, DISEASE_COL]]
        .value_counts()
        .reset_index(name="n_cells")
        .sort_values([CLUSTER_COL, DISEASE_COL])
    )
    counts.to_csv(CELL_COUNTS_CSV, index=False)
    print(f"Saved: {CELL_COUNTS_CSV}")

    # reusable lists
    all_celltypes = adata.obs[CLUSTER_COL].astype(str).unique().tolist()
    all_diseases = adata.obs[DISEASE_COL].astype(str).unique().tolist()

    # -------------------------------------------------
    # Step 2: per cell type gene stats per disease
    # -------------------------------------------------
    print("Step 2: computing per-gene stats per disease per cell type...")
    for ct in all_celltypes:
        ad_ct = adata[adata.obs[CLUSTER_COL].astype(str) == ct]
        if ad_ct.n_obs == 0:
            continue

        # df to collect stats
        stats_df = pd.DataFrame(index=genes)
        for dis in all_diseases:
            ad_sub = ad_ct[ad_ct.obs[DISEASE_COL].astype(str) == dis]
            if ad_sub.n_obs == 0:
                # fill with NaN columns to keep shape
                stats_df[f"{dis}_mean"] = np.nan
                stats_df[f"{dis}_var"] = np.nan
                stats_df[f"{dis}_min"] = np.nan
                stats_df[f"{dis}_max"] = np.nan
                continue

            X = ad_sub.X
            if sp.issparse(X):
                X = X.toarray()
            # axis=0 -> per gene
            stats_df[f"{dis}_mean"] = X.mean(axis=0)
            stats_df[f"{dis}_var"] = X.var(axis=0)
            stats_df[f"{dis}_min"] = X.min(axis=0)
            stats_df[f"{dis}_max"] = X.max(axis=0)

        stats_df.index.name = "gene"
        out_path = os.path.join(STATS_DIR, f"{ct}_stats.csv")
        stats_df.to_csv(out_path)
        print(f"  saved stats for {ct}: {out_path}")

    # -------------------------------------------------
    # Step 3: KW per gene per cell type
    # -------------------------------------------------
    print("Step 3: running Kruskal–Wallis per gene per cell type...")
    for ct in all_celltypes:
        ad_ct = adata[adata.obs[CLUSTER_COL].astype(str) == ct]
        if ad_ct.n_obs == 0:
            continue

        diseases_ct = ad_ct.obs[DISEASE_COL].astype(str)
        uniq_dis = diseases_ct.unique().tolist()
        if len(uniq_dis) < 2:
            continue

        # precompute masks
        masks = {d: (diseases_ct.values == d) for d in uniq_dis}

        kw_results = []
        for g_idx, g in enumerate(genes):
            x = extract_vector(ad_ct.X, g_idx)
            groups = [x[masks[d]] for d in uniq_dis if np.sum(masks[d]) > 0]
            if len(groups) < 2:
                continue
            flat = np.concatenate(groups)
            # if all values are same, skip
            if np.allclose(flat, flat[0]):
                continue
            try:
                stat, pval = kruskal(*groups)
                kw_results.append((g, stat, pval))
            except Exception:
                continue

        if not kw_results:
            continue

        kw_df = pd.DataFrame(kw_results, columns=["gene", "H_statistic", "p_value"])
        # BH FDR
        kw_df = kw_df.sort_values("p_value", kind="mergesort")
        kw_df["FDR"] = benjamini_hochberg(kw_df["p_value"].values)
        out_path = os.path.join(KW_DIR, f"{ct}_KW.csv")
        kw_df.to_csv(out_path, index=False)
        print(f"  saved KW for {ct}: {out_path}")

    # -------------------------------------------------
    # Step 4: effect size per gene per cell type
    #         effect_size = max(mean) - min(mean)
    # -------------------------------------------------
    print("Step 4: computing effect sizes...")
    for ct in all_celltypes:
        stats_path = os.path.join(STATS_DIR, f"{ct}_stats.csv")
        if not os.path.exists(stats_path):
            continue
        stats_df = pd.read_csv(stats_path, index_col=0)

        mean_cols = [c for c in stats_df.columns if c.endswith("_mean")]
        if not mean_cols:
            continue

        eff_df = pd.DataFrame(index=stats_df.index)
        eff_df["effect_size"] = stats_df[mean_cols].max(axis=1) - stats_df[mean_cols].min(axis=1)

        # classify
        conds = [
            (eff_df["effect_size"] < 0.25),
            (eff_df["effect_size"] >= 0.25) & (eff_df["effect_size"] < 0.5),
            (eff_df["effect_size"] >= 0.5),
        ]
        choices = ["low", "medium", "high"]
        eff_df["effect_size_class"] = np.select(conds, choices, default="low")

        eff_df.index.name = "gene"
        out_path = os.path.join(EFFECT_DIR, f"{ct}_effectsize.csv")
        eff_df.to_csv(out_path)
        print(f"  saved effect sizes for {ct}: {out_path}")

    # -------------------------------------------------
    # Step 5: pairwise post-hoc per gene per cell type
    # -------------------------------------------------
    print("Step 5: pairwise post-hoc (Mann–Whitney) per cell type...")
    for ct in all_celltypes:
        ad_ct = adata[adata.obs[CLUSTER_COL].astype(str) == ct]
        if ad_ct.n_obs == 0:
            continue

        diseases_ct = ad_ct.obs[DISEASE_COL].astype(str)
        uniq_dis = diseases_ct.unique().tolist()
        if len(uniq_dis) < 2:
            continue

        masks = {d: (diseases_ct.values == d) for d in uniq_dis}
        pairs = list(combinations(uniq_dis, 2))
        pw_rows = []

        for g_idx, g in enumerate(genes):
            x = extract_vector(ad_ct.X, g_idx)
            for d1, d2 in pairs:
                g1 = x[masks[d1]]
                g2 = x[masks[d2]]
                if g1.size < 2 or g2.size < 2:
                    continue
                try:
                    stat, pval = mannwhitneyu(g1, g2, alternative="two-sided")
                except Exception:
                    continue
                pw_rows.append(
                    {
                        "gene": g,
                        "disease_1": d1,
                        "disease_2": d2,
                        "mw_stat": stat,
                        "mw_p": pval,
                        "n1": int(g1.size),
                        "n2": int(g2.size),
                    }
                )

        if not pw_rows:
            continue

        pw_df = pd.DataFrame(pw_rows)
        pw_df = pw_df.sort_values("mw_p", kind="mergesort")
        pw_df["mw_FDR"] = benjamini_hochberg(pw_df["mw_p"].values)

        out_path = os.path.join(POSTHOC_DIR, f"{ct}_pairwise_posthoc.csv")
        pw_df.to_csv(out_path, index=False)
        print(f"  saved pairwise posthoc for {ct}: {out_path}")

    # -------------------------------------------------
    # Step 6: select genes where UC_inflamed is separated from all others
    # criteria:
    #  - KW_p < STRICT_P and KW_FDR < 0.05
    #  - effect_size > 0.25
    #  - for the gene: all pairs involving UC_inflamed have mw_p < STRICT_P and mw_FDR < 0.05
    # store rows (long format) with KW and effect size info
    # -------------------------------------------------
    print("Step 6: selecting UC_inflamed-separated genes...")
    for ct in all_celltypes:
        kw_path = os.path.join(KW_DIR, f"{ct}_KW.csv")
        eff_path = os.path.join(EFFECT_DIR, f"{ct}_effectsize.csv")
        pw_path = os.path.join(POSTHOC_DIR, f"{ct}_pairwise_posthoc.csv")

        if not (os.path.exists(kw_path) and os.path.exists(eff_path) and os.path.exists(pw_path)):
            continue

        kw_df = pd.read_csv(kw_path)
        eff_df = pd.read_csv(eff_path)
        pw_df = pd.read_csv(pw_path)

        # filter KW
        kw_df_strict = kw_df[(kw_df["p_value"] < STRICT_P) & (kw_df["FDR"] < FDR_THRESH)]
        if kw_df_strict.empty:
            continue

        eff_df = eff_df.rename(columns={"gene": "gene_name"}) if "gene" not in eff_df.columns else eff_df
        # standardize name
        eff_df = eff_df.rename(columns={"gene": "gene"})

        # merge kw + effect
        kw_eff = kw_df_strict.merge(eff_df, on="gene", how="left")

        selected_rows = []
        for _, row in kw_eff.iterrows():
            g = row["gene"]
            eff_val = row["effect_size"]
            if eff_val is None or np.isnan(eff_val) or eff_val <= EFFECT_LOW:
                continue

            # get all pairwise rows for this gene involving UC_inflamed
            pw_g = pw_df[
                (pw_df["gene"] == g)
                & (
                    (pw_df["disease_1"] == UC_DISEASE_NAME)
                    | (pw_df["disease_2"] == UC_DISEASE_NAME)
                )
            ]
            if pw_g.empty:
                continue

            # what other diseases exist for this gene?
            # from KW we know diseases in general, but here we check pairwise completeness
            # to be strict: all such pairs must pass post-hoc strict criteria
            cond_ok = np.all(
                (pw_g["mw_p"] < STRICT_P) & (pw_g["mw_FDR"] < FDR_THRESH)
            )
            if not cond_ok:
                continue

            # if passed, append each UC pair as one row (long format)
            for _, prow in pw_g.iterrows():
                selected_rows.append(
                    {
                        "gene": g,
                        "cell_type": ct,
                        "KW_p": row["p_value"],
                        "KW_FDR": row["FDR"],
                        "effect_size": row["effect_size"],
                        "effect_size_class": row["effect_size_class"],
                        "disease_1": prow["disease_1"],
                        "disease_2": prow["disease_2"],
                        "mw_stat": prow["mw_stat"],
                        "mw_p": prow["mw_p"],
                        "n1": prow["n1"],
                        "n2": prow["n2"],
                        "mw_FDR": prow["mw_FDR"],
                    }
                )

        if selected_rows:
            sel_df = pd.DataFrame(selected_rows)
            out_path = os.path.join(UCSEP_DIR, f"{ct}_UCinflamed.csv")
            sel_df.to_csv(out_path, index=False)
            print(f"  saved UC_inflamed separated genes for {ct}: {out_path}")

    print("Done.")


if __name__ == "__main__":
    main()
