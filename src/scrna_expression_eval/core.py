from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import sparse as sp
from scipy.stats import kruskal, mannwhitneyu


@dataclass(frozen=True)
class AnalysisConfig:
    sample_col: str = "sample_id"
    disease_col: str = "disease"
    cluster_col: str = "major_cluster"
    target_disease: str = "UC_inflamed"
    min_samples_per_group: int = 3
    fdr_threshold: float = 0.05
    effect_threshold: float = 0.25

    def __post_init__(self) -> None:
        if self.min_samples_per_group < 2:
            raise ValueError("min_samples_per_group must be at least 2")
        if not 0.0 < self.fdr_threshold < 1.0:
            raise ValueError("fdr_threshold must be in (0, 1)")
        if self.effect_threshold < 0.0:
            raise ValueError("effect_threshold must be non-negative")


@dataclass(frozen=True)
class PseudobulkResult:
    expression: np.ndarray
    metadata: pd.DataFrame
    gene_names: list[str]


def validate_metadata(
    obs: pd.DataFrame,
    *,
    sample_col: str,
    disease_col: str,
    cluster_col: str,
) -> pd.DataFrame:
    required = {sample_col, disease_col, cluster_col}
    missing = required - set(obs.columns)
    if missing:
        raise ValueError(
            "AnnData obs is missing required columns: "
            + ", ".join(sorted(missing))
        )

    metadata = obs[
        [sample_col, disease_col, cluster_col]
    ].copy()

    for column in [sample_col, disease_col, cluster_col]:
        if metadata[column].isna().any():
            raise ValueError(
                f"Metadata column {column!r} contains missing values."
            )
        metadata[column] = metadata[column].astype(str).str.strip()
        if metadata[column].eq("").any():
            raise ValueError(
                f"Metadata column {column!r} contains empty values."
            )

    disease_counts = (
        metadata.groupby(sample_col, sort=False)[disease_col]
        .nunique()
    )
    inconsistent = disease_counts[disease_counts > 1]
    if len(inconsistent):
        examples = ", ".join(
            map(str, inconsistent.index[:10])
        )
        raise ValueError(
            "Each biological sample must map to exactly one disease. "
            f"Inconsistent sample IDs include: {examples}"
        )

    return metadata


def benjamini_hochberg(
    p_values: Sequence[float],
) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    if values.ndim != 1:
        raise ValueError("p_values must be one-dimensional")
    if values.size == 0:
        return np.asarray([], dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("p_values must be finite")
    if np.any(values < 0.0) or np.any(values > 1.0):
        raise ValueError("p_values must lie in [0, 1]")

    order = np.argsort(values, kind="stable")
    ranked = values[order]
    ranks = np.arange(1, values.size + 1, dtype=float)
    adjusted = ranked * values.size / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)

    output = np.empty_like(adjusted)
    output[order] = adjusted
    return output


def _as_2d_matrix(matrix: object) -> np.ndarray | sp.spmatrix:
    if sp.issparse(matrix):
        if matrix.ndim != 2:
            raise ValueError("Expression matrix must be two-dimensional")
        if not np.isfinite(matrix.data).all():
            raise ValueError("Expression matrix contains non-finite values")
        return matrix

    array = np.asarray(matrix)
    if array.ndim != 2:
        raise ValueError("Expression matrix must be two-dimensional")
    if not np.isfinite(array).all():
        raise ValueError("Expression matrix contains non-finite values")
    return array


def build_pseudobulk_means(
    matrix: object,
    obs: pd.DataFrame,
    gene_names: Sequence[str],
    *,
    sample_col: str,
    disease_col: str,
    cluster_col: str,
    cluster_value: str,
) -> PseudobulkResult:
    expression = _as_2d_matrix(matrix)
    metadata = validate_metadata(
        obs,
        sample_col=sample_col,
        disease_col=disease_col,
        cluster_col=cluster_col,
    )

    if expression.shape[0] != len(metadata):
        raise ValueError(
            "Expression row count must match obs row count."
        )
    if expression.shape[1] != len(gene_names):
        raise ValueError(
            "Expression column count must match gene_names length."
        )

    cluster_mask = (
        metadata[cluster_col].to_numpy()
        == str(cluster_value)
    )
    if not np.any(cluster_mask):
        raise ValueError(
            f"No cells found for cluster {cluster_value!r}."
        )

    cluster_obs = metadata.loc[cluster_mask].copy()
    sample_ids = sorted(
        cluster_obs[sample_col].unique().tolist()
    )
    sample_to_index = {
        sample_id: index
        for index, sample_id in enumerate(sample_ids)
    }
    codes = cluster_obs[sample_col].map(sample_to_index).to_numpy(
        dtype=int
    )
    n_samples = len(sample_ids)
    cell_counts = np.bincount(
        codes,
        minlength=n_samples,
    ).astype(float)

    cluster_expression = expression[cluster_mask]

    if sp.issparse(cluster_expression):
        rows = codes
        columns = np.arange(len(codes), dtype=int)
        aggregator = sp.csr_matrix(
            (
                np.ones(len(codes), dtype=float),
                (rows, columns),
            ),
            shape=(n_samples, len(codes)),
        )
        sums = aggregator @ cluster_expression
        means = sums.multiply(
            (1.0 / cell_counts)[:, None]
        ).toarray()
    else:
        dense = np.asarray(cluster_expression, dtype=float)
        sums = np.zeros(
            (n_samples, dense.shape[1]),
            dtype=float,
        )
        np.add.at(sums, codes, dense)
        means = sums / cell_counts[:, None]

    disease_lookup = (
        metadata[
            [sample_col, disease_col]
        ]
        .drop_duplicates()
        .set_index(sample_col)[disease_col]
        .to_dict()
    )
    result_metadata = pd.DataFrame(
        {
            sample_col: sample_ids,
            disease_col: [
                disease_lookup[sample_id]
                for sample_id in sample_ids
            ],
            cluster_col: str(cluster_value),
            "n_cells": cell_counts.astype(int),
        }
    )

    return PseudobulkResult(
        expression=np.asarray(means, dtype=float),
        metadata=result_metadata,
        gene_names=[str(gene) for gene in gene_names],
    )


def eligible_diseases(
    metadata: pd.DataFrame,
    *,
    disease_col: str,
    min_samples_per_group: int,
) -> tuple[list[str], dict[str, int]]:
    counts = (
        metadata[disease_col]
        .astype(str)
        .value_counts()
        .sort_index()
    )
    counts_dict = {
        str(name): int(count)
        for name, count in counts.items()
    }
    eligible = [
        disease
        for disease, count in counts_dict.items()
        if count >= min_samples_per_group
    ]
    return eligible, counts_dict


def summarize_expression(
    pseudobulk: PseudobulkResult,
    *,
    disease_col: str,
) -> pd.DataFrame:
    expression = pseudobulk.expression
    metadata = pseudobulk.metadata
    diseases = sorted(
        metadata[disease_col]
        .astype(str)
        .unique()
        .tolist()
    )

    output = pd.DataFrame(
        {"gene": pseudobulk.gene_names}
    )
    for disease in diseases:
        mask = (
            metadata[disease_col].astype(str).to_numpy()
            == disease
        )
        values = expression[mask]
        output[f"{disease}_n_samples"] = int(mask.sum())
        output[f"{disease}_mean"] = values.mean(axis=0)
        output[f"{disease}_var"] = (
            values.var(axis=0, ddof=1)
            if values.shape[0] > 1
            else np.nan
        )
        output[f"{disease}_min"] = values.min(axis=0)
        output[f"{disease}_max"] = values.max(axis=0)

    return output


def _epsilon_squared(
    h_statistic: float,
    n_samples: int,
    n_groups: int,
) -> float:
    denominator = n_samples - n_groups
    if denominator <= 0:
        return float("nan")
    value = (
        h_statistic - n_groups + 1
    ) / denominator
    return float(max(0.0, value))


def kruskal_wallis_table(
    pseudobulk: PseudobulkResult,
    *,
    disease_col: str,
    min_samples_per_group: int,
) -> tuple[pd.DataFrame, dict[str, object]]:
    eligible, counts = eligible_diseases(
        pseudobulk.metadata,
        disease_col=disease_col,
        min_samples_per_group=min_samples_per_group,
    )
    if len(eligible) < 2:
        raise ValueError(
            "At least two disease groups must meet the minimum "
            "biological-sample threshold for Kruskal-Wallis testing."
        )

    disease_values = (
        pseudobulk.metadata[disease_col]
        .astype(str)
        .to_numpy()
    )
    masks = {
        disease: disease_values == disease
        for disease in eligible
    }
    n_samples = int(
        sum(int(mask.sum()) for mask in masks.values())
    )

    rows: list[dict[str, object]] = []
    for gene_index, gene in enumerate(
        pseudobulk.gene_names
    ):
        groups = [
            pseudobulk.expression[mask, gene_index]
            for mask in masks.values()
        ]
        flat = np.concatenate(groups)

        if np.allclose(flat, flat[0]):
            statistic = 0.0
            p_value = 1.0
        else:
            statistic, p_value = kruskal(*groups)

        rows.append(
            {
                "gene": gene,
                "H_statistic": float(statistic),
                "p_value": float(p_value),
                "epsilon_squared": _epsilon_squared(
                    float(statistic),
                    n_samples=n_samples,
                    n_groups=len(groups),
                ),
                "n_groups": int(len(groups)),
                "n_samples_total": n_samples,
            }
        )

    result = pd.DataFrame(rows)
    result["FDR"] = benjamini_hochberg(
        result["p_value"].to_numpy(dtype=float)
    )
    result = result.sort_values(
        ["FDR", "p_value", "gene"],
        kind="stable",
    ).reset_index(drop=True)

    metadata = {
        "sample_counts_by_disease": counts,
        "eligible_diseases": eligible,
        "excluded_diseases": [
            disease
            for disease in sorted(counts)
            if disease not in eligible
        ],
        "min_samples_per_group": int(min_samples_per_group),
        "n_gene_tests": int(len(result)),
        "fdr_scope": "all tested genes within this cell type",
    }
    return result, metadata


def effect_size_table(
    pseudobulk: PseudobulkResult,
    *,
    disease_col: str,
    min_samples_per_group: int,
) -> pd.DataFrame:
    eligible, _ = eligible_diseases(
        pseudobulk.metadata,
        disease_col=disease_col,
        min_samples_per_group=min_samples_per_group,
    )
    if len(eligible) < 2:
        raise ValueError(
            "At least two disease groups must meet the minimum "
            "sample threshold for an across-disease effect range."
        )

    disease_values = (
        pseudobulk.metadata[disease_col]
        .astype(str)
        .to_numpy()
    )
    disease_means = np.vstack(
        [
            pseudobulk.expression[
                disease_values == disease
            ].mean(axis=0)
            for disease in eligible
        ]
    )
    maximum_index = np.argmax(
        disease_means,
        axis=0,
    )
    minimum_index = np.argmin(
        disease_means,
        axis=0,
    )
    effect_range = (
        disease_means.max(axis=0)
        - disease_means.min(axis=0)
    )

    return pd.DataFrame(
        {
            "gene": pseudobulk.gene_names,
            "effect_size": effect_range.astype(float),
            "max_mean_disease": [
                eligible[index]
                for index in maximum_index
            ],
            "min_mean_disease": [
                eligible[index]
                for index in minimum_index
            ],
            "n_eligible_diseases": len(eligible),
        }
    ).sort_values(
        ["effect_size", "gene"],
        ascending=[False, True],
        kind="stable",
    ).reset_index(drop=True)


def _rank_biserial(
    u_statistic: float,
    n_first: int,
    n_second: int,
) -> float:
    denominator = n_first * n_second
    if denominator <= 0:
        return float("nan")
    return float(
        2.0 * u_statistic / denominator
        - 1.0
    )


def pairwise_mannwhitney_table(
    pseudobulk: PseudobulkResult,
    *,
    disease_col: str,
    min_samples_per_group: int,
) -> tuple[pd.DataFrame, dict[str, object]]:
    eligible, counts = eligible_diseases(
        pseudobulk.metadata,
        disease_col=disease_col,
        min_samples_per_group=min_samples_per_group,
    )
    if len(eligible) < 2:
        raise ValueError(
            "At least two disease groups must meet the minimum "
            "sample threshold for pairwise testing."
        )

    disease_values = (
        pseudobulk.metadata[disease_col]
        .astype(str)
        .to_numpy()
    )
    rows: list[dict[str, object]] = []

    for disease_1, disease_2 in combinations(
        eligible,
        2,
    ):
        first_mask = disease_values == disease_1
        second_mask = disease_values == disease_2
        n_first = int(first_mask.sum())
        n_second = int(second_mask.sum())

        first_matrix = pseudobulk.expression[
            first_mask
        ]
        second_matrix = pseudobulk.expression[
            second_mask
        ]

        contrast_rows: list[dict[str, object]] = []
        for gene_index, gene in enumerate(
            pseudobulk.gene_names
        ):
            first = first_matrix[:, gene_index]
            second = second_matrix[:, gene_index]
            combined = np.concatenate([first, second])
            if np.allclose(combined, combined[0]):
                statistic = n_first * n_second / 2.0
                p_value = 1.0
                test_status = "all_values_identical"
            else:
                statistic, p_value = mannwhitneyu(
                    first,
                    second,
                    alternative="two-sided",
                    method="auto",
                )
                if not np.isfinite(p_value):
                    p_value = 1.0
                    test_status = "nonfinite_p_conservative"
                else:
                    test_status = "ok"

            contrast_rows.append(
                {
                    "gene": gene,
                    "disease_1": disease_1,
                    "disease_2": disease_2,
                    "mw_stat": float(statistic),
                    "mw_p": float(p_value),
                    "rank_biserial": _rank_biserial(
                        float(statistic),
                        n_first,
                        n_second,
                    ),
                    "n1_samples": n_first,
                    "n2_samples": n_second,
                    "test_status": test_status,
                }
            )

        contrast = pd.DataFrame(contrast_rows)
        contrast["mw_FDR"] = benjamini_hochberg(
            contrast["mw_p"].to_numpy(dtype=float)
        )
        rows.extend(
            contrast.to_dict(orient="records")
        )

    result = pd.DataFrame(rows)
    result = result.sort_values(
        [
            "disease_1",
            "disease_2",
            "mw_FDR",
            "mw_p",
            "gene",
        ],
        kind="stable",
    ).reset_index(drop=True)

    metadata = {
        "sample_counts_by_disease": counts,
        "eligible_diseases": eligible,
        "excluded_diseases": [
            disease
            for disease in sorted(counts)
            if disease not in eligible
        ],
        "min_samples_per_group": int(min_samples_per_group),
        "fdr_scope": (
            "genes corrected separately within each disease-pair "
            "contrast and cell type"
        ),
        "n_contrasts": int(
            len(list(combinations(eligible, 2)))
        ),
    }
    return result, metadata


def select_target_separated_genes(
    kw_table: pd.DataFrame,
    effect_table: pd.DataFrame,
    pairwise_table: pd.DataFrame,
    *,
    target_disease: str,
    eligible_diseases: Sequence[str],
    fdr_threshold: float,
    effect_threshold: float,
) -> pd.DataFrame:
    eligible = [str(value) for value in eligible_diseases]
    if target_disease not in eligible:
        return pd.DataFrame(
            columns=[
                "gene",
                "KW_p",
                "KW_FDR",
                "epsilon_squared",
                "effect_size",
                "max_mean_disease",
                "min_mean_disease",
                "n_target_contrasts_required",
                "n_target_contrasts_passed",
            ]
        )

    other_diseases = [
        disease
        for disease in eligible
        if disease != target_disease
    ]
    if not other_diseases:
        return pd.DataFrame()

    merged = (
        kw_table[
            [
                "gene",
                "p_value",
                "FDR",
                "epsilon_squared",
            ]
        ]
        .rename(
            columns={
                "p_value": "KW_p",
                "FDR": "KW_FDR",
            }
        )
        .merge(
            effect_table[
                [
                    "gene",
                    "effect_size",
                    "max_mean_disease",
                    "min_mean_disease",
                ]
            ],
            on="gene",
            how="inner",
            validate="one_to_one",
        )
    )
    candidates = merged[
        (merged["KW_FDR"] < fdr_threshold)
        & (merged["effect_size"] > effect_threshold)
    ].copy()

    selected: list[dict[str, object]] = []
    for row in candidates.itertuples(index=False):
        gene_pairs = pairwise_table[
            pairwise_table["gene"].eq(row.gene)
            & (
                pairwise_table["disease_1"].eq(target_disease)
                | pairwise_table["disease_2"].eq(target_disease)
            )
        ].copy()

        passed_diseases: set[str] = set()
        for pair in gene_pairs.itertuples(index=False):
            other = (
                pair.disease_2
                if pair.disease_1 == target_disease
                else pair.disease_1
            )
            if (
                other in other_diseases
                and float(pair.mw_FDR) < fdr_threshold
            ):
                passed_diseases.add(str(other))

        if passed_diseases == set(other_diseases):
            selected.append(
                {
                    "gene": row.gene,
                    "KW_p": float(row.KW_p),
                    "KW_FDR": float(row.KW_FDR),
                    "epsilon_squared": float(row.epsilon_squared),
                    "effect_size": float(row.effect_size),
                    "max_mean_disease": row.max_mean_disease,
                    "min_mean_disease": row.min_mean_disease,
                    "n_target_contrasts_required": len(other_diseases),
                    "n_target_contrasts_passed": len(passed_diseases),
                }
            )

    if not selected:
        return pd.DataFrame(
            columns=[
                "gene",
                "KW_p",
                "KW_FDR",
                "epsilon_squared",
                "effect_size",
                "max_mean_disease",
                "min_mean_disease",
                "n_target_contrasts_required",
                "n_target_contrasts_passed",
            ]
        )

    return pd.DataFrame(selected).sort_values(
        [
            "KW_FDR",
            "effect_size",
            "gene",
        ],
        ascending=[True, False, True],
        kind="stable",
    ).reset_index(drop=True)
