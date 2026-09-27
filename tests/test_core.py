import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from scrna_expression_eval.core import (
    AnalysisConfig,
    benjamini_hochberg,
    build_pseudobulk_means,
    effect_size_table,
    kruskal_wallis_table,
    pairwise_mannwhitney_table,
    select_target_separated_genes,
    summarize_expression,
    validate_metadata,
)


def synthetic_obs():
    rows = []
    for disease, samples in {
        "Healthy": ["H1", "H2", "H3"],
        "UC_inflamed": ["U1", "U2", "U3"],
        "CD_inflamed": ["C1", "C2", "C3"],
    }.items():
        for sample in samples:
            for _ in range(2):
                rows.append(
                    {
                        "sample_id": sample,
                        "disease": disease,
                        "major_cluster": "Epithelial",
                    }
                )
    return pd.DataFrame(rows)


def synthetic_matrix():
    obs = synthetic_obs()
    values = np.zeros((len(obs), 3), dtype=float)

    for row_index, row in obs.iterrows():
        if row["disease"] == "Healthy":
            values[row_index] = [0.0, 1.0, 1.0]
        elif row["disease"] == "UC_inflamed":
            values[row_index] = [5.0, 1.0, 2.0]
        else:
            values[row_index] = [0.5, 1.0, 3.0]
    return values


def test_metadata_rejects_sample_with_multiple_diseases():
    obs = pd.DataFrame(
        {
            "sample_id": ["S1", "S1"],
            "disease": ["Healthy", "UC_inflamed"],
            "major_cluster": ["A", "A"],
        }
    )

    with pytest.raises(ValueError, match="exactly one disease"):
        validate_metadata(
            obs,
            sample_col="sample_id",
            disease_col="disease",
            cluster_col="major_cluster",
        )


def test_bh_known_values():
    adjusted = benjamini_hochberg(
        [0.01, 0.04, 0.03, 0.002]
    )
    np.testing.assert_allclose(
        adjusted,
        [0.02, 0.04, 0.04, 0.008],
    )


@pytest.mark.parametrize("as_sparse", [False, True])
def test_pseudobulk_means_use_biological_samples(as_sparse):
    obs = synthetic_obs()
    matrix = synthetic_matrix()
    if as_sparse:
        matrix = sparse.csr_matrix(matrix)

    result = build_pseudobulk_means(
        matrix,
        obs,
        ["G1", "G2", "G3"],
        sample_col="sample_id",
        disease_col="disease",
        cluster_col="major_cluster",
        cluster_value="Epithelial",
    )

    assert result.expression.shape == (9, 3)
    assert len(result.metadata) == 9
    assert set(result.metadata["n_cells"]) == {2}

    healthy = (
        result.metadata["disease"].to_numpy()
        == "Healthy"
    )
    np.testing.assert_allclose(
        result.expression[healthy, 0],
        0.0,
    )


def test_summary_and_nonparametric_tables_are_sample_level():
    result = build_pseudobulk_means(
        synthetic_matrix(),
        synthetic_obs(),
        ["G1", "G2", "G3"],
        sample_col="sample_id",
        disease_col="disease",
        cluster_col="major_cluster",
        cluster_value="Epithelial",
    )

    stats = summarize_expression(
        result,
        disease_col="disease",
    )
    assert (
        stats.loc[
            stats["gene"] == "G1",
            "UC_inflamed_n_samples",
        ].iloc[0]
        == 3
    )

    kw, metadata = kruskal_wallis_table(
        result,
        disease_col="disease",
        min_samples_per_group=3,
    )
    assert metadata["eligible_diseases"] == [
        "CD_inflamed",
        "Healthy",
        "UC_inflamed",
    ]
    assert set(kw["n_samples_total"]) == {9}

    effect = effect_size_table(
        result,
        disease_col="disease",
        min_samples_per_group=3,
    )
    g1 = effect.loc[
        effect["gene"] == "G1"
    ].iloc[0]
    assert g1["effect_size"] == pytest.approx(5.0)
    assert g1["max_mean_disease"] == "UC_inflamed"

    pairwise, pair_meta = pairwise_mannwhitney_table(
        result,
        disease_col="disease",
        min_samples_per_group=3,
    )
    assert pair_meta["n_contrasts"] == 3
    assert set(pairwise["n1_samples"]) == {3}
    assert set(pairwise["n2_samples"]) == {3}


def test_target_selection_requires_every_eligible_target_contrast():
    kw = pd.DataFrame(
        {
            "gene": ["G1", "G2"],
            "p_value": [1e-4, 1e-4],
            "FDR": [0.01, 0.01],
            "epsilon_squared": [0.8, 0.8],
        }
    )
    effects = pd.DataFrame(
        {
            "gene": ["G1", "G2"],
            "effect_size": [2.0, 2.0],
            "max_mean_disease": ["UC_inflamed", "UC_inflamed"],
            "min_mean_disease": ["Healthy", "Healthy"],
        }
    )
    pairwise = pd.DataFrame(
        [
            {
                "gene": "G1",
                "disease_1": "Healthy",
                "disease_2": "UC_inflamed",
                "mw_FDR": 0.01,
            },
            {
                "gene": "G1",
                "disease_1": "CD_inflamed",
                "disease_2": "UC_inflamed",
                "mw_FDR": 0.02,
            },
            {
                "gene": "G2",
                "disease_1": "Healthy",
                "disease_2": "UC_inflamed",
                "mw_FDR": 0.01,
            },
            {
                "gene": "G2",
                "disease_1": "CD_inflamed",
                "disease_2": "UC_inflamed",
                "mw_FDR": 0.20,
            },
        ]
    )

    selected = select_target_separated_genes(
        kw,
        effects,
        pairwise,
        target_disease="UC_inflamed",
        eligible_diseases=[
            "Healthy",
            "UC_inflamed",
            "CD_inflamed",
        ],
        fdr_threshold=0.05,
        effect_threshold=0.25,
    )

    assert selected["gene"].tolist() == ["G1"]


def test_analysis_config_rejects_too_few_samples_per_group():
    with pytest.raises(ValueError, match="at least 2"):
        AnalysisConfig(
            min_samples_per_group=1
        )
