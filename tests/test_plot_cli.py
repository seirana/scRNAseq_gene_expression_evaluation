import sys

import anndata as ad
import numpy as np
import pandas as pd

from scrna_expression_eval.plot_cli import main


def test_gene_plot_cli_uses_sample_level_values(tmp_path, monkeypatch):
    input_path = tmp_path / "tiny.h5ad"
    outdir = tmp_path / "plots"

    obs = pd.DataFrame(
        {
            "sample_id": [
                "H1",
                "H1",
                "H2",
                "H2",
                "U1",
                "U1",
                "U2",
                "U2",
            ],
            "disease": [
                "Healthy",
                "Healthy",
                "Healthy",
                "Healthy",
                "UC_inflamed",
                "UC_inflamed",
                "UC_inflamed",
                "UC_inflamed",
            ],
            "major_cluster": ["Epithelial"] * 8,
        }
    )
    adata = ad.AnnData(
        X=np.asarray(
            [
                [0.0, 1.0],
                [0.2, 1.0],
                [0.1, 1.0],
                [0.3, 1.0],
                [3.0, 1.0],
                [3.2, 1.0],
                [3.1, 1.0],
                [3.3, 1.0],
            ],
            dtype=np.float32,
        ),
        obs=obs,
        var=pd.DataFrame(
            index=["SPOCD1", "CONST"]
        ),
    )
    adata.write_h5ad(input_path)

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "scrna-expression-plot-gene",
            "--input",
            str(input_path),
            "--gene",
            "SPOCD1",
            "--outdir",
            str(outdir),
        ],
    )

    assert main() == 0
    values = pd.read_csv(
        outdir / "SPOCD1_pseudobulk_values.csv"
    )
    assert len(values) == 4
    assert (
        outdir
        / "SPOCD1_Epithelial_pseudobulk_boxplot.png"
    ).exists()
    assert (
        outdir
        / "SPOCD1_pseudobulk_heatmap.png"
    ).exists()
    assert (
        outdir
        / "SPOCD1_effect_range.png"
    ).exists()
