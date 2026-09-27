import json
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from scrna_expression_eval.cli import main


def make_tiny_h5ad(path):
    rows = []
    values = []

    diseases = {
        "Healthy": ["H1", "H2", "H3"],
        "UC_inflamed": ["U1", "U2", "U3"],
    }
    clusters = ["Epithelial", "Myeloid"]

    for cluster_index, cluster in enumerate(clusters):
        for disease, samples in diseases.items():
            for sample_index, sample in enumerate(samples):
                for cell_index in range(2):
                    rows.append(
                        {
                            "sample_id": sample,
                            "disease": disease,
                            "major_cluster": cluster,
                        }
                    )
                    base = (
                        4.0
                        if disease == "UC_inflamed"
                        else 0.0
                    )
                    values.append(
                        [
                            base + 0.1 * cell_index,
                            float(cluster_index + sample_index),
                            1.0,
                        ]
                    )

    adata = ad.AnnData(
        X=sparse.csr_matrix(
            np.asarray(values, dtype=np.float32)
        ),
        obs=pd.DataFrame(rows),
        var=pd.DataFrame(
            index=["SPOCD1", "GENE2", "CONST"]
        ),
    )
    adata.write_h5ad(path)


def test_cli_writes_sample_level_artifacts(tmp_path, monkeypatch):
    input_path = tmp_path / "tiny.h5ad"
    outdir = tmp_path / "outputs"
    make_tiny_h5ad(input_path)

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "scrna-expression-evaluate",
            "--input",
            str(input_path),
            "--outdir",
            str(outdir),
            "--sample-col",
            "sample_id",
            "--min-samples-per-group",
            "3",
        ],
    )

    assert main() == 0

    expected = {
        "cell_count_per_disease_and_celltype.csv",
        "sample_count_per_disease.csv",
        "cell_type_file_map.csv",
        "analysis_summary.csv",
        "run_metadata.json",
    }
    assert expected.issubset(
        {path.name for path in outdir.iterdir()}
    )

    summary = pd.read_csv(
        outdir / "analysis_summary.csv"
    )
    assert set(summary["status"]) == {"ok"}
    assert set(summary["n_biological_samples"]) == {6}

    metadata = json.loads(
        (outdir / "run_metadata.json").read_text(
            encoding="utf-8"
        )
    )
    assert (
        metadata["method"]["inference_unit"]
        == "biological sample/donor"
    )
    assert metadata["input"]["n_cells"] == 24
    assert metadata["input"]["n_biological_samples"] == 6
    assert metadata["input"]["sha256"]


def test_cli_requires_sample_column(tmp_path, monkeypatch):
    input_path = tmp_path / "tiny.h5ad"
    make_tiny_h5ad(input_path)

    adata = ad.read_h5ad(input_path)
    del adata.obs["sample_id"]
    adata.write_h5ad(input_path)

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "scrna-expression-evaluate",
            "--input",
            str(input_path),
        ],
    )

    import pytest

    with pytest.raises(ValueError, match="sample_id"):
        main()
