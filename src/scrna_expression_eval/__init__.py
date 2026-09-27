"""Leakage-aware, sample-level scRNA-seq expression evaluation utilities."""

from .core import (
    AnalysisConfig,
    PseudobulkResult,
    benjamini_hochberg,
    build_pseudobulk_means,
    effect_size_table,
    kruskal_wallis_table,
    pairwise_mannwhitney_table,
    select_target_separated_genes,
    summarize_expression,
    validate_metadata,
)

__all__ = [
    "AnalysisConfig",
    "PseudobulkResult",
    "benjamini_hochberg",
    "build_pseudobulk_means",
    "effect_size_table",
    "kruskal_wallis_table",
    "pairwise_mannwhitney_table",
    "select_target_separated_genes",
    "summarize_expression",
    "validate_metadata",
]
