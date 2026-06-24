"""Visualization helpers for StructPhenotypes.

This module is intentionally a stub for now. The future implementation should
build an interactive py3Dmol view from the combined report.
"""

from __future__ import annotations

from typing import Any


def visualize(report: dict[str, Any]) -> dict[str, Any]:
    """Prepare visualization metadata for a StructPhenotypes report.

    Args:
        report: Combined report from ``main.build_report``.

    Returns:
        A JSON-serializable placeholder describing the intended viewer.
    """
    gene = report.get("gene", "unknown")

    return {
        "status": "not implemented yet",
        "viewer": "py3Dmol",
        "gene": gene,
        "data_sources": {
            "clinvar": _source_has_data(report, "clinvar"),
            "alpha_fold": _source_has_data(report, "alpha_fold"),
            "alpha_missense": _source_has_data(report, "alpha_missense"),
        },
        "planned_layers": [
            "ClinVar phenotype colors for pathogenic variants",
            "AlphaMissense prediction gradient across amino-acid positions",
            "Interactive 3D protein structure loaded from AlphaFold coordinates",
        ],
    }


def _source_has_data(report: dict[str, Any], source_name: str) -> bool:
    """Return whether a report source currently contains non-empty data."""
    source = report.get(source_name, {})
    if not isinstance(source, dict):
        return False

    data = source.get("data", [])
    return bool(data)
