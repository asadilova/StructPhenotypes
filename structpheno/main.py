"""Command-line entry point for StructPhenotypes."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Callable

try:
    from .get_alpha_fold import AlphaFoldRetriever
    from .get_alphamissense import AlphaMissenseRetriever
    from .get_clinvar import ClinVarRetriever
    from .visualize import visualize
except ImportError:
    from get_alpha_fold import AlphaFoldRetriever
    from get_alphamissense import AlphaMissenseRetriever
    from get_clinvar import ClinVarRetriever
    from visualize import visualize


def _retrieve_source(source_name: str, retrieve: Callable[[], Any]) -> dict[str, Any]:
    """Run one data-source retriever without letting it crash the whole report."""
    try:
        return {
            "status": "ok",
            "data": retrieve(),
            "error": None,
        }
    except Exception as exc:
        return {
            "status": "error",
            "data": [],
            "error": f"{source_name} retrieval failed: {exc}",
        }


def build_report(gene: str) -> dict[str, Any]:
    """Build one combined report for a gene symbol."""
    gene = gene.strip()

    clinvar = _retrieve_source("ClinVar", lambda: ClinVarRetriever(gene).get_clinvar_data())

    alpha_fold = _retrieve_source("AlphaFold", lambda: AlphaFoldRetriever(gene).get_alpha_fold_data())

    alpha_missense = _retrieve_source("AlphaMissense", lambda: AlphaMissenseRetriever(gene).get_alpha_missense_data())

    report = {
        "gene": gene,
        "clinvar": clinvar,
        "alpha_fold": alpha_fold,
        "alpha_missense": alpha_missense,
    }
    report["visualization"] = _retrieve_source("Visualization", lambda: visualize(report))
    return report


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Retrieve StructPhenotypes data for a gene.",
    )
    parser.add_argument(
        "gene",
        help="Gene symbol to query, for example SCN2A.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        help="Optional path to write the JSON report.",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Print indented JSON to stdout.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the StructPhenotypes CLI."""
    args = parse_args(argv)
    report = build_report(args.gene)
    json_indent = 2 if args.pretty else None
    json_text = json.dumps(report, indent=json_indent)

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json_text + "\n", encoding="utf-8")
    else:
        print(json_text)

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
