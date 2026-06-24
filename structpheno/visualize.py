"""Local 3Dmol.js visualization helpers for StructPhenotypes."""

from __future__ import annotations

import json
from html import escape
from pathlib import Path
from typing import Any


EXAMPLE_PDB = """\
ATOM      1  N   MET A   1       0.000   1.200   0.000  1.00 20.00           N
ATOM      2  CA  MET A   1       1.100   0.300   0.400  1.00 20.00           C
ATOM      3  C   MET A   1       2.300   1.100   0.900  1.00 20.00           C
ATOM      4  O   MET A   1       2.300   2.300   0.900  1.00 20.00           O
ATOM      5  N   GLU A   2       3.300   0.400   1.300  1.00 20.00           N
ATOM      6  CA  GLU A   2       4.600   1.000   1.800  1.00 20.00           C
ATOM      7  C   GLU A   2       5.500   0.000   2.500  1.00 20.00           C
ATOM      8  O   GLU A   2       5.200  -1.200   2.500  1.00 20.00           O
ATOM      9  N   LEU A   3       6.600   0.500   3.100  1.00 20.00           N
ATOM     10  CA  LEU A   3       7.600  -0.300   3.800  1.00 20.00           C
ATOM     11  C   LEU A   3       8.900   0.400   3.600  1.00 20.00           C
ATOM     12  O   LEU A   3       9.000   1.600   3.700  1.00 20.00           O
ATOM     13  N   LYS A   4       9.900  -0.400   3.300  1.00 20.00           N
ATOM     14  CA  LYS A   4      11.200   0.100   3.000  1.00 20.00           C
ATOM     15  C   LYS A   4      12.000  -0.900   2.200  1.00 20.00           C
ATOM     16  O   LYS A   4      11.600  -2.000   1.900  1.00 20.00           O
ATOM     17  N   ALA A   5      13.200  -0.400   1.900  1.00 20.00           N
ATOM     18  CA  ALA A   5      14.100  -1.200   1.100  1.00 20.00           C
ATOM     19  C   ALA A   5      15.300  -0.400   0.700  1.00 20.00           C
ATOM     20  O   ALA A   5      15.300   0.800   0.900  1.00 20.00           O
ATOM     21  N   GLY A   6      16.300  -1.100   0.100  1.00 20.00           N
ATOM     22  CA  GLY A   6      17.500  -0.500  -0.400  1.00 20.00           C
ATOM     23  C   GLY A   6      17.700  -0.700  -1.900  1.00 20.00           C
ATOM     24  O   GLY A   6      17.400  -1.700  -2.500  1.00 20.00           O
ATOM     25  N   ARG A   7      18.300   0.300  -2.500  1.00 20.00           N
ATOM     26  CA  ARG A   7      18.600   0.200  -3.900  1.00 20.00           C
ATOM     27  C   ARG A   7      20.100   0.100  -4.100  1.00 20.00           C
ATOM     28  O   ARG A   7      20.900   0.800  -3.500  1.00 20.00           O
ATOM     29  N   VAL A   8      20.400  -0.800  -4.900  1.00 20.00           N
ATOM     30  CA  VAL A   8      21.800  -1.000  -5.300  1.00 20.00           C
ATOM     31  C   VAL A   8      22.000  -0.800  -6.800  1.00 20.00           C
ATOM     32  O   VAL A   8      21.200  -1.200  -7.600  1.00 20.00           O
ATOM     33  N   SER A   9      23.100  -0.100  -7.100  1.00 20.00           N
ATOM     34  CA  SER A   9      23.500   0.100  -8.500  1.00 20.00           C
ATOM     35  C   SER A   9      24.700  -0.700  -8.900  1.00 20.00           C
ATOM     36  O   SER A   9      25.700  -0.700  -8.200  1.00 20.00           O
ATOM     37  N   THR A  10      24.600  -1.400 -10.000  1.00 20.00           N
ATOM     38  CA  THR A  10      25.700  -2.200 -10.500  1.00 20.00           C
ATOM     39  C   THR A  10      26.300  -1.500 -11.700  1.00 20.00           C
ATOM     40  O   THR A  10      25.600  -1.100 -12.600  1.00 20.00           O
ATOM     41  N   ILE A  11      27.600  -1.300 -11.700  1.00 20.00           N
ATOM     42  CA  ILE A  11      28.300  -0.700 -12.800  1.00 20.00           C
ATOM     43  C   ILE A  11      29.400  -1.600 -13.300  1.00 20.00           C
ATOM     44  O   ILE A  11      30.000  -2.300 -12.500  1.00 20.00           O
ATOM     45  N   GLN A  12      29.600  -1.500 -14.600  1.00 20.00           N
ATOM     46  CA  GLN A  12      30.600  -2.200 -15.400  1.00 20.00           C
ATOM     47  C   GLN A  12      31.400  -1.200 -16.200  1.00 20.00           C
ATOM     48  O   GLN A  12      30.900  -0.200 -16.600  1.00 20.00           O
TER
END
"""


EXAMPLE_ANNOTATIONS = {
    "clinvar_variants": [
        {
            "residue": 3,
            "amino_acid": "L",
            "significance": "Pathogenic",
            "phenotype": "Example seizure phenotype",
            "color": "#d73027",
        },
        {
            "residue": 7,
            "amino_acid": "R",
            "significance": "Likely pathogenic",
            "phenotype": "Example developmental phenotype",
            "color": "#fc8d59",
        },
        {
            "residue": 11,
            "amino_acid": "I",
            "significance": "Uncertain significance",
            "phenotype": "Example mixed phenotype",
            "color": "#fee08b",
        },
    ],
    "missense_predictions": [
        {"residue": 1, "score": 0.08},
        {"residue": 2, "score": 0.16},
        {"residue": 3, "score": 0.91},
        {"residue": 4, "score": 0.34},
        {"residue": 5, "score": 0.27},
        {"residue": 6, "score": 0.42},
        {"residue": 7, "score": 0.82},
        {"residue": 8, "score": 0.48},
        {"residue": 9, "score": 0.55},
        {"residue": 10, "score": 0.63},
        {"residue": 11, "score": 0.72},
        {"residue": 12, "score": 0.22},
    ],
}


def visualize(
    report: dict[str, Any],
    output_path: str | Path | None = None,
) -> dict[str, Any]:
    """Prepare or write a local 3Dmol.js viewer for a StructPhenotypes report.

    Args:
        report: Combined report from ``main.build_report``.
        output_path: Optional HTML path. If omitted, no file is written.

    Returns:
        JSON-serializable metadata about the local viewer.
    """
    gene = str(report.get("gene", "unknown"))
    resolved_path = Path(output_path) if output_path else None

    metadata: dict[str, Any] = {
        "status": "stub ready",
        "viewer": "3Dmol.js",
        "gene": gene,
        "html_path": str(resolved_path) if resolved_path else None,
        "html_written": False,
        "example_structure": True,
        "data_sources": {
            "clinvar": _source_has_data(report, "clinvar"),
            "alpha_fold": _source_has_data(report, "alpha_fold"),
            "alpha_missense": _source_has_data(report, "alpha_missense"),
        },
        "controls": [
            "Default protein view",
            "ClinVar phenotype coloring",
            "ClinVar significance coloring",
            "AlphaMissense score gradient",
            "Variant residue focus",
        ],
    }

    if resolved_path:
        resolved_path.parent.mkdir(parents=True, exist_ok=True)
        resolved_path.write_text(_render_html(report), encoding="utf-8")
        metadata["html_path"] = str(resolved_path.resolve())
        metadata["html_written"] = True

    return metadata


def _render_html(report: dict[str, Any]) -> str:
    """Render a standalone HTML viewer around example protein data."""
    gene = str(report.get("gene", "unknown"))
    report_json = _safe_json(report)
    pdb_json = _safe_json(EXAMPLE_PDB)
    annotations_json = _safe_json(EXAMPLE_ANNOTATIONS)

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>StructPhenotypes Viewer - {escape(gene)}</title>
  <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
  <style>
    :root {{
      color-scheme: light;
      font-family: Arial, sans-serif;
      --border: #d8dee8;
      --ink: #1f2937;
      --muted: #607085;
      --panel: #f7f9fc;
      --blue: #2459d6;
    }}
    body {{
      margin: 0;
      color: var(--ink);
      background: #ffffff;
    }}
    header {{
      padding: 18px 22px 12px;
      border-bottom: 1px solid var(--border);
    }}
    h1 {{
      margin: 0 0 6px;
      font-size: 22px;
      font-weight: 700;
    }}
    .subtitle {{
      margin: 0;
      color: var(--muted);
      font-size: 14px;
    }}
    .toolbar {{
      display: flex;
      flex-wrap: wrap;
      gap: 12px;
      align-items: end;
      padding: 14px 22px;
      background: var(--panel);
      border-bottom: 1px solid var(--border);
    }}
    label {{
      display: grid;
      gap: 5px;
      font-size: 12px;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.04em;
      color: var(--muted);
    }}
    select,
    button {{
      min-width: 180px;
      height: 36px;
      border: 1px solid var(--border);
      border-radius: 6px;
      background: #ffffff;
      color: var(--ink);
      font-size: 14px;
      padding: 0 10px;
    }}
    button {{
      min-width: 110px;
      cursor: pointer;
      color: #ffffff;
      background: var(--blue);
      border-color: var(--blue);
      font-weight: 700;
    }}
    #viewer {{
      width: 100%;
      height: 68vh;
      min-height: 430px;
      position: relative;
    }}
    #status {{
      padding: 10px 22px;
      color: var(--muted);
      border-top: 1px solid var(--border);
      border-bottom: 1px solid var(--border);
      font-size: 14px;
    }}
    .details {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(230px, 1fr));
      gap: 14px;
      padding: 16px 22px 24px;
    }}
    .panel {{
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 12px;
    }}
    .panel h2 {{
      margin: 0 0 8px;
      font-size: 15px;
    }}
    .panel p,
    .panel li {{
      color: var(--muted);
      font-size: 14px;
      line-height: 1.4;
    }}
    .panel ul {{
      margin: 0;
      padding-left: 18px;
    }}
    @media (max-width: 700px) {{
      select,
      button {{
        min-width: 100%;
      }}
      .toolbar {{
        align-items: stretch;
      }}
      label {{
        width: 100%;
      }}
    }}
  </style>
</head>
<body>
  <header>
    <h1>StructPhenotypes Viewer: {escape(gene)}</h1>
    <p class="subtitle">Local HTML prototype using 3Dmol.js and an embedded example protein structure.</p>
  </header>

  <div class="toolbar">
    <label>
      Color mode
      <select id="colorMode">
        <option value="default">Default protein view</option>
        <option value="phenotype">ClinVar phenotype</option>
        <option value="significance">ClinVar significance</option>
        <option value="missense">AlphaMissense gradient</option>
      </select>
    </label>
    <label>
      Phenotype filter
      <select id="phenotypeFilter"></select>
    </label>
    <label>
      Focus residue
      <select id="focusResidue"></select>
    </label>
    <button id="resetView" type="button">Reset view</button>
  </div>

  <div id="viewer"></div>
  <div id="status">Loading viewer...</div>

  <section class="details">
    <div class="panel">
      <h2>Current Prototype</h2>
      <p>This page uses a small embedded example protein. Once AlphaFold retrieval is implemented, this same viewer can load the real PDB or mmCIF coordinates.</p>
    </div>
    <div class="panel">
      <h2>Variant Layers</h2>
      <ul id="variantList"></ul>
    </div>
    <div class="panel">
      <h2>Next Data Hooks</h2>
      <p>ClinVar residues will drive phenotype and significance colors. AlphaMissense scores will drive the residue gradient.</p>
    </div>
  </section>

  <script>
    const REPORT = {report_json};
    const PDB_TEXT = {pdb_json};
    const ANNOTATIONS = {annotations_json};

    let viewer = null;

    function init() {{
      populateControls();
      populateVariantList();

      if (!window.$3Dmol) {{
        setStatus("3Dmol.js did not load. Check your internet connection for the viewer library.");
        return;
      }}

      viewer = $3Dmol.createViewer("viewer", {{ backgroundColor: "white" }});
      viewer.addModel(PDB_TEXT, "pdb");
      applyStyle();
      viewer.zoomTo();
      viewer.render();
      setStatus("Example protein loaded. Use the dropdowns to switch annotation layers.");
    }}

    function populateControls() {{
      const phenotypeSelect = document.getElementById("phenotypeFilter");
      const phenotypes = [...new Set(ANNOTATIONS.clinvar_variants.map((variant) => variant.phenotype))];
      phenotypeSelect.innerHTML = '<option value="all">All phenotypes</option>' +
        phenotypes.map((phenotype) => `<option value="${{escapeAttr(phenotype)}}">${{phenotype}}</option>`).join("");

      const residueSelect = document.getElementById("focusResidue");
      residueSelect.innerHTML = '<option value="">No residue focus</option>' +
        ANNOTATIONS.clinvar_variants
          .map((variant) => `<option value="${{variant.residue}}">Residue ${{variant.residue}}: ${{variant.significance}}</option>`)
          .join("");

      document.getElementById("colorMode").addEventListener("change", applyStyle);
      phenotypeSelect.addEventListener("change", applyStyle);
      residueSelect.addEventListener("change", focusResidue);
      document.getElementById("resetView").addEventListener("click", resetView);
    }}

    function populateVariantList() {{
      const list = document.getElementById("variantList");
      list.innerHTML = ANNOTATIONS.clinvar_variants
        .map((variant) => `<li>Residue ${{variant.residue}} (${{variant.amino_acid}}): ${{variant.significance}}, ${{variant.phenotype}}</li>`)
        .join("");
    }}

    function applyStyle() {{
      if (!viewer) {{
        return;
      }}

      const mode = document.getElementById("colorMode").value;
      const phenotype = document.getElementById("phenotypeFilter").value;

      viewer.setStyle({{}}, {{ cartoon: {{ color: "spectrum" }}, stick: {{ radius: 0.08 }} }});

      if (mode === "phenotype") {{
        colorClinVarPhenotypes(phenotype);
        setStatus("Coloring example ClinVar residues by phenotype.");
      }} else if (mode === "significance") {{
        colorClinVarSignificance();
        setStatus("Coloring example ClinVar residues by significance.");
      }} else if (mode === "missense") {{
        colorMissenseGradient();
        setStatus("Coloring example residues by AlphaMissense-like prediction score.");
      }} else {{
        setStatus("Default spectrum protein view.");
      }}

      focusResidue(false);
      viewer.render();
    }}

    function colorClinVarPhenotypes(phenotype) {{
      ANNOTATIONS.clinvar_variants.forEach((variant) => {{
        if (phenotype !== "all" && variant.phenotype !== phenotype) {{
          return;
        }}
        setResidueStyle(variant.residue, variant.color, 0.2);
      }});
    }}

    function colorClinVarSignificance() {{
      ANNOTATIONS.clinvar_variants.forEach((variant) => {{
        setResidueStyle(variant.residue, significanceColor(variant.significance), 0.2);
      }});
    }}

    function colorMissenseGradient() {{
      ANNOTATIONS.missense_predictions.forEach((prediction) => {{
        setResidueStyle(prediction.residue, scoreToColor(prediction.score), 0.12);
      }});
    }}

    function setResidueStyle(residue, color, radius) {{
      viewer.setStyle(
        {{ chain: "A", resi: residue }},
        {{ cartoon: {{ color }}, stick: {{ color, radius }}, sphere: {{ color, radius: radius * 1.8 }} }}
      );
    }}

    function focusResidue(zoom = true) {{
      if (!viewer) {{
        return;
      }}

      const residueValue = document.getElementById("focusResidue").value;
      if (!residueValue) {{
        return;
      }}

      const residue = Number(residueValue);
      viewer.addStyle({{ chain: "A", resi: residue }}, {{ sphere: {{ color: "#2459d6", radius: 0.55 }} }});
      if (zoom) {{
        viewer.zoomTo({{ chain: "A", resi: residue }});
      }}
      viewer.render();
    }}

    function resetView() {{
      if (!viewer) {{
        return;
      }}
      document.getElementById("focusResidue").value = "";
      applyStyle();
      viewer.zoomTo();
      viewer.render();
    }}

    function significanceColor(significance) {{
      const normalized = significance.toLowerCase();
      if (normalized.includes("pathogenic") && !normalized.includes("likely")) {{
        return "#b2182b";
      }}
      if (normalized.includes("likely pathogenic")) {{
        return "#ef8a62";
      }}
      if (normalized.includes("uncertain")) {{
        return "#fddbc7";
      }}
      return "#999999";
    }}

    function scoreToColor(score) {{
      const clamped = Math.max(0, Math.min(1, score));
      const red = Math.round(255 * clamped);
      const green = Math.round(180 * (1 - clamped));
      const blue = Math.round(255 * (1 - clamped));
      return `rgb(${{red}}, ${{green}}, ${{blue}})`;
    }}

    function escapeAttr(value) {{
      return String(value).replaceAll("&", "&amp;").replaceAll('"', "&quot;").replaceAll("<", "&lt;");
    }}

    function setStatus(message) {{
      document.getElementById("status").textContent = message;
    }}

    window.addEventListener("load", init);
  </script>
</body>
</html>
"""


def _safe_json(value: Any) -> str:
    """Serialize data for embedding in a script block."""
    return json.dumps(value).replace("</", "<\\/")


def _source_has_data(report: dict[str, Any], source_name: str) -> bool:
    """Return whether a report source currently contains non-empty data."""
    source = report.get(source_name, {})
    if not isinstance(source, dict):
        return False

    data = source.get("data", [])
    return bool(data)
