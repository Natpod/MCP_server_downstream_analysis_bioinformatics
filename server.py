"""Servidor MCP prototipo para análisis de expresión génica.

Herramientas:
  1. identify_metadata   -> orientación, tipo de ID, escala de valores, grupos, calidad
  2. visualize_expression -> heatmap (clustering, z-score) y boxplot por grupo
  3. enrichment_analysis  -> ORA hipergeométrico contra un GMT local (sin red)

Ejecutar:  python server.py            (stdio, para Claude Desktop / Claude Code)
           python server.py --http     (streamable HTTP en :8000)
"""
from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

try:  # SDK mcp >= 2.0
    from mcp.server.mcpserver import Image, MCPServer as Server
except ImportError:  # SDK mcp 1.x
    from mcp.server.fastmcp import FastMCP as Server, Image

import analysis

BASE = Path(__file__).resolve().parent
DEFAULT_GMT = str(BASE / "data" / "demo_hallmark_subset.gmt")
DEFAULT_OUT = str(BASE / "outputs")

mcp = Server(
    "gene-expression",
    instructions=(
        "Herramientas para matrices de expresión génica. Flujo recomendado: "
        "1) identify_metadata para entender el fichero y los grupos; "
        "2) visualize_expression; 3) enrichment_analysis con una lista de genes "
        "o con case_group/control_group para derivarla por expresión diferencial."
    ),
)


def _run_dir(base: str, tool: str) -> str:
    d = Path(base) / f"{datetime.now():%Y%m%d_%H%M%S_%f}_{tool}"
    d.mkdir(parents=True, exist_ok=True)
    return str(d)


def _as_text(obj: dict) -> str:
    return json.dumps(obj, indent=2, ensure_ascii=False, default=str)


@mcp.tool()
def identify_metadata(expression_path: str, sample_metadata_path: str | None = None) -> str:
    """Inspecciona una matriz de expresión (CSV/TSV/XLSX) y devuelve sus metadatos.

    Detecta orientación (genes x muestras o al revés), tipo de identificador
    (símbolo HGNC/MGI, Ensembl, Entrez, RefSeq, sonda Affymetrix), especie probable,
    escala de los valores (counts, log, CPM/TPM, z-score), problemas de calidad y
    grupos experimentales deducidos de los nombres de muestra o de un fichero de
    metadatos opcional (primera columna = ID de muestra).
    """
    return _as_text(analysis.identify_metadata(expression_path, sample_metadata_path))


@mcp.tool()
def visualize_expression(
    expression_path: str,
    genes: list[str] | None = None,
    top_n_variable: int = 50,
    boxplot_genes: list[str] | None = None,
    groups: dict[str, str] | None = None,
    sample_metadata_path: str | None = None,
    group_column: str | None = None,
    output_dir: str = DEFAULT_OUT,
) -> list:
    """Genera un heatmap con clustering jerárquico y boxplots por grupo.

    - genes: lista de genes para el heatmap; si se omite, usa los top_n_variable más variables.
    - boxplot_genes: genes para el boxplot (por defecto, los 6 primeros del heatmap).
    - Grupos: dict muestra->grupo, o sample_metadata_path + group_column, o deducción automática.
    Aplica log2(x+1) si los datos están en escala lineal y z-score por gen en el heatmap.
    Devuelve resumen JSON y las dos imágenes.
    """
    res = analysis.plot_expression(
        expression_path, _run_dir(output_dir, "viz"), genes, top_n_variable, groups, sample_metadata_path, group_column, boxplot_genes
    )
    return [_as_text(res), Image(path=res["heatmap_path"]), Image(path=res["boxplot_path"])]


@mcp.tool()
def enrichment_analysis(
    gene_list: list[str] | None = None,
    expression_path: str | None = None,
    case_group: str | None = None,
    control_group: str | None = None,
    groups: dict[str, str] | None = None,
    sample_metadata_path: str | None = None,
    group_column: str | None = None,
    gmt_path: str = DEFAULT_GMT,
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 1.0,
    min_set_size: int = 5,
    max_set_size: int = 500,
    output_dir: str = DEFAULT_OUT,
) -> list:
    """Over-representation analysis (hipergeométrico + Benjamini-Hochberg) contra un GMT local.

    Dos modos:
    - gene_list: lista de símbolos a testar. Si además se da expression_path, el fondo
      son los genes medidos en la matriz (recomendado); si no, la unión de los gene sets.
    - expression_path + case_group + control_group (sin gene_list): calcula genes
      diferenciales (Welch t-test sobre log2, BH, umbrales padj_cutoff y |log2FC| >= lfc_cutoff)
      y los usa como entrada.
    gmt_path: cualquier GMT (MSigDB Hallmark, GO, KEGG, Reactome...). Por defecto, un subconjunto demo.
    """
    res = analysis.enrichment(
        gmt_path, _run_dir(output_dir, "enrichment"), gene_list, expression_path, groups, sample_metadata_path, group_column,
        case_group, control_group, padj_cutoff, lfc_cutoff, min_set_size, max_set_size,
    )
    content: list = [_as_text(res)]
    if res.get("dotplot_path"):
        content.append(Image(path=res["dotplot_path"]))
    return content


if __name__ == "__main__":
    if "--http" in sys.argv:
        mcp.run(transport="streamable-http")
    else:
        mcp.run()
