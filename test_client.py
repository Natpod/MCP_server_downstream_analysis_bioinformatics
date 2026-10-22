"""Prueba de extremo a extremo: lanza server.py por stdio y llama a las 3 herramientas."""
import asyncio
import json
import sys
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

BASE = Path(__file__).resolve().parent
DATA = BASE / "data"


def text_of(result):
    return next(c.text for c in result.content if c.type == "text")


async def main():
    params = StdioServerParameters(command=sys.executable, args=[str(BASE / "server.py")])
    async with stdio_client(params) as (r, w), ClientSession(r, w) as s:
        await s.initialize()
        tools = await s.list_tools()
        print("Herramientas:", [t.name for t in tools.tools])

        r1 = await s.call_tool("identify_metadata", {"expression_path": str(DATA / "demo_counts.csv")})
        meta = json.loads(text_of(r1))
        print("\n[1] identify_metadata")
        print({k: meta[k] for k in ("orientation_detected", "n_genes", "n_samples", "gene_id_type", "species_guess")})
        print("valores:", meta["values"]["type"], "| grupos:", meta["inferred_groups"]["method"])
        print("avisos:", meta["warnings"])

        r2 = await s.call_tool("visualize_expression", {
            "expression_path": str(DATA / "demo_counts.csv"),
            "top_n_variable": 40,
            "boxplot_genes": ["VEGFA", "CA9", "MKI67", "ISG15", "IFIT1", "GENE0001"],
        })
        viz = json.loads(text_of(r2))
        print("\n[2] visualize_expression ->", viz["heatmap_path"], viz["boxplot_path"])
        print("imágenes devueltas:", sum(c.type == "image" for c in r2.content))
        print("p-valores boxplot:", {k: round(v, 4) for k, v in viz["boxplot_pvalues"].items()})

        r3 = await s.call_tool("enrichment_analysis", {
            "expression_path": str(DATA / "demo_counts.csv"),
            "sample_metadata_path": str(DATA / "demo_sample_metadata.csv"),
            "group_column": "condition",
            "case_group": "tumor",
            "control_group": "control",
        })
        enr = json.loads(text_of(r3))
        print("\n[3] enrichment_analysis (modo DE)")
        print("DE:", {k: enr["de_summary"][k] for k in ("comparison", "n_significant", "n_up", "n_down")})
        for t in enr["top_terms"]:
            print(f"  {t['term']:35s} {t['overlap']:>6s}  FE={t['fold_enrichment']:<6} padj={t['padj']:.2e}")

        r4 = await s.call_tool("enrichment_analysis", {
            "gene_list": ["ISG15", "IFIT1", "MX1", "OAS1", "STAT1", "IRF7", "RSAD2", "VEGFA"],
        })
        enr2 = json.loads(text_of(r4))
        print("\n[3b] enrichment_analysis (lista directa):", enr2["top_terms"][0]["term"], enr2["top_terms"][0]["overlap"])


asyncio.run(main())
