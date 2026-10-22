# gene-expression MCP (prototipo)

Servidor MCP en Python con tres herramientas sobre una matriz de expresión génica (genes x muestras).

| Herramienta | Qué hace |
|---|---|
| `identify_metadata` | Detecta orientación de la matriz (y la transpone si hace falta), tipo de ID de gen (símbolo HGNC/MGI, Ensembl, Entrez, RefSeq, sonda Affy), especie probable, escala de los valores (counts, log, CPM/TPM, z-score), calidad (NA, genes a cero, duplicados) y grupos experimentales a partir de los nombres de muestra o de un fichero de metadatos. |
| `visualize_expression` | Heatmap con clustering jerárquico (z-score por gen, barra de grupos) de una lista de genes o de los N más variables, y boxplots por grupo con Mann-Whitney (2 grupos) o Kruskal-Wallis (>2). Devuelve las imágenes al cliente y las guarda en disco. |
| `enrichment_analysis` | ORA con test hipergeométrico y corrección BH contra cualquier fichero GMT local. Acepta una lista de genes, o calcula genes diferenciales (Welch t-test sobre log2) entre `case_group` y `control_group`. Usa como fondo los genes medidos en la matriz. |

## Instalación

```bash
pip install -r requirements.txt
python make_demo_data.py      # crea data/demo_counts.csv, metadatos y GMT de demo
python test_client.py         # prueba las 3 herramientas a través del protocolo MCP
```

Compatible con el SDK `mcp` 1.x (`FastMCP`) y 2.x (`MCPServer`).

## Conectar a Claude Desktop / Claude Code

`claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "gene-expression": {
      "command": "python",
      "args": ["/ruta/absoluta/gene_mcp/server.py"]
    }
  }
}
```

Claude Code: `claude mcp add gene-expression -- python /ruta/absoluta/gene_mcp/server.py`

Modo HTTP: `python server.py --http` (streamable HTTP en el puerto 8000).

## Formato de entrada

- Expresión: CSV, TSV o XLSX con la primera columna como identificador. Genes en filas o en columnas (se detecta).
- Metadatos opcionales: primera columna con el ID de muestra, resto de columnas con covariables. Pasa `group_column` para elegir la variable de agrupación.
- Si no hay metadatos, los grupos salen de los nombres de muestra quitando el sufijo de réplica (`tumor_3` → `tumor`, `KO-rep2` → `KO`).

## Ejemplo de prompts

- "Inspecciona /datos/counts.tsv y dime qué tipo de datos son y qué grupos hay"
- "Haz un heatmap de los 60 genes más variables y boxplots de VEGFA, CA9 y MKI67"
- "Haz enrichment de los genes diferenciales tumor vs ctrl usando /refs/h.all.v2024.1.Hs.symbols.gmt"

## Limitaciones del prototipo

- El GMT incluido es un subconjunto ilustrativo, no la versión oficial. Para análisis reales descarga los GMT de MSigDB (Hallmark, GO, Reactome, KEGG) y pásalos en `gmt_path`.
- La expresión diferencial es un Welch t-test sobre log2, útil para explorar. Para counts de RNA-seq con rigor, sustitúyelo por DESeq2 (vía `pydeseq2`) o edgeR/limma-voom.
- Solo ORA. GSEA preranked (por ejemplo con `gseapy.prerank`) sería la siguiente herramienta natural.
- No convierte IDs: si la matriz usa Ensembl o Entrez, convierte a símbolos antes del enrichment (la herramienta 1 lo avisa).
- No corrige efectos de lote; `identify_metadata` sí lista columnas candidatas como `batch`.
- Todo funciona sin red, apto para entornos aislados.
