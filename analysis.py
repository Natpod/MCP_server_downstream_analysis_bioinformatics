"""Lógica de análisis independiente del protocolo MCP.

Entrada esperada: matriz de expresión (genes x muestras o muestras x genes)
en CSV/TSV/XLSX. Opcionalmente, un fichero de metadatos de muestras.
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.cluster import hierarchy

# ---------------------------------------------------------------------------
# Carga y detección
# ---------------------------------------------------------------------------

GENE_ID_PATTERNS = {
    "ensembl_gene": re.compile(r"^ENS[A-Z]*G\d{11}(\.\d+)?$"),
    "ensembl_transcript": re.compile(r"^ENS[A-Z]*T\d{11}(\.\d+)?$"),
    "entrez": re.compile(r"^\d+$"),
    "refseq": re.compile(r"^[NX][MR]_\d+(\.\d+)?$"),
    "affy_probe": re.compile(r"^\d+(_[a-z]+)?_at$"),
    "hgnc_symbol": re.compile(r"^[A-Z][A-Z0-9\-\.]{0,14}$"),
    "mgi_symbol": re.compile(r"^[A-Z][a-z0-9\-\.]{1,14}$"),
}

GROUP_KEYWORDS = {
    "control": ["ctrl", "control", "ctl", "normal", "healthy", "wt", "vehicle", "untreated", "sano"],
    "case": ["tumor", "tumour", "case", "disease", "ko", "treated", "trt", "mut", "patient", "enfermo"],
}


def _read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"No existe el fichero: {path}")
    suffix = "".join(path.suffixes).lower()
    if suffix.endswith((".xlsx", ".xls")):
        return pd.read_excel(path, index_col=0)
    sep = "\t" if suffix.endswith((".tsv", ".txt", ".tsv.gz", ".txt.gz")) else ","
    return pd.read_csv(path, sep=sep, index_col=0)


def _classify_ids(ids: list[str]) -> dict[str, float]:
    sample = [str(i) for i in ids[:2000]]
    out = {}
    for name, pat in GENE_ID_PATTERNS.items():
        out[name] = float(np.mean([bool(pat.match(s)) for s in sample])) if sample else 0.0
    return out


def _best_id_type(ids: list[str]) -> tuple[str, float]:
    scores = _classify_ids(ids)
    # Prioridad: los patrones específicos ganan a los genéricos de símbolo
    order = ["ensembl_gene", "ensembl_transcript", "refseq", "affy_probe", "entrez", "hgnc_symbol", "mgi_symbol"]
    for name in order:
        if scores[name] >= 0.6:
            return name, scores[name]
    best = max(scores, key=scores.get)
    return (best if scores[best] > 0.2 else "unknown"), scores[best]


def load_expression(path: str | Path) -> tuple[pd.DataFrame, str]:
    """Devuelve la matriz orientada como genes x muestras y la orientación detectada."""
    df = _read_table(path)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    idx_type, idx_score = _best_id_type(list(df.index))
    col_type, col_score = _best_id_type(list(df.columns))
    looks_like_genes_in_cols = (
        col_type not in ("unknown", "hgnc_symbol", "mgi_symbol") and col_score > idx_score
    ) or (idx_type == "unknown" and col_type != "unknown" and df.shape[1] > df.shape[0])

    if looks_like_genes_in_cols:
        df = df.T
        orientation = "samples_x_genes (transpuesta)"
    else:
        orientation = "genes_x_samples"

    df = df.apply(pd.to_numeric, errors="coerce")
    return df, orientation


def _detect_value_type(df: pd.DataFrame) -> dict:
    vals = df.to_numpy(dtype=float)
    finite = vals[np.isfinite(vals)]
    if finite.size == 0:
        return {"type": "unknown", "reason": "sin valores numéricos"}
    is_int = bool(np.allclose(finite, np.round(finite)))
    vmin, vmax, med = float(finite.min()), float(finite.max()), float(np.median(finite))
    col_sums = df.sum(axis=0)
    cv_sums = float(col_sums.std() / col_sums.mean()) if col_sums.mean() else np.nan

    if vmin < 0:
        t, reason = "log_centered_or_zscore", "hay valores negativos"
    elif is_int and vmax > 100:
        t, reason = "raw_counts", "enteros no negativos con rango amplio"
    elif vmax <= 25 and med < 15:
        t, reason = "log_transformed", "rango comprimido (máx <= 25)"
    elif not is_int and np.allclose(col_sums, 1e6, rtol=0.05):
        t, reason = "CPM/TPM", "las columnas suman ~1e6"
    else:
        t, reason = "normalized_linear", "valores continuos no logarítmicos"
    return {
        "type": t,
        "reason": reason,
        "min": round(vmin, 3),
        "max": round(vmax, 3),
        "median": round(med, 3),
        "integer_valued": is_int,
        "library_size_cv": None if np.isnan(cv_sums) else round(cv_sums, 3),
        "needs_log": t in ("raw_counts", "CPM/TPM", "normalized_linear"),
    }


def _infer_groups_from_names(samples: list[str]) -> dict:
    """Intenta deducir grupos experimentales a partir de los nombres de muestra."""
    labels = {}
    # Grupo = nombre de muestra sin el sufijo de réplica (p.ej. "tumor_3" -> "tumor", "KO-rep2" -> "KO")
    prefixes = {s: re.sub(r"[_\-\.\s]?(rep|r|s)?\d+$", "", s, flags=re.I) for s in samples}
    group_names = sorted(set(prefixes.values()))
    if not (1 < len(group_names) < len(samples)):
        return {"method": "none", "groups": {s: None for s in samples}, "group_sizes": {}, "roles": {}}

    # Rol sugerido de cada grupo (control / case) por palabras clave
    roles = {}
    for g in group_names:
        low = g.lower()
        roles[g] = next(
            (role for role, kws in GROUP_KEYWORDS.items() if any(re.search(rf"(^|[_\-\.\s]){kw}", low) for kw in kws)),
            "unknown",
        )
    return {
        "method": "name_prefix",
        "groups": prefixes,
        "group_sizes": pd.Series(prefixes).value_counts().to_dict(),
        "roles": roles,
    }


def identify_metadata(expression_path: str, sample_metadata_path: str | None = None) -> dict:
    df, orientation = load_expression(expression_path)
    id_type, id_conf = _best_id_type(list(df.index))
    value_info = _detect_value_type(df)

    n_missing = int(df.isna().sum().sum())
    zero_genes = int((df.fillna(0) == 0).all(axis=1).sum())
    dup_genes = int(df.index.duplicated().sum())

    result = {
        "file": str(expression_path),
        "orientation_detected": orientation,
        "n_genes": int(df.shape[0]),
        "n_samples": int(df.shape[1]),
        "samples": list(df.columns),
        "gene_id_type": id_type,
        "gene_id_confidence": round(id_conf, 2),
        "gene_id_examples": list(df.index[:5]),
        "species_guess": _guess_species(list(df.index), id_type),
        "values": value_info,
        "quality": {
            "missing_values": n_missing,
            "all_zero_genes": zero_genes,
            "duplicated_gene_ids": dup_genes,
        },
    }

    if sample_metadata_path:
        meta = _read_table(sample_metadata_path)
        meta.index = meta.index.astype(str)
        overlap = sorted(set(meta.index) & set(df.columns))
        candidate_cols = [
            c for c in meta.columns if 1 < meta[c].nunique() < len(meta) and meta[c].dtype == object
        ]
        result["sample_metadata"] = {
            "file": str(sample_metadata_path),
            "columns": list(meta.columns),
            "samples_matched": len(overlap),
            "samples_unmatched": sorted(set(df.columns) - set(meta.index)),
            "candidate_group_columns": {c: meta[c].value_counts().to_dict() for c in candidate_cols},
        }
    else:
        result["inferred_groups"] = _infer_groups_from_names(list(df.columns))

    result["warnings"] = _warnings(result)
    return result


def _guess_species(ids: list[str], id_type: str) -> str:
    if id_type.startswith("ensembl"):
        prefixes = pd.Series([re.match(r"^(ENS[A-Z]*)[GT]", i).group(1) for i in ids[:500] if re.match(r"^ENS[A-Z]*[GT]", i)])
        top = prefixes.mode().iloc[0] if not prefixes.empty else ""
        return {"ENS": "Homo sapiens", "ENSMUS": "Mus musculus", "ENSRNO": "Rattus norvegicus", "ENSDAR": "Danio rerio"}.get(top, f"desconocida ({top})")
    if id_type == "hgnc_symbol":
        return "probablemente Homo sapiens (símbolos en mayúsculas)"
    if id_type == "mgi_symbol":
        return "probablemente Mus musculus (símbolos tipo Capitalizado)"
    return "no determinada"


def _warnings(r: dict) -> list[str]:
    w = []
    if r["gene_id_type"] == "unknown":
        w.append("No se reconoce el tipo de identificador de gen; el enrichment puede fallar.")
    if r["gene_id_type"] not in ("hgnc_symbol", "mgi_symbol"):
        w.append("Los gene sets por defecto usan símbolos; convierte los IDs antes del enrichment.")
    if r["values"].get("needs_log"):
        w.append("Valores en escala lineal: las visualizaciones aplicarán log2(x+1).")
    if r["quality"]["duplicated_gene_ids"]:
        w.append("Hay IDs de gen duplicados; se promediarán al visualizar.")
    if r["n_samples"] < 3:
        w.append("Menos de 3 muestras: los tests estadísticos no serán fiables.")
    return w


# ---------------------------------------------------------------------------
# Preparación común
# ---------------------------------------------------------------------------


def _prepare(df: pd.DataFrame) -> pd.DataFrame:
    df = df.groupby(level=0).mean()
    df = df.dropna(how="all").fillna(0)
    if _detect_value_type(df).get("needs_log"):
        df = np.log2(df.clip(lower=0) + 1)
    return df


def resolve_groups(df: pd.DataFrame, groups: dict[str, str] | None, sample_metadata_path: str | None, group_column: str | None) -> pd.Series:
    if groups:
        s = pd.Series(groups)
    elif sample_metadata_path and group_column:
        meta = _read_table(sample_metadata_path)
        meta.index = meta.index.astype(str)
        s = meta[group_column].astype(str)
    else:
        s = pd.Series(_infer_groups_from_names(list(df.columns))["groups"])
    s = s.reindex(df.columns)
    return s.fillna("sin_grupo")


# ---------------------------------------------------------------------------
# Visualización
# ---------------------------------------------------------------------------


def plot_expression(
    expression_path: str,
    output_dir: str,
    genes: list[str] | None = None,
    top_n_variable: int = 50,
    groups: dict[str, str] | None = None,
    sample_metadata_path: str | None = None,
    group_column: str | None = None,
    boxplot_genes: list[str] | None = None,
) -> dict:
    raw, _ = load_expression(expression_path)
    df = _prepare(raw)
    grp = resolve_groups(df, groups, sample_metadata_path, group_column)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    missing = []
    if genes:
        missing = [g for g in genes if g not in df.index]
        sel = [g for g in genes if g in df.index]
        if not sel:
            raise ValueError("Ninguno de los genes indicados está en la matriz.")
    else:
        sel = df.var(axis=1).sort_values(ascending=False).head(top_n_variable).index.tolist()

    sub = df.loc[sel]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, 1), axis=0)

    # Heatmap con clustering jerárquico y barra de grupos
    palette = dict(zip(sorted(grp.unique()), sns.color_palette("Set2", grp.nunique())))
    col_colors = grp.map(palette)
    g = sns.clustermap(
        z,
        cmap="RdBu_r",
        center=0,
        vmin=-3,
        vmax=3,
        col_colors=col_colors,
        yticklabels=len(sel) <= 80,
        figsize=(max(6, 0.35 * z.shape[1] + 3), max(5, 0.18 * z.shape[0] + 3)),
        cbar_kws={"label": "z-score"},
        method="average",
        metric="correlation" if z.shape[0] > 2 else "euclidean",
    )
    for lbl, color in palette.items():
        g.ax_col_dendrogram.bar(0, 0, color=color, label=lbl, linewidth=0)
    g.ax_col_dendrogram.legend(loc="center", ncol=min(4, len(palette)), frameon=False)
    heatmap_path = out / "heatmap.png"
    g.savefig(heatmap_path, dpi=150, bbox_inches="tight")
    plt.close(g.fig)

    # Boxplot: genes indicados o los 6 primeros de la selección
    bx_genes = [b for b in (boxplot_genes or sel[:6]) if b in df.index]
    long = (
        df.loc[bx_genes]
        .T.assign(group=grp.values)
        .melt(id_vars="group", var_name="gene", value_name="expression")
    )
    ncol = min(3, len(bx_genes))
    nrow = int(np.ceil(len(bx_genes) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.4 * nrow), squeeze=False)
    pvals = {}
    order = sorted(grp.unique())
    for ax, gene in zip(axes.flat, bx_genes):
        d = long[long.gene == gene]
        sns.boxplot(data=d, x="group", y="expression", hue="group", order=order, palette=palette, ax=ax, showfliers=False, legend=False)
        sns.stripplot(data=d, x="group", y="expression", order=order, color="black", size=3, alpha=0.6, ax=ax)
        samples_by_group = [d.loc[d.group == k, "expression"].values for k in order]
        if len(order) == 2 and all(len(x) >= 2 for x in samples_by_group):
            p = stats.mannwhitneyu(*samples_by_group).pvalue
        elif len(order) > 2 and all(len(x) >= 2 for x in samples_by_group):
            p = stats.kruskal(*samples_by_group).pvalue
        else:
            p = np.nan
        pvals[gene] = None if np.isnan(p) else float(p)
        ax.set_title(f"{gene}" + ("" if np.isnan(p) else f"  (p={p:.2g})"), fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("log2 expresión")
    for ax in list(axes.flat)[len(bx_genes):]:
        ax.axis("off")
    fig.tight_layout()
    boxplot_path = out / "boxplot.png"
    fig.savefig(boxplot_path, dpi=150)
    plt.close(fig)

    return {
        "heatmap_path": str(heatmap_path),
        "boxplot_path": str(boxplot_path),
        "genes_in_heatmap": sel,
        "genes_not_found": missing,
        "groups_used": grp.value_counts().to_dict(),
        "boxplot_pvalues": pvals,
        "boxplot_test": "Mann-Whitney U" if len(order) == 2 else "Kruskal-Wallis",
    }


# ---------------------------------------------------------------------------
# Enrichment (ORA hipergeométrico, sin dependencia de red)
# ---------------------------------------------------------------------------


def read_gmt(path: str | Path) -> dict[str, set[str]]:
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sets[parts[0]] = {g.upper() for g in parts[2:] if g}
    return sets


def differential_genes(df: pd.DataFrame, grp: pd.Series, case: str, control: str, padj_cutoff: float, lfc_cutoff: float) -> pd.DataFrame:
    a = df.loc[:, grp == case]
    b = df.loc[:, grp == control]
    if a.shape[1] < 2 or b.shape[1] < 2:
        raise ValueError(
            f"Se necesitan >=2 muestras por grupo ({case}: {a.shape[1]}, {control}: {b.shape[1]}). "
            f"Grupos disponibles: {grp.value_counts().to_dict()}"
        )
    lfc = a.mean(axis=1) - b.mean(axis=1)
    _, p = stats.ttest_ind(a, b, axis=1, equal_var=False)
    res = pd.DataFrame({"log2FC": lfc, "pvalue": p}, index=df.index).dropna()
    res["padj"] = _bh(res["pvalue"].values)
    res["direction"] = np.where(res["log2FC"] > 0, "up", "down")
    sig = res[(res.padj < padj_cutoff) & (res.log2FC.abs() >= lfc_cutoff)]
    return sig.sort_values("padj")


def _bh(p: np.ndarray) -> np.ndarray:
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def enrichment(
    gmt_path: str,
    output_dir: str,
    gene_list: list[str] | None = None,
    expression_path: str | None = None,
    groups: dict[str, str] | None = None,
    sample_metadata_path: str | None = None,
    group_column: str | None = None,
    case_group: str | None = None,
    control_group: str | None = None,
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 1.0,
    min_set_size: int = 5,
    max_set_size: int = 500,
    top_n_plot: int = 15,
) -> dict:
    gene_sets = read_gmt(gmt_path)
    background = None
    de_summary = None

    if expression_path:
        raw, _ = load_expression(expression_path)
        df = _prepare(raw)
        background = {g.upper() for g in df.index}
        if not gene_list:
            if not (case_group and control_group):
                raise ValueError("Sin gene_list, indica case_group y control_group para obtener genes DE.")
            grp = resolve_groups(df, groups, sample_metadata_path, group_column)
            de = differential_genes(df, grp, case_group, control_group, padj_cutoff, lfc_cutoff)
            gene_list = de.index.tolist()
            de_summary = {
                "comparison": f"{case_group} vs {control_group}",
                "test": "Welch t-test sobre log2, BH",
                "n_significant": len(de),
                "n_up": int((de.direction == "up").sum()),
                "n_down": int((de.direction == "down").sum()),
                "top_genes": de.head(10).round({"log2FC": 3}).reset_index(names="gene").to_dict(orient="records"),
            }
            de.to_csv(Path(output_dir) / "differential_genes.csv")

    if not gene_list:
        raise ValueError("No hay genes de entrada para el enrichment.")

    query = {g.upper() for g in gene_list}
    if background is None:
        background = set().union(*gene_sets.values())
    query &= background
    N, n = len(background), len(query)

    rows = []
    for name, members in gene_sets.items():
        members = members & background
        K = len(members)
        if K < min_set_size or K > max_set_size:
            continue
        overlap = query & members
        k = len(overlap)
        if k == 0:
            continue
        p = stats.hypergeom.sf(k - 1, N, K, n)
        rows.append({
            "term": name,
            "overlap": f"{k}/{K}",
            "k": k,
            "set_size": K,
            "fold_enrichment": (k / n) / (K / N) if n else np.nan,
            "pvalue": p,
            "genes": ";".join(sorted(overlap)),
        })

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    if not rows:
        return {"n_query_genes": n, "background_size": N, "results": [], "de_summary": de_summary, "message": "Sin solapamiento con ningún gene set."}

    res = pd.DataFrame(rows)
    res["padj"] = _bh(res["pvalue"].values)
    res = res.sort_values(["padj", "pvalue"]).reset_index(drop=True)
    table_path = out / "enrichment_results.csv"
    res.to_csv(table_path, index=False)

    top = res.head(top_n_plot).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 0.4 * len(top) + 1.5))
    sc = ax.scatter(
        top["fold_enrichment"], top["term"], s=top["k"] * 25, c=-np.log10(top["padj"]), cmap="viridis", edgecolor="black", linewidth=0.4
    )
    ax.set_xlabel("Fold enrichment")
    ax.set_title(f"ORA: {n} genes vs fondo de {N}")
    fig.colorbar(sc, ax=ax, label="-log10(padj)")
    fig.tight_layout()
    plot_path = out / "enrichment_dotplot.png"
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)

    return {
        "method": "Over-representation analysis (test hipergeométrico, BH)",
        "gene_sets_file": gmt_path,
        "n_query_genes": n,
        "background_size": N,
        "background": "genes de la matriz" if expression_path else "unión de los gene sets",
        "n_terms_tested": len(res),
        "n_significant": int((res.padj < 0.05).sum()),
        "top_terms": [{**r, "fold_enrichment": round(r["fold_enrichment"], 2), "pvalue": float(f"{r['pvalue']:.3g}"), "padj": float(f"{r['padj']:.3g}")} for r in res.head(10).to_dict(orient="records")],
        "table_path": str(table_path),
        "dotplot_path": str(plot_path),
        "de_summary": de_summary,
    }
