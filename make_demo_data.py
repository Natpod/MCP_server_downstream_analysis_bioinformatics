"""Genera datos de demo: counts sintéticos (tumor vs control) y un GMT reducido.

La señal plantada: genes de hipoxia y de ciclo celular suben en tumor,
genes de respuesta a interferón bajan. El GMT es un subconjunto ilustrativo
(no la versión oficial de MSigDB); para análisis reales descarga el GMT completo.
"""
from pathlib import Path

import numpy as np
import pandas as pd

rng = np.random.default_rng(7)
out = Path(__file__).resolve().parent / "data"
out.mkdir(exist_ok=True)

SETS = {
    "DEMO_HYPOXIA": "VEGFA SLC2A1 PGK1 LDHA CA9 BNIP3 ENO1 ALDOA HK2 PDK1 EGLN3 ADM NDRG1 P4HA1 ANKRD37 ERO1A KDM3A MXI1 PFKFB3 LOX",
    "DEMO_E2F_TARGETS": "MKI67 CCNE1 CDK1 MCM2 MCM3 MCM4 MCM5 MCM6 PCNA TOP2A E2F1 RRM2 TYMS CDC6 CDC20 BUB1 AURKA AURKB PLK1 CCNB1",
    "DEMO_INTERFERON_ALPHA": "ISG15 IFIT1 IFIT3 MX1 OAS1 OAS2 IFI6 IFI44 IFI44L RSAD2 STAT1 IRF7 IFITM1 BST2 XAF1 HERC5 USP18 CMPK2 EPSTI1 LY6E",
    "DEMO_OXIDATIVE_PHOSPHORYLATION": "NDUFA1 NDUFA2 NDUFB3 NDUFS1 SDHA SDHB UQCRC1 UQCRC2 COX4I1 COX5A ATP5F1A ATP5F1B CYCS CS IDH3A MDH2 SUCLG1 ACO2 FH OGDH",
    "DEMO_TNFA_SIGNALING_VIA_NFKB": "NFKBIA TNFAIP3 CXCL8 IL6 IL1B CCL2 JUNB FOS EGR1 IER3 SOD2 BIRC3 RELB NFKB1 PTGS2 ICAM1 CD83 ZFP36 KLF6 DUSP1",
    "DEMO_P53_PATHWAY": "CDKN1A MDM2 BAX GADD45A TP53I3 FAS SESN1 RRM2B DDB2 XPC ZMAT3 TNFRSF10B PLK3 TRIAP1 RPS27L FDXR APAF1 SPATA18 POLH TIGAR",
}
with open(out / "demo_hallmark_subset.gmt", "w") as fh:
    for name, genes in SETS.items():
        fh.write(f"{name}\tdemo subset, no oficial\t" + "\t".join(genes.split()) + "\n")

pathway_genes = sorted({g for s in SETS.values() for g in s.split()})
filler = [f"GENE{i:04d}" for i in range(1, 1881)]
genes = pathway_genes + filler

samples = [f"ctrl_{i}" for i in range(1, 7)] + [f"tumor_{i}" for i in range(1, 7)]
is_tumor = np.array([s.startswith("tumor") for s in samples])

base = rng.lognormal(mean=5, sigma=1.5, size=len(genes))
lfc = np.zeros(len(genes))
idx = {g: i for i, g in enumerate(genes)}
for g in SETS["DEMO_HYPOXIA"].split():
    lfc[idx[g]] = rng.uniform(1.5, 3)
for g in SETS["DEMO_E2F_TARGETS"].split():
    lfc[idx[g]] = rng.uniform(1.2, 2.5)
for g in SETS["DEMO_INTERFERON_ALPHA"].split():
    lfc[idx[g]] = -rng.uniform(1.2, 2.5)

mu = base[:, None] * np.where(is_tumor, 2.0 ** lfc[:, None], 1.0)
size_factors = rng.uniform(0.7, 1.4, size=len(samples))
counts = rng.negative_binomial(n=10, p=10 / (10 + mu * size_factors))
pd.DataFrame(counts, index=pd.Index(genes, name="gene"), columns=samples).to_csv(out / "demo_counts.csv")

meta = pd.DataFrame(
    {"condition": np.where(is_tumor, "tumor", "control"), "batch": ["A", "B"] * 6, "sex": rng.choice(["F", "M"], 12)},
    index=pd.Index(samples, name="sample_id"),
)
meta.to_csv(out / "demo_sample_metadata.csv")
print("Datos de demo en", out)
