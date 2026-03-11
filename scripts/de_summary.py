import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Load results ──────────────────────────────────────────
df = pd.read_csv("results/de_analysis/de_results.csv", index_col=0)
df.columns = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]

# ── Summary statistics ────────────────────────────────────
print("=================================")
print("      DE ANALYSIS SUMMARY")
print("=================================")
print(f"Total transcripts tested:  {len(df):,}")
print(f"Significant (padj<0.05):   {(df['padj'] < 0.05).sum():,}")
print(f"Upregulated (LFC>1):       {((df['padj'] < 0.05) & (df['log2FoldChange'] > 1)).sum():,}")
print(f"Downregulated (LFC<-1):    {((df['padj'] < 0.05) & (df['log2FoldChange'] < -1)).sum():,}")
print()

# ── Top 10 most significant ───────────────────────────────
print("Top 10 most significant transcripts:")
top = df.dropna().sort_values("padj").head(10)
print(top[["baseMean", "log2FoldChange", "padj"]].to_string())
print()

# ── Top upregulated ───────────────────────────────────────
print("Top 5 upregulated (treated vs untreated):")
up = df.dropna()[df["log2FoldChange"] > 1].sort_values("padj").head(5)
print(up[["baseMean", "log2FoldChange", "padj"]].to_string())
print()

# ── Top downregulated ─────────────────────────────────
