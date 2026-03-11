import matplotlib
matplotlib.use("Agg")
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import os

# ── Sample metadata ───────────────────────────────────────
metadata = pd.DataFrame({
    "sample": ["SRR1039508", "SRR1039509", "SRR1039512", "SRR1039513"],
    "condition": ["untreated", "treated", "untreated", "treated"]
}).set_index("sample")

print("=================================")
print("      RNA-seq DE ANALYSIS")
print("=================================")
print("Samples:")
print(metadata)
print()

# ── Load Salmon counts ────────────────────────────────────
def load_quant(sample):
    path = f"results/{sample}/quant.sf"
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df["NumReads"].round().astype(int)

print("Loading Salmon quantification files...")
counts = pd.DataFrame({s: load_quant(s) for s in metadata.index})
print(f"Loaded {len(counts):,} transcripts across {len(counts.columns)} samples")

# ── Filter lowly expressed transcripts ───────────────────
min_counts = 10
counts_filtered = counts[counts.sum(axis=1) >= min_counts]
print(f"After filtering (sum >= {min_counts}): {len(counts_filtered):,} transcripts")
print()

# ── DESeq2 analysis ───────────────────────────────────────
print("Running PyDESeq2...")
dds = DeseqDataSet(
    counts=counts_filtered.T,
    metadata=metadata,
    design_factors="condition"
)
dds.deseq2()

# ── Extract results ───────────────────────────────────────
stat_res = DeseqStats(dds, contrast=["condition", "treated", "untreated"])
stat_res.summary()
results = stat_res.results_df

print(f"\nTotal transcripts tested: {len(results):,}")
print(f"Significant (padj<0.05):  {(results['padj'] < 0.05).sum():,}")
print(f"Upregulated (FC>1):       {((results['padj'] < 0.05) & (results['log2FoldChange'] > 1)).sum():,}")
print(f"Downregulated (FC<-1):    {((results['padj'] < 0.05) & (results['log2FoldChange'] < -1)).sum():,}")

# ── Save results ──────────────────────────────────────────
os.makedirs("results/de_analysis", exist_ok=True)
results.to_csv("results/de_analysis/de_results.csv")
print("\nResults saved to results/de_analysis/de_results.csv")

# ── Plot 1: Volcano plot ──────────────────────────────────
print("Generating plots...")
plt.figure(figsize=(10, 6))

# All transcripts
plt.scatter(results["log2FoldChange"], -np.log10(results["pvalue"]),
            alpha=0.3, color="lightgray", s=10, label="Not significant")

# Significant
sig = results[(results["padj"] < 0.05) & (abs(results["log2FoldChange"]) > 1)]
up  = sig[sig["log2FoldChange"] > 1]
dn  = sig[sig["log2FoldChange"] < -1]

plt.scatter(up["log2FoldChange"], -np.log10(up["pvalue"]),
            alpha=0.7, color="firebrick", s=15, label=f"Upregulated ({len(up)})")
plt.scatter(dn["log2FoldChange"], -np.log10(dn["pvalue"]),
            alpha=0.7, color="steelblue", s=15, label=f"Downregulated ({len(dn)})")

plt.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.8, label="p=0.05")
plt.axvline(1,  color="black", linestyle="--", linewidth=0.8)
plt.axvline(-1, color="black", linestyle="--", linewidth=0.8)
plt.xlabel("log2 Fold Change (Treated / Untreated)")
plt.ylabel("-log10(p-value)")
plt.title("Volcano Plot — Dexamethasone Treatment in Airway Cells")
plt.legend()
plt.tight_layout()
plt.savefig("results/de_analysis/volcano_plot.png", dpi=150)
plt.close()
print("Saved: volcano_plot.png")

# ── Plot 2: PCA ───────────────────────────────────────────
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

log_counts = np.log1p(counts_filtered)
scaler = StandardScaler()
scaled = scaler.fit_transform(log_counts.T)

pca = PCA(n_components=2)
pca_coords = pca.fit_transform(scaled)

plt.figure(figsize=(7, 5))
colors = {"untreated": "steelblue", "treated": "firebrick"}
for sample, condition in metadata["condition"].items():
    idx = list(metadata.index).index(sample)
    plt.scatter(pca_coords[idx, 0], pca_coords[idx, 1],
                color=colors[condition], s=120, zorder=5)
    plt.annotate(sample, (pca_coords[idx, 0], pca_coords[idx, 1]),
                 textcoords="offset points", xytext=(8, 4), fontsize=9)

from matplotlib.patches import Patch
legend_elements = [Patch(facecolor="steelblue", label="Untreated"),
                   Patch(facecolor="firebrick", label="Treated")]
plt.legend(handles=legend_elements)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.title("PCA — Sample Clustering")
plt.tight_layout()
plt.savefig("results/de_analysis/pca_plot.png", dpi=150)
plt.close()
print("Saved: pca_plot.png")

# ── Plot 3: Heatmap of top DE genes ──────────────────────
top_genes = sig.reindex(sig["padj"].abs().sort_values().index).head(30)
top_counts = np.log1p(counts_filtered.loc[top_genes.index])

plt.figure(figsize=(10, 10))
sns.heatmap(top_counts, cmap="RdBu_r", center=0,
            xticklabels=metadata.index,
            yticklabels=[idx[:20] for idx in top_counts.index],
            linewidths=0.3)
plt.title("Top 30 Differentially Expressed Transcripts")
plt.xlabel("Sample")
plt.ylabel("Transcript")
plt.tight_layout()
plt.savefig("results/de_analysis/heatmap.png", dpi=150)
plt.close()
print("Saved: heatmap.png")

print("\nAll done! ✅")
print("=================================")
