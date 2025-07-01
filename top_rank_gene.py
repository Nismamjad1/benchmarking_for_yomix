import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Load the CSV files
cosg_df = pd.read_csv("/home/nisma/top20_cosg_all_labels.csv")
scanpy_df = pd.read_csv("/home/nisma/top20_scanpy_wilcoxon_all_labels.csv")
scanpy_df_2 = pd.read_csv("/home/nisma/top20_scanpy_t-test_all_labels.csv")

# Extract N_Breast genes
cosg_genes = set(cosg_df[cosg_df.iloc[:, 0] == "N_Breast"].iloc[:, 1])
scanpy_genes_wil = set(scanpy_df[scanpy_df.iloc[:, 0] == "N_Breast"].iloc[:, 1])
scanpy_genes_ttest = set(scanpy_df_2[scanpy_df_2.iloc[:, 0] == "N_Breast"].iloc[:, 1])

# Yomix N_Breast gene signature
yomix_genes = {
    "ENSG00000180914.10", "ENSG00000235687.8", "ENSG00000240800.1",
    "ENSG00000234918.1", "ENSG00000161649.12", "ENSG00000148513.17",
    "ENSG00000170323.8", "ENSG00000187288.10", "ENSG00000166819.11",
    "ENSG00000181092.9", "ENSG00000135218.17", "ENSG00000226482.1",
    "ENSG00000184811.3", "ENSG00000079435.9", "ENSG00000259603.1",
    "ENSG00000174697.4", "ENSG00000173432.10", "ENSG00000135917.13",
    "ENSG00000259916.1", "ENSG00000228639.6"
}

# Intersections from Yomix perspective
print("\n Yomix ∩ COSG:")
print(sorted(yomix_genes & cosg_genes))

print("\n Yomix ∩ Scanpy Wilcoxon:")
print(sorted(yomix_genes & scanpy_genes_wil))

print("\n Yomix ∩ Scanpy t-test:")
print(sorted(yomix_genes & scanpy_genes_ttest))

print("\n Yomix ∩ Scanpy Wilcoxon ∩ Scanpy t-test:")
print(sorted(yomix_genes & scanpy_genes_wil & scanpy_genes_ttest))

print("\n Yomix ∩ All (COSG ∩ Scanpy Wilcoxon ∩ t-test):")
print(sorted(yomix_genes & cosg_genes & scanpy_genes_wil & scanpy_genes_ttest))

print("\n Yomix Unique Genes (not in any other method):")
unique_yomix = yomix_genes - (cosg_genes | scanpy_genes_wil | scanpy_genes_ttest)
print(sorted(unique_yomix))

# 1. Venn: Yomix vs COSG
plt.figure(figsize=(5, 5))
venn2([yomix_genes, cosg_genes], set_labels=("Yomix", "COSG"))
plt.title("Yomix vs COSG")
plt.tight_layout()
plt.show()

# 2. Venn: Yomix vs Scanpy Wilcoxon
plt.figure(figsize=(5, 5))
venn2([yomix_genes, scanpy_genes_wil], set_labels=("Yomix", "Scanpy Wilcoxon"))
plt.title("Yomix vs Scanpy Wilcoxon")
plt.tight_layout()
plt.show()

# 3. Venn: Yomix vs Scanpy t-test
plt.figure(figsize=(5, 5))
venn2([yomix_genes, scanpy_genes_ttest], set_labels=("Yomix", "Scanpy t-test"))
plt.title("Yomix vs Scanpy t-test")
plt.tight_layout()
plt.show()