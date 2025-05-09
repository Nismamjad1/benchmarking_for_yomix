import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# Load the CSVs
scanpy_df = pd.read_csv("/home/nisma/new_yomix/yomix/project/output/TCGA/benchmark_mcc_scores_TCGA_scanpy_wilcoxon_one-vs-rest.csv")
cosg_df = pd.read_csv("/home/nisma/new_yomix/yomix/project/output/TCGA/benchmark_mcc_scores_TCGA_cosg.csv")
yomix_df = pd.read_csv("/home/nisma/new_yomix/yomix/project/output/TCGA/yomix _TCGA- Sheet1.csv")
scan_py_2 = pd.read_csv("/home/nisma/new_yomix/yomix/project/output/TCGA/benchmark_mcc_scores_TCGA_scanpy_t-test_one-vs-rest.csv")

# Add method labels
scanpy_df["Method"] = "Scanpy_wilcoxon"
cosg_df["Method"] = "COSG"
yomix_df["Method"] = "Yomix"
scan_py_2["Method"] = "Scanpy_t_test"

# Normalize benchmark names
for df in [scanpy_df, cosg_df, yomix_df, scan_py_2]:
    df["Benchmark"] = df["Benchmark"].str.strip().str.lower()

# Keep only necessary columns: Benchmark, Method, and 20_features_mcc
dfs = []
for df in [scanpy_df, cosg_df, yomix_df, scan_py_2]:
    temp_df = df[["Benchmark", "Method", "1_features_mcc"]].copy()
    temp_df.rename(columns={"1_features_mcc": "MCC_Score"}, inplace=True)
    dfs.append(temp_df)

# Combine
mcc_df = pd.concat(dfs)

# Pivot for heatmap
heatmap_df = mcc_df.pivot(index="Benchmark", columns="Method", values="MCC_Score")
heatmap_df = heatmap_df.fillna(0)  # optional: fill NaNs with 0

# Plot
plt.figure(figsize=(14, 24))
sns.heatmap(heatmap_df, annot=True, fmt=".2f", cmap="coolwarm", linewidths=.5, linecolor='black', annot_kws={"size": 5})
plt.title("MCC Score Heatmap (1 Features Only)")
plt.xlabel("Method")
plt.ylabel("Benchmark")
plt.xticks(rotation=45, fontsize=10)
plt.yticks(fontsize=8)
plt.tight_layout()
plt.show()
