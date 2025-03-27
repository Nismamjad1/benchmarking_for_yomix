# Re-import necessary libraries and reload data after code reset
import pandas as pd
import matplotlib.pyplot as plt
import os
import zipfile
import pandas as pd
import matplotlib.pyplot as plt

# Load both Scanpy and COSG MCC files
scanpy_df = pd.read_csv("/home/nisma/new_yomix/yomix/project/output/benchmark_mcc_scores_scanpy_wilcoxon_one-vs-rest.csv")
cosg_df = pd.read_csv("/home/nisma/new_yomix/yomix/project/output/benchmark_mcc_scores_cosg.csv")

# Add method labels
scanpy_df["Method"] = "Scanpy"
cosg_df["Method"] = "COSG"

# Ensure COSG has the same columns (fill missing ones)
for col in scanpy_df.columns:
    if col not in cosg_df.columns:
        cosg_df[col] = pd.NA

# Combine both
combined_df = pd.concat([scanpy_df, cosg_df], ignore_index=True)

# Extract only MCC columns
mcc_cols = [col for col in combined_df.columns if "_features_mcc" in col]
id_vars = ["Benchmark", "Method"]

# Melt for easier plotting
mcc_df = combined_df.melt(id_vars=id_vars, value_vars=mcc_cols,
                          var_name="Feature_Size", value_name="MCC_Score")

# Extract number of features
mcc_df["Feature_Size"] = mcc_df["Feature_Size"].str.extract(r"(\d+)_features_mcc").astype(int)

# Plot per label: Scanpy & COSG together
display_limit = 10
displayed = 0

for benchmark in mcc_df["Benchmark"].unique():
    plt.figure(figsize=(8, 6))
    for method in mcc_df["Method"].unique():
        subset = mcc_df[(mcc_df["Benchmark"] == benchmark) & (mcc_df["Method"] == method)]
        subset = subset.sort_values("Feature_Size")
        plt.plot(subset["Feature_Size"], subset["MCC_Score"], marker='o', label=method)

    plt.title(f"{benchmark} - MCC Score vs Feature Size (COSG & Scanpy)")
    plt.xlabel("Number of Features")
    plt.ylabel("MCC Score")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    displayed += 1
    if displayed >= display_limit:
        break
