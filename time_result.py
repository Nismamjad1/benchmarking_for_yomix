import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load CSVs !!!!!!!! Absolute paths, doesn't work on other people's computer
# scanpy_df = pd.read_csv("/home/nisma/new_yomix/yomix/project/output/benchmark_mcc_scores_TCGA_scanpy_wilcoxon_one-vs-rest.csv")
# cosg_df = pd.read_csv("/home/nisma/new_yomix/yomix/project/output/benchmark_mcc_scores_TCGA_cosg.csv")
# yomix_df = pd.read_csv("/home/nisma/new_yomix/yomix/project/output/yomix _TCGA- Sheet1.csv")
# scan_py_2= pd.read_csv("/home/nisma/new_yomix/yomix/project/output/benchmark_mcc_scores_TCGA_scanpy_t-test_one-vs-rest.csv")

# Relative paths
scanpy_df = pd.read_csv("output/TCGA/benchmark_mcc_scores_TCGA_scanpy_wilcoxon_one-vs-rest.csv")
cosg_df = pd.read_csv("output/TCGA/benchmark_mcc_scores_TCGA_cosg.csv")
yomix_df = pd.read_csv("output/TCGA/yomix-Sheet1.csv")
scan_py_2= pd.read_csv("output/TCGA/benchmark_mcc_scores_TCGA_scanpy_t-test_one-vs-rest.csv")


# Add method label
scanpy_df["Method"] = "Scanpy_wilcoxon"
cosg_df["Method"] = "COSG"
yomix_df["Method"]="Yomix"
scan_py_2["Method"]="Scanpy_t-test"

# Combine time columns into a single DataFrame
time_df = pd.concat([
    scanpy_df[["Benchmark", "DE_Time_Taken", "Method"]],
    # cosg_df[["Benchmark", "DE_Time_Taken", "Method"]], !!!! DE_Time_Taken column doesn't exist
    cosg_df[["Benchmark", "DE_Time_Taken", "Method"]],
    yomix_df[["Benchmark", "DE_Time_Taken", "Method"]],
    scan_py_2[["Benchmark", "DE_Time_Taken", "Method"]]
])
# Standardize Benchmark labels (Remove extra spaces and convert to lowercase for consistency)
time_df["Benchmark"] = time_df["Benchmark"].str.strip().str.lower()

# Remove duplicates and keep the first occurrence
time_df.drop_duplicates(subset=["Benchmark", "Method"], keep='first', inplace=True)


# Pivot the DataFrame for plotting
pivot = time_df.pivot_table(index="Benchmark", columns="Method", values="DE_Time_Taken", aggfunc="mean").fillna(0)
pivot = pivot.sort_index()

# Define positions and width for the bars
methods = pivot.columns
x = np.arange(len(pivot.index))
width = 0.2  # Smaller width to fit all bars without overlap

fig, ax = plt.subplots(figsize=(16, 6))

# Plot each method with a separate bar
bars = []
for i, method in enumerate(methods):
    bar = ax.bar(x + i * width - (len(methods) - 1) * width / 2, pivot[method], width, label=method)
    bars.append(bar)

# Adjust x-ticks and labels
ax.set_xticks(x)
ax.set_xticklabels(pivot.index, rotation=90, ha="center")

# Add titles and labels
ax.set_title("Time Taken per Benchmark TCGA (Scanpy vs COSG vs Yomix)")
ax.set_ylabel("Time (seconds)")
ax.set_xlabel("Benchmark")
ax.legend()

plt.tight_layout()
plt.show()
