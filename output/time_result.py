import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load CSVs
scanpy_df = pd.read_csv("/home/nisma/new_yomix/yomix/project/output/benchmark_mcc_scores_scanpy_wilcoxon_one-vs-rest.csv")
cosg_df = pd.read_csv("/home/nisma/new_yomix/yomix/project/output/benchmark_mcc_scores_cosg.csv")

# Add method label
scanpy_df["Method"] = "Scanpy"
cosg_df["Method"] = "COSG"

# Combine time columns
time_df = pd.concat([
    scanpy_df[["Benchmark", "DE_Time_Taken", "Method"]],
    cosg_df[["Benchmark", "DE_Time_Taken", "Method"]]
])

# Pivot table: rows=Benchmark, columns=Method
pivot = time_df.pivot(index="Benchmark", columns="Method", values="DE_Time_Taken")

# Sort by Scanpy time or alphabetical
pivot = pivot.sort_index()

# Bar plot settings
x = np.arange(len(pivot.index))  # number of benchmarks
width = 0.35  # bar width

fig, ax = plt.subplots(figsize=(16, 6))
bars1 = ax.bar(x - width/2, pivot["Scanpy"], width, label="Scanpy")
bars2 = ax.bar(x + width/2, pivot["COSG"], width, label="COSG")

# One label per benchmark
ax.set_xticks(x)
ax.set_xticklabels(pivot.index, rotation=90, ha="center")

ax.set_title("Time Taken per Benchmark (Scanpy vs COSG)")
ax.set_ylabel("Time (seconds)")
ax.set_xlabel("Benchmark")
ax.legend()
plt.tight_layout()
plt.show()
