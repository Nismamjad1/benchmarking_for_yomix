import pandas as pd
import matplotlib.pyplot as plt

# Load the CSVs
scanpy_df = pd.read_csv("/home/nisma/benchmark/benchmarking_for_yomix/output/TCGA/benchmark_mcc_scores_TCGA_scanpy_wilcoxon_one-vs-rest.csv")
cosg_df = pd.read_csv("/home/nisma/benchmark/benchmarking_for_yomix/output/TCGA/benchmark_mcc_scores_TCGA_cosg.csv")
yomix_df = pd.read_csv("/home/nisma/benchmark/benchmarking_for_yomix/output/TCGA/yomix-Sheet1.csv")
scan_py_2 = pd.read_csv("/home/nisma/benchmark/benchmarking_for_yomix/output/TCGA/benchmark_mcc_scores_TCGA_scanpy_t-test_one-vs-rest.csv")

# Extract and standardize the time columns
df_ttest = scan_py_2[["Benchmark", "DE_Time_Taken"]].copy()
df_ttest["Method"] = "Scanpy_t-test"

df_wilcoxon = scanpy_df[["Benchmark", "DE_Time_Taken"]].copy()
df_wilcoxon["Method"] = "Scanpy_wilcoxon"

df_yomix = yomix_df[["Benchmark", "DE_Time_Taken"]].copy()
df_yomix["Method"] = "Yomix"

df_cosg = cosg_df[["Benchmark", "Time_Taken"]].rename(columns={"Time_Taken": "DE_Time_Taken"}).copy()
df_cosg["Method"] = "COSG"

# Combine all timing data into one DataFrame
all_times = pd.concat([df_ttest, df_wilcoxon, df_yomix, df_cosg])

# Calculate the average and standard deviation per method
summary = all_times.groupby("Method")["DE_Time_Taken"].agg(["mean", "std"]).reset_index()

# Plotting
fig, ax = plt.subplots(figsize=(8, 5))
colors = ['green', 'yellow', 'orange', 'red']
bars = ax.bar(summary["Method"], summary["mean"], yerr=summary["std"], capsize=5, color=colors, edgecolor='black')

# Annotate bars with average values
for bar, avg in zip(bars, summary["mean"]):
    height = bar.get_height()
    ax.annotate(f'{avg:.2f}s',
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3), textcoords="offset points",
                ha='center', va='bottom', fontsize=9)

# Final formatting
ax.set_ylabel("Average Time (seconds)")
ax.set_title("Average Gene Ranking Time per Method (TCGA)")
plt.xticks(rotation=15)
plt.tight_layout()
plt.show()
