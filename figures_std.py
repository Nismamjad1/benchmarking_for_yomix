import pandas as pd
import matplotlib.pyplot as plt

# Load the CSVs
scanpy_wilcoxon = pd.read_csv("output/benchmark/benchmarking_for_yomix/output/TCGA/benchmark_mcc_scores_TCGA_scanpy_wilcoxon_one-vs-rest.csv")
cosg = pd.read_csv("output/benchmark/benchmarking_for_yomix/output/TCGA/benchmark_mcc_scores_TCGA_cosg.csv")
yomix = pd.read_csv("output/benchmark/benchmarking_for_yomix/output/TCGA/yomix - Sheet1.csv")
scanpy_ttest = pd.read_csv("output/benchmark/benchmarking_for_yomix/output/TCGA/benchmark_mcc_scores_TCGA_scanpy_t-test_one-vs-rest.csv")
# Add method labels
yomix["Method"] = "Yomix"
scanpy_wilcoxon["Method"] = "Scanpy Wilcoxon"
scanpy_ttest["Method"] = "Scanpy T-Test"
cosg["Method"] = "COSG"

merged = pd.concat([yomix, scanpy_wilcoxon, scanpy_ttest, cosg], ignore_index=True)

for col in merged.columns:
    if any(x in col for x in ["mcc", "precision", "recall", "f1"]):
        merged[col] = pd.to_numeric(merged[col], errors='coerce')

target_labels = [
    "N_Liver_vs_rest","T_OV_vs_Rest", "T_CESC_vs_Rest", "T_SARC_vs_Rest", "T_DLBC_vs_Rest", "T_CHOL_vs_Rest",
    "N_Adipose_Tissue_vs_Rest", "N_Heart_vs_Rest", "N_Vagina_vs_Rest", "N_Cervix_Uteri_vs_Rest"
]

merged["Benchmark"] = merged["Benchmark"].str.lower().str.strip()
target_labels = [x.lower() for x in target_labels]
filtered = merged[merged["Benchmark"].isin(target_labels)]


feature_sizes = [1, 3, 10, 20]

for label in target_labels:
    plt.figure(figsize=(10, 6))
    label_df = filtered[filtered["Benchmark"] == label]

    for method in label_df["Method"].unique():
        method_df = label_df[label_df["Method"] == method]

        mcc_means = []
        mcc_stds = []

        for size in feature_sizes:
            mean_col = f"{size}_features_mcc"
            std_col = f"{size}_features_mcc_std"
            mcc_means.append(method_df[mean_col].mean())
            mcc_stds.append(method_df[std_col].mean())

        plt.errorbar(feature_sizes, mcc_means, yerr=mcc_stds,
                     label=method, capsize=4, marker='o', linestyle='-')

    plt.title(f"MCC Scores for {label}")
    plt.xlabel("Number of Features")
    plt.ylabel("Mean MCC Score")
    plt.xticks(feature_sizes)
    plt.grid(True)
    plt.legend(title="Method")
    plt.tight_layout()
    plt.show()  