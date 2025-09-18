import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path
import os

def heatmap_improved(file_path):
    """
    Generates a single figure with 4 heatmap subplots for a single dataset file.
    """
    metrics = ["mcc", "f1_score", "precision", "recall"]
    nbr_of_features = 1 

    try:
        df = pd.read_csv(file_path, index_col=0)
    except FileNotFoundError:
        print(f"Error: The file was not found at {file_path}")
        return

    
    df_filtered = (
        df[df["nb_genes"] == nbr_of_features]
        .groupby(["method", "label_vs_rest"], as_index=False)
        .mean(numeric_only=True)
    )

    
    sns.set_theme(style="white", context="paper", font_scale=1.2)
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()

    dataset_name = Path(file_path).stem

    for i, metric in enumerate(metrics):
        df_pivot = df_filtered.pivot_table(index="label_vs_rest", columns="method", values=metric)

        sns.heatmap(
            df_pivot,
            ax=axes[i],
            annot=True,
            fmt=".2f",
            cmap="vlag",
            linewidths=0.5,
            linecolor="black",
        )
        axes[i].set_title(f"{metric.replace('_', ' ').title()}", fontsize=16, weight='bold')
        axes[i].set_xlabel("Method", fontsize=12)
        axes[i].set_ylabel("Class Comparison", fontsize=12)
        axes[i].tick_params(axis='x', labelsize=10, rotation=45)
        axes[i].tick_params(axis='y', labelsize=8)

    fig.suptitle(f"Performance Metrics for {dataset_name.title()} ({nbr_of_features} Feature)", fontsize=20, weight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    
    output_filename = f"heatmap_{dataset_name}_final.png"
    plt.savefig(output_filename, dpi=300)
    print(f"Final heatmap plot saved as '{output_filename}'")


def heatmap_average_improved():
    """
    Generates a single figure with 4 heatmap subplots for an average across all datasets.
    """
    # Use the specified path for the result files
    result_dir = "/home/nisma/bench_new/benchmarking_for_yomix/result"
    
    if not os.path.isdir(result_dir):
        print(f"Error: The specified directory does not exist: {result_dir}")
        return

    datasets = ["citeseq", "pbmc", "lawlor", "meth", "tcga"]
    
    
    dfs = []
    for dataset in datasets:
        file_path = os.path.join(result_dir, f"{dataset}.csv")
        if os.path.exists(file_path):
            dfs.append(pd.read_csv(file_path, index_col=0).query("nb_genes==20").assign(dataset=dataset))
        else:
            print(f"Warning: Data file not found for dataset '{dataset}' and will be skipped.")
    
    if not dfs:
        print("No data files were found to process.")
        return

    df_all = pd.concat(dfs, ignore_index=True)
    df_avg = df_all.groupby(["dataset", "method"])[["mcc", "f1_score", "precision", "recall"]].mean().reset_index()

    metrics = ["mcc", "f1_score", "precision", "recall"]
    
    
    sns.set_theme(style="white", context="paper", font_scale=1.2)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for i, metric in enumerate(metrics):
        df_pivot = df_avg.pivot(index="method", columns="dataset", values=metric)
        
        sns.heatmap(
            df_pivot,
            ax=axes[i],
            annot=True,
            fmt=".2f",
            cmap="vlag",
            linewidths=0.5,
            linecolor="black"
        )
        axes[i].set_title(f"Average {metric.replace('_', ' ').title()}", fontsize=15, weight='bold')
        axes[i].set_xlabel("Dataset", fontsize=13, weight='bold')
        axes[i].set_ylabel("Method", fontsize=13, weight='bold')
        axes[i].tick_params(axis='x', labelsize=11, rotation=45)
        axes[i].tick_params(axis='y', labelsize=11, rotation=0)

    fig.suptitle("Average Performance Metrics Across Datasets (20 Features)", fontsize=20, weight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    
    plt.savefig('heatmap_average_comparison_final.png', dpi=300)
    print("Final heatmap average comparison plot saved as 'heatmap_average_comparison_final.png'")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate heatmap plots for performance metrics."
    )
    parser.add_argument(
        "file", type=str, nargs="?", default=None, 
        help="Path to a single CSV file to generate a detailed heatmap. If not provided, an average heatmap across all datasets will be generated."
    )

    args = parser.parse_args()

    if args.file is None:
        print("No file argument provided. Generating average heatmap across all datasets...")
        heatmap_average_improved()
    else:
        filearg = Path(args.file)
        print(f"File argument provided. Generating detailed heatmap for {filearg.name}...")
        heatmap_improved(filearg.absolute())