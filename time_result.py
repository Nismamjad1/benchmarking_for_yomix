import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import argparse
from pathlib import Path

def compare_time_per_method_improved():
    """
    Creates a publication-quality grouped barplot with sample/feature counts on the x-axis,
    reading data from a specified path.
    """
    
    dataset_info = {
        'citeseq': {'samples': 8617, 'features': 2000},
        'lawlor':  {'samples': 603, 'features': 2003},
        'pbmc':    {'samples': 2638, 'features': 13714},
        'meth': {'samples': 1077,  'features': 428230},
        'tcga':    {'samples': 10541,'features': 8000}
    }
    
    method_colors = {
         'yomix': '#0173B2',
        'cosg': '#DE8F05',
        'scanpy_wilcoxon': '#029E73',
        'scanpy_t-test': '#D55E00'
    }

    result_dir = "/home/nisma/bench_new/benchmarking_for_yomix/result" # Update this path as needed
    
    if not os.path.isdir(result_dir):
        print(f"Error: The specified directory does not exist: {result_dir}")
        return

    file_names = [fn for fn in os.listdir(result_dir) if fn.endswith("_runtime.csv")]
    if not file_names:
        print(f"No '_runtime.csv' files found in '{result_dir}'.")
        return

    all_dfs = []
    for file_name in file_names:
        df_tmp = pd.read_csv(os.path.join(result_dir, file_name), index_col=0)
        df_long = df_tmp.reset_index().melt(id_vars='index', var_name='run', value_name='time')
        df_long = df_long.rename(columns={'index': 'method'})
        df_long["dataset"] = file_name.replace("_runtime.csv", "")
        all_dfs.append(df_long)
        
    df_all_runtimes = pd.concat(all_dfs, ignore_index=True)

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.7)
    plt.figure(figsize=(14, 8)) 

    # --- CHANGE: REORDER THE LIST TO PLACE 'tcga' AT THE END ---
    dataset_order = ['meth', 'pbmc', 'tcga', 'citeseq', 'lawlor']

    ax = sns.barplot(
        data=df_all_runtimes,
        x="dataset",
        y="time",
        hue="method",
        palette=method_colors,
        errorbar="sd",
        order=dataset_order  # This line uses the new order
    )

    ax.set_xlabel("Dataset", fontsize=16, weight='bold')
    ax.set_ylabel("Average Runtime (s)", fontsize=16, weight='bold')
    ax.set_title("Comparison of Method Runtimes Across Datasets", fontsize=18, weight='bold')

    ax.set_yscale("log")

    positions = ax.get_xticks()
    new_labels = []
    current_labels = [label.get_text() for label in ax.get_xticklabels()]
    for dataset_name in current_labels:
        info = dataset_info.get(dataset_name, {'samples': '?', 'features': '?'})
        display_name = 'sarcomas' if dataset_name == 'meth' else dataset_name
        new_labels.append(
            f"{display_name}\nS: {info['samples']}\nF: {info['features']}"
        )

    ax.set_xticks(positions)
    ax.set_xticklabels(new_labels)
        
    plt.xticks(rotation=0, ha='center', fontsize=14)
    plt.yticks(fontsize=14)
    
    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(
        handles, labels,
        title='Method',
        loc='upper right',
        fontsize=14,
        title_fontsize=15,
        frameon=True
    )

    plt.tight_layout(rect=[0, 0, 0.9, 1])
    # I've changed the output filename to reflect the new order
    plt.savefig('runtime_comparison_tcga_last.svg', dpi=300)
    print("Final reordered runtime comparison plot saved as 'runtime_comparison_tcga_last.svg'")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "file", type=str, nargs="?", default=None, help="the _runtime.csv file to open"
    )
    args = parser.parse_args()
    
    if args.file is None:
        print(
            "WARNING : no file passed as parameter, "
            'compute the plot on all "_runtime.csv" files in the result folder'
        )
        compare_time_per_method_improved()
    else:
        print(f"INFO: file argument '{args.file}' ignored.")
        print('Computing the plot on all "_runtime.csv" files in the result folder.')
        compare_time_per_method_improved()