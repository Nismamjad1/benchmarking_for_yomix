import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import shutil
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
    

    result_dir = "/home/nisma/bench_new/benchmarking_for_yomix/result"
    
    if not os.path.isdir(result_dir):
        print(f"Error: The specified directory does not exist: {result_dir}")
        return

    file_names = [fn for fn in os.listdir(result_dir) if fn.endswith("_runtime.csv")]
    if not file_names:
        print(f"No '_runtime.csv' files found in '{result_dir}'.")
        return

    dfs = []
    for file_name in file_names:
        df_tmp = pd.read_csv(os.path.join(result_dir, file_name), index_col=0)
        df_long = df_tmp.reset_index().melt(id_vars='index', var_name='run', value_name='time')
        df_long = df_long.rename(columns={'index': 'method'})
        
        df_stats = df_long.groupby('method')['time'].agg(['mean', 'std']).reset_index()
        df_stats["dataset"] = file_name.replace("_runtime.csv", "")
        dfs.append(df_stats)
        
    df_runtime = pd.concat(dfs).set_index(['dataset', 'method'])

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.7)
    plt.figure(figsize=(14, 8)) 

    ax = sns.barplot(
        data=df_runtime.reset_index(),
        x="dataset",
        y="mean",
        hue="method",
        palette="colorblind"
    )

    
    for bar in ax.patches:
        x_center = bar.get_x() + bar.get_width() / 2.
        height = bar.get_height()
        dataset = ax.get_xticklabels()[np.argmin([abs(x_center - tick.get_position()[0]) for tick in ax.get_xticklabels()])].get_text()
        try:
            std_dev = df_runtime.xs(dataset).loc[np.isclose(df_runtime.xs(dataset)['mean'], height)]['std'].values[0]
            ax.errorbar(x=x_center, y=height, yerr=std_dev, fmt='none', c='black', capsize=3)
        except (IndexError, KeyError):
            pass

    
    ax.set_xlabel("Dataset", fontsize=16, weight='bold')
    ax.set_ylabel("Average Runtime (s)", fontsize=16, weight='bold')
    ax.set_title("Comparison of Method Runtimes Across Datasets", fontsize=18, weight='bold')

    ax.set_yscale("log")

    
    positions = ax.get_xticks()

    
    new_labels = []
    current_labels = [label.get_text() for label in ax.get_xticklabels()]
    for dataset_name in current_labels:
        info = dataset_info.get(dataset_name, {'samples': '?', 'features': '?'})
        new_labels.append(
            f"{dataset_name}\nS: {info['samples']}\nF: {info['features']}"
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
        fontsize=9,
        title_fontsize=2,
        markerscale=0.8,
        handlelength=0.5,      
        handletextpad=0.00001,     
        labelspacing=0.00001,      
        borderaxespad=0.00001,     
        borderpad=0.00001,         
        frameon=True
    )


    plt.setp(legend.get_texts(), fontsize='14')
    plt.setp(legend.get_title(), fontsize='15')

    plt.tight_layout(rect=[0, 0, 0.9, 1])
    plt.savefig('runtime_comparison_final.png', dpi=300)
    print("Final runtime comparison plot saved as 'runtime_comparison_final.png'")


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