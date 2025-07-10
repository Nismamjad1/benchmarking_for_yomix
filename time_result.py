import pandas as pd
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
import os
import numpy as np


def compare_time_per_label(file):

    runtime = pd.read_csv(file, index_col=0)
    ax = runtime.T.plot.bar()
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    # Add titles and labels
    ax.set_title("Time Taken per method on TCGA (Scanpy vs COSG vs Yomix)")
    ax.set_ylabel("Time (seconds)")
    ax.set_xlabel("Benchmark")
    ax.legend()

    plt.tight_layout()
    plt.show()


def compare_time_per_method():

    file_names = [
        fn for fn in os.listdir("result") if fn.split("_")[-1] == "runtime.csv"
    ]
    dfs = []
    for dataset in file_names:
        df_tmp = pd.read_csv(f"result/{dataset}", index_col=0)
        df_tmp["mean"] = df_tmp.mean(axis=1)
        df_tmp["std"] = df_tmp.std(axis=1)
        df_tmp["dataset"] = dataset.split("_")[0]

        dfs.append(df_tmp[["mean", "std", "dataset"]])
    df_runtime = pd.concat(dfs)
    df_runtime["method"] = df_runtime.index
    df_runtime = df_runtime[df_runtime["method"].isin(["cosg"])]
    datasets = df_runtime["dataset"].unique()
    methods = df_runtime["method"].unique()
    n_methods = len(methods)
    bar_width = 0.8 / n_methods
    for i, method in enumerate(methods):
        x = np.arange(len(datasets)) + (i - n_methods / 2) * bar_width + bar_width / 2
        method_df = df_runtime[df_runtime["method"] == method]
        plt.bar(x, method_df["mean"], width=bar_width, label=method)
        plt.errorbar(
            x,
            method_df["mean"],
            yerr=method_df["std"],
            fmt="none",
            capsize=5,
            color="black",
        )

    plt.xticks(ticks=np.arange(len(datasets)), labels=datasets)
    plt.xlabel("Dataset")
    plt.ylabel("Average runtime (s)")
    plt.title("Feature selection methods' average runtime across different datasets")
    plt.legend(title="Method")
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.savefig("plots/time_comparison.png")
    plt.show()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "file", type=str, nargs="?", default=None, help="the _runtime.csv file to open"
    )

    args = parser.parse_args()

    if args.file is None:
        compare_time_per_method()
    else:
        filearg = Path(args.file)
        compare_time_per_label(filearg.absolute())
