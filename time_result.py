import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path


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
    # TODO get all *_runtime.csv files from result folders instead of hardcoding paths
    datasets = ["tcga", "meth", "pbmc"]
    dfs = []
    for dataset in datasets:
        df_tmp = pd.read_csv(f"result/{dataset}_runtime.csv", index_col=0)
        df_tmp["mean"] = df_tmp.mean(axis=1)
        df_tmp["std"] = df_tmp.std(axis=1)
        df_tmp["dataset"] = dataset

        dfs.append(df_tmp[["mean", "std", "dataset"]])
    df_runtime = pd.concat(dfs)
    df_runtime["method"] = df_runtime.index

    sns.barplot(data=df_runtime, x="dataset", y="mean", hue="method", errorbar=None)

    a = 0
    for i, row in df_runtime.iterrows():
        # TODO change formula for errorbars to make it flexible
        plt.errorbar(
            x=a // 4 + (-0.3 + 0.2 * (a % 4)),
            y=row["mean"],
            yerr=row["std"],
            fmt="none",
            capsize=5,
            color="black",
        )
        a += 1
    plt.xlabel("Dataset")
    plt.ylabel("Average runtime (s)")
    plt.title("Feature selection methods' average runtime across different datasets")
    plt.legend(title="Method")
    plt.grid(axis="y", linestyle="--", alpha=0.7)

    plt.savefig("foo.png")
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
