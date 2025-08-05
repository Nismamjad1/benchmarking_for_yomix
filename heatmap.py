import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path


def heatmap(file):
    metrics = ["mcc", "precision", "f1_score", "recall"]

    df = pd.read_csv(file, index_col=0)
    nbr_of_features = 1
    df = (
        df[df["nb_genes"] == nbr_of_features]
        .groupby(["method", "label_vs_rest"], as_index=False)
        .mean(numeric_only=True)
    )
    for m in metrics:

        df_tmp = df.pivot_table(index="label_vs_rest", columns="method", values=m)

        ax = sns.heatmap(
            df_tmp,
            # annot=True,
            fmt=".2f",
            cmap="coolwarm",
            linewidths=0.5,
            linecolor="black",
            # cbar_kws={'label': 'MCC'}
        )

        ax.collections[0].colorbar.set_label(m)

        plt.title(f"{m} Score Heatmap ({nbr_of_features} Features Only)")
        plt.xlabel("Method")
        plt.ylabel("Class")
        plt.xticks(rotation=45, fontsize=10)
        plt.yticks(fontsize=8)
        plt.tight_layout()
        plt.show()

def heatmap_average():

    datasets=["citeseq","pbmc","lawlor","meth", "tcga"]
    dfs = [
        pd.read_csv(f"result/{dataset}.csv", index_col=0).query("nb_genes==20").assign(dataset=dataset)
        for dataset in datasets
    ]
    
    df_all = pd.concat(dfs, ignore_index=True)
    df_avg = df_all.groupby(["dataset", "method"])[["mcc", "precision", "f1_score", "recall"]].mean().reset_index()

    metrics = ["mcc", "precision", "f1_score", "recall"]
    for metric in metrics:
        df_pivot = df_avg.pivot(index="dataset", columns="method", values=metric)
        fig, ax = plt.subplots(figsize=(8, 4))
        ax = sns.heatmap(
            df_pivot,
            annot=True,
            fmt=".2f",
            cmap="coolwarm",
            linewidths=0.5,
            linecolor="black"
        )
        ax.set_title(f"Average {metric.lower()} per method across datasets for 20 features")
        ax.set_xlabel("Method", labelpad=15)
        ax.set_ylabel("Dataset", labelpad=15)
        ax.tick_params(axis='x', rotation=45)

        fig.tight_layout()
        plt.show()



if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "file", type=str, nargs="?", default=None, help="the _runtime.csv file to open"
    )

    args = parser.parse_args()
    if args.file is None:
        heatmap_average()
    else:
        filearg = Path(args.file)
        heatmap(filearg.absolute())
