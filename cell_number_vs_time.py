import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path
import anndata as ad
import scanpy as sc

def cellnum_vs_time():
    datasets = ["citeseq", "pbmc", "lawlor"]
    cell_counts = {}
    for dataset in datasets:
        adata = sc.read_h5ad(f"data/{dataset}.h5ad")
        cell_counts[dataset] = adata.n_obs

    runtime_dfs = []
    for dataset in datasets:
        df = pd.read_csv(f"result/{dataset}_runtime.csv", index_col=0)
        df["mean_time"] = df.mean(axis=1) 
        df["dataset"] = dataset
        runtime_dfs.append(df[["mean_time", "dataset"]])

    df_runtime = pd.concat(runtime_dfs)
    df_runtime["method"] = df_runtime.index
    df_runtime = df_runtime.reset_index(drop=True)
    df_runtime["cell_count"] = df_runtime["dataset"].map(cell_counts)


    plt.figure(figsize=(8, 6))
    for method in df_runtime["method"].unique():
        subset = df_runtime[df_runtime["method"] == method].sort_values("cell_count")
        plt.plot(subset["cell_count"], subset["mean_time"], marker='o', label=method)

    plt.title("Runtime vs Number of Cells per Method")
    plt.xlabel("Number of Cells")
    plt.ylabel("Average Runtime (seconds)")
    plt.legend(title="Method")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__=="__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("file", type=str, nargs="?", default=None, help="the _runtime.csv file to open")
    args = parser.parse_args()
    if args.file is None:
        cellnum_vs_time() 