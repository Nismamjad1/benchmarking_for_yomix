import scanpy as sc
import numpy as np
import pandas as pd
from method.cosg_file import run_cosg
from method.scan_py import run_benchmark
from result.result import save_results
from pathlib import Path
import argparse

# Configurations  ("T_LUAD", "T_LUSC"),
   # ("T_STAD", "T_PAAD"),
   # ("T_GBM", "T_LGG"),
   # ("T_LIHC", "T_CHOL"),
   # ("T_KIRP", "T_KIRC"),
   # ("T_UCEC", "T_UCS"),
   # ("T_CESC", "T_ESCA"),
   # ("T_THYM", "T_HNSC"),
signature_sizes = [1,3,10,20] #we change signatutre sizes to [1,3,10,20]


benchmark_problems = [
    ("CD8 T", "B"), #("T_BRCA", "T_BLCA")
    # ("T_SKCM", "T_UVM"),
    
   
]
"""
   you can add more benchmark problems here
   switch between one-vs-rest and pairwise
   switch between cosg and scanpy(t-test, logreg, wilcoxon)
"""

# Main benchmarking clearly defined
def main(
    adata,
    marker_method,        # can switch to "scanpy" vs "cosg"
    comparison_mode,   # can switch to "one-vs-rest" vs "pairwise"
    classifier_method="svm"
):
    ranked_genes = {}
    if comparison_mode == "one-vs-rest":
        labels = adata.obs["label"].unique()
        #define rest here what is insdie rest
        labels = [label for label in labels if label != "rest"]
        benchmarks = [(label, "rest") for label in labels]
    else:
        benchmarks = benchmark_problems

    for label_a, label_b in benchmarks:
        print(f"\n Comparing: {label_a} vs {label_b}")

<<<<<<< HEAD
        # # Prepare labels clearly
        # if label_b == "rest":
        #     adata.obs["binary_labels"] = np.where(adata.obs["label"] == label_a, label_a, "rest")
        #     groups = [label_a]
        #     reference = "rest"
        # else:
        #     adata.obs['binary_labels'] = adata.obs['label'].replace({label_a: 'Cluster1', label_b: 'Cluster2'})
        #     groups = ["Cluster1"]
        #     reference = "Cluster2"
=======
        # Prepare labels clearly
        if cancer_b == "rest":
            adata.obs["binary_labels"] = np.where(adata.obs["label"] == cancer_a, cancer_a, "rest")
            groupby = "binary_labels"
            groups = [cancer_a]
            reference = "rest"
        else:
            adata.obs['binary_labels'] = adata.obs['label'].replace({cancer_a: 'Cluster1', cancer_b: 'Cluster2'})
            groupby = "binary_labels"
            groups = ["Cluster1"]
            reference = "Cluster2"
>>>>>>> 512b67e (add things)

        # Run chosen gene-ranking method clearly
        if marker_method == "cosg":
            marker_genes_df = run_cosg(adata, signature_sizes, label_a, groupby="label")
            key = f"{label_a}_vs_Rest"
            if key in marker_genes_df:
                ranked_genes[key] = marker_genes_df[key]
            else:
                print(f" Key {key} not found in COSG result.")

        elif marker_method=="scanpy":  # scanpy method
<<<<<<< HEAD
            method="wilcoxon"  # Can switch between "t-test", "logreg", "wilcoxon"
            print(f"Calling run_benchmark with groups={label_a}, reference={label_b}, groupby=binary_labels")
            ranked_genes.update( run_benchmark(adata, signature_sizes, "label", label_a, label_b, method=method,  classifier="svm"))
=======
            method="t-test"  # Can switch between "t-test", "logreg", "wilcoxon"
            print(f"Calling run_benchmark with groups={groups}, reference={reference}, groupby={groupby}")
            ranked_genes.update( run_benchmark(adata, groups, signature_sizes, comparison_mode, reference, groupby, cancer_a, cancer_b, method=method,  classifier="svm"))
>>>>>>> 512b67e (add things)
            print("\nFinal Benchmark Results:")
            print(ranked_genes)
    return ranked_genes

if __name__ == "__main__":
<<<<<<< HEAD
    parser = argparse.ArgumentParser(description="Yomix command-line tool")

    parser.add_argument(
        "file", type=str, nargs="?", default=None, help="the .ha5d file to open"
    )
        
    args = parser.parse_args()

    argument = args.file

    if argument:
        assert (
            args.file is not None
        ), "yomix: error: the following arguments are required: file"
        filearg = Path(args.file)
    else:
        filearg = Path(__file__).parent / "data" / "pbmc.h5ad"

    xd = sc.read_h5ad(filearg.absolute())
    #xd = sc.read_h5ad("/home/nisma/new_yomix/yomix/xd_tcga_labels_umap.h5ad")

=======
    datasets = {
        "TCGA": "/home/nisma/new_yomix/yomix/xd_tcga_labels_umap.h5ad",
        "PBMC": "/home/nisma/new_yomix/yomix/pbmc.h5ad",
        "Lawlor": "/home/nisma/new_yomix/yomix/lawlor.h5ad"
    }

   
    selected_dataset = "PBMC"  

    match selected_dataset:
        case "TCGA":
            dataset_path = datasets["TCGA"]
            print(" Running benchmark on TCGA dataset")
        case "PBMC":
            dataset_path = datasets["PBMC"]
            print(" Running benchmark on PBMC dataset")
        case "Lawlor":
            dataset_path = datasets["Lawlor"]
            print(" Running benchmark on Lawlor dataset")
        case _:
            print(" Invalid dataset! Defaulting to TCGA.")
            dataset_path = datasets["TCGA"]
    adata=sc.read_h5ad(dataset_path)
    #adata = sc.read_h5ad("/home/nisma/new_yomix/yomix/xd_tcga_labels_umap.h5ad")
    marker_method = "cosg"#  Can switch between "cosg" and "scanpy"
>>>>>>> 512b67e (add things)
    comparison_mode = "one-vs-rest"  #  Can switch between "one-vs-rest" and "pairwise"
    classifier_method = "svm"  #  Can switch between "svm", "logistic", "tree", "forest", "boosting"
    # method = "wilcoxon"  #  Can switch between "t-test", "logreg", "wilcoxon"

    results=main(
        xd,
        marker_method="cosg",
        comparison_mode=comparison_mode,
        classifier_method="svm"     # logistic, tree, forest, boosting
    )

    # save_results(results, comparison_mode, output_dir="results") #add directory where you want to save 