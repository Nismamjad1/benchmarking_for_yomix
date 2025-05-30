import scanpy as sc
import numpy as np
import pandas as pd
from method.cosg_file import run_cosg
from method.scan_py import run_benchmark
from result.result import save_results

# Configurations  ("T_LUAD", "T_LUSC"),
   # ("T_STAD", "T_PAAD"),
   # ("T_GBM", "T_LGG"),
   # ("T_LIHC", "T_CHOL"),
   # ("T_KIRP", "T_KIRC"),
   # ("T_UCEC", "T_UCS"),
   # ("T_CESC", "T_ESCA"),
   # ("T_THYM", "T_HNSC"),
signature_sizes = [1] #we change signatutre sizes to [1,3,10,20]


benchmark_problems = [
    ("T_BRCA", "T_BLCA"),
    ("T_SKCM", "T_UVM"),
    
   
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
    results = {}
    ranked_genes = {}
    adata_temp=adata.copy()
    if comparison_mode == "one-vs-rest":
        labels = adata_temp.obs["labels"].unique()
        #define rest here what is insdie rest
        labels = [label for label in labels if label != "rest"]
        benchmarks = [(label, "rest") for label in labels]
    else:
        benchmarks = benchmark_problems

    for cancer_a, cancer_b in benchmarks:
        if cancer_a == "N_Regional_Lymph_Node":
            print(f"Skipping {cancer_a} as it is N_Regional_Lymph_Node")
            continue  # Skip to the next iteration
        print(f"\n Comparing: {cancer_a} vs {cancer_b}")

        # Prepare labels clearly
        if cancer_b == "rest":
            adata_temp.obs["binary_labels"] = np.where(adata_temp.obs["labels"] == cancer_a, cancer_a, "rest")
            groupby = "binary_labels"
            groups = [cancer_a]
            reference = "rest"
        else:
            adata.obs['binary_labels'] = adata_temp.obs['label'].replace({cancer_a: 'Cluster1', cancer_b: 'Cluster2'})
            groupby = "binary_labels"
            groups = ["Cluster1"]
            reference = "Cluster2"

        # Run chosen gene-ranking method clearly
        if marker_method == "cosg":
            marker_genes_df = run_cosg(adata_temp, signature_sizes, groupby=groupby)
            key = f"{cancer_a}_vs_Rest"
            if key in marker_genes_df:
                ranked_genes[key] = marker_genes_df[key]
            else:
                print(f" Key {key} not found in COSG result.")

        elif marker_method=="scanpy":  # scanpy method
            method="wilcoxon"  # Can switch between "t-test", "logreg", "wilcoxon"
            print(f"Calling run_benchmark with groups={groups}, reference={reference}, groupby={groupby}")
            ranked_genes.update( run_benchmark(adata_temp, groups, signature_sizes, comparison_mode, reference, groupby, cancer_a, cancer_b, method=method,  classifier="svm"))
            print("\nFinal Benchmark Results:")
            print(ranked_genes)
    return ranked_genes

if __name__ == "__main__":
    datasets = {
        "TCGA": "/home/nisma/new_yomix/yomix/xd_tcga_labels_umap.h5ad",
        "PMBC": "/home/nisma/new_yomix/yomix/pmbc_data.h5ad",
        "Lawlor": "/home/nisma/new_yomix/yomix/lawlor.h5ad"
    }

   
    selected_dataset = "TCGA"  

    match selected_dataset:
        case "TCGA":
            dataset_path = datasets["TCGA"]
            print(" Running benchmark on TCGA dataset")
        case "PMBC":
            dataset_path = datasets["PMBC"]
            print(" Running benchmark on PMBC dataset")
        case "Lawlor":
            dataset_path = datasets["Lawlor"]
            print(" Running benchmark on Lawlor dataset")
        case _:
            print(" Invalid dataset! Defaulting to TCGA.")
            dataset_path = datasets["TCGA"]
    adata=sc.read_h5ad(dataset_path)
    #adata = sc.read_h5ad("/home/nisma/new_yomix/yomix/xd_tcga_labels_umap.h5ad")
    marker_method = "scanpy"#  Can switch between "cosg" and "scanpy"
    comparison_mode = "one-vs-rest"  #  Can switch between "one-vs-rest" and "pairwise"
    classifier_method = "svm"  #  Can switch between "svm", "logistic", "tree", "forest", "boosting"
    method = "wilcoxon"  #  Can switch between "t-test", "logreg", "wilcoxon"
    
    results=main(
        adata,
        marker_method=marker_method,
        comparison_mode=comparison_mode,
        classifier_method="svm"     # logistic, tree, forest, boosting
    )
    
    save_results(results, marker_method, comparison_mode, method, selected_dataset, output_dir="/home/nisma/new_yomix/yomix/project/output/PBMC") #add directory where you want to save 