import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
from yomix_signature import compute_signature
import argparse
from scipy.sparse import issparse
from sklearn.model_selection import train_test_split
from sklearn.metrics import matthews_corrcoef
from sklearn.svm import SVC
import time
from sklearn.metrics import precision_score, recall_score, f1_score
import cosg

# Configurations  ("T_LUAD", "T_LUSC"),
# ("T_STAD", "T_PAAD"),
# ("T_GBM", "T_LGG"),
# ("T_LIHC", "T_CHOL"),
# ("T_KIRP", "T_KIRC"),
# ("T_UCEC", "T_UCS"),
# ("T_CESC", "T_ESCA"),
# ("T_THYM", "T_HNSC"),

benchmark_problems = [
    ("CD8 T", "B"),  # ("T_BRCA", "T_BLCA")
    # ("T_SKCM", "T_UVM"),
]


# Main benchmarking clearly defined
def main(
    xd,
    comparison_mode,  # can switch to "one-vs-rest" vs "pairwise"
    output_filename,  # name of outptu
    label_column,
    classifier_method="svm",
    nb_clf_runs=10,  # Number of runs for the classifier
    signatures_size=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20],
):
    """
    Runs benchmarking of feature selection methods.
    This function compares different marker selection methods
    (yomix, cosg, scanpy_wilcoxon, scanpy_t-test)
    and evaluates their performance using a classifier
    (default: SVM) across various gene signature sizes.
    It supports both "one-vs-rest" and "pairwise" comparison modes.
    Parameters
    ----------
    xd : AnnData
        Annotated data matrix (e.g., from Scanpy) containing gene
        expression data and cell metadata.
    marker_method : str
        Marker selection method to use ("scanpy" or "cosg").
    comparison_mode : str
        Mode of comparison, either "one-vs-rest" or "pairwise".
    output_filename : str
        Name of the output CSV file (without extension) to save
        benchmarking results.
    label_column : str
        Column in `xd.obs` containing cell type or class labels.
    classifier_method : str, optional
        Classifier to use for evaluation (default is "svm").
    nb_clf_runs : int, optional
        Number of classifier runs with different random seeds (default is 10).
    Returns
    -------
    ranked_genes : dict
        Dictionary containing ranked gene lists for each method and comparison."""

    all_methods = ["yomix", "cosg", "scanpy_wilcoxon", "scanpy_t-test"]
    runtime = {}
    results = []
    ranked_genes = {}
    ranked_genes = {}

    if comparison_mode == "one-vs-rest":
        labels = xd.obs[label_column].unique()
        # define rest here what is insdie rest
        labels = [label for label in labels if label != "rest"]
        benchmarks = [(label, "rest") for label in labels]
    else:
        benchmarks = benchmark_problems

    for label_a, label_b in benchmarks:
        signature_key = label_a + "_vs_" + label_b
        runtime[signature_key] = {}
        ranked_genes[signature_key] = {}
        print(f"\n Comparing: {label_a} vs {label_b}")
        # Prepare labels clearly
        if label_b == "rest":
            xd.obs["binary_labels"] = np.where(
                xd.obs[label_column] == label_a, label_a, "rest"
            )
        else:
            xd.obs["binary_labels"] = xd.obs[label_column].replace(
                {label_a: label_a, label_b: label_b}
            )
        groupby = "binary_labels"
        start_time = time.time()
        cosg.cosg(
            xd,
            key_added="cosg",
            use_raw=False,
            mu=100,
            expressed_pct=0.05,
            remove_lowly_expressed=True,
            n_genes_user=20,
            groupby=groupby,
        )
        runtime[signature_key]["cosg"] = time.time() - start_time
        ranked_genes[signature_key]["cosg"] = pd.DataFrame(
            xd.uns["cosg"]["names"], columns=xd.uns["cosg"]["names"].dtype.names
        )[label_a].values
        methods_scanpy = ["wilcoxon", "t-test"]

        for method_sc in methods_scanpy:
            start_time = time.time()
            sc.tl.rank_genes_groups(
                xd,
                groupby=groupby,
                groups=[label_a],
                reference=label_b,
                method=method_sc,
            )
            runtime[signature_key]["scanpy_" + method_sc] = time.time() - start_time
            ranked_genes[signature_key]["scanpy_" + method_sc] = xd.uns[
                "rank_genes_groups"
            ]["names"][label_a]

        xd.obs["id"] = [i for i in range(xd.obs.shape[0])]
        indices_label = xd[xd.obs[label_column] == label_a, :].obs["id"].to_list()
        start_time = time.time()
        genes, _, _ = compute_signature(
            adata=xd,
            means=xd.var["mean_values"],
            stds=xd.var["standard_deviations"],
            obs_indices_A=indices_label,
        )
        runtime[signature_key]["yomix"] = time.time() - start_time
        ranked_genes[signature_key]["yomix"] = xd.var.iloc[genes].index.tolist()

        for method in all_methods:
            for size in signatures_size:
                selected_genes = ranked_genes[signature_key][method][:size]
                X_subset = xd[:, selected_genes].X
                X_subset = (
                    X_subset.toarray() if hasattr(X_subset, "toarray") else X_subset
                )
                y_binary = np.where(xd.obs.binary_labels == label_a, 1, 0)

                for run in range(nb_clf_runs):
                    X_train, X_test, y_train, y_test = train_test_split(
                        X_subset,
                        y_binary,
                        test_size=0.3,
                        stratify=y_binary,
                        random_state=run,
                    )
                    clf = SVC(
                        kernel="linear", class_weight="balanced", random_state=run
                    )
                    # clf = KNeighborsClassifier()
                    clf.fit(X_train, y_train)

                    y_pred = clf.predict(X_test)
                    results.append(
                        {
                            "method": method,
                            "mcc": matthews_corrcoef(y_test, y_pred),
                            "precision": precision_score(y_test, y_pred),
                            "f1_score": f1_score(y_test, y_pred),
                            "recall": recall_score(y_test, y_pred),
                            "label_vs_rest": label_a,
                            "nb_genes": size,
                            "model": "svm",
                        }
                    )
    runtime_df = pd.DataFrame(runtime)
    runtime_df.to_csv(f"result/{output_filename}_runtime.csv")
    res_df = pd.DataFrame(results)
    res_df.to_csv(f"result/{output_filename}.csv")
    return ranked_genes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Yomix benchmark")
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

    def _to_dense(x):
        if issparse(x):
            return x.todense()
        else:
            return x

    xd.X = np.asarray(_to_dense(xd.X))
    min_norm = np.min(xd.X, axis=0)
    max_norm = np.max(xd.X, axis=0)
    xd.X = np.divide(xd.X - min_norm, max_norm - min_norm + 1e-8)

    def var_mean_values(adata) -> np.ndarray:
        return np.squeeze(np.asarray(np.mean(adata.X, axis=0)))

    def var_standard_deviations(adata) -> np.ndarray:
        return np.squeeze(np.asarray(np.std(adata.X, axis=0)))

    if "mean_values" not in xd.var:
        xd.var["mean_values"] = var_mean_values(xd)
        xd.var["standard_deviations"] = var_standard_deviations(xd)

    comparison_mode = "one-vs-rest"  # Can switch between "one-vs-rest" and "pairwise"
    classifier_method = (
        "svm"  # Can switch between "svm", "logistic", "tree", "forest", "boosting"
    )
    results = main(
        xd,
        comparison_mode=comparison_mode,
        output_filename="test_pbmc",
        label_column="label",
        classifier_method="svm",
        signatures_size=[1, 10, 20],
        nb_clf_runs=3,
    )

    # save_results(results, comparison_mode, output_dir="results")
    # add directory where you want to save
