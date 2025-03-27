import cosg
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import matthews_corrcoef
from sklearn.svm import SVC
import time
from sklearn.metrics import precision_score, recall_score, f1_score


def run_cosg(adata, n_genes=10, groupby="labels"):
    """
    Runs COSG analysis on the AnnData object (adata), retrieves top marker genes per cluster,
    and returns benchmarking metrics (MCC, Precision, Recall, F1) for selected genes.

    Parameters:
    - adata: AnnData object with single-cell data
    - n_genes: Number of top marker genes per cluster
    - groupby: Column name in adata.obs to group cells

    Returns:
    - Dictionary with structure matching scanpy, e.g., {"T_BRCA_vs_Rest": {...}, ...}
    """
    start_time = time.time()

    cosg.cosg(
        adata,
        key_added='cosg',
        use_raw=False,
        mu=100,
        expressed_pct=0.1,
        remove_lowly_expressed=True,
        n_genes_user=n_genes,
        groupby=groupby
    )

    end_time = time.time()
    cosg_time = end_time - start_time
    print(f"COSG Analysis completed in {cosg_time:.2f} seconds")

    marker_genes_df = pd.DataFrame(
        adata.uns['cosg']['names'],
        columns=adata.uns['cosg']['names'].dtype.names
    ).iloc[:n_genes]

    marker_gene_scores_df = pd.DataFrame(
        adata.uns['cosg']['scores'],
        columns=marker_genes_df.columns
    ).iloc[:n_genes, :]

    print("Available marker gene columns:", marker_genes_df.columns.tolist())

    signature_sizes = [1, 3, 10, 20]
    n_runs = 10

    labels_all = adata.obs[groupby].values
    final_results = {}

    # Remove 'rest' column if present
    marker_genes_df = marker_genes_df[[col for col in marker_genes_df.columns if col.lower() != "rest"]]

    for size in signature_sizes:
        print(f"\n=== Signature Size: {size} ===")

        for label in marker_genes_df.columns:
            print(f"\nEvaluating: {label} vs Rest")

            selected_genes = marker_genes_df[label].iloc[:size].tolist()
            selected_genes = list(pd.unique(selected_genes))
            gene_indices = [np.where(adata.var_names == gene)[0][0] for gene in selected_genes]

            X_subset = adata.X[:, gene_indices]
            X_subset = X_subset.toarray() if hasattr(X_subset, 'toarray') else X_subset

            y_binary = np.where(labels_all == label, 1, 0)

            if sum(y_binary) < 2 or sum(y_binary == 0) < 2:
                print(f"Skipped '{label}': Not enough samples.")
                continue

            mcc_scores = []
            precision_scores = []
            recall_scores = []
            f1_scores = []

            for run in range(n_runs):
                X_train, X_test, y_train, y_test = train_test_split(
                    X_subset, y_binary,
                    test_size=0.3, stratify=y_binary, random_state=run
                )

                clf = SVC(kernel='linear', class_weight='balanced', random_state=run)
                clf.fit(X_train, y_train)

                y_pred = clf.predict(X_test)
                mcc_scores.append(matthews_corrcoef(y_test, y_pred))
                precision_scores.append(precision_score(y_test, y_pred, average="binary"))
                recall_scores.append(recall_score(y_test, y_pred, average="binary"))
                f1_scores.append(f1_score(y_test, y_pred, average="binary"))

            mean_mcc = np.mean(mcc_scores)
            mean_precision = np.mean(precision_scores)
            mean_recall = np.mean(recall_scores)
            mean_f1 = np.mean(f1_scores)

            key = f"{label}_vs_Rest"

            # Initialize per-label result if first time
            if key not in final_results:
                final_results[key] = {}

            final_results[key]["Time_Taken"] = cosg_time
            final_results[key][f"{size}_features_mcc"] = mean_mcc
            final_results[key][f"{size}_features_precision"] = mean_precision
            final_results[key][f"{size}_features_recall"] = mean_recall
            final_results[key][f"{size}_features_f1"] = mean_f1

    return final_results
