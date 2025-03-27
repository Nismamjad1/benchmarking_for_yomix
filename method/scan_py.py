import numpy as np
import scanpy as sc
from sklearn.model_selection import train_test_split
from sklearn.metrics import matthews_corrcoef
from sklearn.svm import SVC
import time
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import classification_report, precision_score, recall_score, f1_score
import pdb
# Define benchmark problems (pairwise)
benchmark_problems = [
    ("T_BRCA", "T_BLCA"),
    ("T_SKCM", "T_UVM"),
   
]

signature_sizes = [1]
num_runs = 10

def run_benchmark(adata, groups, comparison_mode, reference, groupby, cancer_a, cancer_b, method,  classifier="svm"):
    
    results = {}
    
   
    group_counts = adata.obs[groupby].value_counts()
    print(group_counts)  # Check how many samples exist per group

    # Filter out groups with fewer than 2 cells
    valid_groups = group_counts[group_counts >= 2].index.tolist()
    '''
    if len(valid_groups) < 2:
        pdb.set_trace()
        print(f"Skipping {cancer_a} vs {cancer_b}: Not enough samples.")
        return'
    '''
    # Differential expression analysis clearly
    start_time = time.time()
    sc.tl.rank_genes_groups(adata, groupby=groupby, groups=groups, reference=reference, method=method)
    end_time = time.time()
    de_time_taken = end_time - start_time
    print(f" Time taken for ranking genes ({cancer_a} vs {cancer_b}): {de_time_taken:.2f} seconds")

    # Access ranked genes
    ranked_genes = adata.uns['rank_genes_groups']['names'][groups[0]]

    results[f"{cancer_a}_vs_{cancer_b}"] = {"DE_Time_Taken": de_time_taken}

    # Access clearly ranked genes
    ranked_genes = adata.uns['rank_genes_groups']['names'][groups[0]]

    results[f"{cancer_a}_vs_{cancer_b}"] = {}

    for size in [1, 3, 10, 20]:
        top_genes = ranked_genes[:size]

        mcc_scores = []
        precision_scores=[]
        recall_scores=[]
        f1_scores=[]
        

        for run in range(10):
            if cancer_b == "rest":
                idx_A = adata.obs.index[adata.obs[groupby] == cancer_a]
                idx_B = adata.obs.index[adata.obs[groupby] != cancer_a]
            else:
                idx_A = adata.obs.index[adata.obs[groupby] == "Cluster1"]
                idx_B = adata.obs.index[adata.obs[groupby] == "Cluster2"]

            # Extract clearly expression data
            subset_A = adata[idx_A, top_genes].X.toarray() if hasattr(adata[idx_A, top_genes].X, "toarray") else adata[idx_A, top_genes].X
            subset_B = adata[idx_B, top_genes].X.toarray() if hasattr(adata[idx_B, top_genes].X, "toarray") else adata[idx_B, top_genes].X

            # Combine clearly subsets and labels
            X = np.vstack((subset_A, subset_B))
            y = np.concatenate((np.ones(len(idx_A)), np.zeros(len(idx_B))))

            # Train-test split
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, stratify=y, random_state=np.random.randint(10000))

            # Classifier clearly chosen
            if classifier == "logistic":
                clf = LogisticRegression(class_weight='balanced', max_iter=1000)
            elif classifier == "boosting":
                clf = HistGradientBoostingClassifier()
            elif classifier == "forest":
                clf = RandomForestClassifier()
            elif classifier == "tree":
                clf = DecisionTreeClassifier()
            else:  # Default clearly SVM
                clf = SVC(kernel='linear', class_weight='balanced')

            clf.fit(X_train, y_train)

            # Clearly evaluate predictions
            y_pred = clf.predict(X_test)
            mcc = matthews_corrcoef(y_test, y_pred)
            mcc_scores.append(mcc)
            precision = precision_score(y_test, y_pred, average="binary")
            recall = recall_score(y_test, y_pred, average="binary")
            f1 = f1_score(y_test, y_pred, average="binary")

            precision_scores.append(precision)
            recall_scores.append(recall)
            f1_scores.append(f1)

        mean_mcc = np.mean(mcc_scores)
        # Compute mean scores over multiple runs
        mean_mcc = np.mean(mcc_scores)
        mean_precision = np.mean(precision_scores)
        mean_recall = np.mean(recall_scores)
        mean_f1 = np.mean(f1_scores)

        print(f"Signature size: {size}, Mean MCC (over {len(mcc_scores)} runs): {mean_mcc:.4f}")
        results[f"{cancer_a}_vs_{cancer_b}"]["DE_Time_Taken"] = de_time_taken
        results[f"{cancer_a}_vs_{cancer_b}"][f"{size}_features_mcc"] = mean_mcc
        results[f"{cancer_a}_vs_{cancer_b}"][f"{size}_features_precision"] = mean_precision
        results[f"{cancer_a}_vs_{cancer_b}"][f"{size}_features_recall"] = mean_recall
        results[f"{cancer_a}_vs_{cancer_b}"][f"{size}_features_f1"] = mean_f1    
    

    

    return results