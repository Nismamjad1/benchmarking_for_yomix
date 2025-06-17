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

num_runs = 10

def run_benchmark(adata, signature_sizes , comparison_mode, groupby, label_a, label_b, method,  classifier="svm"):
    
    results = {}
    
   
    group_counts = adata.obs[groupby].value_counts()
    print(group_counts)  # Check how many samples exist per group

    # Filter out groups with fewer than 2 cells
    valid_groups = group_counts[group_counts >= 2].index.tolist()
    '''
    if len(valid_groups) < 2:
        pdb.set_trace()
        print(f"Skipping {label_a} vs {label_b}: Not enough samples.")
        return'
    '''
    # Differential expression analysis clearly
    start_time = time.time()
    sc.tl.rank_genes_groups(adata, groupby=groupby, groups=[label_a], reference=label_b, method=method)
    end_time = time.time()
    de_time_taken = end_time - start_time
    print(f" Time taken for ranking genes ({label_a} vs {label_b}): {de_time_taken:.2f} seconds")

    # Access ranked genes
    ranked_genes = adata.uns['rank_genes_groups']['names'][label_a]

    results[f"{label_a}_vs_{label_b}"] = {"DE_Time_Taken": de_time_taken}

    # Access clearly ranked genes
    ranked_genes = adata.uns['rank_genes_groups']['names'][label_a]

    results[f"{label_a}_vs_{label_b}"] = {}

    for size in signature_sizes:
        top_genes = ranked_genes[:size]

        mcc_scores = []
        precision_scores=[]
        recall_scores=[]
        f1_scores=[]
        

        for run in range(10):
            idx_A = adata.obs.index[adata.obs[groupby] == label_a]
            if label_b == "rest":
                idx_B = adata.obs.index[adata.obs[groupby] != label_a]
            else:
                idx_B = adata.obs.index[adata.obs[groupby] == label_b]

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
        std_precision = np.std(precision_scores)

        mean_precision = np.mean(precision_scores)
        std_recall = np.std(recall_scores)
        mean_recall = np.mean(recall_scores)
        std_f1 = np.std(f1_scores)
        mean_f1 = np.mean(f1_scores)
        std_mcc = np.std(mcc_scores)

        print(f"Signature size: {size}, Mean MCC (over {len(mcc_scores)} runs): {mean_mcc:.4f}")
        results[f"{label_a}_vs_{label_b}"]["DE_Time_Taken"] = de_time_taken
        
        results[f"{label_a}_vs_{label_b}"][f"{size}_features_mcc"] = mean_mcc
        results[f"{label_a}_vs_{label_b}"][f"{size}_features_mcc_std"] = std_mcc
        
        
        
        results[f"{label_a}_vs_{label_b}"][f"{size}_features_precision"] = mean_precision
        results[f"{label_a}_vs_{label_b}"][f"{size}_features_precision_std"] = std_precision

        results[f"{label_a}_vs_{label_b}"][f"{size}_features_recall"] = mean_recall
        results[f"{label_a}_vs_{label_b}"][f"{size}_features_recall_std"] = std_recall

        results[f"{label_a}_vs_{label_b}"][f"{size}_features_f1"] = mean_f1
        results[f"{label_a}_vs_{label_b}"][f"{size}_features_f1_std"] = std_f1   

    return results