import scanpy as sc
import anndata
import pandas as pd
import time
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, matthews_corrcoef, confusion_matrix, classification_report
#from scipy.stats import wasserstein_distance
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from imblearn.over_sampling import SMOTE  # 📊 For class imbalance handling i impoterd this and downled this
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
#import scvelo as scv
'''

results_file = "/home/nisma/write/pbmc3k.h5ad"
adata = anndata.read_h5ad(results_file)
print(adata)
sc.pl.umap(adata, color=["leiden"])


# Start the timer
start_time = time.time()

# Compare cluster "0" to all other clusters
sc.tl.rank_genes_groups(adata, groupby="leiden", groups=["2"], reference="rest", method="wilcoxon")

# End the timer
end_time = time.time()

# Print elapsed time
print(f"Time taken for comparing cluster 0 to all other clusters: {end_time - start_time:.2f} seconds")

# Plot the top 20 ranked genes for cluster "0"
sc.pl.rank_genes_groups(adata, groups=["2"], n_genes=20, sharey=False)

# Access and print the top-ranked genes for cluster "0"
top_genes_cluster_0 = adata.uns["rank_genes_groups"]["names"]["2"]
print(f"\nTop 20 genes for cluster 0 compared to all other clusters:")
print(top_genes_cluster_0[:20])


'''
'''
# Specify the cluster you want to analyze
cluster_of_interest = "2"  # Change this to the cluster you're interested in

# Filter cells belonging to the specific cluster
adata_cluster = adata[adata.obs['leiden'] == cluster_of_interest]
# Start the timer
start_time = time.time()


# Perform rank_genes_groups
sc.tl.rank_genes_groups(adata_cluster, "leiden", method="t-test")

# End the timer
end_time = time.time()
# Print the elapsed time
print(f"Time taken for cluster {cluster_of_interest}: {end_time - start_time:.2f} seconds")

# Access and print top-ranked genes for the cluster
top_genes = adata_cluster.uns['rank_genes_groups']['names']
print(f"\nTop genes for cluster {cluster_of_interest}:")
print(top_genes['0'][:10])  



# Print the elapsed time
print(f"Time taken: {end_time - start_time:.2f} seconds")

# Access and print top-ranked genes
top_genes = adata.uns['rank_genes_groups']['names']  # Access the ranked gene names
for group in top_genes.dtype.names:  # Iterate through the groups (e.g., clusters)
    print(f"\nTop genes for group {group}:")
    print(top_genes[group][:10])

'''


#################################################################for the tcga data#################################
import scanpy as sc
import anndata
import pandas as pd
import time
import numpy as np
#from scipy.stats import wasserstein_distance
# Load the data
'''
results_file = "/home/nisma/write/xd_tcga.h5ad"


adata = anndata.read_h5ad(results_file)

adata.var_names_make_unique()
sc.pl.highest_expr_genes(adata, n_top=20)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
# annotate the group of mitochondrial genes as "mt"
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)


# Initial data cleanup: Filter cells based on QC metrics
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)
sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")

# Filter cells with high mitochondrial content or extreme gene counts
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

# Normalize total counts and handle potential data issues
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


# Identify highly variable genes
sc.pl.highly_variable_genes(adata)
#adata.raw = adata

# Filter the dataset to include only highly variable genes
adata = adata[:, adata.var.highly_variable]

# Regress out unwanted sources of variation
#sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])


# Check and handle NaN values in adata.X
print("Checking for NaN values in adata.X before PCA...")
if np.any(np.isnan(adata.X)):
    print("NaN values detected in adata.X. Replacing with zeros.")
    adata.X = np.nan_to_num(adata.X)  # Replace NaNs with 0

# Re-check for infinite values
if np.any(np.isinf(adata.X)):
    print("Infinite values detected in adata.X. Replacing with maximum finite values.")
    adata.X[np.isinf(adata.X)] = np.nanmax(adata.X[np.isfinite(adata.X)])

# Scale the data to unit variance and limit max values
sc.pp.scale(adata, max_value=10)

# Perform PCA for dimensionality reduction
sc.tl.pca(adata, svd_solver="arpack")
#sc.pl.pca(adata, color="CST3")
sc.pl.pca_variance_ratio(adata, log=True)

# Save intermediate results
results_file = "/home/nisma/write/results.h5ad"
adata.write(results_file)
# Build the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Perform Leiden clustering to generate community labels
sc.tl.leiden(
    adata,
    resolution=0.9,
    random_state=0,
    n_iterations=2,
    directed=False,
)


# PAGA graph abstraction
sc.tl.paga(adata)  # Now Leiden labels are available
sc.pl.paga(adata, plot=False)  # Set plot=True to visualize

# Compute UMAP using PAGA initialization
sc.tl.umap(adata, init_pos="paga")
sc.pl.umap(adata, color=["leiden"])

# Save results
results_file = "/home/nisma/write/xd_tcga_huge.h5ad"
adata.write(results_file)
adata = sc.read_h5ad(results_file)
print(adata)

start_time = time.time()
sc.tl.rank_genes_groups(adata, groupby="leiden", groups=["26"], reference="rest", method="wilcoxon")
end_time = time.time()

# Print elapsed time
print(f"Time taken for rank_genes_groups: {end_time - start_time:.2f} seconds")
expression_summaries = {}
wasserstein_results = {}


top_genes = adata.uns["rank_genes_groups"]["names"]
cluster_of_interest = "26"
top_20_genes = top_genes[cluster_of_interest][:20]
print(f"\n🔬 Top 20 genes for Cluster {cluster_of_interest}:\n{top_20_genes}")

# Get the indices of samples belonging to cluster 26
cluster_26_cells = adata.obs_names[adata.obs["leiden"] == "26"]

# Print the sample names
print("Sample names in cluster 26:")
print(cluster_26_cells.tolist())
# Extract data for cluster 26
cluster_26_data = adata[adata.obs["leiden"] == "26", :]
print("Shape of cluster 26 subset:", cluster_26_data.shape)

# Subset Data for Cluster 26 and Rest
cluster_26_cells = adata.obs[adata.obs["leiden"] == cluster_of_interest].index
rest_cells = adata.obs[adata.obs["leiden"] != cluster_of_interest].index

# Subset the AnnData object for the rest cells
rest_data = adata[rest_cells, :]

# Print the shape of the rest data
print(f"Shape of Rest Data: {rest_data.shape}")  # (num_samples, num_genes)

cluster_26_data = adata[cluster_26_cells, top_20_genes].X.toarray() if hasattr(adata[cluster_26_cells, top_20_genes].X, "toarray") else adata[cluster_26_cells, top_20_genes].X
rest_data = adata[rest_cells, top_20_genes].X.toarray() if hasattr(adata[rest_cells, top_20_genes].X, "toarray") else adata[rest_cells, top_20_genes].X

# ✅ Combine Data and Labels
X = np.vstack((cluster_26_data, rest_data))  # Combined data
y = np.concatenate((np.ones(cluster_26_data.shape[0]), np.zeros(rest_data.shape[0])))  # 1 for Cluster 27, 0 for Rest

print(f"\n📊 Shape of Combined Data     (X): {X.shape}")
print(f"📊 Label Distribution in y: {np.bincount(y.astype(int))}")

# ✅ Train-Test Split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, stratify=y, random_state=42)

# ✅ Apply SMOTE to Balance Classes
smote = SMOTE(random_state=42)
X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)

print(f"Before SMOTE → Class Distribution: {np.bincount(y_train.astype(int))}")
print(f"After SMOTE  → Class Distribution: {np.bincount(y_train_resampled.astype(int))}")

# ✅ Train Logistic Regression Classifier
clf = LogisticRegression(max_iter=1000, class_weight='balanced')
clf.fit(X_train_resampled, y_train_resampled)

# ✅ Classifier Evaluation
y_pred = clf.predict(X_test)

accuracy = accuracy_score(y_test, y_pred)
mcc = matthews_corrcoef(y_test, y_pred)
cm = confusion_matrix(y_test, y_pred)

print(f"\n Classifier Performance with SMOTE:")
print(f"🔹 Accuracy: {accuracy:.4f}")
print(f"🔹 MCC Score: {mcc:.4f}")
print(f"🔹 Confusion Matrix:\n{cm}")
print("\n Classification Report:\n", classification_report(y_test, y_pred))


for group in top_genes.dtype.names:  # Iterate through groups
    print(f"\nTop genes for group {group}:")
    top_20_genes = top_genes[group][:20]
    print(top_20_genes)
    '''
'''
    # Calculate the expression levels for the top 20 genes of this group
    adata_top_genes = adata[:, top_20_genes]  # Subset the data to include only the top 20 genes
    expression_sum = adata_top_genes.X.sum(axis=1)  # Sum expression levels across all cells
    
    # Store the results in a DataFrame for this group
    expression_summaries[group] = pd.DataFrame(
        data={k
            "Cell": adata.obs_names,
            "Expression_Sum": expression_sum,
        }
    )
    # Print the first few rows of the summary for this group
    print(expression_summaries[group].head())

   ####################################################################################
    # 📏 Compute Wasserstein Distance Between Cluster 0 and Rest for the Top 20 print("clusteer_0_data", cluster_0_cells )Genes
    ####################################################################################

    print(f"\n🔬 Wasserstein Distances for Group {group}:")
    
    # Split data into Cluster 0 and the rest
    cluster_0_cells = adata.obs[adata.obs["leiden"] == "27"].index
    print("clusteer_0_data_index", cluster_0_cells )
    rest_cells = adata.obs[adata.obs["leiden"] != "27"].index

    cluster_0_data = adata[cluster_0_cells, top_20_genes]
    print("clusteer_0_data", cluster_0_cells )
    rest_data = adata[rest_cells, top_20_genes]

    # Initialize list to store distances for each gene
    gene_wasserstein_distances = []

    for gene in top_20_genes:
        # Correct indexing by converting gene name to index
        gene_idx = cluster_0_data.var_names.get_loc(gene)
        # Extract gene expression data for Cluster 0 and Rest
        expr_cluster_0 = cluster_0_data[:, gene_idx].X.toarray().flatten() if hasattr(cluster_0_data[:, gene_idx].X, "toarray") else cluster_0_data[:, gene_idx].X.flatten()
        print("expr_cluster",expr_cluster_0)
        expr_rest = rest_data[:, gene_idx].X.toarray().flatten() if hasattr(rest_data[:, gene_idx].X, "toarray") else rest_data[:, gene_idx].X.flatten()
        print("expr_rest",expr_rest)

        # Compute Wasserstein Distance
        distance = wasserstein_distance(expr_cluster_0, expr_rest)
        gene_wasserstein_distances.append((gene, distance))

        print(f"Gene: {gene}, Wasserstein Distance: {distance:.4f}")

    # Store Wasserstein results in DataFrame
    wasserstein_results[group] = pd.DataFrame(gene_wasserstein_distances, columns=["Gene", "Wasserstein_Distance"])

    # Calculate and print the average Wasserstein Distance for this group
    avg_distance = np.mean([dist for _, dist in gene_wasserstein_distances])
    print(f"\n📏 Average Wasserstein Distance for Group {group}: {avg_distance:.4f}\n")


#from scipy.stats import wasserstein_distance
# Load the data
# Save intermediate results
results_file = "/home/nisma/write/xd_tcga_labels_umap.h5ad"
adata = anndata.read_h5ad(results_file)
adata.write(results_file)
# Build the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Perform Leiden clustering to generate community labels
sc.tl.leiden(
    adata,
    resolution=0.9,
    random_state=0,
    n_iterations=2,
    directed=False,
)

# PAGA graph abstraction
sc.tl.paga(adata)  # Now Leiden labels are available
sc.pl.paga(adata, plot=False)  # Set plot=True to visualize

# Compute UMAP using PAGA initialization
sc.tl.umap(adata, init_pos="paga")
sc.pl.umap(adata, color=["leiden"])

# Save results
results_file = "/home/nisma/write/xd_tcga_labels_umap.h5ad"
adata.write(results_file)
adata = sc.read_h5ad(results_file)
print(adata)

start_time = time.time()
sc.tl.rank_genes_groups(adata, groupby="leiden", groups=["26"], reference="rest", method="wilcoxon")
end_time = time.time()

top_genes = adata.uns["rank_genes_groups"]["names"]
cluster_of_interest = "26"
top_20_genes = top_genes[cluster_of_interest][:20]
print(f"\n🔬 Top 20 genes for Cluster {cluster_of_interest}:\n{top_20_genes}")
'''
import scanpy as sc
import anndata
import time

# Load the data
results_file = "/home/nisma/write/xd_tcga_labels_umap.h5ad"
adata = anndata.read_h5ad(results_file)

'''

# Build the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Perform Leiden clustering to generate community labels
sc.tl.leiden(
    adata,
    resolution=0.9,
    random_state=0,
    n_iterations=2,
    directed=False,
)

# PAGA graph abstraction
sc.tl.paga(adata)  # Now Leiden labels are available
sc.pl.paga(adata, plot=False)  # Set plot=True to visualize

# Compute UMAP using PAGA initialization
sc.tl.umap(adata, init_pos="paga")
sc.pl.umap(adata, color=["leiden"])

# Save results
adata.write(results_file)

# Reload saved data
adata = sc.read_h5ad(results_file)
print(adata)

# Perform differential expression analysis using 'labels' instead of 'leiden'
start_time = time.time()
sc.tl.rank_genes_groups(adata, groupby="labels", groups=["T_BRCA", "T_BLCA"], method="wilcoxon")
end_time = time.time()



# Create a new grouping in the metadata
adata.obs['new_labels'] = adata.obs['labels'].replace({
    'T_BRCA': 'T_CANCER', 
    'T_BLCA': 'T_CANCER'
})

# Run the differential expression analysis
start_time = time.time()
sc.tl.rank_genes_groups(adata, groupby="new_labels", groups=["T_CANCER"], method="wilcoxon")
end_time = time.time()

print(f"Time taken for gene ranking analysis: {end_time - start_time:.2f} seconds")

# Identify cell indices for T_CANCER and the rest
cancer_cells = adata.obs.index[adata.obs["new_labels"] == "T_CANCER"]
rest_cells = adata.obs.index[adata.obs["new_labels"] != "T_CANCER"]
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
# Access ranked genes
result = adata.uns['rank_genes_groups']
top_20_genes= result['names']['T_CANCER'][:20]

# Print the top 20 differentially expressed genes
print("Top 20 differentially expressed genes for T_CANCER vs. rest:")
for i, gene in enumerate(top_20_genes, 1):
    print(f"{i}. {gene}")
  # Get top genes for T_CANCER vs. rest


# Extract expression data for T_CANCER and the rest
cancer_data = adata[cancer_cells, top_20_genes].X.toarray() if hasattr(adata[cancer_cells, top_20_genes].X, "toarray") else adata[cancer_cells, top_20_genes].X
rest_data = adata[rest_cells, top_20_genes].X.toarray() if hasattr(adata[rest_cells, top_20_genes].X, "toarray") else adata[rest_cells, top_20_genes].X

# Combine data and create labels (1 for T_CANCER, 0 for rest)
X = np.vstack((cancer_data, rest_data))
y = np.concatenate((np.ones(cancer_data.shape[0]), np.zeros(rest_data.shape[0])))

print(f"\n Shape of Combined Data (X): {X.shape}")
print(f" Label Distribution in y: {np.bincount(y.astype(int))}")

# Train-Test Split (70% train, 30% test)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, stratify=y, random_state=42)

# Train Logistic Regression Classifier without SMOTE
clf = LogisticRegression(max_iter=1000, class_weight='balanced')
clf.fit(X_train, y_train)

# Classifier Evaluation
y_pred = clf.predict(X_test)

# Calculate Accuracy and MCC Score
accuracy = accuracy_score(y_test, y_pred)
mcc = matthews_corrcoef(y_test, y_pred)

# Generate Confusion Matrix
cm = confusion_matrix(y_test, y_pred)

# Print evaluation metrics
print(f"\n🔹 Accuracy: {accuracy:.4f}")
print(f"🔹 MCC Score: {mcc:.4f}")
print(f"🔹 Confusion Matrix:\n{cm}")
print("\n Classification Report:\n", classification_report(y_test, y_pred))


'''

'''
# Train-Test Split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, stratify=y, random_state=42)

# Train Logistic Regression Classifier (without SMOTE)
clf = LogisticRegression(max_iter=1000, class_weight='balanced')
clf.fit(X_train, y_train)

# Classifier Evaluation
y_pred = clf.predict(X_test)

accuracy = accuracy_score(y_test, y_pred)
mcc = matthews_corrcoef(y_test, y_pred)
cm = confusion_matrix(y_test, y_pred)

print(f"\n Accuracy: {accuracy:.4f}")
print(f" MCC Score: {mcc:.4f}")
print(f" Confusion Matrix:\n{cm}")
print("\n Classification Report:\n", classification_report(y_test, y_pred))
'''
'''
# Define benchmark problems
benchmark_problems = [
    ("T_BRCA", "T_BLCA"),
    ("T_SKCM", "T_UVM"),
    ("T_LUAD", "T_LUSC"),
    ("T_STAD", "T_PAAD"),
    ("T_GBM", "T_LGG"),
    ("T_LIHC", "T_CHOL"),
    ("T_KIRP", "T_KIRC"),
    ("T_UCEC", "T_UCS"),
    ("T_CESC", "T_ESCA"),
    ("T_THYM", "T_HNSC"),
]

# Signature sizes to test
signature_sizes = [1, 3, 10, 20, 25]

# Function to run the benchmark
def run_benchmark(adata, method="wilcoxon"):
    results = {}
    
    for cancer_1, cancer_2 in benchmark_problems:
        print(f"\nRunning benchmark for: {cancer_1} vs {cancer_2}")
        
        # Create new labels
        adata.obs['new_labels'] = adata.obs['labels'].replace({cancer_1: 'Cluster1', cancer_2: 'Cluster2'})

        # Differential expression analysis
        start_time = time.time()
        sc.tl.rank_genes_groups(adata, groupby="new_labels", groups=["Cluster1"], reference="Cluster2", method=method)
        end_time = time.time()
        print(f"Time taken for gene ranking analysis: {end_time - start_time:.2f} seconds")

        # Extract top genes and evaluate at different signature sizes
        ranked_genes = adata.uns['rank_genes_groups']['names']['Cluster1']
        results[f"{cancer_1}_vs_{cancer_2}"] = {}
        
        for size in signature_sizes:
            top_genes = ranked_genes[:size]
            
            # Extract expression data
            cluster1_cells = adata.obs.index[adata.obs['new_labels'] == 'Cluster1']
            cluster2_cells = adata.obs.index[adata.obs['new_labels'] == 'Cluster2']
            import numpy as np

            cluster1_data = adata[cluster1_cells, top_genes].X.toarray() if hasattr(adata[cluster1_cells, top_genes].X, "toarray") else adata[cluster1_cells, top_genes].X
            cluster2_data = adata[cluster2_cells, top_genes].X.toarray() if hasattr(adata[cluster2_cells, top_genes].X, "toarray") else adata[cluster2_cells, top_genes].X
            
            # Combine data and create labels
            X = np.vstack((cluster1_data, cluster2_data))
            y = np.concatenate((np.ones(cluster1_data.shape[0]), np.zeros(cluster2_data.shape[0])))
            
            # Train-test split
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, stratify=y, random_state=42)
            
            # Train classifier
            clf = LogisticRegression(max_iter=1000, class_weight='balanced')
            clf.fit(X_train, y_train)
            
            # Evaluate classifier
            y_pred = clf.predict(X_test)
            accuracy = accuracy_score(y_test, y_pred)
            mcc = matthews_corrcoef(y_test, y_pred)
            cm = confusion_matrix(y_test, y_pred)
            
            results[f"{cancer_1}_vs_{cancer_2}"][f"{size}_features"] = {
                'accuracy': accuracy,
                'mcc': mcc,
                'confusion_matrix': cm.tolist(),
                'classification_report': classification_report(y_test, y_pred, output_dict=True)
            }
            
            print(f"\n🔹 Results for {size} features:")
            print(f"Accuracy: {accuracy:.4f}")
            print(f"MCC Score: {mcc:.4f}")
            print(f"Confusion Matrix: {cm}")
            
            mcc = matthews_corrcoef(y_test, y_pred)
    return mcc


# Number of runs for averaging
num_runs = 2
signature_sizes = [1, 3, 10, 20]  # Only these sizes are considered

# Initialize a dictionary to store MCC values
mcc_results = {size: [] for size in signature_sizes}

# Run the benchmark multiple times
for _ in range(num_runs):
    for size in signature_sizes:
        mcc = run_benchmark(adata, method="wilcoxon")  # Run benchmark for each size
        mcc_results[size].append(mcc)  # Store MCC value

# Compute mean MCC for each signature size
mean_mcc_results = {size: np.mean(mcc_values) for size, mcc_values in mcc_results.items()}

# Print results
print("\nMean MCC for Each Feature Set:")
for size, mean_mcc in mean_mcc_results.items():
    print(f"  {size} features: Mean MCC = {mean_mcc:.4f}")
'''

import numpy as np
import time
import scanpy as sc
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, matthews_corrcoef, confusion_matrix, classification_report

# Define benchmark problems
benchmark_problems = [
    ("T_BRCA", "T_BLCA"),
    ("T_SKCM", "T_UVM"),
    ("T_LUAD", "T_LUSC"),
    ("T_STAD", "T_PAAD"),
    ("T_GBM", "T_LGG"),
    ("T_LIHC", "T_CHOL"),
    ("T_KIRP", "T_KIRC"),
    ("T_UCEC", "T_UCS"),
    ("T_CESC", "T_ESCA"),
    ("T_THYM", "T_HNSC"),
]

# Signature sizes to test
signature_sizes = [1, 3, 10, 20]
num_runs = 100  # Run each signature size 10 times

def run_benchmark(adata, method_n="wilcoxon", method="tree"):
    results = {}

    for cancer_1, cancer_2 in benchmark_problems:
        print(f"\nRunning benchmark for: {cancer_1} vs {cancer_2}")
        
        # Create new labels
        adata.obs['new_labels'] = adata.obs['labels'].replace({cancer_1: 'Cluster1', cancer_2: 'Cluster2'})

        # Differential expression analysis
        start_time = time.time()
        sc.tl.rank_genes_groups(adata, groupby="new_labels", groups=["Cluster1"], reference="Cluster2", method=method_n)
        end_time = time.time()
        print(f"Time taken for gene ranking analysis: {end_time - start_time:.2f} seconds")

        # Extract top genes and evaluate at different signature sizes
        ranked_genes = adata.uns['rank_genes_groups']['names']['Cluster1']
        results[f"{cancer_1}_vs_{cancer_2}"] = {}

        for size in signature_sizes:
            top_genes = ranked_genes[:size]
            mcc_scores = []

            for _ in range(num_runs):  # Run each size 10 times
                # Extract expression data
                cluster1_cells = adata.obs.index[adata.obs['new_labels'] == 'Cluster1']
                cluster2_cells = adata.obs.index[adata.obs['new_labels'] == 'Cluster2']

                cluster1_data = adata[cluster1_cells, top_genes].X.toarray() if hasattr(adata[cluster1_cells, top_genes].X, "toarray") else adata[cluster1_cells, top_genes].X
                cluster2_data = adata[cluster2_cells, top_genes].X.toarray() if hasattr(adata[cluster2_cells, top_genes].X, "toarray") else adata[cluster2_cells, top_genes].X

                # Combine data and create labels
                X = np.vstack((cluster1_data, cluster2_data))
                y = np.concatenate((np.ones(cluster1_data.shape[0]), np.zeros(cluster2_data.shape[0])))

                # Train-test split
                X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, stratify=y, random_state=np.random.randint(10000))

                match method:
                    case "logistic":
                        clf = LogisticRegression(max_iter=1000, class_weight='balanced')
                    case "boosting":
                        clf = HistGradientBoostingClassifier()
                    case "svm":
                        clf = LinearSVC(max_iter=1000, class_weight='balanced')
                    case "forest":
                        clf = RandomForestClassifier()
                    case "tree":
                        clf = DecisionTreeClassifier()

                clf.fit(X_train, y_train)

                # Evaluate classifier
                y_pred = clf.predict(X_test)
                mcc = matthews_corrcoef(y_test, y_pred)
                mcc_scores.append(mcc)  # Store MCC for this run

            # Compute mean MCC after 10 runs
            mean_mcc = np.mean(mcc_scores)

            results[f"{cancer_1}_vs_{cancer_2}"][f"{size}_features"] = {
                'mean_mcc': mean_mcc,
                'mcc_scores': mcc_scores  # Store all MCC scores for reference
            }

            print(f"\n🔹 Results for {size} features (Mean over {num_runs} runs):")
            print(f"Mean MCC Score: {mean_mcc:.4f}")
    
    return results

# Run the benchmark
final_results = run_benchmark(adata, method_n="wilcoxon", method="tree")
