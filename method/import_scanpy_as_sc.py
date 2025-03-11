import scanpy as sc
import anndata
import pandas as pd
import time
import numpy as np
from scipy.stats import wasserstein_distance
import os

'''
# Load the data
filepath="/home/nisma/new_yomix/yomix"
results_file = "/home/nisma/new_yomix/yomix/meth_kid.h5ad"
adata = anndata.read_h5ad(results_file)
print(adata)
# Preprocessing: Normalize and log-transform data
sc.pp.normalize_total(adata, target_sum=1e4)  # Normalize to 10,000 reads per cell
sc.pp.log1p(adata)  # Log transformation (natural log)

# Perform PCA
print("Computing PCA...")
sc.pp.pca(adata)

# Compute UMAP
print("Computing neighborhood graph for UMAP...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)  # Adjust parameters as needed
print("Computing UMAP...")
sc.tl.umap(adata)

# Print results to verify
print("PCA and UMAP computed successfully!")
print(adata)


output_dir = os.path.dirname(filepath)  # Get the directory of the input file
output_file = os.path.join(output_dir, "meth_kid_processed.h5ad")  # Create output path
adata.write(output_file)
print(f"Processed file saved as: {output_file}")
'''
# Load the data
filepath="/home/nisma/new_yomix/yomix"
results_file = "/home/nisma/new_yomix/yomix/aldinger20.processed.h5ad"
adata = anndata.read_h5ad(results_file)
print(adata)