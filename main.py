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

import scanpy as sc
import anndata
import time
from method.scanpymine import scanpymine
# Load the data
results_file = "/home/nisma/write/xd_tcga_labels_umap.h5ad"
adata = anndata.read_h5ad(results_file)

