# Add this at the top of your script
from rpy2.robjects import r, pandas2ri
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import tempfile
import os

pandas2ri.activate()
base = importr("base")
utils = importr("utils")
seurat = importr("Seurat")

def run_seurat_markers(adata, groupby, group_1, group_2):
    
    # Save AnnData to temporary file
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        temp_path = tmp.name
        adata.write(temp_path)

    # R script to read and analyze
    r_code = f"""
    library(Seurat)
    library(SeuratDisk)
    Convert("{temp_path}", dest = "h5seurat", overwrite = TRUE)
    seurat_obj <- LoadH5Seurat("{temp_path.replace('.h5ad', '.h5seurat')}")
    seurat_obj$binary_labels <- seurat_obj${groupby}
    markers <- FindMarkers(seurat_obj, ident.1 = "{group_1}", ident.2 = "{group_2}", group.by = "binary_labels", test.use = "wilcox")
    markers$gene <- rownames(markers)
    markers
    """
    result_df = r(r_code)
    os.remove(temp_path)
    return pandas2ri.rpy2py(result_df)
