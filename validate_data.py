# This piece verifies the data inside the h5ad file
# First install required packages
!pip -q install scanpy pandas numpy matplotlib seaborn scipy

# Import required libraries
import scanpy as sc
import pandas as pd
import numpy as np
import os
import scipy.sparse

def inspect_h5ad_file(file_path):
    """
    Comprehensive inspection of h5ad file structure and contents
    """
    print("\n" + "="*50)
    print("DETAILED H5AD FILE INSPECTION")
    print("="*50)

    # 1. Basic File Check
    print("\n1. FILE CHECK")
    print("-"*20)
    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' not found!")
        return

    # 2. Load Data
    try:
        print("Loading data...")
        adata = sc.read_h5ad(file_path)
        print("✓ File loaded successfully")
    except Exception as e:
        print(f"Error loading file: {str(e)}")
        return

    # 3. Basic Information
    print("\n2. BASIC INFORMATION")
    print("-"*20)
    print(f"• Number of cells (observations): {adata.n_obs}")
    print(f"• Number of genes (variables): {adata.n_vars}")
    print(f"• Size in memory: {adata.X.data.nbytes / 1024**2:.2f} MB" if scipy.sparse.issparse(adata.X)
          else f"• Size in memory: {adata.X.nbytes / 1024**2:.2f} MB")

    # 4. Expression Matrix Information
    print("\n3. EXPRESSION MATRIX (adata.X) INFORMATION")
    print("-"*20)
    print(f"• Type: {type(adata.X)}")
    print(f"• Shape: {adata.X.shape}")
    print("• Data statistics:")
    X_sample = adata.X.data if scipy.sparse.issparse(adata.X) else adata.X.flatten()
    print(f"  - Min value: {np.min(X_sample):.4f}")
    print(f"  - Max value: {np.max(X_sample):.4f}")
    print(f"  - Mean value: {np.mean(X_sample):.4f}")
    print(f"  - Median value: {np.median(X_sample):.4f}")
    print(f"  - Sparsity: {(X_sample == 0).sum() / len(X_sample):.2%}")

    # 5. Cell Metadata
    print("\n4. CELL METADATA (adata.obs)")
    print("-"*20)
    print("Available columns:")
    for col in adata.obs.columns:
        n_unique = adata.obs[col].nunique()
        print(f"\n• {col}:")
        print(f"  - Data type: {adata.obs[col].dtype}")
        print(f"  - Number of unique values: {n_unique}")
        if n_unique < 10:  # If few unique values, show them
            print(f"  - Unique values: {sorted(adata.obs[col].unique())}")
        print(f"  - Missing values: {adata.obs[col].isnull().sum()}")

    # 6. Gene Metadata
    print("\n5. GENE METADATA (adata.var)")
    print("-"*20)
    print("Available columns:")
    for col in adata.var.columns:
        print(f"\n• {col}:")
        print(f"  - Data type: {adata.var[col].dtype}")
        print(f"  - Number of unique values: {adata.var[col].nunique()}")
        print(f"  - Missing values: {adata.var[col].isnull().sum()}")

    # 7. Check for Existing Results
    print("\n6. EXISTING ANALYSES")
    print("-"*20)
    if adata.uns:
        print("Available items in uns:")
        for key in adata.uns.keys():
            print(f"• {key}")
    else:
        print("No existing analyses found in uns")

    if adata.obsm:
        print("\nAvailable items in obsm:")
        for key in adata.obsm.keys():
            print(f"• {key}")
    else:
        print("\nNo dimensional reductions found in obsm")

    # 8. Sample Gene Names
    print("\n7. SAMPLE GENE NAMES")
    print("-"*20)
    print("First 10 genes:")
    print(list(adata.var_names[:10]))

    return adata

# Run the inspection
print("Current working directory:", os.getcwd())
print("\nFiles in directory:")
!ls

file_path = "alzheimer_single_cell_data.h5ad"
adata = inspect_h5ad_file(file_path)

# Optional: Save inspection results to a text file
def save_inspection_results(adata, output_file="data_inspection_results.txt"):
    """
    Save key information about the dataset to a text file
    """
    import sys
    import datetime

    original_stdout = sys.stdout
    with open(output_file, 'w') as f:
        sys.stdout = f
        print(f"Data Inspection Results - {datetime.datetime.now()}")
        inspect_h5ad_file(file_path)
        sys.stdout = original_stdout

    print(f"\nInspection results saved to {output_file}")

# Uncomment to save results:
# save_inspection_results(adata)
