import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse
import os

def comprehensive_alzheimer_analysis(adata, output_dir="analysis_results"):
    """
    Comprehensive analysis pipeline combining detailed visualization with clear explanations
    """
    print("\n=== STEP 1: Dataset Overview and Quality Control ===")
    print("This dataset contains:")
    print(f"• {adata.n_obs} cells")
    print(f"• {adata.n_vars} genes")
    
    # Store raw counts
    adata.raw = adata
    
    # Cell type distribution
    print("\nCell Types in the dataset:")
    for cell_type in adata.obs['Cell.Types'].unique():
        count = sum(adata.obs['Cell.Types'] == cell_type)
        percentage = (count/len(adata.obs))*100
        print(f"• {cell_type}: {count} cells ({percentage:.1f}%)")

    # Quality metrics visualization
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # RNA metrics violin plots
    sc.pl.violin(adata, ['nCount_RNA', 'nFeature_RNA', 'percent.mt'],
                groupby='disease', rotation=45, ax=axes[0, 0])
    axes[0, 0].set_title('Quality Metrics by Disease Status')
    
    # Cell type composition pie chart
    cell_counts = adata.obs['Cell.Types'].value_counts()
    axes[0, 1].pie(cell_counts, labels=cell_counts.index, autopct='%1.1f%%')
    axes[0, 1].set_title('Cell Type Distribution')
    
    # Disease state distribution
    disease_counts = adata.obs['disease'].value_counts()
    axes[1, 0].bar(disease_counts.index, disease_counts.values)
    axes[1, 0].set_title('Disease State Distribution')
    axes[1, 0].set_ylabel('Number of Cells')
    
    # Braak stage distribution
    braak_counts = adata.obs['Braak'].value_counts()
    axes[1, 1].bar(braak_counts.index, braak_counts.values)
    axes[1, 1].set_title('Braak Stage Distribution')
    axes[1, 1].set_ylabel('Number of Cells')
    
    plt.tight_layout()
    plt.show()

    print("\n=== STEP 2: Disease and Cell Type Analysis ===")
    # Cross-tabulation of cell types and disease status
    cross_tab = pd.crosstab(adata.obs['Cell.Types'], adata.obs['disease'])
    print("\nNumber of cells for each cell type in AD vs Normal:")
    print(cross_tab)
    
    # Stacked bar plot
    plt.figure(figsize=(12, 6))
    cross_tab.plot(kind='bar', stacked=True)
    plt.title("Cell Type Distribution by Disease Status")
    plt.xlabel("Cell Type")
    plt.ylabel("Number of Cells")
    plt.legend(title="Disease Status")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

    print("\n=== STEP 3: Braak Stage Analysis ===")
    print("\nDistribution of Braak stages:")
    for stage in braak_counts.index:
        percentage = (braak_counts[stage]/len(adata.obs))*100
        print(f"• Stage {stage}: {braak_counts[stage]} cells ({percentage:.1f}%)")
    
    # Analyze gene expression changes across Braak stages
    sc.tl.rank_genes_groups(adata, 'Braak', method='wilcoxon')
    print("\nTop differentially expressed genes across Braak stages:")
    sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)

    print("\n=== STEP 4: Cell Type-Specific Analysis ===")
    cell_types = adata.obs['Cell.Types'].unique()
    
    for cell_type in cell_types:
        print(f"\nAnalyzing {cell_type}")
        cell_type_data = adata[adata.obs['Cell.Types'] == cell_type]
        
        # Compare gene expression between AD and normal
        sc.tl.rank_genes_groups(cell_type_data, 'disease', 
                              groups=['Alzheimer disease'], 
                              reference='normal',
                              method='wilcoxon')
        
        print(f"\nTop differentially expressed genes in {cell_type}:")
        sc.pl.rank_genes_groups(cell_type_data, n_genes=10, sharey=False)

    print("\n=== STEP 5: Demographic Analysis ===")
    # Age distribution
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=adata.obs, x='disease', y='Age')
    plt.title('Age Distribution by Disease Status')
    plt.show()
    
    # Sex distribution
    sex_disease_cross = pd.crosstab(adata.obs['sex'], adata.obs['disease'])
    plt.figure(figsize=(8, 6))
    sex_disease_cross.plot(kind='bar')
    plt.title('Sex Distribution by Disease Status')
    plt.xlabel('Sex')
    plt.ylabel('Number of Cells')
    plt.legend(title='Disease Status')
    plt.tight_layout()
    plt.show()

    print("\n=== STEP 6: Dimensional Reduction Visualization ===")
    # Compute UMAP if not already computed
    if 'X_umap' not in adata.obsm_keys():
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pp.pca(adata, svd_solver='arpack')
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
    
    # Create UMAP visualizations
    fig, axes = plt.subplots(2, 2, figsize=(15, 15))
    axes = axes.flatten()
    
    sc.pl.umap(adata, color='disease', ax=axes[0], show=False, title='Disease Status')
    sc.pl.umap(adata, color='Cell.Types', ax=axes[1], show=False, title='Cell Types')
    sc.pl.umap(adata, color='Braak', ax=axes[2], show=False, title='Braak Stage')
    sc.pl.umap(adata, color='sex', ax=axes[3], show=False, title='Sex')
    
    plt.tight_layout()
    plt.show()

    # Save results if output directory is specified
    if output_dir:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Save the processed AnnData object
        adata.write(f"{output_dir}/processed_data.h5ad")
        
        # Save cell type markers
        if 'rank_genes_groups' in adata.uns:
            pd.DataFrame(adata.uns['rank_genes_groups']['names']).to_csv(
                f"{output_dir}/cell_type_markers.csv")
        
        print(f"\nResults saved to {output_dir}/")
    
    return adata

# Example usage:
if __name__ == "__main__":
    print("Loading data...")
    adata = sc.read_h5ad("alzheimer_single_cell_data.h5ad")
    print("Running analysis...")
    adata = comprehensive_alzheimer_analysis(adata)