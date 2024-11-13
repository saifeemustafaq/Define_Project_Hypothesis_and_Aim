# Detailed Code Explanation: Single-cell RNA Analysis Pipeline

## 1. Initial Setup and Imports

```python
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
```

This section imports all necessary libraries:
- `scanpy`: A library specifically designed for single-cell RNA analysis
- `seaborn` and `matplotlib`: For creating visualizations
- `pandas` and `numpy`: For data manipulation
- `scipy.stats`: For statistical analysis

## 2. Data Loading and Gene Definition

```python
adata = sc.read_h5ad("dataset.h5ad")

genes_of_interest = ['ENSG00000176697', 'ENSG00000183454', 'ENSG00000186868',
                    'ENSG00000136244', 'ENSG00000131095']
gene_name_map = {
    'ENSG00000176697': 'BDNF',
    'ENSG00000183454': 'GRIN2A',
    'ENSG00000186868': 'MAPT',
    'ENSG00000136244': 'IL6',
    'ENSG00000131095': 'GFAP'
}
```

Here, the code:
1. Loads the dataset from an H5AD file (a format specifically designed for single-cell data)
2. Defines a list of genes using their Ensembl IDs
3. Creates a mapping dictionary to convert Ensembl IDs to readable gene names

## 3. Data Filtering and Matrix Creation

```python
adata_subset = adata[:, [gene for gene in genes_of_interest if gene in adata.var_names]]

expression_matrix = pd.DataFrame(
    adata_subset.X.toarray(),
    columns=[gene_name_map.get(g, g) for g in adata_subset.var_names],
    index=adata_subset.obs.index
)
```

This section:
1. Filters the dataset to include only the genes of interest
2. Creates a pandas DataFrame containing the expression values
3. Converts gene IDs to their readable names using the mapping
4. Preserves the original sample indices

## 4. Metadata Integration

```python
expression_matrix['disease'] = adata_subset.obs['disease'].values
expression_matrix['sex'] = adata_subset.obs['sex'].values
expression_matrix['brain_region'] = adata_subset.obs['Brain.Region'].values
```

The code adds three important metadata columns to the expression matrix:
- Disease status
- Sex
- Brain region
This allows for grouped analysis and comparison across these variables.

## 5. Visualization: Gene Expression by Disease

```python
fig, axes = plt.subplots(2, 3, figsize=(20, 12))
axes = axes.flatten()

for idx, (gene_id, gene_name) in enumerate(gene_name_map.items()):
    if idx < len(axes):
        sns.violinplot(data=expression_matrix, x='disease', y=gene_name,
                      ax=axes[idx], inner='box', palette='Set2')
        axes[idx].set_title(f'{gene_name} Expression by Disease')
        axes[idx].set_xlabel('')
        axes[idx].tick_params(axis='x', rotation=45)
```

This visualization section:
1. Creates a 2x3 grid of subplots
2. Generates violin plots for each gene showing:
   - Distribution of expression values across disease groups
   - Box plots inside the violins showing quartiles
   - Color-coded by disease group
3. Customizes the appearance with rotated labels and titles

## 6. Statistical Analysis

```python
print("\nSummary Statistics:")
for gene_name in gene_name_map.values():
    print(f"\n{gene_name} statistics by disease:")
    summary = expression_matrix.groupby('disease')[gene_name].agg(['mean', 'std', 'count'])
```

The statistical analysis includes:
1. Summary statistics for each gene:
   - Mean expression
   - Standard deviation
   - Sample count
   - Grouped by disease status

## 7. Correlation Analysis

```python
correlation_matrix = expression_matrix[list(gene_name_map.values())].corr()
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm',
            center=0, fmt='.2f', square=True)
```

This section:
1. Calculates correlations between all genes
2. Creates a heatmap showing:
   - Correlation strength (color-coded)
   - Numerical correlation values
   - Symmetrical display

## 8. Statistical Testing

```python
for gene_name in gene_name_map.values():
    groups = [group[gene_name].values for name, group in expression_matrix.groupby('disease')]
    h_stat, p_val = stats.kruskal(*groups)
    eta_squared = (h_stat - len(groups) + 1) / (n - len(groups))
```

The statistical testing:
1. Performs Kruskal-Wallis tests for each gene
2. Calculates:
   - H-statistic (test statistic)
   - p-value (statistical significance)
   - Effect size (η²) to measure the strength of the relationship

## 9. Brain Region Analysis

```python
for idx, (gene_id, gene_name) in enumerate(gene_name_map.items()):
    sns.violinplot(data=expression_matrix, x='brain_region', y=gene_name,
                  ax=axes[idx], inner='box', palette='Set3')
```

The brain region analysis:
1. Creates violin plots similar to the disease analysis
2. Shows expression patterns across different brain regions
3. Calculates summary statistics for each brain region

## Key Features of the Code

1. **Modularity**: Each analysis component is separate and can be modified independently

2. **Error Handling**:
   - Checks for gene presence in the dataset
   - Uses get() method for safe dictionary access
   - Handles missing data appropriately

3. **Visualization Best Practices**:
   - Consistent styling
   - Clear labeling
   - Appropriate color schemes
   - Efficient use of space

4. **Statistical Rigor**:
   - Multiple levels of analysis
   - Appropriate statistical tests
   - Effect size calculations
   - Comprehensive summary statistics

## Data Flow

1. Raw Data → Filtered Dataset → Expression Matrix
2. Expression Matrix + Metadata → Analysis Components
3. Analysis Results → Visualizations and Statistics

This structure allows for:
- Easy modification of analysis parameters
- Addition of new analysis components
- Scalability to different datasets
- Clear data lineage tracking
