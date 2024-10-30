# Define_Project_Hypothesis_and_Aim

## Validate Data

1. You can run the validate data file first to understand the kind of information this dataset will provide

```text
Files in directory:
alzheimer_single_cell_data.h5ad  checkdata.ipynb  sample_data

==================================================
DETAILED H5AD FILE INSPECTION
==================================================

1. FILE CHECK
--------------------
Loading data...
✓ File loaded successfully

2. BASIC INFORMATION
--------------------
• Number of cells (observations): 23197
• Number of genes (variables): 33091
• Size in memory: 185.39 MB

3. EXPRESSION MATRIX (adata.X) INFORMATION
--------------------
• Type: <class 'scipy.sparse._csr.csr_matrix'>
• Shape: (23197, 33091)
• Data statistics:
  - Min value: 0.2920
  - Max value: 7.8774
  - Mean value: 1.4397
  - Median value: 1.3271
  - Sparsity: 0.00%

4. CELL METADATA (adata.obs)
--------------------
Available columns:

• nCount_RNA:
  - Data type: float64
  - Number of unique values: 8705
  - Missing values: 0

• nFeature_RNA:
  - Data type: int32
  - Number of unique values: 4172
  - Missing values: 0

• percent.mt:
  - Data type: float64
  - Number of unique values: 22064
  - Missing values: 0

• SORT:
  - Data type: category
  - Number of unique values: 3
  - Unique values: ['AT8', 'MAP2', 'MAP2control']
  - Missing values: 0

• Amyloid:
  - Data type: category
  - Number of unique values: 3
  - Unique values: ['C3', 'DP (C0)', 'No']
  - Missing values: 0

• Age:
  - Data type: category
  - Number of unique values: 14
  - Missing values: 0

• RIN:
  - Data type: category
  - Number of unique values: 10
  - Missing values: 0

• nCount_SCT:
  - Data type: float64
  - Number of unique values: 4470
  - Missing values: 0

• nFeature_SCT:
  - Data type: int32
  - Number of unique values: 3513
  - Missing values: 0

• nCount_Exon:
  - Data type: float64
  - Number of unique values: 4665
  - Missing values: 0

• nFeature_Exon:
  - Data type: int32
  - Number of unique values: 2831
  - Missing values: 0

• PMI:
  - Data type: category
  - Number of unique values: 15
  - Missing values: 0

• Braak:
  - Data type: category
  - Number of unique values: 4
  - Unique values: ['0', 'I', 'II', 'VI']
  - Missing values: 0

• Sample.ID:
  - Data type: category
  - Number of unique values: 24
  - Missing values: 0

• Cell.Types:
  - Data type: category
  - Number of unique values: 7
  - Unique values: ['In1_LHX6-PVALB', 'In2_LHX6-PVALB-Chandelier', 'In3_LHX6-SST', 'In4_LHX6-SST-NPY', 'In5_LHX6-ADARB2-LAMP5', 'In6_ADARB2-LAMP5', 'In7_ADARB2-CALB2']
  - Missing values: 0

• tissue_ontology_term_id:
  - Data type: category
  - Number of unique values: 1
  - Unique values: ['UBERON:0000451']
  - Missing values: 0

• assay_ontology_term_id:
  - Data type: category
  - Number of unique values: 2
  - Unique values: ['EFO:0009899', 'EFO:0009922']
  - Missing values: 0

• disease_ontology_term_id:
  - Data type: category
  - Number of unique values: 2
  - Unique values: ['MONDO:0004975', 'PATO:0000461']
  - Missing values: 0

• cell_type_ontology_term_id:
  - Data type: category
  - Number of unique values: 1
  - Unique values: ['CL:0000498']
  - Missing values: 0

• development_stage_ontology_term_id:
  - Data type: category
  - Number of unique values: 14
  - Missing values: 0

• self_reported_ethnicity_ontology_term_id:
  - Data type: category
  - Number of unique values: 4
  - Unique values: ['HANCESTRO:0005', 'HANCESTRO:0014', 'HANCESTRO:0568', 'unknown']
  - Missing values: 0

• sex_ontology_term_id:
  - Data type: category
  - Number of unique values: 2
  - Unique values: ['PATO:0000383', 'PATO:0000384']
  - Missing values: 0

• is_primary_data:
  - Data type: bool
  - Number of unique values: 1
  - Unique values: [True]
  - Missing values: 0

• organism_ontology_term_id:
  - Data type: category
  - Number of unique values: 1
  - Unique values: ['NCBITaxon:9606']
  - Missing values: 0

• donor_id:
  - Data type: category
  - Number of unique values: 16
  - Missing values: 0

• suspension_type:
  - Data type: category
  - Number of unique values: 1
  - Unique values: ['cell']
  - Missing values: 0

• tissue_type:
  - Data type: category
  - Number of unique values: 1
  - Unique values: ['tissue']
  - Missing values: 0

• cell_type:
  - Data type: category
  - Number of unique values: 1
  - Unique values: ['inhibitory interneuron']
  - Missing values: 0

• assay:
  - Data type: category
  - Number of unique values: 2
  - Unique values: ["10x 3' v2", "10x 3' v3"]
  - Missing values: 0

• disease:
  - Data type: category
  - Number of unique values: 2
  - Unique values: ['Alzheimer disease', 'normal']
  - Missing values: 0

• organism:
  - Data type: category
  - Number of unique values: 1
  - Unique values: ['Homo sapiens']
  - Missing values: 0

• sex:
  - Data type: category
  - Number of unique values: 2
  - Unique values: ['female', 'male']
  - Missing values: 0

• tissue:
  - Data type: category
  - Number of unique values: 1
  - Unique values: ['prefrontal cortex']
  - Missing values: 0

• self_reported_ethnicity:
  - Data type: category
  - Number of unique values: 4
  - Unique values: ['African American', 'European', 'Hispanic or Latin American', 'unknown']
  - Missing values: 0

• development_stage:
  - Data type: category
  - Number of unique values: 14
  - Missing values: 0

• observation_joinid:
  - Data type: object
  - Number of unique values: 23197
  - Missing values: 0

5. GENE METADATA (adata.var)
--------------------
Available columns:

• feature_is_filtered:
  - Data type: bool
  - Number of unique values: 1
  - Missing values: 0

• feature_name:
  - Data type: category
  - Number of unique values: 33091
  - Missing values: 0

• feature_reference:
  - Data type: category
  - Number of unique values: 1
  - Missing values: 0

• feature_biotype:
  - Data type: category
  - Number of unique values: 1
  - Missing values: 0

• feature_length:
  - Data type: category
  - Number of unique values: 5063
  - Missing values: 0

• feature_type:
  - Data type: category
  - Number of unique values: 22
  - Missing values: 0

6. EXISTING ANALYSES
--------------------
Available items in uns:
• citation
• schema_reference
• schema_version
• title

Available items in obsm:
• X_pca
• X_umap

7. SAMPLE GENE NAMES
--------------------
First 10 genes:
['ENSG00000278915', 'ENSG00000168454', 'ENSG00000139180', 'ENSG00000229177', 'ENSG00000204564', 'ENSG00000116717', 'ENSG00000254418', 'ENSG00000114654', 'ENSG00000257894', 'ENSG00000198398']
```

## 
