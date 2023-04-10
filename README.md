# Differential expression analysis with limma
This script performs differential expression analysis on a gene expression dataset using the limma package in R. The dataset is loaded from the Gene Expression Omnibus (GEO) database using the getGEO() function from the GEOquery R package. The dataset being loaded is identified by the accession code "GSE21962".

# Requirements
- GEOquery
- limma
- umap

# Steps
1. Load the packages (make sure to install them via console before loading them)
2. Load the dataset using getGEO() function
3. Modify feature variable labels of a GEOquery object gset
4. Handle missing values using the KNN impute
5. Perform log2 transformation on the filtered expression matrix
6. Normalize the data using quantile-based normalization
7. Assign samples to groups and set up design matrix
8. Fit linear model using lmFit() function
9. Set up contrasts of interest and recalculate model coefficients
10. Compute statistics and table of top significant genes
11. Calculate how many genes are there that satisfy the conditions of p.value and fold change value range
12. Copy significant DEGs to another table
13. Subset relevant columns from the DEGs table
14. Write the DEGs table to a csv file
15. Build histogram of P-values for all genes
16. Summarize test results as "up", "down" or "not expressed"
17. Create a Venn diagram of results
18. Create a Q-Q plot for t-statistic
19. Create a volcano plot (log P-value vs log fold change)
