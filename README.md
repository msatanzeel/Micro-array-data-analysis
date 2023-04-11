# Differential expression analysis with limma (Main.r)
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


# Preprocessing of Data (Data preparation.ipynb)

The Jupyter Notebook is divided into multiple cells, each containing a different section of the code. The purpose of this notebook is to preprocess gene expression data and prepare it for machine learning analysis. Here is a brief description of what each cell does:

1. **Loading the data:** In this cell, the complete gene expression data is loaded into a Pandas DataFrame using the read_csv() function from the pandas library. The index column is set to 0, and the top columns are displayed using the head() function.

2. **Removing genes with missing values:** This cell calculates the threshold for the number of non-null values required and removes genes with less than the threshold number of non-null values. The copy() function is used to store a copy of the original dataset. The dimensions of the new dataset are displayed using the shape attribute.

3. **Loading Differentially Expressed Genes Id's:** In this cell, the IDs of differentially expressed genes are loaded from a CSV file using the read_csv() function from the pandas library. A new DataFrame is created containing only the 'ID' column. The to_frame() method is used to convert it to a DataFrame.

4. **Selecting the DEG's subset:** This cell creates a Boolean mask to select rows with matching IDs and then selects matching rows from all_data_copy using the loc[] accessor. The number of rows that match the IDs is displayed using the sum() function.

5. **Transposing the data and handling missing values:** In this cell, the rows and columns of the dataframe are transposed. The missing values in the dataset are handled using KNN imputation, and the column names are set as row names.

6. **Preparing a results dataframe:** In this cell, a string is created representing the status of the samples as either infected or control. A new DataFrame is created with a column called 'status', and the string is split into individual characters and assigned to the 'status' column.

7. **Exporting the processed data:** In this cell, the processed data is exported to CSV files for further analysis using the to_csv() function from the pandas library.

This Jupyter Notebook can be used as a template for preprocessing gene expression data and preparing it for machine learning analysis.
