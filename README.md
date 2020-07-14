# abcCellmap

#### R package for the prediction of blood cells at single-cell level

abcCellmap is a R package to predict the blood cell types for large amounts of data, which extends the function of cell map from ABC website(http://scrna.sklehabc.com/). 

Users can achieve the prediction by using two approaches (Scmap and Seurat). Single cells in our ABC reference are labeled by 43 different RNA clusters according to unsupervised clustering of single-cell transcriptional profiles, and also labeled by 32 immunophenotypic cell types, involving HSPC, B cell, T cell, NK cell, Neutrophil, Monocyte and Erythrocyte population.

### Installation

```R
#This require devtools  
install.packages('devtools')
library(devtools)
# abcCellmap requires dplyr/SingleCellExperiment/SummarizedExperiment/Seurat/scmap/reshape2
install_github("pzhulab/abcCellmap")
```

### Quick Start

```R
library(abcCellmap)
# query.exp is an example of expression matrix from ABC project(Atlas of Blood Cells).
> data(query.exp)
> dim(query.exp)
[1] 500   300
# 500 genes expression value of 300 cells

# abcCellmap has one parameter: queryData, a tab-delimited txt format file containing the expression matrix. 
# Each row should be a gene symbol and each column should be a single cell. 
# We prefer the expression data is presented by UMIs (Unique Molecular Identifiers) per gene in each single cell. 
> mapresult <- abcCellmap(queryData = query.exp)

```
 The return of abcCellmap() is a dataframe which contains 12 columns.
*  queryCell: the cell information in the query data
*  Seurat.RNACluster:  the RNA cluster predicted by Seurat
*  Seurat.RNACluster.score:  the prediction score of RNA cluster by Seurat
*  Seurat.Immunophenotype:  the immunophenotypic cell type predicted by Seurat
*  Seurat.Immunophenotype.score:  the prediction score of immunophenotypic cell type by Seurat
*  scmap.RNACluster:  the RNA cluster predicted by scmap
*  scmap.RNACluster.score:  the prediction score of RNA cluster by scmap
*  scmap.Immunophenotype:  the immunophenotypic cell type predicted by scmap
*  scmap.Immunophenotype.score:  the prediction score of immunophenotypic cell type by scmap
*  scmap.Cell:  the nearest single cell in our ABC reference predicted by scmap
*  scmap.Cell.score:  the prediction score of the nearest cell by scmap
*  pertype:  the percentages of top 2 immunophenotypic cell types in corresponding Seurat.RNACluster result

