#' Human Blood Cell Landscape mapping
#'
#' Human Blood Cell Landscape R package for single cell mapping.
#'
#' @param queryData a space-delimited txt format file containing the expression matrix. Each row should be a gene symbol and each column should be a single cell. We prefer the expression data is presented by UMIs (Unique Molecular Identifiers) per gene in each single cell.
#'
#' @return Users can predict the cell types of hematopoietic cells by implementing two approaches (Scmap and Seurat). Cells in our ABC are labeled by 43 different RNA clusters according to unsupervised clustering of single-cell transcriptional profiles, and also labeled by 32 immunophenotypic cell types, involving HSPC, B cell, T cell, NK cell, Neutrophil, Monocyte and Erythrocyte population. The format of result are as follows:
#'
#' queryCell, Seurat.RNACluster, Seurat.RNACluster.score, Seurat.Immunophenotype, Seurat.Immunophenotype.score, scmap.RNACluster, scmap.RNACluster.score, scmap.Immunophenotype, scmap.Immunophenotype.score, scmap.Cell, scmap.Cell.score, pertype.
#'
#' Of which, "queryCell" is the cell information in the query data,
#'
#' "Seurat.RNACluster" is the RNA cluster predicted by Seurat,
#'
#' "Seurat.RNACluster.score" is the prediction score of RNA cluster by Seurat,
#'
#' "Seurat.Immunophenotype" is the immunophenotypic cell type predicted by Seurat,
#'
#' "Seurat.Immunophenotype.score" is the prediction score of immunophenotypic cell type by Seurat,
#'
#' "scmap.RNACluster" is the RNA cluster predicted by scmap,
#'
#' "scmap.RNACluster.score" is the prediction score of RNA cluster by scmap,
#'
#' "scmap.Immunophenotype" is the immunophenotypic cell type predicted by scmap,
#'
#' "scmap.Immunophenotype.score" is the prediction score of immunophenotypic cell type by scmap,
#'
#' "scmap.Cell" is the nearest single cell in our ABC reference predicted by scmap,
#'
#' "scmap.Cell.score" is the prediction score of the nearest cell by scmap,
#'
#' "pertype" means the percentages of top 2 immunophenotypic cell types in corresponding Seurat.RNACluster result.
#'
#' @import dplyr SingleCellExperiment SummarizedExperiment
#'
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures FindTransferAnchors TransferData
#'
#' @importFrom scmap scmapCluster scmapCell
#'
#' @importFrom reshape2 melt
#'
#' @export
#'
#' @examples mapresult<-abcCellmap(queryData=query.exp)

abcCellmap <- function(queryData){
  seurat_data <- as(as.matrix(queryData), "dgCMatrix")
  seuratdata <- CreateSeuratObject(counts=seurat_data, project="abc")
  seuratdata <- NormalizeData(seuratdata, verbose=FALSE)
  seuratdata <- FindVariableFeatures(seuratdata, selection.method = "vst", nfeatures = 5)
  seuratref_data <- as(as.matrix(abc.exp), "dgCMatrix")
  seuratrefdata <- CreateSeuratObject(counts=seuratref_data, project="abc")
  seuratrefdata <- NormalizeData(seuratrefdata, verbose=FALSE)
  seuratrefdata <- FindVariableFeatures(seuratrefdata, selection.method = "vst", nfeatures = 5)
  seuratrefdata@meta.data$cluster=seuratcluster$cluster
  seuratrefdata@meta.data$immunophenotype=seuratcluster$immunophenotype
  seuratdata@assays$RNA@var.features=seuratmarker
  seuratrefdata@assays$RNA@var.features=seuratmarker
  rna_seuratanchors <- FindTransferAnchors(reference = seuratrefdata, query = seuratdata, dims = 1:30)
  rna_seuratpredictions <- TransferData(anchorset = rna_seuratanchors, refdata = seuratrefdata$cluster,dims = 1:30)
  rna_seuratresult<-data.frame(rownames(rna_seuratpredictions),rna_seuratpredictions[c("predicted.id","prediction.score.max")])
  names(rna_seuratresult)<-c("queryCell","Seurat.RNACluster","Seurat.RNACluster.score")
  rna_seuratresult$queryCell<-as.character(rna_seuratresult$queryCell)
  seuratdata@assays$RNA@var.features=seuratimmumarker
  seuratrefdata@assays$RNA@var.features=seuratimmumarker
  imm_seuratanchors <- FindTransferAnchors(reference = seuratrefdata, query = seuratdata, dims = 1:30)
  imm_seuratpredictions <- TransferData(anchorset = imm_seuratanchors, refdata = seuratrefdata$immunophenotype,dims = 1:30)
  imm_seuratresult<-data.frame(rownames(imm_seuratpredictions),imm_seuratpredictions[c("predicted.id","prediction.score.max")])
  names(imm_seuratresult)<-c("queryCell","Seurat.Immunophenotype","Seurat.Immunophenotype.score")
  imm_seuratresult$queryCell<-as.character(imm_seuratresult$queryCell)
  samp<-SingleCellExperiment(assays=list(normcounts=as.matrix(queryData)))
  logcounts(samp)<-log2(normcounts(samp)+1)
  rowData(samp)$feature_symbol<-rownames(samp)
  samp<-samp[!duplicated(rownames(samp)), ]
  rnascmap=scmapCluster(projection=samp, index_list=list(predicted_cluster=scmapcluindex), threshold=0.0)
  rnascmapresult=data.frame(colnames(samp), rnascmap$scmap_cluster_labs, rnascmap$scmap_cluster_siml)
  colnames(rnascmapresult)<-c("queryCell","scmap.RNACluster","scmap.RNACluster.score")
  rnascmapresult$queryCell<-as.character(rnascmapresult$queryCell)
  immscmap=scmapCluster(projection=samp, index_list=list(predicted_cluster=scmapimmindex), threshold=0.0)
  immscmapresult=data.frame(colnames(samp), immscmap$scmap_cluster_labs, immscmap$scmap_cluster_siml)
  colnames(immscmapresult)<-c("queryCell","scmap.Immunophenotype","scmap.Immunophenotype.score")
  immscmapresult$queryCell<-as.character(immscmapresult$queryCell)
  scmapCell_results=scmapCell(samp, list(predicted_cell=scmapcellindex))
  sc1<-melt(as.data.frame(scmapCell_results$predicted_cell$cells))
  sc2<-melt(as.data.frame(scmapCell_results$predicted_cell$similarities))
  scall<-cbind(sc1,sc2)
  scall<-scall[-3]
  colnames(scall)<-c("queryCell","CellID","scmap.Cell.score")
  scall <- scall %>% left_join(scmapcellID, by='CellID')
  scmapcellresult<-scall %>% group_by(queryCell) %>% filter(scmap.Cell.score == max(scmap.Cell.score)) %>% ungroup() %>% select(c(queryCell,scmap.Cell,scmap.Cell.score))
  scmapcellresult$queryCell<-as.character(scmapcellresult$queryCell)
  cluimmu$Seurat.RNACluster<-as.character(cluimmu$Seurat.RNACluster)
  mapresults<-rna_seuratresult %>% full_join(imm_seuratresult,by = "queryCell") %>% full_join(rnascmapresult,by = "queryCell") %>% full_join(immscmapresult,by = "queryCell") %>% full_join(scmapcellresult,by = "queryCell")
  mapresults<-mapresults%>%left_join(cluimmu,by="Seurat.RNACluster")
  return(mapresults)
}
