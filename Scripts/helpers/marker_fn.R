library(Seurat)

marker_fn <- function(data, identities, weights, cluster_id1, cluster_id2=NULL, marker_args=NULL){
  # cannot test cluster against itself
  # --> check for cluster collision
  if((!is.null(cluster_id2)) && (cluster_id1 == cluster_id2)){
    p_values <- rep(1., nrow(data))
    names(p_values) <- rownames(data)
    return(p_values)
  } else {
    # prepare seurat object
    seurat_obj <- CreateSeuratObject(counts=data)
    seurat_obj <- NormalizeData(seurat_obj, verbose=FALSE)
    Idents(seurat_obj) <- identities
    
    markers <- FindMarkers(
      seurat_obj, ident.1=cluster_id1, ident.2=cluster_id2, test.use="wilcox",
      logfc.threshold=0.0, min.pct=0., # no pre-selection of genes
      verbose=FALSE
    )
    
    # format results
    p_values <- markers$p_val_adj
    names(p_values) <- rownames(markers)
  }
  
  return(p_values)
}
