library(Seurat)

# define function which will find N clusters by iteratively adapting resolution
FindNClusters <- function(obj, n_clusters=NULL, initial_resolution=0.8, max_iter=100, resolution_lower_bound=0.0, resolution_upper_bound=Inf){
  if(is.null(n_clusters)){
    return(FindClusters(obj, verbose=FALSE))
  }
  
  resolution <- initial_resolution
  for(i in 1:max_iter){
    temp_seurat_obj <- FindClusters(obj, resolution=resolution, verbose=FALSE)
    temp_idents <- Idents(temp_seurat_obj)
    n_idents <- length(unique(temp_idents))
    
    if(n_idents == n_clusters){
      break
    } else if(n_idents < n_clusters){
      # increase resolution
      resolution_lower_bound <- resolution
      resolution <- ifelse(
        is.infinite(resolution_upper_bound),
        resolution * 2,
        (resolution + resolution_upper_bound) / 2
      )
    } else {
      # decrease resolution
      resolution_upper_bound <- resolution
      resolution <- resolution <- ifelse(
        resolution_lower_bound == 0.,
        resolution / 2,
        (resolution + resolution_lower_bound) / 2
      )
    }
  }
  
  if(i == max_iter){
    warn(paste0(
      "Could not find ", n_clusters, " in ", max_iter, " iterations.\n",
      "Current resolution = ", resolution, " with bounds (", resolution_lower_bound, ", ", resolution_upper_bound, ").\n",
      "Current number of clusters is ", n_idents, "."
    ))
  }
  
  return(temp_seurat_obj)
}

cluster_fn <- function(data, subsample_indices, cluster_args=NULL){
  if(is.null(cluster_args)){
    cluster_args <- list()
  }
  
  seurat_obj <- CreateSeuratObject(counts=data)
  seurat_obj <- NormalizeData(seurat_obj, verbose=FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose=FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose=FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose=FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, verbose=FALSE)
  seurat_obj <- FindNClusters(seurat_obj, n_clusters=cluster_args$n_clusters)
  identities <- Idents(seurat_obj)
  return(identities)
}