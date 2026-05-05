library(Matrix)
library(Seurat)
library(future)
library(dplyr)

options(future.globals.maxSize=10 * 1024^3) # Increase max allowed size (10 GB) for SCTransform

# Load data (Download from GEO Database, Accession No. GSE175634)
data_dir <- "data/cardiomyocyte_differentiation/"
results_dir <- "results/cardiomyocyte_differentiation/"

count_data <- readMM(gzfile(paste0(data_dir, "GSE175634_cell_counts.mtx.gz")))
barcodes <- read.delim(paste0(data_dir, "GSE175634_cell_indices.tsv.gz"), header=TRUE)
features <- read.delim(paste0(data_dir, "GSE175634_gene_indices_counts.tsv.gz"), header=TRUE)

colnames(count_data) <- barcodes$cell_name
rownames(count_data) <- features$gene_name


metadata <- read.delim(paste0(data_dir, "GSE175634_cell_metadata.tsv.gz"), header=TRUE)
metadata$diffday <- as.factor(as.numeric(gsub("day", "", metadata$diffday)))


# Quality Control described in paper has already been performed!


### Paper Pipeline ###
# see github:
# https://github.com/jmp448/sc-dynamic-eqtl/blob/master/code/
#
# files:
# - "normalize_seurat.R"
# - "process_seurat.R"
#
# paper describes switching from seurat to scanpy in the process, but clustering
# is also possible with seurat alone


seurat_obj <- CreateSeuratObject(counts=count_data, project="cardiomyocyte_differentiation")
seurat_obj$diffday <- metadata$diffday
seurat_obj$article_celltype <- metadata$type

# cell cycle regression
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- CellCycleScoring(
  seurat_obj,
  s.features=cc.genes.updated.2019$s.genes,
  g2m.features=cc.genes.updated.2019$g2m.genes,
  set.ident=TRUE
)
seurat_obj$cc_difference <- seurat_obj$S.Score - seurat_obj$G2M.Score

seurat_obj <- SCTransform(seurat_obj, vars.to.regress="cc_difference", variable.features.n=1000)

# clustering
seurat_obj <- RunPCA(seurat_obj, npcs=100, features=VariableFeatures(seurat_obj)[1:1000])
seurat_obj <- FindNeighbors(seurat_obj, dims=1:50, nn.eps=0.1)

seurat_obj <- FindClusters(seurat_obj, resolution=0.15, n.start=10)
seurat_obj$low_res_clusters <- seurat_obj$seurat_clusters
seurat_obj <- FindClusters(seurat_obj, algorithm=4, resolution=0.35) # this is mentioned in the article, but gives other results
seurat_obj$med_res_clusters <- seurat_obj$seurat_clusters
seurat_obj <- FindClusters(seurat_obj, resolution=0.65, n.start=10)
seurat_obj$high_res_clusters <- seurat_obj$seurat_clusters



# visualization
seurat_obj <- RunUMAP(seurat_obj, dims=1:50, min.dist=0.15, n.neighbors=30)
# p <- {(DimPlot(seurat_obj, group.by="article_celltype") | DimPlot(seurat_obj, group.by="diffday")) /
#   (DimPlot(seurat_obj, group.by="low_res_clusters") | DimPlot(seurat_obj, group.by="med_res_clusters")) /
#   (DimPlot(seurat_obj, group.by="high_res_clusters") | ggplot())}
# ggsave(paste0(data_dir, "umap_different_clusterings.jpg"), p, width=6000, height=7500, units="px")


# test for markers
Idents(seurat_obj) <- seurat_obj$low_res_clusters # use low res louvain clusters (closest to paper results)
seurat_cluster_ids <- levels(Idents(seurat_obj))
df_list <- list()
for(i in seq_along(seurat_cluster_ids)){
  for(j in seq_along(seurat_cluster_ids)){

    test_id1 <- seurat_cluster_ids[i]

    if(j < i){
      next()
    } else if(i == j){
      test_id2 <- NULL
    } else{
      test_id2 <- seurat_cluster_ids[j]
    }

    markers <- FindMarkers(
      seurat_obj, ident.1=test_id1, ident.2=test_id2, test.use="wilcox",
      logfc.threshold=0.0, min.pct=0. # no pre-selection of genes
    )

    markers$cluster1 <- test_id1
    markers$cluster2 <- if(is.null(test_id2)) NA else test_id2
    markers$marker <- rownames(markers)
    rownames(markers) <- NULL

    df_list[[length(df_list)+1]] <- markers
  }
}

double_dipping_markers <- bind_rows(df_list)
umap <- Embeddings(seurat_obj, reduction="umap")
idents <- seurat_obj@meta.data[ , c("diffday", "article_celltype", "low_res_clusters", "med_res_clusters", "high_res_clusters")]

ddip_results <- list(
  idents=idents,
  markers=double_dipping_markers,
  umap=umap
)
save(ddip_results, file=paste0(results_dir, "ddip_results.Rdata"))
