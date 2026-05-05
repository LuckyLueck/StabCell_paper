library(Matrix)
library(Seurat)
library(dplyr)
library(future)
library(doFuture)
library(RhpcBLASctl)
library(countsplit)

set.seed(18032026)

options(future.globals.maxSize=50 * 1024^3) # Increase max allowed size (20 GB) for NormalizeData & SCTransform
blas_set_num_threads(64)
plan(multicore, workers=4)
# plan(multisession, workers=50)
registerDoFuture()

data_dir <- "data/cardiomyocyte_differentiation/"
results_dir <- "results/cardiomyocyte_differentiation/"

# Load data (Download from GEO Database, Accession No. GSE175634)
count_data <- readMM(gzfile(paste0(data_dir, "GSE175634_cell_counts.mtx.gz")))
barcodes <- read.delim(paste0(data_dir, "GSE175634_cell_indices.tsv.gz"), header=TRUE)
features <- read.delim(paste0(data_dir, "GSE175634_gene_indices_counts.tsv.gz"), header=TRUE)

count_data <- as(count_data, "CsparseMatrix")
colnames(count_data) <- barcodes$cell_name
rownames(count_data) <- features$gene_name

metadata <- read.delim(paste0(data_dir, "GSE175634_cell_metadata.tsv.gz"), header=TRUE)
metadata$diffday <- as.factor(as.numeric(gsub("day", "", metadata$diffday)))


# split data for clustering and marker detection
# split <- countsplit( ################# for testing purposes
#   count_data[
#     unique(c(
#       1:1000,
#       match(c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes), rownames(count_data))
#     )),
#     1:10000
#   ],
#   epsilon=c(0.5, 0.5)
# )
split <- countsplit(count_data, epsilon=c(0.5, 0.5))
data_clust <- split[[1]]
data_marker <- split[[2]]


# Create seurat object with both splits as different assays
seurat_cluster <- CreateSeuratObject(
  counts=data_clust,
  assay="CLUSTER",
  project="countsplit"
)
seurat_cluster$diffday <- metadata$diffday
seurat_cluster$article_celltype <- metadata$type


# cell cycle regression
DefaultAssay(seurat_cluster) <- "CLUSTER"
seurat_cluster <- NormalizeData(seurat_cluster)
seurat_cluster <- CellCycleScoring(
  seurat_cluster,
  s.features=cc.genes.updated.2019$s.genes,
  g2m.features=cc.genes.updated.2019$g2m.genes
)
seurat_cluster$cc_difference <- seurat_cluster$S.Score - seurat_cluster$G2M.Score
seurat_cluster$s_score <- seurat_cluster$S.Score
seurat_cluster$g2m_score <- seurat_cluster$G2M.Score
seurat_cluster$phase <- seurat_cluster$Phase
seurat_cluster <- SCTransform(
  seurat_cluster,
  vars.to.regress="cc_difference",
  assay="CLUSTER",
  variable.features.n=1000
)


# clustering
seurat_cluster <- RunPCA(seurat_cluster, npcs=100, features=VariableFeatures(seurat_cluster))
seurat_cluster <- FindNeighbors(seurat_cluster, dims=1:50, nn.eps=0.1)
# seurat_cluster <- FindClusters(seurat_cluster, resolution=0.15, n.start=10)
# seurat_cluster$countsplit_idents_low <- Idents(seurat_cluster)
# seurat_cluster <- FindClusters(seurat_cluster, resolution=0.35, n.start=10)
# seurat_cluster$countsplit_idents_med <- Idents(seurat_cluster)
# seurat_cluster <- FindClusters(seurat_cluster, resolution=0.65, n.start=10)
# seurat_cluster$countsplit_idents_high <- Idents(seurat_cluster)
seurat_cluster <- FindClusters(seurat_cluster, resolution=0.2, n.start=10)
seurat_cluster$countsplit_idents <- Idents(seurat_cluster)

countsplit_idents <- Idents(seurat_cluster)

# visualization
seurat_cluster <- RunUMAP(seurat_cluster, dims=1:50, min.dist=0.15, n.neighbors=30)
# ggplot2::ggsave(
#   filename=paste0(data_dir, "countsplit_clusters.jpg"),
#   plot=DimPlot(seurat_cluster, group.by="countsplit_idents"),
#   width=4000, height=4000, units="px",
# )


# prepare marker detection
seurat_marker <- CreateSeuratObject(
  counts=data_marker,
  assay="MARKER",
  project="count_split"
)
seurat_marker$diffday <- metadata$diffday
seurat_marker$article_celltype <- metadata$type
seurat_marker$countsplit_idents <- countsplit_idents
Idents(seurat_marker) <- "countsplit_idents"
seurat_marker <- NormalizeData(seurat_marker)


# DE analysis
seurat_cluster_ids <- levels(Idents(seurat_marker))
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
      seurat_marker, ident.1=test_id1, ident.2=test_id2, test.use="wilcox",
      logfc.threshold=0.0, min.pct=0., # no pre-selection of genes
      assay="MARKER", slot="data"
    )

    markers$cluster1 <- test_id1
    markers$cluster2 <- if(is.null(test_id2)) NA else test_id2
    markers$marker <- rownames(markers)
    rownames(markers) <- NULL

    df_list[[length(df_list)+1]] <- markers
  }
}


# save results
countsplit_markers <- bind_rows(df_list)
umap <- Embeddings(seurat_cluster, reduction="umap")

countsplit_results <- list(
  idents=countsplit_idents,
  markers=countsplit_markers,
  umap=umap
)
save(countsplit_results, file=paste0(results_dir, "countsplit_results.Rdata"))
