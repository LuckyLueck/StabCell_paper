library(Matrix)
library(Seurat)
library(dplyr)
library(RhpcBLASctl)
library(doFuture)
library(progressr)

library(StabCell)


options(future.globals.maxSize=10 * 1024^3)

set.seed(26082025)

data_dir <- "data/cardiomyocyte_differentiation/"
results_dir <- "results/cardiomyocyte_differentiation/"


# Load data (Download from GEO Database, Accession No. GSE175634)
cat("load count data...\n")

count_data <- readMM(gzfile(paste0(data_dir, "GSE175634_cell_counts.mtx.gz")))
barcodes <- read.delim(paste0(data_dir, "GSE175634_cell_indices.tsv.gz"), header=TRUE)
features <- read.delim(paste0(data_dir, "GSE175634_gene_indices_counts.tsv.gz"), header=TRUE)

colnames(count_data) <- barcodes$cell_name
rownames(count_data) <- features$gene_name

# Load Double Dipping idents as initial clustering
cat("load initial clustering...\n")
load(paste0(results_dir, "ddip_results.Rdata"))
initial_idents <- ddip_results$idents$low_res_clusters

# subsampling of dataset for speed test?
subsampling_factor <- 0.
if(subsampling_factor != 0.){
  cat("subsample data for speed tests...\n")
  n <- length(initial_idents)
  rnd_idx <- sample(
    1:n,
    as.integer(n / subsampling_factor),
    replace=FALSE
  )
} else {
  rnd_idx <- TRUE
}
initial_idents <- initial_idents[rnd_idx]
count_data <- count_data[ , rnd_idx]

# convert to column sparse matrix for seurat
count_data <- as(count_data, "CsparseMatrix")

cat("num cells:", ncol(count_data), "\n")
cat("num genes:", nrow(count_data), "\n")


# prepare StabCell
cat("prepare StabCell\n")

PFER_max <- c(1, 5, 10, 25, 50, 100)
steps <- NULL
p_value_threshold <- 0.05


results_file <- paste0(results_dir, "stabcell_results__subsampling_factor=", subsampling_factor, ".Rdata")
checkpoint_file <- NULL # paste0(results_dir, "stabcell_checkpoint__subsampling_factor=", subsampling_factor, ".Rdata")

# different worker setup for different clustering vs other steps
# (steps must be strictly increasing)
step_settings <- list(
  list(steps=1,    num_workers=4,  blas_threads=48),
  list(steps=2:10, num_workers=20,  blas_threads=1)
)

cluster_fn <- function(data, subsample_indices, cluster_args){
  seurat_obj <- CreateSeuratObject(counts=data)

  # normalization & cell cycle regression
  seurat_obj <- NormalizeData(seurat_obj, verbose=FALSE)
  seurat_obj <- CellCycleScoring(
    seurat_obj,
    s.features=cc.genes.updated.2019$s.genes,
    g2m.features=cc.genes.updated.2019$g2m.genes,
    set.ident=TRUE
  )
  seurat_obj$cc_difference <- seurat_obj$S.Score - seurat_obj$G2M.Score
  seurat_obj <- SCTransform(seurat_obj, vars.to.regress="cc_difference", variable.features.n=1000, verbose=FALSE)

  # actual clustering
  seurat_obj <- RunPCA(seurat_obj, npcs=100, features=VariableFeatures(seurat_obj)[1:1000], verbose=FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims=1:50, nn.eps=0.1, verbose=FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution=0.15, n.start=10, verbose=FALSE)

  identities <- Idents(seurat_obj)
  return(identities)
}

test_fn <- function(data, identities, weights, test_id1, test_id2=NULL, test_args=NULL){
  # cannot test cluster against itself
  # --> check for cluster collision
  if((!is.null(test_id2)) && (test_id1 == test_id2)){
    p_values <- rep(1., nrow(data))
    names(p_values) <- rownames(data)
    return(p_values)
  } else {
    # prepare seurat object
    seurat_obj <- CreateSeuratObject(counts=data)
    seurat_obj <- NormalizeData(seurat_obj, verbose=FALSE)
    Idents(seurat_obj) <- identities

    markers <- FindMarkers(
      seurat_obj, ident.1=test_id1, ident.2=test_id2, test.use="wilcox",
      logfc.threshold=0.0, min.pct=0., # no pre-selection of genes
      verbose=FALSE
    )

    # format results
    p_values <- markers$p_val_adj
    names(p_values) <- rownames(markers)
  }

  return(p_values)
}

multicore_setup <- function(num_workers, blas_threads){
  cat(paste0("Change multicore setup: workers=", num_workers, ", threads=", blas_threads, "\n"))
  blas_set_num_threads(blas_threads)
  plan(multicore, workers=num_workers)
  registerDoFuture()
}

StabCell_steps <- function(steps, checkpoint){
  StabCell(
    data=count_data,
    reference_identities=initial_idents,
    cluster_fn=cluster_fn,
    test_fn=test_fn,
    PFER_max=PFER_max,
    p_value_threshold=p_value_threshold,
    n_subsample_markers=n_subsample_markers,
    cluster_args=NULL,
    test_args=NULL,
    test_global=TRUE,
    test_local=TRUE,
    subsampling_iterations=50,
    subsampling_fraction=0.5,
    complementary_pairs=TRUE,
    cell_stability_metric="gini",
    assumption=c("none", "unimodal", "r_concavity"),
    ranking="p_val",
    theta=NULL,
    cells_in_columns=TRUE,
    packages=c("Seurat"),    # packages that are needed for clustering or testing
    export=NULL,             # additional custom objects needed for clustering or testing
    steps_to_run=steps,
    checkpoint=checkpoint,
    checkpoint_path=checkpoint_file,
    verbose=TRUE
  )
}

# run StabCell
cat("run StabCell...\n")
if(is.null(steps) || 1 %in% steps){
  checkpoint <- NULL
} else if(!is.null(checkpoint_file)) {
  cat("load checkpoint...\n")
  load(checkpoint_file)
} else {
  stop("Could not load checkpoint!")
}

time_start <- Sys.time()
for(i in 1:length(step_settings)){

  if(!is.null(steps)){
    run_steps <- intersect(step_settings[[i]]$steps, steps)
  } else {
    run_steps <- step_settings[[i]]$steps
  }

  if(length(run_steps) == 0) next

  # create multicore setup for specific step(s)
  multicore_setup(
    num_workers=step_settings[[i]]$num_workers,
    blas_threads=step_settings[[i]]$blas_threads
  )

  # run step(s)
  cat("run steps", paste(run_steps, collapse=", "), "\n")
  with_progress({
    checkpoint <- StabCell_steps(steps=run_steps, checkpoint=checkpoint)
  })
}
time_end <- Sys.time()

stabcell_results <- checkpoint

time <- list(
  time_start=time_start,
  time_end=time_end,
  duration=time_end - time_start
)

save(
  stabcell_results,
  rnd_idx,
  steps,
  step_settings,
  time,
  file=results_file
)


print(summary(warnings()))




