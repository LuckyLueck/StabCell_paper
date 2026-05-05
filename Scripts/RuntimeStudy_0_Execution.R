library(RhpcBLASctl)
library(doRNG)
library(doFuture)
library(progressr)

library(dplyr)

library(splatter)
library(Seurat)
library(countsplit)
library(StabCell)

################
#### config ####
################
blas_threads <- 4
num_workers_stabcell <- 25

options(future.globals.maxSize=12 * 1024^3)

blas_set_num_threads(blas_threads)

seed <- 22032026

n_cells <- c(500, 1000, 5000, 10000, 20000, 35000, 50000)
n_celltypes <- 10

PFER_max <- 25
assumption <- "r_concavity"

repetitions <- 25


cluster_fn <- function(data, subsample_indices, cluster_args=NULL){
  seurat_obj <- CreateSeuratObject(counts=data)
  seurat_obj <- NormalizeData(seurat_obj, verbose=FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose=FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose=FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose=FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, verbose=FALSE)
  seurat_obj <- FindClusters(seurat_obj, verbose=FALSE)
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




sim_grid <- expand.grid(
  n_cells=n_cells,
  iter=1:repetitions
) %>%
  arrange(n_cells, iter) %>%
  mutate(setting=paste(
    paste0("cells=", n_cells),
    paste0("iter=", iter),
    sep="__"
  ))

cat("total number of scenario settings:", nrow(sim_grid), "\n")

# remove settings that have been run before
sim_grid <- sim_grid %>%
  filter(!file.exists(paste0(
    "results/runtime_study__n_celltypes=", n_celltypes,
    "/results_runtime__", setting, ".Rdata"
  )))

n_settings <- nrow(sim_grid)

cat("num settings left to run:", n_settings, "\n")


with_progress({

  p <- progressor(steps=n_settings)

  for(i in 1:n_settings){
    plan(sequential)
    registerDoSEQ()

    t_start <- Sys.time()

    n_cells <- sim_grid[i, "n_cells"]
    iter <- sim_grid[i, "iter"]
    setting <- sim_grid[i, "setting"]

    results_file <- paste0(
      "results/runtime_study__n_celltypes=", n_celltypes,
      "/results_runtime__", setting, ".Rdata"
    )


    set.seed(seed + iter)

    #######################
    #### simulate data ####
    #######################
    {
      params <- newSplatParams(
        # see: https://bioconductor.org/packages/release/bioc/vignettes/splatter/inst/doc/splat_params.html

        seed=seed + iter,

        nGenes=2500,
        batchCells=n_cells,

        # cell types
        group.prob=rep(1 / n_celltypes, n_celltypes),

        # library size (= counts per cell)
        lib.loc=5.,
        lib.scale=0.2,

        # dropout probabilities
        dropout.type="experiment", # options: "none", "experiment", "batch, "group", "cell"
        dropout.mid=50,
        dropout.shape=5e-2, # smaller -> more dropouts ?

        # differential expression
        de.prob=0.05,  # probability of differentially expressed gene, can be set differently for each group
        de.downProb=0.5,  # probability of downregulation vs upregulation
        de.facLoc=15.,  # location of logNorm distributed DE factors
        de.facScale=0.4,  # scale of logNorm distributed DE factors

      )

      sim <- splatSimulate(
        params,
        method="groups",
        verbose=FALSE
      )
      count_data <- counts(sim)
      count_data <- as(count_data, "CsparseMatrix")
    }
    message(paste0("setting ", i, "/", n_settings, ": finished data simulation"))

    ####################
    #### Run Seurat ####
    ####################
    {
      t_start_seurat <- Sys.time()

      seurat_obj <- CreateSeuratObject(counts=count_data)
      seurat_obj <- NormalizeData(seurat_obj, verbose=FALSE)
      seurat_obj <- FindVariableFeatures(seurat_obj, verbose=FALSE)
      seurat_obj <- ScaleData(seurat_obj, verbose=FALSE)
      seurat_obj <- RunPCA(seurat_obj, verbose=FALSE)
      seurat_obj <- FindNeighbors(seurat_obj, verbose=FALSE)
      seurat_obj <- FindClusters(seurat_obj, verbose=FALSE)
      seurat_idents <- Idents(seurat_obj)

      t_clust_seurat <- Sys.time()

      markers_ddip <- FindAllMarkers(
        seurat_obj, test.use="wilcox",
        logfc.threshold=0.0, min.pct=0., # no pre-selection of genes
        verbose=FALSE
      )

      t_marker_seurat <- Sys.time()

      duration_clust_seurat <- t_clust_seurat - t_start_seurat
      duration_marker_seurat <- t_marker_seurat - t_clust_seurat
      duration_total_seurat <- t_marker_seurat - t_start_seurat
    }
    message(paste0("setting ", i, "/", n_settings, ": finished double dipping"))

    ################################################################################
    ################################################################################
    ####                                                                        ####
    ####                                CountSplit                              ####
    ####                                                                        ####
    ################################################################################
    ################################################################################
    {
      t_start_countsplit <- Sys.time()

      # split dataset
      split <- countsplit(
        count_data,
        epsilon=c(0.5, 0.5),
        overdisps=rep(Inf, ncol(count_data)) # Poisson splitting
      )
      data_clust <- split[[1]]
      data_marker <- split[[2]]

      t_split_countsplit <- Sys.time()

      # clustering on one split
      seurat_cluster <- CreateSeuratObject(counts=data_clust, assay="CLUSTER")
      seurat_cluster <- NormalizeData(seurat_cluster, verbose=FALSE)
      seurat_cluster <- FindVariableFeatures(seurat_cluster, verbose=FALSE)
      seurat_cluster <- ScaleData(seurat_cluster, verbose=FALSE)
      seurat_cluster <- RunPCA(seurat_cluster, verbose=FALSE)
      seurat_cluster <- FindNeighbors(seurat_cluster, verbose=FALSE)
      seurat_cluster <- FindClusters(seurat_cluster, verbose=FALSE)
      countsplit_idents <- Idents(seurat_cluster)

      t_clust_countsplit <- Sys.time()

      # marker detection on other split
      seurat_marker <- CreateSeuratObject(counts=data_marker, assay="MARKER")
      seurat_marker$countsplit_idents <- countsplit_idents
      Idents(seurat_marker) <- "countsplit_idents"
      seurat_marker <- NormalizeData(seurat_marker, verbose=FALSE)
      markers_ddip <- FindAllMarkers(
        seurat_obj, test.use="wilcox",
        logfc.threshold=0.0, min.pct=0., # no pre-selection of genes
        verbose=FALSE
      )

      t_marker_countsplit <- Sys.time()

      duration_split_countsplit <- t_split_countsplit - t_start_countsplit
      duration_clust_countsplit <- t_clust_countsplit - t_split_countsplit
      duration_marker_countsplit <- t_marker_countsplit - t_clust_countsplit
      duration_total_countsplit <- t_marker_countsplit - t_start_countsplit
    }
    message(paste0("setting ", i, "/", n_settings, ": finished count splitting"))

    ################################################################################
    ################################################################################
    ####                                                                        ####
    ####                                StabCell                                ####
    ####                                                                        ####
    ################################################################################
    ################################################################################
    plan(multicore, workers=num_workers_stabcell)
    registerDoFuture()
    {
      t_start_stabcell <- Sys.time()

      stabcell_results <- StabCell(
        data=count_data,
        reference_identities=seurat_idents,
        cluster_fn=cluster_fn,
        test_fn=test_fn,
        cluster_args=NULL,
        test_args=NULL,
        PFER_max=PFER_max,
        assumption=assumption,
        packages=c("Seurat"),    # packages that are needed for clustering or testing
        verbose=FALSE
      )
      t_end_stabcell <- Sys.time()

      duration_refclust_stabcell <- duration_clust_seurat
      duration_subsamples_stabcell <- t_end_stabcell - t_start_stabcell
      duration_total_stabcell <- duration_refclust_stabcell + duration_subsamples_stabcell
    }
    message(paste0("setting ", i, "/", n_settings, ": finished StabCell"))

    ######################
    #### save results ####
    ######################
    {
      runtimes <- list(
        n_cells=n_cells,
        iter=iter,
        ddip=list(
          duration_clust=duration_clust_seurat,
          duration_marker=duration_marker_seurat,
          duration_total=duration_total_seurat
        ),
        countsplit=list(
          duration_split=duration_split_countsplit,
          duration_clust=duration_clust_countsplit,
          duration_marker=duration_marker_countsplit,
          duration_total=duration_total_countsplit
        ),
        stabcell=list(
          duration_refclust=duration_refclust_stabcell,
          duration_subsamples=duration_subsamples_stabcell,
          duration_total=duration_total_stabcell
        )
      )

      save(runtimes, file=results_file)
    }
    message(paste0("setting ", i, "/", n_settings, ": finished saving"))
    message(paste0("setting ", i, "/", n_settings, ": total time = ", Sys.time() - t_start))
    p()

  }
})


cat(paste0("finished at ", strftime(Sys.time())))

print(summary(warnings()))


