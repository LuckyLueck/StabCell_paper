setwd("/mnt/mikasa_data/nlueck/StabCell_paper")

library(RhpcBLASctl)
library(doRNG)
library(doFuture)
library(progressr)

library(dplyr)

library(splatter)
library(Seurat)
library(countsplit)
library(StabCell)
source("Scripts/helpers/cluster_fn.R")
source("Scripts/helpers/marker_fn.R")

################
#### config ####
################
blas_threads <- 8
num_workers <- 25

options(future.globals.maxSize=12 * 1024^3)

blas_set_num_threads(blas_threads)

seed <- 01042026

n_celltypes <- 10

num_nearest_neighbors <- 10  # when transferring cluster labels

B <- 50
PFER_max <- 25
assumption <- c("none", "unimodal", "r_concavity")

repetitions <- 100


sim_grid <- expand.grid(iter=1:repetitions) %>%
  mutate(setting=paste0("iter=", iter))

cat("total number of scenario settings:", nrow(sim_grid), "\n")

# remove settings that have been run before
sim_grid <- sim_grid %>%
  filter(!file.exists(paste0("results/independent_reference/results_indepRef__", setting, ".Rdata")))

n_settings <- nrow(sim_grid)

cat("num settings left to run:", n_settings, "\n")


plan(multicore, workers=num_workers)
registerDoFuture()

with_progress({

  p <- progressor(steps=n_settings)

  x <- foreach(i=1:n_settings) %dorng% {
    plan(sequential)
    registerDoSEQ()

    iter <- sim_grid[i, "iter"]
    setting <- sim_grid[i, "setting"]

    results_file <- paste0("results/independent_reference/results_indepRef__", setting, ".Rdata")


    set.seed(seed + iter)

    #######################
    #### simulate data ####
    #######################
    {
      params <- newSplatParams(
        # see: https://bioconductor.org/packages/release/bioc/vignettes/splatter/inst/doc/splat_params.html

        seed=seed + iter,

        nGenes=10000,
        batchCells=5000,

        # cell types
        group.prob=rep(1 / n_celltypes, n_celltypes),

        # library size (= counts per cell)
        lib.loc=11.,
        lib.scale=0.2,

        # dropout probabilities
        dropout.type="experiment", # options: "none", "experiment", "batch, "group", "cell"
        dropout.mid=50,
        dropout.shape=5e-2, # smaller -> more dropouts ?

        # differential expression
        de.prob=0.05,  # probability of differentially expressed gene, can be set differently for each group
        de.downProb=0.5,  # probability of downregulation vs upregulation
        de.facLoc=0.1,  # location of logNorm distributed DE factors
        de.facScale=0.4,  # scale of logNorm distributed DE factors

      )

      sim <- splatSimulate(
        params,
        method="groups",
        verbose=FALSE
      )
      count_data <- as(counts(sim), "CsparseMatrix")
      true_idents <- as.numeric(substr(colData(sim)$Group, 6, 10))
      true_fold_changes <- as.data.frame(rowData(sim)[startsWith(colnames(rowData(sim)), "DEFacGroup")])
    }



    ################################################################################
    ################################################################################
    ####                                                                        ####
    ####                                StabCell                                ####
    ####                                                                        ####
    ################################################################################
    ################################################################################

    # without splitting
    {
      reference_idents_nosplit <- cluster_fn(count_data, NULL)
      without_progress({
        stabcell_results_no_split <- StabCell(
          data=count_data,
          cluster_fn=cluster_fn,
          PFER_max=PFER_max,
          reference_identities=reference_idents_nosplit,
          marker_fn=marker_fn,
          p_value_threshold=0.05,
          cluster_args=NULL,
          marker_args=NULL,
          marker_global=TRUE,
          marker_local=FALSE,
          subsampling_iterations=B,
          subsampling_fraction=0.5,
          complementary_pairs=TRUE,
          assumption=assumption,
          cells_in_columns=TRUE,
          packages=c("Seurat"),
          verbose=FALSE
        )
      })
    }

    # with splitting for clustering and marker detection
    {
      # split data
      n_cells <- ncol(count_data)
      idx_split <- sample.int(n_cells, n_cells %/% 2, replace=FALSE)
      count_data_clust <- count_data[ , idx_split]
      count_data_marker <- count_data[ , -idx_split]

      # perform clustering
      seurat_obj_clust <- CreateSeuratObject(counts=count_data_clust)
      seurat_obj_clust <- NormalizeData(seurat_obj_clust, verbose=FALSE)
      clustering_norm_data <- GetAssayData(seurat_obj_clust, layer="data")
      seurat_obj_clust <- FindVariableFeatures(seurat_obj_clust, verbose=FALSE)
      seurat_obj_clust <- ScaleData(seurat_obj_clust, verbose=FALSE)
      seurat_obj_clust <- RunPCA(seurat_obj_clust, verbose=FALSE)
      seurat_obj_clust <- FindNeighbors(seurat_obj_clust, verbose=FALSE)
      seurat_obj_clust <- FindClusters(seurat_obj_clust, verbose=FALSE)
      clustering_idents <- as.numeric(Idents(seurat_obj_clust))

      # project marker data in clustering data pca space
      # - scaling mean and sd must be calculated manually
      clustering_embeddings <- Embeddings(seurat_obj_clust, reduction="pca")
      clustering_loadings <- Loadings(seurat_obj_clust, reduction="pca")
      pca_features <- rownames(clustering_loadings)

      scale_means <- rowMeans(clustering_norm_data[pca_features, ])
      scale_sds <- apply(clustering_norm_data[pca_features, ], 1, sd)
      scale_sds[scale_sds == 0] <- 0.01 # same as seurat

      seurat_obj_marker <- CreateSeuratObject(counts=count_data_marker)
      seurat_obj_marker <- NormalizeData(seurat_obj_marker, verbose=FALSE)
      marker_norm_data <- GetAssayData(seurat_obj_marker, layer="data")
      marker_scale_data <- (marker_norm_data[pca_features, ] - scale_means) / scale_sds
      clipping_value <- 10 # same as seurat
      marker_scale_data[marker_scale_data > clipping_value] <- clipping_value
      marker_embeddings <- t(marker_scale_data) %*% clustering_loadings

      # find marker data cluster labels with knn in pca space
      marker_idents <- class::knn(
        train=clustering_embeddings,
        test=marker_embeddings,
        cl=clustering_idents,
        k=num_nearest_neighbors
      )

      # generate subsamples and their identities
      n_cells_marker <- ncol(count_data_marker)
      n_subsample <- n_cells_marker %/% 2
      subsampling_identities <- matrix(data=NaN, nrow=2*B, ncol=n_cells_marker)
      for(b in 1:B){
        idx_subsample <- sample.int(n_cells_marker, n_subsample, replace=FALSE)
        for(complementary_factor in c(1, -1)){
          if((complementary_factor == -1) && (n_subsample %% 2 == 1)){
            idx_subsample <- c(idx_subsample, sample(1:n_cells_marker[-idx_subsample], 1))
          }
          idx_subsample <- complementary_factor * idx_subsample
          row_idx <- if(complementary_factor == 1) 2 * b - 1 else 2 * b
          subsampling_identities[row_idx, idx_subsample] <- marker_idents[idx_subsample]
        }
      }

      # perform marker detection in subsamples
      without_progress({
        marker_alignment <- MarkerStabilityAlignment(
          reference_identities=marker_idents,
          subsampling_identities=subsampling_identities
        )
        subsampling_p_values <- SubsampleMarkers(
          data=count_data_marker,
          reference_identities=marker_idents,
          subsampling_identities=subsampling_identities,
          marker_alignment=marker_alignment,
          cell_stability_score=NULL,
          marker_fn=marker_fn,
          marker_args=NULL,
          packages="Seurat",
          verbose=FALSE
        )
  
        # find stable markers
        marker_selection_frequencies <- CalculateMarkerSelectionFrequency(
          subsampling_p_values=subsampling_p_values,
          p_value_threshold=0.05,
          n_subsample_markers=NULL
        )
        stable_markers_global <- FindAllGlobalStableMarkers(
          marker_selection_frequencies=marker_selection_frequencies,
          PFER_max=PFER_max,
          assumption=assumption,
          B=B,
          theta=NULL
        )
      })

      stabcell_results_with_split <- list(
        idx_split=idx_split,
        clustering_idents=clustering_idents,
        marker_idents=marker_idents,
        subsample_idents=subsampling_identities,
        subsampling_p_values=subsampling_p_values,
        marker_selection_frequencies=marker_selection_frequencies,
        stable_markers_global=stable_markers_global
      )
    }


    ######################
    #### save results ####
    ######################
    {
      results <- list(
        iter=iter,
        count_data=count_data,
        true_fold_changes=true_fold_changes,
        true_idents=true_idents,
        no_split=stabcell_results_no_split,
        with_split=stabcell_results_with_split
      )

      save(results, file=results_file)
    }

    p()
    return(NULL)
  }
})


cat(paste0("finished at ", strftime(Sys.time())))

print(summary(warnings()))


