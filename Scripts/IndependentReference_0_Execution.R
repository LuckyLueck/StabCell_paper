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

cluster_fn <- function(data, subsample_indices, cluster_args=NULL){
  seurat_obj <- CreateSeuratObject(counts=data)
  seurat_obj <- NormalizeData(seurat_obj, verbose=FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose=FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose=FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose=FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, verbose=FALSE)
  seurat_obj <- FindClusters(seurat_obj, verbose=FALSE)
  identities <- as.numeric(Idents(seurat_obj))
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
      count_data <- counts(sim)
      count_data <- as(count_data, "CsparseMatrix")
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
      stabcell_results_no_split <- StabCell(
        data=count_data,
        reference_identities=reference_idents_nosplit,
        cluster_fn=cluster_fn,
        test_fn=test_fn,
        cluster_args=NULL,
        test_args=NULL,
        p_value_threshold=0.05,
        subsampling_iterations=B,
        PFER_max=PFER_max,
        assumption=assumption,
        packages=c("Seurat"),    # packages that are needed for clustering or testing
        verbose=FALSE
      )
      stabcell_results_no_split[["reference_idents"]] <- reference_idents_nosplit
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
      marker_alignment <- MarkerStabilityAlignment(
        reference_identities=marker_idents,
        subsampling_identities=subsampling_identities
      )
      subsampling_p_values <- SubsampleTesting(
        data=count_data_marker,
        reference_identities=marker_idents,
        subsampling_identities=subsampling_identities,
        marker_alignment=marker_alignment,
        cell_stability_score=NULL,
        test_fn=test_fn,
        packages="Seurat",
        verbose=FALSE
      )

      # find stable markers
      selection_frequencies <- CalculateSelectionFrequency(
        subsampling_p_values=subsampling_p_values,
        p_value_threshold=0.05,
        n_subsample_markers=NULL
      )
      stable_markers_global <- FindAllGlobalStableMarkers(
        selection_frequencies=selection_frequencies,
        PFER_max=PFER_max,
        assumption=assumption,
        B=B,
        theta=NULL
      )

      stabcell_results_with_split <- list(
        idx_split=idx_split,
        clustering_idents=clustering_idents,
        marker_idents=marker_idents,
        subsample_idents=subsampling_identities,
        subsampling_p_values=subsampling_p_values,
        selection_frequencies=selection_frequencies,
        stable_markers_global=stable_markers_global
      )
    }

    #### diagnostics ####

    # split_idents <- rep(NA, ncol(count_data))
    # split_idents[idx_split] <- clustering_idents
    # split_idents[-idx_split] <- marker_idents
    #
    # table(reference_idents_nosplit, true_idents)
    # table(stabcell_results_no_split$clustering_mode, true_idents)
    # table(split_idents, true_idents)
    # {
    #   # select ids corresponding to true_ident == 1
    #   seurat_id <- 5
    #   stabcell_id <- 5
    #   split_id <- 2
    #
    #   # select correct DEFacGroup!!!
    #   true_markers <- rownames(true_fold_changes %>% select(DEFacGroup1) %>% filter(DEFacGroup1 != 1))
    #   false_markers <- setdiff(rownames(count_data), true_markers)
    #
    #   seurat_marker <- test_fn(count_data, reference_idents_nosplit, test_id1=seurat_id)
    #   seurat_marker <- names(seurat_marker[seurat_marker <= 0.05])
    #   fp_markers_seurat <- length(intersect(false_markers, seurat_marker))
    #   tp_markers_seurat <- length(intersect(true_markers, seurat_marker))
    #
    #   stabcell_marker <- (stabcell_results_no_split$stable_markers_global %>% filter(id == !!stabcell_id))$marker
    #   fp_markers_stabcell <- length(intersect(false_markers, stabcell_marker))
    #   tp_markers_stabcell <- length(intersect(true_markers, stabcell_marker))
    #
    #   split_marker <- (stabcell_results_with_split$stable_markers_global %>% filter(id == !!split_id))$marker
    #   fp_markers_split <- length(intersect(false_markers, split_marker))
    #   tp_markers_split <- length(intersect(true_markers, split_marker))
    #
    #
    #   c(
    #     FP_seurat=fp_markers_seurat,
    #     TP_seurat=tp_markers_seurat,
    #     FP_stabcell=fp_markers_stabcell,
    #     TP_stabcell=tp_markers_stabcell,
    #     FP_split=fp_markers_split,
    #     TP_split=tp_markers_split
    #   )
    # }
    # plotfn <- function(coords, idents){
    #    ggplot(
    #          data=as.data.frame(coords),
    #          mapping=aes(x=PC_1, y=PC_2, color=idents)
    #      ) +
    #          geom_point(size=0.1)
    # }
    # plotfn(clustering_embeddings, as.factor(clustering_idents)) + plotfn(marker_embeddings, marker_idents)

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


