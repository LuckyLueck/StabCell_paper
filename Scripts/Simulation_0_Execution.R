library(RhpcBLASctl)
library(doFuture)
library(doRNG)
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
num_workers <- 40

options(future.globals.maxSize=16 * 1024^3)

blas_set_num_threads(blas_threads)
plan(multicore, workers=num_workers)
registerDoFuture()

seed <- 08052025

scenarios <- c("2Cluster_Simulation", "NCluster_Simulation")
set_scenario_settings <- function(scenario){
  if(scenario == "2Cluster_Simulation"){
    n_celltypes <<- c(1, 2)
    n_clusters <<- 2  # always find 2 clusters when clustering
    test_local <<- FALSE  # local & global are the same for 2 cell types!
  } else if(scenario == "NCluster_Simulation") {
    n_celltypes <<- c(1, 10)
    n_clusters <<- NULL # no restriction on clusters to be found
    test_local <<- TRUE
  } else {
    stop("Unknown scenario!")
  }
}

repetitions <- 100
fixed_num_subsample_markers <- FALSE

subsampling_fraction <- 0.5
subsampling_iterations <- 50
complementary_pairs <- TRUE
PFER_max <- c(1, 5, 10, 25, 50, 100)
cell_stability_metric <- c("gini", "entropy", "L2norm")[1]

de.probs <- c(0.01, 0.05)   # larger values --> less sparse
de.facLocs <- c(0.1, 0.4)   # larger value --> stronger differences between cell types
lib.locs <- c(9.3, 10.3, 11.) # ~10k/25k/50k reads per cell, default: lib.loc=11
# hist(colSums(counts(splatSimulate(newSplatParams(batchCells=5000, nGenes=100, lib.loc=11)))))
# exp(c(9.3, 10.3, 11) + 0.2 / 2)
# log(c(10000, 25000, 50000)) - 0.2 / 2


if(fixed_num_subsample_markers){
  n_subsample_markers <- c(10, 25, 50, 100, 250, 500)
  p_value_threshold <- NULL
} else {
  n_subsample_markers <- NULL
  p_value_threshold <- c(0.01, 0.05, 0.1)
}




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
  seurat_obj <- CreateSeuratObject(counts=data)
  seurat_obj <- NormalizeData(seurat_obj, verbose=FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose=FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose=FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose=FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, verbose=FALSE)
  seurat_obj <- FindNClusters(seurat_obj, n_clusters=n_clusters)
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






with_progress({
  for(scenario in scenarios){

    cat(
      "\n",
      paste(rep("#", 14 + nchar(scenario)), collapse=""), "\n",
      "#### run", scenario, "####\n",
      paste(rep("#", 14 + nchar(scenario)), collapse=""), "\n"
    )
    cat(paste0("start at ", strftime(Sys.time())), "\n")

    set_scenario_settings(scenario)

    # construct grid with simulation parameters
    sim_grid <- expand.grid(
      n_celltypes=n_celltypes,
      lib.loc=lib.locs,
      de.prob=de.probs,
      de.facLoc=de.facLocs,
      iter=1:repetitions
    ) %>%
      filter(
        # for one cell type differential expression is meaningless
        # --> only 1 parameter set required
        !(n_celltypes == 1 & (de.prob != de.probs[1] | de.facLoc != de.facLocs[1]))
      ) %>%
      arrange(n_celltypes, lib.loc, de.facLoc, de.prob, iter) %>%
      mutate(setting=paste(
        paste0("celltypes=", n_celltypes),
        ifelse(complementary_pairs, "complementary", "subsampling"),
        paste0("iter=", subsampling_iterations),
        paste0("lib.loc=", lib.loc),
        paste0("de.prob=", de.prob),
        paste0("de.facLoc=", de.facLoc),
        paste0("subsampling=", subsampling_fraction),
        paste0("iter=", iter),
        sep="__"
      ))

    cat("total number of scenario settings:", nrow(sim_grid), "\n")

    # remove settings that have been run before
    sim_grid <- sim_grid %>%
      filter(!file.exists(paste0("results/", scenario, "/results_simulation__", setting, ".Rdata")))

    cat("num scenarios left to run:", nrow(sim_grid), "\n")

    if(nrow(sim_grid) == 0) next


    p <- progressor(steps=nrow(sim_grid))

    foreach(
      i=1:nrow(sim_grid),
      .combine="bind_rows",
      .export="StabCell"
    ) %dorng% {

      # no parallelization within the simulation setting (but resore old plan)
      on.exit(plan(plan()), add=TRUE)
      plan(sequential)


      de.prob <- sim_grid[i, "de.prob"]
      de.facLoc <- sim_grid[i, "de.facLoc"]
      n_celltypes <- sim_grid[i, "n_celltypes"]
      lib.loc <- sim_grid[i, "lib.loc"]
      iter <- sim_grid[i, "iter"]

      setting <- sim_grid[i, "setting"]

      results_file <- paste0("results/", scenario, "/results_simulation__", setting, ".Rdata")


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
          group.prob=rep(1. / n_celltypes, n_celltypes),

          # library size (= counts per cell)
          lib.loc=lib.loc,
          lib.scale=0.2,

          # dropout probabilities
          dropout.type="experiment", # options: "none", "experiment", "batch, "group", "cell"
          dropout.mid=50,
          dropout.shape=5e-2, # smaller -> more dropouts ?

          # differential expression
          de.prob=de.prob,  # probability of differentially expressed gene, can be set differently for each group
          de.downProb=0.5,  # probability of downregulation vs upregulation
          de.facLoc=de.facLoc,  # location of logNorm distributed DE factors
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
        if(length(true_idents) == 0){
          true_idents <- rep(1, ncol(count_data))
        }
        true_fold_changes <- as.data.frame(rowData(sim)[startsWith(colnames(rowData(sim)), "DEFacGroup")])
      }
      # cat("simulation finished")
      # print(Sys.time())

      ################################################################################
      ################################################################################
      ####                                                                        ####
      ####                              Double Dipping                            ####
      ####                                                                        ####
      ################################################################################
      ################################################################################
      {
        seurat_obj <- CreateSeuratObject(counts=count_data)
        seurat_obj <- NormalizeData(seurat_obj, verbose=FALSE)
        seurat_obj <- FindVariableFeatures(seurat_obj, verbose=FALSE)
        seurat_obj <- ScaleData(seurat_obj, verbose=FALSE)
        seurat_obj <- RunPCA(seurat_obj, verbose=FALSE)
        seurat_obj <- FindNeighbors(seurat_obj, verbose=FALSE)
        seurat_obj <- FindNClusters(seurat_obj, n_clusters=n_clusters)

        seurat_cluster_ids <- levels(Idents(seurat_obj))
        df_list <- list()
        for(i in seq_along(seurat_cluster_ids)){
          for(j in seq_along(seurat_cluster_ids)){

            test_id1 <- seurat_cluster_ids[i]

            if(j < i){
              next
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
        seurat_idents <- Idents(seurat_obj)

      }
      # cat("ddip finished")
      # print(Sys.time())

      ################################################################################
      ################################################################################
      ####                                                                        ####
      ####                                CountSplit                              ####
      ####                                                                        ####
      ################################################################################
      ################################################################################
      {
        # split dataset
        split <- countsplit(
          count_data,
          epsilon=c(0.5, 0.5),
          overdisps=rep(Inf, ncol(count_data)) # for Poisson
        )
        data_clust <- split[[1]]
        data_marker <- split[[2]]

        # clustering on one split
        seurat_cluster <- CreateSeuratObject(counts=data_clust, assay="CLUSTER")
        seurat_cluster <- NormalizeData(seurat_cluster, verbose=FALSE)
        seurat_cluster <- FindVariableFeatures(seurat_cluster, verbose=FALSE)
        seurat_cluster <- ScaleData(seurat_cluster, verbose=FALSE)
        seurat_cluster <- RunPCA(seurat_cluster, verbose=FALSE)
        seurat_cluster <- FindNeighbors(seurat_cluster, verbose=FALSE)
        seurat_cluster <- FindNClusters(seurat_cluster, n_clusters=n_clusters)
        countsplit_idents <- Idents(seurat_cluster)

        # marker detection on other split
        seurat_marker <- CreateSeuratObject(counts=data_marker, assay="MARKER")
        seurat_marker$countsplit_idents <- countsplit_idents
        Idents(seurat_marker) <- "countsplit_idents"
        seurat_marker <- NormalizeData(seurat_marker, verbose=FALSE)

        countsplit_cluster_ids <- levels(countsplit_idents)
        df_list <- list()
        for(i in seq_along(countsplit_cluster_ids)){
          for(j in seq_along(countsplit_cluster_ids)){

            test_id1 <- countsplit_cluster_ids[i]

            if(j < i){
              next
            } else if(i == j){
              test_id2 <- NULL
            } else{
              test_id2 <- countsplit_cluster_ids[j]
            }

            markers <- FindMarkers(
              seurat_marker, ident.1=test_id1, ident.2=test_id2, test.use="wilcox",
              logfc.threshold=0.0, min.pct=0., # no pre-selection of genes
              verbose=FALSE
            )

            markers$cluster1 <- test_id1
            markers$cluster2 <- if(is.null(test_id2)) NA else test_id2
            markers$marker <- rownames(markers)
            rownames(markers) <- NULL

            df_list[[length(df_list)+1]] <- markers
          }
        }

        countsplit_markers <- bind_rows(df_list)
      }
      # cat("csplit finished")
      # print(Sys.time())

      ################################################################################
      ################################################################################
      ####                                                                        ####
      ####                                StabCell                                ####
      ####                                                                        ####
      ################################################################################
      ################################################################################

      stabcell_results <- StabCell(
        data=count_data,
        reference_identities=seurat_idents,
        cluster_fn=cluster_fn,
        test_fn=test_fn,
        cluster_args=NULL,
        test_args=NULL,
        p_value_threshold=p_value_threshold,
        n_subsample_markers=n_subsample_markers,
        PFER_max=PFER_max,
        test_global=TRUE,
        test_local=test_local,
        cell_stability_metric=cell_stability_metric,
        subsampling_iterations=subsampling_iterations,
        subsampling_fraction=subsampling_fraction,
        complementary_pairs=complementary_pairs,
        ranking="p_val",
        assumption=c("none", "unimodal", "r_concavity"),
        theta=NULL,
        steps_to_run=NULL,
        cells_in_columns=TRUE,
        checkpoint_path=NULL,
        checkpoint=NULL,
        packages=c("Seurat"),    # packages that are needed for clustering or testing
        export=c("FindNClusters", "n_clusters"),   # additional custom objects needed for clustering or testing
        verbose=FALSE
      )
      # cat("stabcell finished")
      # print(Sys.time())

      ######################
      #### save results ####
      ######################

      results <- list(
        config=list(
          n_celltypes=n_celltypes,
          lib.loc=lib.loc,
          de.prob=de.prob,
          de.facLoc=de.facLoc,
          subsampling_fraction=subsampling_fraction,
          complementary_pairs=complementary_pairs,
          n_iterations=subsampling_iterations
        ),
        count_data=count_data,
        true_fold_changes=true_fold_changes,
        true_idents=true_idents,
        double_dipping=list(
          markers=double_dipping_markers,
          identities=seurat_idents
        ),
        countsplit=list(
          markers=countsplit_markers,
          identities=countsplit_idents
        ),
        stabcell=stabcell_results
      )

      save(results, file=results_file)

      p()
      return(NULL)
    }

    print(summary(warnings()))
    cat(paste0("finished at ", strftime(Sys.time())), "\n")
  }
})





