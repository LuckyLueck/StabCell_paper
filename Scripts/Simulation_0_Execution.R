setwd("/mnt/mikasa_data/nlueck/StabCell_paper")
Sys.setenv(
  RETICULATE_ENABLE_PYTHON_FINALIZER="yes",
  MPLBACKEND="Agg"
)

library(RhpcBLASctl)
library(doFuture)
library(doRNG)
library(progressr)

library(dplyr)

library(splatter)
library(Seurat)
library(countsplit)
library(ClusterDE)
# library(ZINBMM) # slightly modified version to also return beta coefficient matrix and beta_global vector

library(StabCell)
source("Scripts/helpers/cluster_fn.R")
source("Scripts/helpers/marker_fn.R")
source("Scripts/helpers/utils.R")


################
#### config ####
################
options(
  future.globals.maxSize=16 * 1024^3,
  nwarnings=500000
)

# set up python venv if not exists
if(!reticulate::virtualenv_exists("./.venv")){
  reticulate::virtualenv_create(
    "./.venv",
    version="3.12",
    packages=c("numpy", "pandas", "scipy", "matplotlib", "scikit-learn", "joblib", "truncated_normal")
  )
}

blas_threads <- 4
num_workers <- 50

plan(multisession, workers=num_workers)
registerDoFuture()

seed <- 08052025

scenarios <- c("2Cluster_Simulation", "NCluster_Simulation")
set_scenario_settings <- function(scenario){
  if(scenario == "2Cluster_Simulation"){
    n_celltypes <<- c(1, 2)
    n_clusters <<- 2  # always find 2 clusters when clustering
  } else if(scenario == "NCluster_Simulation") {
    n_celltypes <<- c(1, 10)
    n_clusters <<- NULL # no restriction on clusters to be found
  } else {
    stop("Unknown scenario!")
  }
}

n_cells <- 5000
n_genes <- 10000

repetitions <- 100
fixed_num_subsample_markers <- FALSE

subsampling_iterations <- 50
PFER_max <- c(1, 5, 10, 25, 50, 100)
cell_stability_metric <- c("gini", "entropy", "L2norm")[1]

de.probs <- c(0.01, 0.05)   # larger values --> less sparse
de.facLocs <- c(0.1, 0.4)   # larger value --> stronger differences between cell types
lib.scale <- 0.2
lib.size <- c(10000, 25000, 50000)



marker_args <- NULL

if(fixed_num_subsample_markers){
  n_subsample_markers <- c(10, 25, 50, 100, 250, 500)
  p_value_threshold <- NULL
} else {
  n_subsample_markers <- NULL
  p_value_threshold <- c(0.01, 0.05, 0.1)
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
      lib.size=lib.size,
      de.prob=de.probs,
      de.facLoc=de.facLocs,
      iter=1:repetitions
    ) %>%
      filter(
        # for one cell type differential expression is meaningless
        # --> only 1 parameter set required
        !(n_celltypes == 1 & (de.prob != de.probs[1] | de.facLoc != de.facLocs[1]))
      ) %>%
      arrange(desc(n_celltypes), lib.size, de.facLoc, de.prob, iter) %>%
      mutate(setting=paste(
        paste0("celltypes=", n_celltypes),
        paste0("lib.size=", lib.size),
        paste0("de.prob=", de.prob),
        paste0("de.facLoc=", de.facLoc),
        paste0("iter=", iter),
        sep="__"
      ))
    
    cat("total number of scenario settings:", nrow(sim_grid), "\n")

    # remove settings that have been run before
    sim_grid <- sim_grid %>%
      filter(!file.exists(paste0("results/", scenario, "/results_simulation__", setting, ".Rdata")))

    cat("num scenarios left to run:", nrow(sim_grid), "\n")

    if(nrow(sim_grid) == 0) next

    if(reticulate::py_available(initialize=FALSE)){
      stop("Python started before parallelization!")
    }
    
    p <- progressor(steps=nrow(sim_grid))
    foreach(
      i=1:nrow(sim_grid),
      .combine="bind_rows",
      .export="StabCell",
      .options.future=list(chunk.size=1) # one task at a time for better workload balancing
    ) %dorng% {
      
      blas_set_num_threads(blas_threads)
      omp_set_num_threads(blas_threads)
      
      # no parallelization within the simulation setting (but resore old plan)
      old_plan <- plan()
      on.exit(plan(old_plan), add=TRUE)
      plan(sequential)


      de.prob <- sim_grid[i, "de.prob"]
      de.facLoc <- sim_grid[i, "de.facLoc"]
      n_celltypes <- sim_grid[i, "n_celltypes"]
      lib.size <- sim_grid[i, "lib.size"]
      iter <- sim_grid[i, "iter"]
      
      cluster_args <- list(n_clusters=n_clusters)

      setting <- sim_grid[i, "setting"]

      results_file <- paste0("results/", scenario, "/results_simulation__", setting, ".Rdata")


      set.seed(seed + iter)

      #######################
      #### simulate data ####
      #######################
      {
        t_start <- Sys.time()
        params <- newSplatParams(
          # see: https://bioconductor.org/packages/release/bioc/vignettes/splatter/inst/doc/splat_params.html

          seed=seed + iter,

          nGenes=n_genes,
          batchCells=n_cells,

          # cell types
          group.prob=rep(1. / n_celltypes, n_celltypes),

          # library size (= counts per cell)
          lib.loc=log(lib.size) - lib.scale^2 / 2,
          lib.scale=lib.scale,

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
        count_data <- as(counts(sim), "CsparseMatrix")
        # assay(sim, "Dropout")

        true_idents <- as.numeric(substr(colData(sim)$Group, 6, 10))
        if(length(true_idents) == 0){
          true_idents <- rep(1, ncol(count_data))
        }
        true_fold_changes <- as.data.frame(rowData(sim)[startsWith(colnames(rowData(sim)), "DEFacGroup")])
        # cat("Simulation:", format(Sys.time() - t_start), "\n")
      }

      
      
      ################################################################################
      ################################################################################
      ####                                                                        ####
      ####                              Oracle Markers                            ####
      ####                                                                        ####
      ################################################################################
      ################################################################################
      {
        t_start <- Sys.time()
        
        # marker detection with true idents
        unique_true_idents <- unique(true_idents)
        if(length(unique_true_idents) > 1){
          df_list <- list()
          for(oracle_id in unique_true_idents){
            oracle_cluster_markers <- marker_fn(
              data=count_data,
              identities=true_idents,
              weights=NULL,
              cluster_id1=oracle_id,
              cluster_id2=NULL,
              marker_args=marker_args
            )
            df_list[[length(df_list)+1]] <- data.frame(
              cluster1=oracle_id,
              cluster2=NA,
              marker=names(oracle_cluster_markers),
              pval=oracle_cluster_markers
            )
          }

          oracle_markers <- bind_rows(df_list)
          rownames(oracle_markers) <- NULL
        } else {
          oracle_markers <- data.frame(
            cluster1=unique_true_idents,
            cluster2=NA,
            marker=rownames(count_data),
            pval=1
          )
        }
        # cat("Oracle:", format(Sys.time() - t_start), "\n")
      }



      ################################################################################
      ################################################################################
      ####                                                                        ####
      ####                              Double Dipping                            ####
      ####                                                                        ####
      ################################################################################
      ################################################################################
      {
        t_start <- Sys.time()
        
        # clustering
        ddip_idents <- cluster_fn(data=count_data, subsample_indices=NULL, cluster_args=cluster_args)

        # marker detection
        df_list <- list()
        for(ddip_id in levels(ddip_idents)){
          ddip_cluster_markers <- marker_fn(
            data=count_data,
            identities=ddip_idents,
            weights=NULL,
            cluster_id1=ddip_id,
            cluster_id2=NULL,
            marker_args=marker_args
          )
          df_list[[length(df_list)+1]] <- data.frame(
            cluster1=ddip_id,
            cluster2=NA,
            marker=names(ddip_cluster_markers),
            pval=ddip_cluster_markers
          )
        }

        ddip_markers <- bind_rows(df_list)
        rownames(ddip_markers) <- NULL
        # cat("DDip:", format(Sys.time() - t_start), "\n")
      }



      ################################################################################
      ################################################################################
      ####                                                                        ####
      ####                                CountSplit                              ####
      ####                                                                        ####
      ################################################################################
      ################################################################################
      {
        t_start <- Sys.time()
        
        # split dataset
        split <- countsplit(
          t(count_data),
          epsilon=c(0.5, 0.5),
          overdisps=rep(Inf, nrow(count_data)) # for Poisson
        )
        data_clust <- t(split[[1]])
        data_marker <- t(split[[2]])

        # clustering on one split
        countsplit_idents <- cluster_fn(data=data_clust, subsample_indices=NULL, cluster_args=cluster_args)

        # marker detection on other split
        df_list <- list()
        for(countsplit_id in levels(countsplit_idents)){
          countsplit_cluster_markers <- marker_fn(
            data=data_marker,
            identities=countsplit_idents,
            weights=NULL,
            cluster_id1=countsplit_id,
            cluster_id2=NULL,
            marker_args=marker_args
          )
          df_list[[length(df_list)+1]] <- data.frame(
            cluster1=countsplit_id,
            cluster2=NA,
            marker=names(countsplit_cluster_markers),
            pval=countsplit_cluster_markers
          )
        }

        countsplit_markers <- bind_rows(df_list)
        rownames(countsplit_markers) <- NULL
        # cat("Countsplit full", format(Sys.time() - t_start), "\n")
      }


      
      ################################################################################
      ################################################################################
      ####                                                                        ####
      ####                               ClusterDE                                ####
      ####                                                                        ####
      ################################################################################
      ################################################################################
      {
        t_start <- Sys.time()
        
        # create null data and cluster it
        # --> assumes single cluster under H0, may not hold for 1-vs-rest testing
        null_data <- suppress_pois_check_warnings(suppressMessages(constructNull(
          obj=CreateSeuratObject(counts=count_data),
          family="poisson",
          n_cores=blas_threads,
          seed=seed+iter
        )))
        null_idents <- cluster_fn(data=null_data, subsample_indices=NULL, cluster_args=cluster_args)
        
        # real-null-cluster-mapping based on relative frequency
        # matches power for real and null scenario
        ddip_freq <- prop.table(table(ddip_idents))
        null_freq <- prop.table(table(null_idents))
        closest_idx <- vapply(
          ddip_freq,
          function(freq) which.min(abs(null_freq - freq)),
          integer(1)
        )
        ddip_null_mapping <- data.frame(
          ddip_id=transfer_dtype(x=names(ddip_freq), y=ddip_idents),
          null_id=transfer_dtype(x=names(null_freq)[closest_idx], y=null_idents)
        )
        
        # marker detection
        df_list <- list()
        for(ddip_id in levels(ddip_idents)){
          
          # testing on null data
          null_id <- ddip_null_mapping %>% filter(ddip_id == !!ddip_id) %>% pull(null_id)
          null_markers <- marker_fn(
            data=null_data,
            identities=null_idents,
            weights=NULL,
            cluster_id1=null_id,
            cluster_id2=NULL,
            marker_args=marker_args
          )
          
          # error control step
          withr::with_options(list(progress_enabled=FALSE), {
                              
            df_list[[length(df_list)+1]] <- callDE(
              targetScores=ddip_markers %>%
                filter(cluster1 == ddip_id) %>%
                pull(pval, name=marker),
              nullScores=null_markers
            ) %>%
              mutate(
                cluster1=ddip_id,
                cluster2=NA,
                marker=gene,
                pval=(-1 * (record >= 0.5)) + 1, # turn binary marker indicator into "pval"
                .keep="none"
              )
            
          })
        }
        
        clusterde_markers <- bind_rows(df_list)
        rownames(clusterde_markers) <- NULL
        # cat("ClusterDE:", format(Sys.time() - t_start), "\n")
      }
      


      ################################################################################
      ################################################################################
      ####                                                                        ####
      ####                                TN test                                 ####
      ####                                                                        ####
      ################################################################################
      ################################################################################
      {
        t_start <- Sys.time()
        
        # intialize Python
        reticulate::use_virtualenv("./.venv", required=TRUE)
        source("Scripts/helpers/python.R")
        python_blas(1L) # parallelize via "n_cores" below
        numpy <- reticulate::import("numpy", convert=TRUE)
        numpy$random$seed(as.integer(seed + iter))

        # split data
        idx_split <- sample.int(n=n_cells, size=n_cells %/% 2, replace=FALSE)
        counts_clust <- count_data[ ,  idx_split]
        counts_test <- count_data[ , -idx_split]

        # cluster on clustering split
        tntest_idents_clust <- cluster_fn(data=counts_clust, subsample_indices=NULL, cluster_args=cluster_args)
        tntest_ids <- unique(as.character(tntest_idents_clust))

        tntest_idents <- rep(NA_character_, n_cells)
        tntest_idents[idx_split] <- as.character(tntest_idents_clust)

        # normalize data (required for Gaussian assumption)
        X_clust <- CreateSeuratObject(counts=counts_clust) %>%
          NormalizeData(normalization.method="LogNormalize", scale.factor=10000, verbose=FALSE) %>%
          LayerData(layer="data") %>%
          as.matrix() %>%
          t()
        X_test <- CreateSeuratObject(counts=counts_test) %>%
          NormalizeData(normalization.method="LogNormalize", scale.factor=10000, verbose=FALSE) %>%
          LayerData(layer="data") %>%
          as.matrix() %>%
          t()

        # TN test 1-vs-rest (with reticulate for Python code)
        # approximation since pooled "rest" is most likely not single Gaussian, but multimodal
        success <- FALSE
        df_list <- list()
        tryCatch(
          {
            for(tntest_id in tntest_ids){
              # construct binary label and class weights
              idx_target_clust <- tntest_idents_clust == tntest_id
              labels_clust <- factor(ifelse(idx_target_clust, "target", "rest"), levels=c("rest", "target"))
    
              n_total <- length(labels_clust)
              n_target <- sum(idx_target_clust)
              n_rest <- n_total - n_target
              class_weights <- c(
                rest=n_total / (2 * n_rest),
                target=n_total / (2 * n_target)
              )
    
              # fit single class svm and recover hyperplane
              # decision(x)=X %*% w + intercept
              singleclass_svm <- e1071::svm(
                x=X_clust,
                y=labels_clust,
                type="C-classification",
                kernel="linear",
                cost=100,
                scale=FALSE,
                class.weights=class_weights
              )
              w <- drop(t(singleclass_svm$coefs) %*% singleclass_svm$SV)
              intercept <- -as.numeric(singleclass_svm$rho)
    
              # hyper plane orientation
              pred_target <- as.character(predict(singleclass_svm, newdata=X_clust)) == "target"
              score_clust <- X_clust %*% w + intercept
              agreement <- mean((score_clust > 0 ) == pred_target)
              agreement_flipped <- mean((score_clust <= 0 ) == pred_target)
              if(agreement_flipped > agreement){
                w <- -1 * w
                intercept <- -1 * intercept
              }
              score_test <- drop(X_test %*% w + intercept)
              idx_target_test <- score_test > 0
    
              # perform TN test in Python
              pval <- tn_python(inner_threads=1L)$tn_test(
                  y=X_test[!idx_target_test, , drop=FALSE],
                  z=X_test[idx_target_test, , drop=FALSE],
                  a=numpy$asarray(w),
                  b=intercept,
                  verbose=FALSE,
                  num_iters=100000L,
                  num_cores=as.integer(blas_threads)
                )
    
              df_list[[length(df_list) + 1]] <- data.frame(
                cluster1=tntest_id,
                cluster2=NA,
                marker=colnames(X_test),
                pval=p.adjust(pval, method="bonferroni")
              )
            }
            success <- TRUE
          },
          finally=try(shutdown_joblib(kill_workers=!success), silent=TRUE)
        )
        tntest_markers <- bind_rows(df_list)
        # cat("TN Test:"; format(Sys.time() - t_start), "\n")
      }



      ################################################################################
      ################################################################################
      ####                                                                        ####
      ####                                StabCell                                ####
      ####                                                                        ####
      ################################################################################
      ################################################################################

      without_progress({
        t_start <- Sys.time()
        stabcell_results <- StabCell(
          data=count_data,
          cluster_fn=cluster_fn,
          PFER_max=PFER_max,
          reference_identities=ddip_idents,
          marker_fn=marker_fn,
          p_value_threshold=p_value_threshold,
          n_subsample_markers=n_subsample_markers,
          cluster_args=cluster_args,
          marker_args=marker_args,
          marker_global=TRUE,
          marker_local=FALSE,
          subsampling_iterations=subsampling_iterations,
          subsampling_fraction=0.5,
          complementary_pairs=TRUE,
          cell_stability_metric=cell_stability_metric,
          assumption=c("none", "unimodal", "r_concavity"),
          ranking="p_val",
          theta=NULL,
          cells_in_columns=TRUE,
          packages=c("Seurat"),        # packages that are needed for clustering or testing
          export=c("FindNClusters"),   # additional custom objects needed for clustering or testing
          steps_to_run=NULL,
          checkpoint_path=NULL,
          checkpoint=NULL,
          verbose=FALSE
        )
        # cat("StabCell:", format(Sys.time() - t_start), "\n")
      })

      ######################
      #### save results ####
      ######################

      results <- list(
        config=list(
          n_celltypes=n_celltypes,
          lib.size=lib.size,
          de.prob=de.prob,
          de.facLoc=de.facLoc,
          n_iterations=subsampling_iterations
        ),
        sim=list(
          count_data=count_data,
          true_idents=true_idents,
          true_fold_changes=true_fold_changes
        ),
        oracle=list(
          markers=oracle_markers
        ),
        double_dipping=list(
          markers=ddip_markers,
          identities=ddip_idents
        ),
        countsplit=list(
          markers=countsplit_markers,
          identities=countsplit_idents
        ),
        clusterde=list(
          markers=clusterde_markers
        ),
        tntest=list(
          markers=tntest_markers,
          identities=tntest_idents
        ),
        stabcell=stabcell_results
      )

      save(results, file=results_file)

      p()
      return(NULL)
    }

    cat(paste0("finished at ", strftime(Sys.time())), "\n")
  }
})

print(summary(warnings()))
