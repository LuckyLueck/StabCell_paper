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

library(StabCell)
source("Scripts/helpers/cluster_fn.R")
source("Scripts/helpers/marker_fn.R")
source("Scripts/helpers/utils.R")

################
#### config ####
################
blas_threads <- 4L
num_workers <- 25L

options(future.globals.maxSize=12 * 1024^3)

blas_set_num_threads(blas_threads)
omp_set_num_threads(blas_threads)


seed <- 22032026

n_cells <- c(500, 1000, 5000, 10000, 20000, 35000, 50000)
n_celltypes <- 2
cluster_args <- list(n_clusters=2)

PFER_max <- 25
assumption <- "r_concavity"

repetitions <- 25


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
  filter(!file.exists(paste0("results/runtime_study/results_runtime__", setting, ".Rdata")))

n_settings <- nrow(sim_grid)

cat("num settings left to run:", n_settings, "\n")


# intialize Python
reticulate::use_virtualenv("./.venv", required=TRUE)
source("Scripts/helpers/python.R")
python_blas(blas_threads)
numpy <- reticulate::import("numpy", convert=TRUE)
tn_test <- tn_python(inner_threads=blas_threads)$tn_test

with_progress({

  p <- progressor(steps=n_settings)

  for(i in 1:n_settings){

    t_start <- Sys.time()
    
    plan(sequential)
    registerDoSEQ()

    n_cells <- sim_grid[i, "n_cells"]
    iter <- sim_grid[i, "iter"]
    setting <- sim_grid[i, "setting"]

    results_file <- paste0(
      "results/runtime_study/results_runtime__", setting, ".Rdata"
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
        lib.loc=9.23, # lib.size ~= 10'000
        lib.scale=0.2,

        # dropout probabilities
        dropout.type="experiment", # options: "none", "experiment", "batch, "group", "cell"
        dropout.mid=50,
        dropout.shape=5e-2, # smaller -> more dropouts ?

        # differential expression
        de.prob=0.05,  # probability of differentially expressed gene, can be set differently for each group
        de.downProb=0.5,  # probability of downregulation vs upregulation
        de.facLoc=0.4,    # location of logNorm distributed DE factors
        de.facScale=0.4,  # scale of logNorm distributed DE factors

      )

      sim <- splatSimulate(
        params,
        method="groups",
        verbose=FALSE
      )
      count_data <- as(counts(sim), "CsparseMatrix")
    }
    message(paste0("setting ", i, "/", n_settings, ": finished data simulation"))

    ##################
    #### Run DDip ####
    ##################
    {
      t_start_ddip <- Sys.time()

      ddip_idents <- cluster_fn(data=count_data, subsample_indices=NULL, cluster_args=cluster_args)
      
      t_clust_ddip <- Sys.time()
      
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
      
      markers_ddip <- bind_rows(df_list)
      rownames(markers_ddip) <- NULL

      t_marker_ddip <- Sys.time()

      duration_clust_ddip <- t_clust_ddip - t_start_ddip
      duration_marker_ddip <- t_marker_ddip - t_clust_ddip
      duration_total_ddip <- t_marker_ddip - t_start_ddip
    }
    message(paste0("setting ", i, "/", n_settings, ": finished double dipping (", round(as.numeric(duration_total_ddip, units="mins"), 2), "mins)"))

    
    
    ################################################################################
    ################################################################################
    ####                                                                        ####
    ####                                CountSplit                              ####
    ####                                                                        ####
    ################################################################################
    ################################################################################
    {
      t_start_countsplit <- Sys.time()
      
      theta <- rep(Inf, nrow(count_data))  # Poisson assumption
      
      split <- countsplit(
        t(count_data),
        epsilon=c(0.5, 0.5),
        overdisps=theta
      )
      data_clust <- t(split[[1]])
      data_marker <- t(split[[2]])
      
      t_split_countsplit <- Sys.time()
      
      # clustering on one split
      countsplit_idents <- cluster_fn(data=data_clust, subsample_indices=NULL, cluster_args=cluster_args)
      
      t_clust_countsplit <- Sys.time()
      
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
      
      markers_countsplit <- bind_rows(df_list)
      rownames(markers_countsplit) <- NULL

      t_marker_countsplit <- Sys.time()

      duration_split_countsplit <- t_split_countsplit - t_start_countsplit
      duration_clust_countsplit <- t_clust_countsplit - t_split_countsplit
      duration_marker_countsplit <- t_marker_countsplit - t_clust_countsplit
      duration_total_countsplit <- t_marker_countsplit - t_start_countsplit
    }
    message(paste0("setting ", i, "/", n_settings, ": finished count splitting (", round(as.numeric(duration_total_countsplit, units="mins"), 2), "mins)"))

    
    
    ################################################################################
    ################################################################################
    ####                                                                        ####
    ####                               ClusterDE                                ####
    ####                                                                        ####
    ################################################################################
    ################################################################################
    {
      t_start_clusterde <- Sys.time()
      
      # create null data and cluster it
      null_data <- suppress_pois_check_warnings(suppressMessages(constructNull(
        obj=CreateSeuratObject(counts=count_data),
        family="poisson",
        n_cores=num_workers,
        seed=seed+iter
      )))
      if(any(!is.finite(null_data))) stop("non finite value in `null_data`!")
      null_idents <- cluster_fn(data=null_data, subsample_indices=NULL, cluster_args=cluster_args)
      
      t_clust_clusterde <- Sys.time()
      
      # real-null-cluster-mapping based on relative frequency
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
      
      df_list <- list()
      for(ddip_id in levels(ddip_idents)){
        null_id <- ddip_null_mapping %>% filter(ddip_id == !!ddip_id) %>% pull(null_id)
        null_markers <- marker_fn(
          data=null_data,
          identities=null_idents,
          weights=NULL,
          cluster_id1=null_id,
          cluster_id2=NULL,
          marker_args=marker_args
        )
        
        withr::with_options(list(progress_enabled=FALSE), {
          
          df_list[[length(df_list)+1]] <- callDE(
            targetScores=markers_ddip %>%
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
      
      markers_clusterde <- bind_rows(df_list)
      rownames(markers_clusterde) <- NULL
      
      t_marker_clusterde <- Sys.time()
      
      duration_clust_clusterde <- duration_clust_ddip + (t_clust_clusterde - t_start_clusterde)
      duration_marker_clusterde <- duration_marker_ddip + (t_marker_clusterde - t_clust_clusterde)
      duration_total_clusterde <- duration_clust_clusterde + duration_marker_clusterde
    }
    message(paste0("setting ", i, "/", n_settings, ": finished ClusterDE (", round(as.numeric(duration_total_clusterde, units="mins"), 2), "mins)"))
    
    ################################################################################
    ################################################################################
    ####                                                                        ####
    ####                                TN test                                 ####
    ####                                                                        ####
    ################################################################################
    ################################################################################
    {
      t_start_tntest <- Sys.time()
      
      # split data
      numpy$random$seed(as.integer(seed + iter))
      idx_split <- sample.int(n=n_cells, size=n_cells %/% 2, replace=FALSE)
      counts_clust <- count_data[ ,  idx_split]
      counts_test <- count_data[ , -idx_split]
      
      t_split_tntest <- Sys.time()
      
      # cluster on clustering split
      tntest_idents_clust <- cluster_fn(data=counts_clust, subsample_indices=NULL, cluster_args=cluster_args)
      tntest_ids <- unique(as.character(tntest_idents_clust))
      
      tntest_idents <- rep(NA_character_, n_cells)
      tntest_idents[idx_split] <- as.character(tntest_idents_clust)
      
      t_clust_tntest <- Sys.time()
      
      # normalize data
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
      # approximation since pooled "rest" is most likely not Gaussian, but multimodal
      df_list <- list()
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
        pval <- tn_test(
          y=X_test[!idx_target_test, , drop=FALSE],
          z=X_test[idx_target_test, , drop=FALSE],
          a=numpy$asarray(w),
          b=intercept,
          verbose=FALSE,
          num_iters=100000L,
          num_cores=num_workers
        )
        
        df_list[[length(df_list) + 1]] <- data.frame(
          cluster1=tntest_id,
          cluster2=NA,
          marker=colnames(X_test),
          pval=p.adjust(pval, method="bonferroni")
        )
      }
      markers_tntest <- bind_rows(df_list)
      row.names(markers_tntest) <- NULL

      t_marker_tntest <- Sys.time()
      
      duration_split_tntest <- t_split_tntest - t_start_tntest
      duration_clust_tntest <- t_clust_tntest - t_split_tntest
      duration_marker_tntest <- t_marker_tntest - t_clust_tntest
      duration_total_tntest <- t_marker_tntest - t_start_tntest
    }
    message(paste0("setting ", i, "/", n_settings, ": finished TN test (", round(as.numeric(duration_total_tntest, units="mins"), 2), "mins)"))
    
    
    
    ################################################################################
    ################################################################################
    ####                                                                        ####
    ####                                StabCell                                ####
    ####                                                                        ####
    ################################################################################
    ################################################################################
    plan(multicore, workers=num_workers)
    registerDoFuture()
    {
      t_start_stabcell <- Sys.time()
      
      without_progress({
        stabcell_results <- StabCell(
          data=count_data,
          cluster_fn=cluster_fn,
          PFER_max=PFER_max,
          reference_identities=ddip_idents,
          marker_fn=marker_fn,
          p_value_threshold=0.05,
          n_subsample_markers=NULL,
          cluster_args=cluster_args,
          marker_args=marker_args,
          marker_global=TRUE,
          marker_local=FALSE,
          subsampling_iterations=50,
          subsampling_fraction=0.5,
          complementary_pairs=TRUE,
          cell_stability_metric="gini",
          assumption=assumption,
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
      })
      
      t_end_stabcell <- Sys.time()

      duration_refclust_stabcell <- duration_clust_ddip
      duration_subsamples_stabcell <- t_end_stabcell - t_start_stabcell
      duration_total_stabcell <- duration_refclust_stabcell + duration_subsamples_stabcell
    }
    message(paste0("setting ", i, "/", n_settings, ": finished StabCell (", round(as.numeric(duration_total_stabcell, units="mins"), 2), "mins)"))

    ######################
    #### save results ####
    ######################
    {
      runtimes <- list(
        n_cells=n_cells,
        iter=iter,
        ddip=list(
          duration_clust=duration_clust_ddip,
          duration_marker=duration_marker_ddip,
          duration_total=duration_total_ddip
        ),
        countsplit=list(
          duration_split=duration_split_countsplit,
          duration_clust=duration_clust_countsplit,
          duration_marker=duration_marker_countsplit,
          duration_total=duration_total_countsplit
        ),
        clusterde=list(
          duration_clust=duration_clust_clusterde,
          duration_marker=duration_marker_clusterde,
          duration_total=duration_total_clusterde
        ),
        tntest=list(
          duration_split=duration_split_tntest,
          duration_clust=duration_clust_tntest,
          duration_marker=duration_marker_tntest,
          duration_total=duration_total_tntest
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
    message(paste0("setting ", i, "/", n_settings, ": total time = ", round(as.numeric(Sys.time() - t_start, units="mins"), 2), "mins"))
    p()

  }
})


cat(paste0("finished at ", strftime(Sys.time()), "\n"))

print(summary(warnings()))


