setwd("/mnt/mikasa_data/nlueck/StabCell_paper")

library(doFuture)
library(dplyr)
library(RhpcBLASctl)
library(mclust)
library(progressr)

source("Scripts/helpers/metrics.R")

blas_threads <- 2
num_workers <- 100

options(future.globals.maxSize=4 * 1024^3)

blas_set_num_threads(blas_threads)
plan(multicore, workers=num_workers)
registerDoFuture()

p_val_thresh_alt <- 0.05 # threshold to select markers in alternative approaches (double dipping, countsplit)

scenarios <- c(
  "2Cluster_Simulation",
  "NCluster_Simulation"
)[1]


find_closest_cluster <- StabCell:::find_closest_cluster

for(scenario in scenarios){

  cat(
    "\n",
    paste(rep("#", 14 + nchar(scenario)), collapse=""), "\n",
    "#### run", scenario, "####\n",
    paste(rep("#", 14 + nchar(scenario)), collapse=""), "\n"
  )
  cat(paste0("start at ", strftime(Sys.time())), "\n")

  files <- list.files(paste0("results/", scenario), pattern="^results_simulation__.*\\.Rdata$")

  cat("num files: ", length(files), "\n")


  with_progress({

    p <- progressor(along=files)

    df_list <- foreach(
      x=1:length(files),
      .packages=c("dplyr", "foreach", "tidyr", "Matrix")
    ) %dopar% {

      ###################
      #### load data ####
      ###################

      file <- files[x]
      out_file <- paste0("results/", scenario, "/metrics/", gsub("results_simulation_", "results_metrics_", file))

      # load metrics if already calculated
      if(file.exists(out_file)){

        load(out_file)

        p()
        return( # exit foreach iteration early
          list(
            marker_metrics=marker_metrics,
            clust=df_clust,
            counts=df_counts,
            auc=auc,
            skip=x,
            error=NA
          )
        )
      }


      # try to load simulation results, exit foreach iter if it fails
      results <- tryCatch(
        expr={
          load(paste0("results/", scenario, "/", file))
          results
        },
        error=function(e){
          message(paste("Failed to load file: ", x, "with error:"))
          message(conditionMessage(e))
          NULL
        }
      )
      if(is.null(results)){
        p()
        return(
          list(
            marker_metrics=list(
              ddip=list(
                global=data.frame(),
                truth_alignment=data.frame(),
                stabcell_alignment=data.frame()
              ),
              countsplit=list(
                global=data.frame(),
                truth_alignment=data.frame(),
                stabcell_alignment=data.frame()
              ),
              clusterde=list(
                global=data.frame(),
                truth_alignment=data.frame(),
                stabcell_alignment=data.frame()
              ),
              tntest=list(
                global=data.frame(),
                truth_alignment=data.frame(),
                stabcell_alignment=data.frame()
              ),
              oracle=list(
                global=data.frame(),
                truth_alignment=data.frame(),
                stabcell_alignment=data.frame()
              ),
              stabcell=list(
                global=data.frame(),
                truth_alignment=data.frame(),
                stabcell_alignment=data.frame()
              )
            ),
            clust=data.frame(),
            counts=data.frame(),
            auc=data.frame(),
            skip=NA,
            error=x
          )
        )
      }

      # unpack config data
      n_celltypes <- results$config$n_celltypes
      lib.size <- results$config$lib.size
      de.prob <- results$config$de.prob
      de.facLoc <- results$config$de.facLoc
      subsampling_iterations <- results$config$n_iterations
      iter <- as.numeric(sub(".*__iter=([0-9]*)\\.Rdata", "\\1", file))

      # unpack simulation results
      count_data <- results$sim$count_data
      true_idents <- results$sim$true_idents %>% as.character()
      true_fold_changes <- results$sim$true_fold_changes
      
      # oracle results
      oracle_markers <- results$oracle$markers %>% mutate(cluster1=as.character(cluster1))
      
      # double dipping results
      ddip_idents <- results$double_dipping$identities %>% as.character()
      ddip_markers <- results$double_dipping$markers %>% mutate(cluster1=as.character(cluster1))

      # countsplit results
      countsplit_idents <- results$countsplit$identities %>% as.character()
      countsplit_markers <- results$countsplit$markers %>% mutate(cluster1=as.character(cluster1))
      
      # clusterde results
      clusterde_markers <- results$clusterde$markers %>% mutate(cluster1=as.character(cluster1))
      
      # tntest results
      tntest_idents <- results$tntest$identities %>% as.character()
      tntest_markers <- results$tntest$markers %>% mutate(cluster1=as.character(cluster1))

      # stabcell results
      subsampling_identities <- results$stabcell$subsampling_identities
      clustering_alignment <- results$stabcell$clustering_alignment
      stabcell_idents <- results$stabcell$stable_clustering %>% as.character()
      cell_stability_score <- if(ncol(results$stabcell$assignment_probabilities) == 1) rep(1, nrow(results$stabcell$assignment_probabilities)) else results$stabcell$cell_stability_score
      marker_alignment <- results$stabcell$marker_alignment
      stable_markers_global <- results$stabcell$stable_markers_global %>% mutate(id=as.character(id))
      

      ##########################
      #### evaluate markers ####
      ##########################

      # find true marker genes
      if(n_celltypes == 1){

        # no true marker genes
        true_markers_global <- data.frame(id=numeric(), marker=character(), lfc=numeric())

      } else if(n_celltypes == 2){

        # with only 2 clusters, LFC != 0 from the other cluster is
        # automatically also marker for the cluster under investigation
        marker_idx <- (true_fold_changes[, "DEFacGroup1"] != true_fold_changes[, "DEFacGroup2"])
        lfc <- log2(true_fold_changes[marker_idx, "DEFacGroup1"] / true_fold_changes[marker_idx, "DEFacGroup2"])

        true_markers <- data.frame(
          id1="1",
          id2="2",
          marker=rownames(true_fold_changes)[marker_idx],
          lfc=lfc
        )
        true_markers_global <- bind_rows(
          true_markers %>% select(-id2) %>% rename(id=id1),
          true_markers %>% select(-id1) %>% rename(id=id2) %>% mutate(lfc=-1 * lfc)
        )

      } else { # n_celltypes > 2
        
        # with more than 2 clusters, LFC != 0 from the other clusters is not treated
        # as a marker as it will be "dilluted" by pooling the other clusters
        true_idents_unique <- unique(true_idents)
        n_true_idents <- length(true_idents_unique)

        # find true global markers
        df_list <- list()
        for(i in 1:n_true_idents){
          id <- true_idents_unique[i]

          marker_idx <- true_fold_changes[ , paste0("DEFacGroup", id)] != 1.
          lfc <- log2(true_fold_changes[marker_idx, paste0("DEFacGroup", id)])

          df_list[[length(df_list) + 1]] <- data.frame(
            id=id,
            marker=rownames(true_fold_changes)[marker_idx],
            lfc=lfc
          )
        }
        true_markers_global <- bind_rows(df_list)

      }


      marker_metrics <- list()
      for(method in c("oracle", "ddip", "countsplit", "clusterde", "tntest", "stabcell")){

        if(method == "oracle"){
          method_idents <- true_idents
          method_markers <- oracle_markers
        } else if(method == "ddip"){
          method_idents <- ddip_idents
          method_markers <- ddip_markers
        } else if(method == "countsplit"){
          method_idents <- countsplit_idents
          method_markers <- countsplit_markers
        } else if(method == "clusterde"){
          method_idents <- ddip_idents  # use ddip identities here!
          method_markers <- clusterde_markers
        } else if(method == "tntest"){
          method_idents <- tntest_idents  # only half of cells available!
          method_markers <- tntest_markers
        } else if(method == "stabcell"){
          method_idents <- stabcell_idents
        } else {
          stop("Unknown method!")
        }

        # align method clusters to true clusters
        valid_idents <- !is.na(method_idents)
        method_idents_unique <- unique(method_idents[valid_idents])
        simulation_ids <- rep(NA_character_, length(method_idents_unique))

        for(i in seq_along(method_idents_unique)){
          simulation_ids[i] <- find_closest_cluster(
            reference_allocations=method_idents[valid_idents],
            reference_cluster_id=method_idents_unique[i],
            comparison_allocations=true_idents[valid_idents]
          )
        }
        df_align_simulation <- data.frame(
          method_id=method_idents_unique,
          simulation_id_aligned=simulation_ids
        )

        # calculate marker metrics
        if(method == "stabcell"){

          # global marker metrics
          settings_global <- stable_markers_global %>%
            group_by(p_val_thresh, q, PFER_max, assumption, selection_freq_thresh, exceeds_one, id) %>%
            summarize(.groups="drop")
          
          if(nrow(settings_global) == 0) return(data.frame())

          df_metrics_global <- foreach(
            i=1:nrow(settings_global),
            .combine="bind_rows",
            .packages=c("dplyr")
          ) %do% {
            q <- settings_global$q[i]
            p_val_thresh <- settings_global$p_val_thresh[i]
            PFER_max <- settings_global$PFER_max[i]
            assumption <- settings_global$assumption[i]
            selection_freq_thresh <- settings_global$selection_freq_thresh[i]
            exceeds_one <- settings_global$exceeds_one[i]
            stabcell_id <- settings_global$id[i]
            simulation_id <- df_align_simulation[df_align_simulation$method_id == stabcell_id, "simulation_id_aligned"]

            selected_markers_df <- (
              stable_markers_global %>%
                filter(
                  q == !!q,
                  p_val_thresh == !!p_val_thresh,
                  PFER_max == !!PFER_max,
                  assumption == !!assumption,
                  selection_freq_thresh == !!selection_freq_thresh,
                  id == !!stabcell_id,
                  !is.na(marker)
                )
            )
            unselected_markers <- setdiff(rownames(true_fold_changes), selected_markers_df$marker)

            true_markers_df <- (true_markers_global %>% filter(id == !!simulation_id))
            false_markers <- setdiff(rownames(true_fold_changes), true_markers_df$marker)

            df_temp <- cbind(
              n_celltypes=n_celltypes,
              n_iter=subsampling_iterations,
              lib.size=lib.size,
              de.prob=de.prob,
              de.facLoc=de.facLoc,
              iter=iter,

              q=q,
              p_val_thresh=p_val_thresh,
              PFER_max=PFER_max,
              assumption=assumption,
              selection_freq_thresh=selection_freq_thresh,
              exceeds_one=exceeds_one,

              id=stabcell_id,
              truth_id_aligned=simulation_id,

              calculate_metrics(
                selected_markers=selected_markers_df$marker,
                unselected_markers=unselected_markers,
                true_markers=true_markers_df$marker,
                false_markers=false_markers
              ),
              mean_cell_stability_score=mean(cell_stability_score[stabcell_idents == stabcell_id]),

              top1_agreement=top_k_marker_agreement(
                true_log_fold_changes=true_markers_df,
                method_marker_df=selected_markers_df,
                k=1
              ),
              top5_agreement=top_k_marker_agreement(
                true_log_fold_changes=true_markers_df,
                method_marker_df=selected_markers_df,
                k=5
              ),
              top10_agreement=top_k_marker_agreement(
                true_log_fold_changes=true_markers_df,
                method_marker_df=selected_markers_df,
                k=10
              )
            )


            return(df_temp)

          }

          df_metrics_global <- df_metrics_global %>% mutate(
            n_iter=as.factor(n_iter),
            lib.size=as.factor(lib.size),
            de.prob=as.factor(de.prob),
            de.facLoc=as.factor(de.facLoc),
            PFER_max=as.factor(PFER_max),
            assumption=as.factor(assumption)
          )

        } else {  # other methods

          # find alignment to stabcell clustering
          stabcell_idents_unique <- unique(stabcell_idents)
          method_ids <- rep(NA_character_, length(stabcell_idents_unique))
          for(i in seq_along(stabcell_idents_unique)){
            method_ids[i] <- find_closest_cluster(
              reference_allocations=stabcell_idents[valid_idents],
              reference_cluster_id=stabcell_idents_unique[i],
              comparison_allocations=method_idents[valid_idents]
            )
          }
          df_align_stabcell <- data.frame(
            stabcell_id=stabcell_idents_unique,
            method_id_aligned=method_ids
          )

          # global marker metrics
          df_metrics_global <- foreach(
            i=1:nrow(df_align_simulation),
            .combine="bind_rows",
            .packages=c("dplyr")
          ) %do% {
            method_id <- df_align_simulation[i, "method_id"]
            simulation_id <- df_align_simulation[i, "simulation_id_aligned"]

            selected_markers_df <- (
              method_markers %>%
                filter(
                  cluster1 == !!method_id,
                  is.na(cluster2),
                  pval <= !!p_val_thresh_alt
                )
            )
            unselected_markers <- setdiff(rownames(true_fold_changes), selected_markers_df$marker)

            true_markers_df <- (true_markers_global %>% filter(id == !!simulation_id))
            false_markers <- setdiff(rownames(true_fold_changes), true_markers_df$marker)

            df_temp <- cbind(
              n_celltypes=n_celltypes,
              n_iter=subsampling_iterations,
              lib.size=lib.size,
              de.prob=de.prob,
              de.facLoc=de.facLoc,
              iter=iter,

              id=method_id,
              truth_id_aligned=simulation_id,

              calculate_metrics(
                selected_markers=selected_markers_df$marker,
                unselected_markers=unselected_markers,
                true_markers=true_markers_df$marker,
                false_markers=false_markers
              ),

              top1_agreement=top_k_marker_agreement(
                true_log_fold_changes=true_markers_df,
                method_marker_df=selected_markers_df,
                k=1
              ),
              top5_agreement=top_k_marker_agreement(
                true_log_fold_changes=true_markers_df,
                method_marker_df=selected_markers_df,
                k=5
              ),
              top10_agreement=top_k_marker_agreement(
                true_log_fold_changes=true_markers_df,
                method_marker_df=selected_markers_df,
                k=10
              )
            )

            df_temp

          }

          df_metrics_global <- df_metrics_global %>%
            mutate(
              n_iter=as.factor(n_iter),
              lib.size=as.factor(lib.size),
              de.prob=as.factor(de.prob),
              de.facLoc=as.factor(de.facLoc)
            )

          # transfer stabcell ids
          df_metrics_global <- df_metrics_global %>%
            full_join(df_align_stabcell, by=c("id"="method_id_aligned")) %>%
            rename(
              method_id=id,
              id=stabcell_id
            )

        }

        marker_metrics[[method]] <- list(
          global=df_metrics_global,
          truth_alignment=df_align_simulation,
          stabcell_alignment=if(method != "stabcell") df_align_stabcell else NULL
        )

      }
      
      #######################
      #### AUROC / AUPRC ####
      #######################
      
      if(scenario == "2Cluster_Simulation" & n_celltypes > 1){
        scores <- results$sim$true_fold_changes %>%
          mutate(is_marker=as.integer(DEFacGroup1 != DEFacGroup2)) %>%
          tibble::rownames_to_column(var="gene") %>%
          select(-DEFacGroup1, -DEFacGroup2) %>%
          full_join(
            results$double_dipping$markers %>%
              filter(cluster1 == "0") %>%
              mutate(pval=-1 * pval + 1) %>%
              select(marker, pval) %>%
              dplyr::rename(`Double Dipping`=pval),
            by=c("gene"="marker")
          ) %>%
          full_join(
            results$countsplit$markers %>%
              filter(cluster1 == "0") %>%
              mutate(pval=-1 * pval + 1) %>%
              select(marker, pval) %>%
              dplyr::rename(Countsplit=pval),
            by=c("gene"="marker")
          ) %>%
          full_join(
            results$tntest$markers %>%
              filter(cluster1 == "0") %>%
              mutate(pval=-1 * pval + 1) %>%
              select(marker, pval) %>%
              dplyr::rename(`TN Test`=pval),
            by=c("gene"="marker")
          ) %>%
          full_join(
            results$stabcell$marker_selection_frequencies %>%
              filter(
                id1 == "0",
                p_val_thresh == 0.05,
              ) %>%
              select(marker, selection_freq) %>%
              dplyr::rename(StabCell=selection_freq),
            by=c("gene"="marker")
          )
        
        curves <- precrec::evalmod(
          scores=scores %>% select(-is_marker, -gene),
          labels=scores$is_marker,
          modnames=c("Double Dipping", "Countsplit", "TN Test", "StabCell"),
          posclass=1L
        )
        
        auc <- precrec::auc(curves) %>%
          select(-dsids) %>%
          mutate(
            n_celltypes=n_celltypes,
            n_iter=subsampling_iterations,
            lib.size=lib.size,
            de.prob=de.prob,
            de.facLoc=de.facLoc,
            iter=iter
          )
        
      } else {
        auc <- data.frame()
      }
      
      

      ############################
      #### clustering metrics ####
      ############################
      valid_tntest_idx <- !is.na(tntest_idents)
      
      df_clust <- data.frame(
        n_celltypes=n_celltypes,
        n_iter=subsampling_iterations,
        lib.size=lib.size,
        de.prob=de.prob,
        de.facLoc=de.facLoc,
        iter=iter,
        n_clusters_true=length(unique(true_idents)),
        n_clusters_ddip=length(unique(ddip_idents)),
        n_clusters_countsplit=length(unique(countsplit_idents)),
        n_clusters_tntest=length(unique(tntest_idents[valid_tntest_idx])),
        n_clusters_stabcell=length(unique(stabcell_idents)),
        ari_ddip=adjustedRandIndex(true_idents, ddip_idents),
        ari_countsplit=adjustedRandIndex(true_idents, countsplit_idents),
        ari_tntest=adjustedRandIndex(true_idents[valid_tntest_idx], tntest_idents[valid_tntest_idx]),
        ari_stabcell=adjustedRandIndex(true_idents, stabcell_idents),
        mean_cell_stability=mean(cell_stability_score)
      )


      ################
      #### counts ####
      ################
      cell_counts <- colSums(count_data)

      df_counts <- data.frame(
        n_celltypes=n_celltypes,
        n_iter=subsampling_iterations,
        lib.size=lib.size,
        de.prob=de.prob,
        de.facLoc=de.facLoc,
        iter=iter,
        cell_counts_mean=mean(cell_counts),
        cell_counts_sd=sd(cell_counts)
      )


      ######################
      #### save results ####
      ######################

      save(
        marker_metrics,
        df_clust,
        df_counts,
        auc,

        file=out_file
      )

      p()
      return(
        list(
          marker_metrics=marker_metrics,
          clust=df_clust,
          counts=df_counts,
          auc=auc,
          skip=NA,
          error=NA
        )
      )
    }
  })

  # combine iteration metrics into one data frame
  skipped_files <- numeric()
  error_files <- numeric()
  list_stabcell_global <- list()
  list_ddip_global <- list()
  list_countsplit_global <- list()
  list_clusterde_global <- list()
  list_tntest_global <- list()
  list_oracle_global <- list()
  list_clust <- list()
  list_counts <- list()
  list_auc <- list()

  for(x in 1:length(df_list)){
    list_stabcell_global[[length(list_stabcell_global) + 1]] <- df_list[[x]]$marker_metrics$stabcell$global
    list_ddip_global[[length(list_ddip_global) + 1]] <- df_list[[x]]$marker_metrics$ddip$global
    list_countsplit_global[[length(list_countsplit_global) + 1]] <- df_list[[x]]$marker_metrics$countsplit$global
    list_clusterde_global[[length(list_countsplit_global) + 1]] <- df_list[[x]]$marker_metrics$clusterde$global
    list_tntest_global[[length(list_countsplit_global) + 1]] <- df_list[[x]]$marker_metrics$tntest$global
    list_oracle_global[[length(list_countsplit_global) + 1]] <- df_list[[x]]$marker_metrics$oracle$global
    list_clust[[length(list_clust) + 1]] <- df_list[[x]]$clust
    list_counts[[length(list_counts) + 1]] <- df_list[[x]]$counts
    list_auc[[length(list_auc) + 1]] <- df_list[[x]]$auc

    if(!is.na(df_list[[x]]$skip)) skipped_files <- c(skipped_files, x)
    if(!is.na(df_list[[x]]$error)) error_files <- c(error_files, x)
  }

  df_stabcell_global <- bind_rows(list_stabcell_global)
  df_ddip_global <- bind_rows(list_ddip_global)
  df_countsplit_global <- bind_rows(list_countsplit_global)
  df_clusterde_global <- bind_rows(list_clusterde_global)
  df_tntest_global <- bind_rows(list_tntest_global)
  df_oracle_global <- bind_rows(list_oracle_global)
  df_clust <- bind_rows(list_clust)
  df_counts <- bind_rows(list_counts)
  df_auc <- bind_rows(list_auc)

  save(
    df_stabcell_global,
    df_ddip_global,
    df_countsplit_global,
    df_clusterde_global,
    df_tntest_global,
    df_oracle_global,
    df_clust,
    df_counts,
    df_auc,

    file=paste0("results/", scenario, "/metrics/results_metrics_total.Rdata")
  )


  if(length(skipped_files) > 0) cat("\n\nSkipped the following files:\n", skipped_files, "\n")
  if(length(error_files) > 0) cat("\n\nCould not load the following files:", paste(error_files, "-", files[error_files]), "\n", sep="\n")

  cat("finished processing at ", strftime(Sys.time()), "\n")

  summary(warnings())
}


