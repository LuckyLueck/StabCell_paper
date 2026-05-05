library(doFuture)
library(dplyr)
library(RhpcBLASctl)
library(mclust)
library(progressr)

blas_threads <- 2
num_workers <- 60

options(future.globals.maxSize=4 * 1024^3)

blas_set_num_threads(blas_threads)
plan(multicore, workers=num_workers)
registerDoFuture()

p_val_thresh_alt <- 0.05 # threshold to select markers in alternative approaches (double dipping, countsplit)

scenarios <- c(
  "2Cluster_Simulation",
  "NCluster_Simulation"
)

calculate_metrics <- function(
  selected_markers,
  unselected_markers,
  true_markers,
  false_markers
){
  n_markers <- length(selected_markers)
  tp_markers <- length(intersect(selected_markers, true_markers))
  fp_markers <- length(intersect(selected_markers, false_markers))
  tn_markers <- length(intersect(unselected_markers, false_markers))
  fn_markers <- length(intersect(unselected_markers, true_markers))

  accuracy <- (tp_markers + tn_markers) / (tp_markers + tn_markers + fp_markers + fn_markers)
  precision <- tp_markers / (tp_markers + fp_markers)
  recall <- tp_markers / (tp_markers + fn_markers)
  true_negative_rate <- tn_markers / (tn_markers + fp_markers)
  F1_score <- precision * recall / (precision + recall)

  results <- data.frame(
    n_markers=n_markers,
    tp_markers=tp_markers,
    fp_markers=fp_markers,
    tn_markers=tn_markers,
    fn_markers=fn_markers,
    accuracy=accuracy,
    precision=precision,
    recall=recall,
    true_negative_rate=true_negative_rate,
    F1_score=F1_score
  )

  return(results)
}

top_k_marker_agreement <- function(
  true_log_fold_changes,
  log_fold_changes,
  k,
  seurat=FALSE
){

  if(nrow(true_log_fold_changes) == 0) return(NA)

  top_k_truth <- (true_log_fold_changes %>% arrange(desc(abs(lfc))))[1:k, "marker"]
  if(seurat){
    top_k_stable <- (log_fold_changes %>% arrange(desc(abs(avg_log2FC))))[1:k, "marker"]
  } else {
    top_k_stable <- (log_fold_changes %>% arrange(desc(abs(lfc_mean))))[1:k, "marker"]
  }
  agreement <- length(intersect(top_k_truth, top_k_stable)) / k
  return(agreement)
}

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
                local=data.frame(),
                truth_alignment=data.frame(),
                stabcell_alignment=data.frame()
              ),
              countsplit=list(
                global=data.frame(),
                local=data.frame(),
                truth_alignment=data.frame(),
                stabcell_alignment=data.frame()
              ),
              stabcell=list(
                global=data.frame(),
                local=data.frame(),
                truth_alignment=data.frame(),
                stabcell_alignment=data.frame()
              )
            ),
            clust=data.frame(),
            counts=data.frame(),
            skip=NA,
            error=x
          )
        )
      }

      # unpack config data
      n_celltypes <- results$config$n_celltypes
      lib.loc <- results$config$lib.loc
      de.prob <- results$config$de.prob
      de.facLoc <- results$config$de.facLoc
      subsampling_iterations <- results$config$n_iterations
      complementary_pairs <- results$config$complementary_pairs
      subsampling_fraction <- results$config$subsampling_fraction
      iter <- as.numeric(sub(".*__iter=([0-9]*)\\.Rdata", "\\1", file))

      # unpack simulation results
      count_data <- results$count_data
      true_idents <- results$true_idents
      true_fold_changes <- results$true_fold_changes
      umap <- results$umap

      # double dipping results
      ddip_idents <- results$double_dipping$identities
      ddip_markers <- results$double_dipping$markers

      # countsplit results
      countsplit_idents <- results$countsplit$identities
      countsplit_markers <- results$countsplit$markers

      # stabcell results
      subsampling_identities <- results$stabcell$subsampling_identities
      clustering_alignment <- results$stabcell$clustering_alignment
      stabcell_idents <- results$stabcell$clustering_mode
      cell_stability_score <- if(ncol(results$stabcell$assignment_probabilities) == 1) rep(1, nrow(results$stabcell$assignment_probabilities)) else results$stabcell$cell_stability_score
      marker_alignment <- results$stabcell$marker_alignment
      stable_markers_global <- results$stabcell$stable_markers_global
      stable_markers_local <- results$stabcell$stable_markers_local


      ##########################
      #### evaluate markers ####
      ##########################

      # find true marker genes
      if(n_celltypes == 1){

        # no true marker genes
        true_markers_local <- data.frame(id1=numeric(), id2=numeric(), marker=character(), lfc=numeric())
        true_markers_global <- data.frame(id=numeric(), marker=character(), lfc=numeric())

      } else if(n_celltypes == 2){

        # local and global markers are identical
        marker_idx <- (true_fold_changes[, "DEFacGroup1"] != true_fold_changes[, "DEFacGroup2"])
        lfc <- log2(true_fold_changes[marker_idx, "DEFacGroup1"] / true_fold_changes[marker_idx, "DEFacGroup2"])

        true_markers_local <- data.frame(
          id1=1,
          id2=2,
          marker=rownames(true_fold_changes)[marker_idx],
          lfc=lfc
        )
        true_markers_global <- bind_rows( # local & global marker markers are same when n_celltypes == 2
          true_markers_local %>% select(-id2) %>% rename(id=id1),
          true_markers_local %>% select(-id1) %>% rename(id=id2) %>% mutate(lfc=-1 * lfc)
        )

      } else { # n_celltypes > 2

        true_idents_unique <- unique(true_idents)
        n_true_idents <- length(true_idents_unique)

        # find true local markers
        df_list = list()
        for(i in 1:(n_true_idents - 1)){
          for(j in (i+1):n_true_idents){
            id1 <- true_idents_unique[i]
            id2 <- true_idents_unique[j]

            marker_idx <- true_fold_changes[ , paste0("DEFacGroup", id1)] != true_fold_changes[ , paste0("DEFacGroup", id2)]
            lfc <- log2(
              true_fold_changes[marker_idx, paste0("DEFacGroup", id1)]
              / true_fold_changes[marker_idx, paste0("DEFacGroup", id2)]
            )

            df_list[[length(df_list) + 1]] <- data.frame(
              id1=id1,
              id2=id2,
              marker=rownames(true_fold_changes)[marker_idx],
              lfc=lfc
            )

          }
        }
        true_markers_local <- bind_rows(df_list)


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
      for(method in c("ddip", "countsplit", "stabcell")){

        if(method == "ddip"){
          method_idents <- ddip_idents
          method_markers <- ddip_markers
        } else if(method == "countsplit"){
          method_idents <- countsplit_idents
          method_markers <- countsplit_markers
        } else if(method == "stabcell"){
          method_idents <- stabcell_idents
        } else {
          stop("Unknown method!")
        }

        # align method clusters to true clusters
        method_idents_unique <- unique(method_idents)
        simulation_ids <- rep(NA, length(method_idents_unique))

        for(i in seq_along(method_idents_unique)){
          simulation_ids[i] <- find_closest_cluster(
            reference_allocations=method_idents,
            reference_cluster_id=method_idents_unique[i],
            comparison_allocations=true_idents
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
              complementary_pairs=complementary_pairs,
              n_celltypes=n_celltypes,
              n_iter=subsampling_iterations,
              lib.loc=lib.loc,
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
                log_fold_changes=selected_markers_df,
                k=1
              ),
              top5_agreement=top_k_marker_agreement(
                true_log_fold_changes=true_markers_df,
                log_fold_changes=selected_markers_df,
                k=5
              ),
              top10_agreement=top_k_marker_agreement(
                true_log_fold_changes=true_markers_df,
                log_fold_changes=selected_markers_df,
                k=10
              )
            )


            return(df_temp)

          }

          df_metrics_global <- df_metrics_global %>% mutate(
            complementary_pairs=as.factor(complementary_pairs),
            n_iter=as.factor(n_iter),
            lib.loc=as.factor(lib.loc),
            de.prob=as.factor(de.prob),
            de.facLoc=as.factor(de.facLoc),
            PFER_max=as.factor(PFER_max),
            assumption=as.factor(assumption)
          )


          # local marker metrics
          if(n_celltypes > 2){

            settings_local <- stable_markers_local %>%
              group_by(q, p_val_thresh, PFER_max, assumption, selection_freq_thresh, exceeds_one, id1, id2) %>%
              summarize(.groups="drop")

            df_metrics_local <- foreach(
              i=1:nrow(settings_local),
              .combine="bind_rows",
              .packages=c("dplyr")
            ) %do% {
              q <- settings_local$q[i]
              p_val_thresh <- settings_local$p_val_thresh[i]
              PFER_max <- settings_local$PFER_max[i]
              assumption <- settings_local$assumption[i]
              selection_freq_thresh <- settings_local$selection_freq_thresh[i]
              exceeds_one <- settings_local$exceeds_one[i]
              stabcell_id1 <- settings_local$id1[i]
              stabcell_id2 <- settings_local$id2[i]
              simulation_id1 <- df_align_simulation[df_align_simulation$method_id == stabcell_id1, "simulation_id_aligned"]
              simulation_id2 <- df_align_simulation[df_align_simulation$method_id == stabcell_id2, "simulation_id_aligned"]

              selected_markers_df <- (
                stable_markers_local %>%
                  filter(
                    q == !!q,
                    p_val_thresh == !!p_val_thresh,
                    PFER_max == !!PFER_max,
                    assumption == !!assumption,
                    selection_freq_thresh == !!selection_freq_thresh,
                    ((id1 == !!stabcell_id1) & (id2 == !!stabcell_id2)) | ((id1 == !!stabcell_id2) & (id2 == !!stabcell_id1)),
                    !is.na(marker)
                  )
              )
              unselected_markers <- setdiff(rownames(true_fold_changes), selected_markers_df$marker)

              if(simulation_id1 == simulation_id2){
                # both clusters are aligned to the same **true** cluster
                df_temp <- data.frame(
                  complementary_pairs=complementary_pairs,
                  n_celltypes=n_celltypes,
                  n_iter=subsampling_iterations,
                  lib.loc=lib.loc,
                  de.prob=de.prob,
                  de.facLoc=de.facLoc,
                  iter=iter,

                  q=q,
                  p_val_thresh=p_val_thresh,
                  PFER_max=PFER_max,
                  assumption=assumption,
                  selection_freq_thresh=selection_freq_thresh,
                  exceeds_one=exceeds_one,

                  id1=stabcell_id1,
                  id2=stabcell_id2,
                  truth_id1_aligned=simulation_id1,
                  truth_id2_aligned=simulation_id2,

                  n_markers=length(selected_markers_df$marker),
                  tp_markers=NA,
                  fp_markers=NA,
                  tn_markers=NA,
                  fn_markers=NA,
                  accuracy=NA,
                  precision=NA,
                  recall=NA,
                  true_negative_rate=NA,
                  F1_score=NA,

                  top1_agreement=NA,
                  top5_agreement=NA,
                  top10_agreement=NA
                )
                return(df_temp)  # exit foreach iteration early
              }

              true_markers_df <- (
                true_markers_local %>%
                  filter(((id1 == !!simulation_id1) & (id2 == !!simulation_id2)) | ((id1 == !!simulation_id2) & (id2 == !!simulation_id1)))
              )
              false_markers <- setdiff(rownames(true_fold_changes), true_markers_df$marker)

              df_temp <- cbind(
                complementary_pairs=complementary_pairs,
                n_celltypes=n_celltypes,
                n_iter=subsampling_iterations,
                lib.loc=lib.loc,
                de.prob=de.prob,
                de.facLoc=de.facLoc,
                iter=iter,

                q=q,
                p_val_thresh=p_val_thresh,
                PFER_max=PFER_max,
                assumption=assumption,
                selection_freq_thresh=selection_freq_thresh,
                exceeds_one=exceeds_one,

                id1=stabcell_id1,
                id2=stabcell_id2,
                truth_id1_aligned=simulation_id1,
                truth_id2_aligned=simulation_id2,

                calculate_metrics(
                  selected_markers=selected_markers_df$marker,
                  unselected_markers=unselected_markers,
                  true_markers=true_markers_df$marker,
                  false_markers=false_markers
                ),

                top1_agreement=top_k_marker_agreement(
                  true_log_fold_changes=true_markers_df,
                  log_fold_changes=selected_markers_df,
                  k=1
                ),
                top5_agreement=top_k_marker_agreement(
                  true_log_fold_changes=true_markers_df,
                  log_fold_changes=selected_markers_df,
                  k=5
                ),
                top10_agreement=top_k_marker_agreement(
                  true_log_fold_changes=true_markers_df,
                  log_fold_changes=selected_markers_df,
                  k=10
                )
              )

              df_temp

            }

            df_metrics_local <- df_metrics_local %>% mutate(
              complementary_pairs=as.factor(complementary_pairs),
              n_iter=as.factor(n_iter),
              lib.loc=as.factor(lib.loc),
              de.prob=as.factor(de.prob),
              de.facLoc=as.factor(de.facLoc),
              PFER_max=as.factor(PFER_max),
              assumption=as.factor(assumption)
            )
          } else {
            df_metrics_local <- data.frame()
          }

        } else {  # ddip & countsplit

          # find alignment to stabcell clustering
          stabcell_idents_unique <- unique(stabcell_idents)
          method_ids <- rep(NA, length(stabcell_idents_unique))
          for(i in seq_along(stabcell_idents_unique)){
            method_ids[i] <- find_closest_cluster(
              reference_allocations=stabcell_idents,
              reference_cluster_id=stabcell_idents_unique[i],
              comparison_allocations=method_idents
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
                  p_val_adj <= !!p_val_thresh_alt
                )
            )
            unselected_markers <- setdiff(rownames(true_fold_changes), selected_markers_df$marker)

            true_markers_df <- (true_markers_global %>% filter(id == !!simulation_id))
            false_markers <- setdiff(rownames(true_fold_changes), true_markers_df$marker)

            df_temp <- cbind(
              complementary_pairs=complementary_pairs,
              n_celltypes=n_celltypes,
              n_iter=subsampling_iterations,
              lib.loc=lib.loc,
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
                log_fold_changes=selected_markers_df,
                k=1,
                seurat=TRUE
              ),
              top5_agreement=top_k_marker_agreement(
                true_log_fold_changes=true_markers_df,
                log_fold_changes=selected_markers_df,
                k=5,
                seurat=TRUE
              ),
              top10_agreement=top_k_marker_agreement(
                true_log_fold_changes=true_markers_df,
                log_fold_changes=selected_markers_df,
                k=10,
                seurat=TRUE
              )
            )

            df_temp

          }

          df_metrics_global <- df_metrics_global %>%
            mutate(
              complementary_pairs=as.factor(complementary_pairs),
              n_iter=as.factor(n_iter),
              lib.loc=as.factor(lib.loc),
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

          # calculate local metrics
          if(n_celltypes > 2){
            df_metrics_local <- foreach(
              i=1:(length(method_idents_unique) - 1),
              .combine="bind_rows",
              .packages=c("dplyr")
            ) %:%
              foreach(
                j=(i + 1):length(method_idents_unique),
                .combine="bind_rows"
              ) %do% {
                method_id1  <- df_align_simulation[i, "method_id"]
                method_id2  <- df_align_simulation[j, "method_id"]
                simulation_id1 <- df_align_simulation[i, "simulation_id_aligned"]
                simulation_id2 <- df_align_simulation[j, "simulation_id_aligned"]

                selected_markers_df <- (
                  method_markers %>%
                    filter(
                      ((cluster1 == !!method_id1) & (cluster2 == !!method_id2)) | ((cluster1 == !!method_id2) & (cluster2 == !!method_id1)),
                      p_val_adj <= !!p_val_thresh_alt
                    )
                )
                unselected_markers <- setdiff(rownames(true_fold_changes), selected_markers_df$marker)

                if(simulation_id1 == simulation_id2){
                  # both clusters are aligned to the same **true** cluster
                  df_temp <- data.frame(
                    complementary_pairs=complementary_pairs,
                    n_celltypes=n_celltypes,
                    n_iter=subsampling_iterations,
                    de.prob=de.prob,
                    de.facLoc=de.facLoc,
                    iter=iter,

                    id1=method_id1,
                    id2=method_id2,
                    truth_id1_aligned=simulation_id1,
                    truth_id2_aligned=simulation_id2,

                    n_markers=length(selected_markers_df$marker),
                    tp_markers=NA,
                    fp_markers=NA,
                    tn_markers=NA,
                    fn_markers=NA,
                    accuracy=NA,
                    precision=NA,
                    recall=NA,
                    true_negative_rate=NA,
                    F1_score=NA,

                    top1_agreement=NA,
                    top5_agreement=NA,
                    top10_agreement=NA
                  )
                } else {

                  true_markers_df <- (
                    true_markers_local %>%
                      filter(((id1 == !!simulation_id1) & (id2 == !!simulation_id2)) | ((id1 == !!simulation_id2) & (id2 == !!simulation_id1)))
                  )
                  false_markers <- setdiff(rownames(true_fold_changes), true_markers_df$marker)


                  df_temp <- data.frame(
                    complementary_pairs=complementary_pairs,
                    n_celltypes=n_celltypes,
                    n_iter=subsampling_iterations,
                    lib.loc=lib.loc,
                    de.prob=de.prob,
                    de.facLoc=de.facLoc,
                    iter=iter,

                    id1=method_id1,
                    id2=method_id2,
                    truth_id1_aligned=simulation_id1,
                    truth_id2_aligned=simulation_id2,

                    calculate_metrics(
                      selected_markers=selected_markers_df$marker,
                      unselected_markers=unselected_markers,
                      true_markers=true_markers_df$marker,
                      false_markers=false_markers
                    ),

                    top1_agreement=top_k_marker_agreement(
                      true_log_fold_changes=true_markers_df,
                      log_fold_changes=selected_markers_df,
                      k=1,
                      seurat=TRUE
                    ),
                    top5_agreement=top_k_marker_agreement(
                      true_log_fold_changes=true_markers_df,
                      log_fold_changes=selected_markers_df,
                      k=5,
                      seurat=TRUE
                    ),
                    top10_agreement=top_k_marker_agreement(
                      true_log_fold_changes=true_markers_df,
                      log_fold_changes=selected_markers_df,
                      k=10,
                      seurat=TRUE
                    )
                  )
                }

                df_temp

              }

            df_metrics_local <- df_metrics_local %>%
              mutate(
                complementary_pairs=as.factor(complementary_pairs),
                n_iter=as.factor(n_iter),
                lib.loc=as.factor(lib.loc),
                de.prob=as.factor(de.prob),
                de.facLoc=as.factor(de.facLoc)
              )

            # transfer stabcell ids
            df_metrics_local <- df_metrics_local %>%
              full_join(df_align_stabcell, by=c("id1"="method_id_aligned"), relationship="many-to-many") %>%
              rename(
                method_id1=id1,
                id1=stabcell_id
              ) %>%
              full_join(df_align_stabcell, by=c("id2"="method_id_aligned"), relationship="many-to-many") %>%
              rename(
                method_id2=id2,
                id2=stabcell_id
              )

          } else {
            df_metrics_local <- data.frame()
          }
        }

        marker_metrics[[method]] <- list(
          global=df_metrics_global,
          local=df_metrics_local,
          truth_alignment=df_align_simulation,
          stabcell_alignment=if(method != "stabcell") df_align_stabcell else NULL
        )

      }

      ############################
      #### clustering metrics ####
      ############################
      df_clust <- data.frame(
        complementary_pairs=complementary_pairs,
        n_celltypes=n_celltypes,
        n_iter=subsampling_iterations,
        lib.loc=lib.loc,
        de.prob=de.prob,
        de.facLoc=de.facLoc,
        iter=iter,
        n_clusters_true=length(unique(true_idents)),
        n_clusters_ddip=length(unique(ddip_idents)),
        n_clusters_countsplit=length(unique(countsplit_idents)),
        n_clusters_stabcell=length(unique(stabcell_idents)),
        ari_ddip=adjustedRandIndex(true_idents, ddip_idents),
        ari_countsplit=adjustedRandIndex(true_idents, countsplit_idents),
        ari_stabcell=adjustedRandIndex(true_idents, stabcell_idents),
        mean_cell_stability=mean(cell_stability_score)
      )


      ################
      #### counts ####
      ################
      cell_counts <- colSums(count_data)

      df_counts <- data.frame(
        complementary_pairs=complementary_pairs,
        n_celltypes=n_celltypes,
        n_iter=subsampling_iterations,
        lib.loc=lib.loc,
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

        file=out_file
      )

      p()
      return(
        list(
          marker_metrics=marker_metrics,
          clust=df_clust,
          counts=df_counts,
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
  list_stabcell_local <- list()
  list_ddip_global <- list()
  list_ddip_local <- list()
  list_countsplit_global <- list()
  list_countsplit_local <- list()
  list_clust <- list()
  list_counts <- list()

  for(x in 1:length(df_list)){
    list_stabcell_global[[length(list_stabcell_global) + 1]] <- df_list[[x]]$marker_metrics$stabcell$global
    list_stabcell_local[[length(list_stabcell_local) + 1]] <- df_list[[x]]$marker_metrics$stabcell$local
    list_ddip_global[[length(list_ddip_global) + 1]] <- df_list[[x]]$marker_metrics$ddip$global
    list_ddip_local[[length(list_ddip_local) + 1]] <- df_list[[x]]$marker_metrics$ddip$local
    list_countsplit_global[[length(list_countsplit_global) + 1]] <- df_list[[x]]$marker_metrics$countsplit$global
    list_countsplit_local[[length(list_countsplit_local) + 1]] <- df_list[[x]]$marker_metrics$countsplit$local
    list_clust[[length(list_clust) + 1]] <- df_list[[x]]$clust
    list_counts[[length(list_counts) + 1]] <- df_list[[x]]$counts

    if(!is.na(df_list[[x]]$skip)) skipped_files <- c(skipped_files, x)
    if(!is.na(df_list[[x]]$error)) error_files <- c(error_files, x)
  }

  df_stabcell_global <- bind_rows(list_stabcell_global)
  df_stabcell_local <- bind_rows(list_stabcell_local)
  df_ddip_local <- bind_rows(list_ddip_local)
  df_ddip_global <- bind_rows(list_ddip_global)
  df_countsplit_local <- bind_rows(list_countsplit_local)
  df_countsplit_global <- bind_rows(list_countsplit_global)
  df_clust <- bind_rows(list_clust)
  df_counts <- bind_rows(list_counts)

  save(
    df_stabcell_global,
    df_stabcell_local,
    df_ddip_global,
    df_ddip_local,
    df_countsplit_global,
    df_countsplit_local,
    df_clust,
    df_counts,

    file=paste0("results/", scenario, "/metrics/results_metrics_total.Rdata")
  )


  if(length(skipped_files) > 0) cat("\n\nSkipped the following files:\n", skipped_files, "\n")
  if(length(error_files) > 0) cat("\n\nCould not load the following files:", paste(error_files, "-", files[error_files]), "\n", sep="\n")

  cat("finished processing at ", strftime(Sys.time()), "\n")

  summary(warnings())
}


