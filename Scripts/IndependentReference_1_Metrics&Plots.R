library(RhpcBLASctl)
library(doFuture)
library(progressr)
library(mclust)
library(dplyr)
library(ggplot2)
library(patchwork)

blas_threads <- 8
num_workers <- 25

my_colors <- c(
  "#e33e50",
  "#F2B705",
  "#274b29",
  "#2E7D32",
  "#6bb16f",
  "#04756F",
  "#6cceba",
  "#8c7f9d",
  "#563d2c",
  "#D95100"
)


blas_set_num_threads(blas_threads)
options(future.globals.maxSize=16 * 1024^3)
plan(multicore, workers=num_workers)
registerDoFuture()


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

results_dir <- "results/independent_reference"
results_file <- paste0(results_dir, "/metrics_indepRef.Rdata")
plot_dir <- "plots/independent_reference"

if(file.exists(results_file)){
  load(results_file)
} else {
  files <- list.files(results_dir, pattern="^results_indepRef__.*\\.Rdata$")

  cat("num files: ", length(files), "\n")


  with_progress({

    p <- progressor(along=files)

    results_list <- foreach(
      x=1:length(files),
      .packages=c("dplyr", "foreach", "tidyr", "Matrix")
    ) %dopar% {

      plan(sequential)
      registerDoSEQ()

      ###################
      #### load data ####
      ###################

      file <- files[x]

      # try to load simulation results, exit foreach iter if it fails
      results <- tryCatch(
        expr={
          load(paste0(results_dir, "/", file))
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
              split=list(
                global=data.frame(),
                truth_alignment=data.frame()
              ),
              nosplit=list(
                global=data.frame(),
                truth_alignment=data.frame()
              )
            ),
            clust=data.frame(),
            error=x
          )
        )
      }

      iter <- results$iter

      # unpack simulation results
      count_data <- results$count_data
      true_idents <- results$true_idents
      true_fold_changes <- results$true_fold_changes

      # "with split" results
      split_idx <- results$with_split$idx_split
      split_idents_cluster <- results$with_split$clustering_idents
      split_idents_marker <- results$with_split$marker_idents
      split_idents <- rep(NA, ncol(count_data))
      split_idents[split_idx] <- split_idents_cluster
      split_idents[-split_idx] <- split_idents_marker
      split_markers <- results$with_split$stable_markers_global

      # "without split" results
      nosplit_idents <- results$no_split$clustering_mode
      nosplit_markers <- results$no_split$stable_markers_global



      ##########################
      #### evaluate markers ####
      ##########################

      # find true marker genes
      true_idents_unique <- unique(true_idents)
      n_true_idents <- length(true_idents_unique)

      temp_list <- list()
      for(i in 1:n_true_idents){
        id <- true_idents_unique[i]

        marker_idx <- true_fold_changes[ , paste0("DEFacGroup", id)] != 1.
        lfc <- log2(true_fold_changes[marker_idx, paste0("DEFacGroup", id)])

        temp_list[[length(temp_list) + 1]] <- data.frame(
          id=id,
          marker=rownames(true_fold_changes)[marker_idx],
          lfc=lfc
        )
      }
      true_markers <- bind_rows(temp_list)

      # alignment for comparison
      # -> for every cluster in default StabCell (nosplit) find closest
      #    cluster in split approach
      nosplit_ids <- unique(nosplit_idents)
      split_ids <- rep("", length(nosplit_ids))

      for(i in seq_along(nosplit_ids)){
        split_ids[i] <- find_closest_cluster(
          reference_allocations=split_idents,
          reference_cluster_id=nosplit_ids[i],
          comparison_allocations=nosplit_idents
        )
      }
      method_alignment_df <- data.frame(
        nosplit_id=nosplit_ids,
        split_id=split_ids
      )

      # calculate metrics
      marker_metrics <- list()
      for(method in c("split", "nosplit")){

        if(method == "split"){
          # combine clustering and knn marker idents
          method_idents <- as.character(split_idents)
          method_markers <- split_markers %>% mutate(id=as.character(id))
        } else if(method == "nosplit"){
          method_idents <- nosplit_idents
          method_markers <- nosplit_markers
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


        # global marker metrics
        df_metrics <- foreach(
          i=1:nrow(df_align_simulation),
          .combine="bind_rows",
          .packages=c("dplyr")
        ) %do% {
          method_id <- df_align_simulation$method_id[i]
          simulation_id <- df_align_simulation$simulation_id_aligned[i]

          true_markers_df <- (true_markers %>% filter(id == !!simulation_id))
          false_markers <- setdiff(rownames(true_fold_changes), true_markers_df$marker)

          df_temp <- foreach(
            assumption=c("none", "unimodal", "r_concavity"),
            .combine="bind_rows",
            .packages=c("dplyr")
          ) %do% {
            selected_markers_df <- (
              method_markers %>%
                mutate(id=as.character(id)) %>%
                filter(
                  assumption == !!assumption,
                  id == as.character(!!method_id),
                  !is.na(marker)
                )
            )
            unselected_markers <- setdiff(rownames(true_fold_changes), selected_markers_df$marker)

            df_temp <- cbind(
              iter=iter,

              id=method_id,
              truth_id_aligned=simulation_id,
              assumption=assumption,

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

            return(df_temp)
          }
          return(df_temp)
        }

        # align cluster ids between methods
        if(method == "split"){
          df_metrics <- df_metrics %>%
            left_join(method_alignment_df, by=c("id"="split_id")) %>%
            rename(split_id=id, id=nosplit_id)
        } else {
          df_metrics <- df_metrics %>%
            left_join(method_alignment_df, by=c("id"="nosplit_id"))
        }


        marker_metrics[[method]] <- list(
          global=df_metrics,
          truth_alignment=df_align_simulation
        )
      }


      ############################
      #### clustering metrics ####
      ############################
      df_clust <- data.frame(
        iter=iter,
        n_clusters_true=length(unique(true_idents)),
        n_clusters_split=length(unique(split_idents)),
        n_clusters_nosplit=length(unique(nosplit_idents)),
        ari_split=adjustedRandIndex(true_idents, split_idents),
        ari_nosplit=adjustedRandIndex(true_idents, nosplit_idents)
      )

      p()
      return(
        list(
          marker_metrics=marker_metrics,
          clust=df_clust,
          error=NA
        )
      )
    }
  })

  # combine iteration metrics into one data frame
  error_files <- numeric()
  list_split <- list()
  list_nosplit <- list()
  list_clust <- list()

  for(x in 1:length(results_list)){
    list_split[[length(list_split) + 1]] <- results_list[[x]]$marker_metrics$split$global
    list_nosplit[[length(list_nosplit) + 1]] <- results_list[[x]]$marker_metrics$nosplit$global
    list_clust[[length(list_clust) + 1]] <- results_list[[x]]$clust

    if(!is.na(results_list[[x]]$error)) error_files <- c(error_files, x)
  }

  df_split <- bind_rows(list_split)
  df_nosplit <- bind_rows(list_nosplit)
  df_clust <- bind_rows(list_clust)

  save(
    df_split,
    df_nosplit,
    df_clust,

    file=results_file
  )


  if(length(error_files) > 0) cat("\n\nCould not load the following files:", paste(error_files, "-", files[error_files]), "\n", sep="\n")
}







id <- 1
plot_metrics <- c("fp_markers", "tp_markers")
{
plot_list <- list()

# marker metrics
for(plot_metric in plot_metrics){
  p <- ggplot(
    data=bind_rows(list(
      df_split %>%
        filter(id == !!id) %>%
        select(!!plot_metric, assumption) %>%
        mutate(method="Independent Reference"),
      df_nosplit %>%
        filter(id == !!id) %>%
        select(!!plot_metric, assumption) %>%
        mutate(method="StabCell")
    )) %>% mutate(
      method=factor(method, levels=c("StabCell", "Independent Reference")),
      assumption=factor(assumption, levels=c("none", "unimodal", "r_concavity")),
    ),
    mapping=aes(x=method, y=!!sym(plot_metric), color=assumption)
  ) +
    geom_boxplot() +
    labs(
      x=element_blank(),
      y=switch(
        plot_metric,
        "fp_markers"="False Positives",
        "tp_markers"="True Positives",
        "fn_markers"="False Negatives",
        "tn_markers"="True Negatives",
        "true_negative_rate"="True Negative Rate",
        "F1_score"="F1 Score",
        stringr::str_to_title(plot_metric)
      ),
      color=element_blank()
    ) +
    theme(
      legend.position="top",
      panel.background=element_blank(),
      panel.grid=element_blank(),
      panel.grid.major.y=element_line(color="lightgrey"),
      # axis.ticks.x=element_blank(),
      # axis.text.x=element_blank()
    ) +
    scale_color_manual(values=my_colors[c(1, 2, 4)])

  if(plot_metric == "fp_markers"){
    p <- p + geom_hline(mapping=aes(yintercept=25), color="black", linetype="dashed")
  }
  if(!stringr::str_ends(plot_metric, "markers")){
    p <- p + coord_cartesian(ylim=0:1)
  }

  plot_list[[length(plot_list) + 1]] <- p
}

# # clustering metrics
# plot_list[[length(plot_list) + 1]] <- ggplot(
#   data=bind_rows(list(
#     data.frame(
#       method="Independent Reference",
#       value=df_clust$ari_split
#     ),
#     data.frame(
#       method="StabCell",
#       value=df_clust$ari_nosplit
#     )
#   )) %>% mutate(method=factor(method, levels=c("StabCell", "Independent Reference"))),
#   mapping=aes(x=method, y=value, color=method)
# ) +
#   geom_boxplot() +
#   labs(
#     x=element_blank(),
#     y="Adjusted Rand Index",
#     color=element_blank()
#   ) +
#   theme(
#     legend.position="none",
#     panel.background=element_blank(),
#     panel.grid=element_blank(),
#     panel.grid.major.y=element_line(color="lightgrey"),
#     # axis.ticks.x=element_blank(),
#     # axis.text.x=element_blank()
#   ) +
#   coord_cartesian(ylim=0:1)

# combine
p <- wrap_plots(plot_list) +
  plot_layout(guides = "collect") &
  theme(legend.position="top")
plot(p)

ggsave(
  plot=p,
  filename="IndependentReference.eps",
  path=plot_dir,
  width=17, height=7, unit="cm", dpi=350
)
ggsave(
  plot=p,
  filename="IndependentReference.tif",
  path=plot_dir,
  width=17, height=7, unit="cm", dpi=350
)
}


