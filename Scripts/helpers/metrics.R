library(dplyr)

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
  F1_score <- 2 * precision * recall / (precision + recall)
  
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
    method_marker_df,
    k
){
  # rank true markers by LFC
  # rank method markers by p value
  
  if(nrow(true_log_fold_changes) == 0) return(NA)
  
  if("p_val_mean" %in% colnames(method_marker_df)){
    method_marker_df <- method_marker_df %>% dplyr::rename(pval=p_val_mean)
  }
  
  top_k_truth <- true_log_fold_changes %>% 
    arrange(desc(abs(lfc))) %>%
    head(k) %>%
    pull(marker)
  top_k_method <- method_marker_df %>%
    arrange(pval) %>%
    head(k) %>%
    pull(marker)
  
  agreement <- length(intersect(top_k_truth, top_k_method)) / k
  return(agreement)
}
