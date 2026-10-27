library(ggplot2)
library(patchwork)
library(ggh4x)
library(purrr)
library(ggnewscale)


my_colors <- c(
  "Double Dipping"="#e33e50",
  "Double Dipping\nClusterDE"="#e33e50",
  "Countsplit"="#F2B705",
  "StabCell"="#2E7D32",
  "ClusterDE"="#8c7f9d",
  "TN Test"="#563d2c",
  "Oracle"="#D95100",

  "none"="#274b29",
  "unimodal"="#2E7D32",
  "r-concavity"="#6bb16f",

  "1 Celltype"="#04756F",
  "2 Celltypes"="#6cceba",
  "5 Celltypes"="#6cceba",
  "10 Celltypes"="#6cceba"
)


dynamic_device <- function(format){
  if(format %in% c("pdf", "svg", "eps")) return(cairo_ps)
  return(NULL)
}


pretty_name <- function(x){
  switch (x,
          "lib.size"="Library Size",
          "de.prob"="DE Probability",
          "de.facLoc"="DE Strength",
          "p_val_thresh"=expression(lambda),
          "n_subsample_markers"="q",
          "PFER_max"=expression(PFER[max]),
          "n_celltypes"="Cell Types",
          "fp_markers"="False Positives",
          "tp_markers"="True Positives",
          x
  )
}

custom_labeller <- function(labels, sep=": "){
  output <- list()
  labels <- as.data.frame(labels)
  for(i in 1:ncol(labels)){
    var_name <- gsub(" ", ' * " " * ', pretty_name(colnames(labels)[i]))
    var_outputs <- list()
    for(j in 1:nrow(labels)){
      var_value <- gsub(" ", ' * " " * ', labels[j, i])
      txt <- if(is.na(var_value)) "" else paste0(var_name, " * '", sep, "' * ", var_value)
      var_outputs[[j]] <- parse(text=txt)
    }
    output[[i]] <- var_outputs
  }
  output
}


extract_present_factor_levels <- function(x){
  levels(x)[levels(x) %in% x]
}

plot_marker_metric <- function(
    data,
    yvar="fp_markers",
    rows=c("de.prob", "de.facLoc"),
    cols="lib.size",
    plot_mean=TRUE
) {
  pfer_levels <- data %>%
    filter(method == "StabCell", !is.na(PFER_max)) %>%
    pull(PFER_max) %>%
    extract_present_factor_levels()
  
  assumption_levels <- data %>%
    filter(method == "StabCell") %>%
    pull(assumption) %>%
    extract_present_factor_levels()
  
  method_levels <- data %>%
    filter(method != "StabCell") %>%
    pull(method) %>%
    extract_present_factor_levels()
  
  
  plot_data <- data %>%
    mutate(
      x_plot=factor(
        if_else(method == "StabCell", as.character(PFER_max), ""),
        levels=c(as.character(pfer_levels), "")
      ),
      
      # determines which boxplots are dodged within each x value
      series=factor(
        if_else(method == "StabCell", as.character(assumption), as.character(method)),
        levels=c(assumption_levels, method_levels)
      )
    )
  
  
  
  # ------------------------------------------------------------------
  # Base plot
  # ------------------------------------------------------------------
  
  my_dodge <- position_dodge2(width=0.75)
  
  p <- ggplot(
    plot_data,
    aes(x=x_plot, y=.data[[yvar]], group=interaction(x_plot, series))
  ) +
    
    # StabCell asusmptions
    geom_boxplot(
      data=plot_data %>% filter(method == "StabCell"),
      mapping=aes(color=series),
      position=my_dodge
    ) +
    scale_color_manual(
      name="StabCell Assumption", values=my_colors, 
      breaks=assumption_levels, guide=guide_legend(nrow=1, order=1)
    ) +
    
    # other methods
    new_scale_color() +
    geom_boxplot(
      data=plot_data %>% filter(method != "StabCell"),
      mapping=aes(color=series),
      position=my_dodge
    ) +
    scale_color_manual(
      name="Method", values=my_colors, 
      breaks=method_levels, guide=guide_legend(nrow=1, order=2)
    )
  
  # ------------------------------------------------------------------
  # Means
  # ------------------------------------------------------------------
  
  if (plot_mean) {
    
    mean_data <- plot_data %>%
      group_by(across(all_of(unique(c("x_plot", "series", "method", "PFER_max", rows, cols))))) %>%
      summarise(mean_y=mean(.data[[yvar]], na.rm=TRUE), .groups="drop") %>%
      mutate(
        mean_gt_PFERmax=yvar == "fp_markers" & method == "StabCell" & !is.na(PFER_max) & mean_y > as.numeric(as.character(PFER_max)),
        mean_color=if_else(mean_gt_PFERmax, "red", "black")
      )
    
    
    p <- p +
      new_scale_color() +
      geom_point(
        data=mean_data,
        aes(x=x_plot, y=mean_y, group=series, color=I(mean_color)),
        inherit.aes=FALSE,
        position=my_dodge,
        shape=4,
        show.legend=FALSE
      ) +
      scale_color_identity(guide="none")
  }
  
  
  # ------------------------------------------------------------------
  # dashed lines as reference level
  # ------------------------------------------------------------------
  
  if (yvar == "fp_markers") {
    
    pfer_lines <- plot_data %>%
      filter(method == "StabCell", !is.na(PFER_max)) %>%
      distinct(across(all_of(unique(c("x_plot", "PFER_max", rows, cols))))) %>%
      mutate(
        # Position of this category on the discrete x scale.
        # This is ONLY used for drawing the segment.
        x_center=match(as.character(x_plot), levels(plot_data$x_plot)),
        PFER_max_value=as.numeric(as.character(PFER_max))
      )
    
    p <- p +
      geom_segment(
        data=pfer_lines,
        aes(x=x_center - 0.45, xend=x_center + 0.45, y=PFER_max_value, yend=PFER_max_value),
        inherit.aes=FALSE,
        linetype="dashed",
        linewidth=0.5,
        color="black"
      )
  }
  
  
  
  # ------------------------------------------------------------------
  # Facets + styling
  # ------------------------------------------------------------------
  
  p +
    facet_nested(
      rows=vars(!!!rlang::syms(rows)),
      cols=vars(!!!rlang::syms(cols)),
      labeller=custom_labeller,
      scales="free_y"
    ) +
    labs(
      x=pretty_name("PFER_max"),
      y=pretty_name(yvar)
    ) +
    scale_x_discrete(breaks=pfer_levels) +
    theme(
      panel.background=element_blank(),
      panel.grid.major.x=element_blank(),
      panel.grid.major.y=element_line(color="lightgrey"),
      
      legend.position="top",
      legend.direction="horizontal",
      legend.box="vertical",
      legend.box.just="center",
      legend.title.position="left",
      legend.margin=margin(t=-5, b=-2)
    ) +
    coord_cartesian(ylim=c(0, NA))
}



