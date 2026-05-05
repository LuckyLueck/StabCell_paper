library(ggplot2)
library(ggh4x)
library(patchwork)
library(dplyr)


# NOTE:
# - Precision (& therefore F1 score) are undefined when no marker is selected (divide 0 by 0)
# - Recall is undefined if there are no simulated markers, i.e., when n_celltypes == 1

use_simulation <- c(
  "2Cluster_Simulation",
  "NCluster_Simulation"
)[1]

load(paste0("results/", use_simulation, "/metrics/results_metrics_total.Rdata"))

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

# see which settings do not have 100 sets of iterations
# -> this can happen if clusters get merged and are no longer present
print(
  df_stabcell_global %>%
    group_by(complementary_pairs, n_celltypes, n_iter, lib.loc, de.prob, de.facLoc, PFER_max, assumption, p_val_thresh, id) %>%
    summarize(n=n(), .groups="drop") %>%
    filter(n != 100) %>%
    group_by(complementary_pairs, n_celltypes, n_iter, lib.loc, de.prob, de.facLoc, PFER_max, id) %>%
    summarize(n=mean(n), .groups="drop"),
  n=500
)


increasing_value_factor <- function(x){
  if(is.factor(x)){
    z <- factor(x, levels=sort(as.numeric(levels(x))))
  } else if(is.numeric(x)){
    z <- factor(x, levels=sort(unique(x)))
  } else {
    warning("Did not transform to increasing value factor.")
    z <- x
  }
  return(z)
}

# pretty print for variable names
pprint_variable <- function(var){
  switch (var,
          "lib.loc"="Library Size",
          "de.prob"="DE Probability",
          "de.facLoc"="DE Strength",
          "p_val_thresh"=expression(lambda),
          "n_subsample_markers"="q",
          "PFER_max"=expression(PFER[max]),
          "n_celltypes"="Cell Types",
          var
  )
}

# Nice label of variable name and value in plots
custom_labeller <- function(labels, sep=": "){
  output <- list()
  labels <- as.data.frame(labels)
  for(i in 1:ncol(labels)){
    var_name <- gsub(" ", ' * " " * ', pprint_variable(colnames(labels)[i]))
    var_outputs <- list()
    for(j in 1:nrow(labels)){
      var_label <- gsub(" ", ' * " " * ', labels[j, i])
      txt <- paste0(var_name, " * '", sep, "' * ", var_label)
      var_outputs[[j]] <- parse(text=txt)
      output[[i]] <- var_outputs
    }
  }
  output
}
class(custom_labeller) <- "labeller"

# mutate variables to have nicer names
mutate_variables <- function(df, with_metric=TRUE){
  df <- df %>%
    mutate(
      de.prob=increasing_value_factor(de.prob),
      de.facLoc=increasing_value_factor(de.facLoc),
      lib.loc=factor(
        case_when(lib.loc == 9.3 ~ "low", lib.loc == 10.3 ~ "medium",lib.loc == 11 ~ "large"),
        levels=c("low", "medium", "large")
      )
    )
  if(with_metric){
    df <- df %>%
      mutate(
        Metric=case_when(
          Metric == "n_markers" ~ "Selected",
          Metric == "tp_markers" ~ "True Positives",
          Metric == "fp_markers" ~ "False Positives",
          Metric == "tn_markers" ~ "True Negatives",
          Metric == "fn_markers" ~ "False Negatives",
          Metric == "accuracy" ~ "Accuracy",
          Metric == "precision" ~ "Precision",
          Metric == "recall" ~ "Recall",
          Metric == "true_negative_rate" ~ "True Negative Rate",
          Metric == "F1_score" ~ "F1 Score",
          Metric == "mean_cell_stability_score" ~ "Mean Cell Stability",
          Metric == "top1_agreement" ~ "Top 1 Agreement",
          Metric == "top5_agreement" ~ "Top 5 Agreement",
          Metric == "top10_agreement" ~ "Top 10 Agreement",
          TRUE ~ Metric
        )
      )
  }
  return(df)
}

# assert that id1 < id2
sort_ids <- function(df, global){
  if(global) return(df)
  df %>%
    mutate(
      temp_id=id1,
      id1=min(id1, id2),
      id2=max(temp_id, id2)
    ) %>%
    select(-temp_id)
}

# join metrics of different approaches and get df in log format for plotting
get_plot_df <- function(stabcell_metrics, ddip_metrics, countsplit_metrics, global=TRUE){
  df_temp_stabcell <- stabcell_metrics %>%
    select(-c(
      "selection_freq_thresh",
      if(global) "truth_id_aligned" else c("truth_id1_aligned", "truth_id2_aligned")
    )) %>%
    reshape2::melt(
      id.vars=c(
        if(global) "id" else c("id1", "id2"),
        "complementary_pairs", "n_celltypes", "n_iter", "lib.loc", "de.prob",
        "de.facLoc", "q", "PFER_max", "exceeds_one", "assumption", "p_val_thresh", "iter"
      ),
      variable.name="Metric",
      value.name="StabCell"
    ) %>%
    mutate_variables

  df_temp_ddip <- ddip_metrics %>%
    select(-c(
      if(global) c("method_id", "truth_id_aligned") else c("method_id1", "method_id2", "truth_id1_aligned", "truth_id2_aligned")
    )) %>%
    reshape2::melt(
      id.vars=c(
        if(global) "id" else c("id1", "id2"),
        "complementary_pairs", "n_celltypes", "n_iter", "lib.loc", "de.prob",
        "de.facLoc", "iter"
      ),
      variable.name="Metric",
      value.name="Double Dipping"
    ) %>%
    mutate_variables

  df_temp_countsplit <- countsplit_metrics %>%
    select(-c(
      if(global) c("method_id", "truth_id_aligned") else c("method_id1", "method_id2", "truth_id1_aligned", "truth_id2_aligned")
    )) %>%
    reshape2::melt(
      id.vars=c(
        if(global) "id" else c("id1", "id2"),
        "complementary_pairs", "n_celltypes", "n_iter", "lib.loc", "de.prob",
        "de.facLoc", "iter"
      ),
      variable.name="Metric",
      value.name="Countsplit"
    ) %>%
    mutate_variables

  df_plot <- df_temp_stabcell %>%
    left_join(
      df_temp_ddip,
      by=c(
        if(global) "id" else c("id1", "id2"),
        "complementary_pairs", "n_celltypes", "n_iter", "lib.loc", "de.prob", "de.facLoc", "iter", "Metric"
      )
    ) %>%
    left_join(
      df_temp_countsplit,
      by=c(
        if(global) "id" else c("id1", "id2"),
        "complementary_pairs", "n_celltypes", "n_iter", "lib.loc", "de.prob", "de.facLoc", "iter", "Metric"
      )
    ) %>%
    reshape2::melt(
      id.vars=c(
        if(global) "id" else c("id1", "id2"),
        "complementary_pairs", "n_celltypes", "n_iter", "lib.loc", "de.prob",
        "de.facLoc", "q", "PFER_max", "exceeds_one", "assumption", "p_val_thresh", "Metric", "iter"
      ),
      variable.name="Method",
      value.name="value"
    ) %>%
    filter(!(Method %in% c("Double Dipping", "Countsplit") & Metric == "Mean Cell Stability")) %>%
    mutate(
      assumption=factor(
        case_when(assumption == "r_concavity" ~ "r-concavity",  TRUE ~ assumption),
        levels=c("none", "unimodal", "r-concavity")
      )
    )

  return(df_plot)
}

df_plot_global <- get_plot_df(df_stabcell_global, df_ddip_global, df_countsplit_global, global=TRUE)
df_clust <- df_clust %>% mutate(
  de.prob=increasing_value_factor(de.prob),
  de.facLoc=increasing_value_factor(de.facLoc),
  lib.loc=factor(
    case_when(lib.loc == 9.3 ~ "low", lib.loc == 10.3 ~ "medium",lib.loc == 11 ~ "large"),
    levels=c("low", "medium", "large")
  )
)

if(use_simulation == "NCluster_Simulation"){
  # all selected markers are FPs if n_celltypes == 1
  df_stabcell_local <- df_stabcell_local %>% mutate(fp_markers=ifelse(n_celltypes == 1, n_markers, fp_markers))
  df_ddip_local <- df_ddip_local %>% mutate(fp_markers=ifelse(n_celltypes == 1, n_markers, fp_markers))
  df_plot_local  <- get_plot_df(df_stabcell_local,  df_ddip_local,  df_countsplit_local,  global=FALSE)
}

# view settings with PFER_max breach
df_plot_global %>%
  filter(
    p_val_thresh == 0.05,
    Method == "StabCell",
    Metric == "False Positives"
  ) %>%
  group_by(n_celltypes, assumption, PFER_max, de.prob, de.facLoc, lib.loc) %>%
  summarize(m = mean(value), .groups="drop") %>%
  filter(m > as.numeric(as.character(PFER_max))) %>%
  arrange(n_celltypes, PFER_max, de.prob, de.facLoc, lib.loc, assumption)

#########################################
#### investigate cluster performance ####
#########################################

# number of clusters per scenario
{
  ggplot(
    data=df_clust,
    mapping=aes(x=lib.loc, y=n_clusters_stabcell)
  ) +
    geom_boxplot() +
    facet_nested(de.facLoc + de.prob ~ n_celltypes, labeller=custom_labeller)
}

{
df_temp <- df_clust %>%
  select(-n_clusters_true, -ari_ddip, -ari_stabcell, -mean_cell_stability) %>%
  tidyr::pivot_longer(
    cols=c("n_clusters_ddip", "n_clusters_countsplit", "n_clusters_stabcell"),
    names_to="Method",
    values_to="n_clusters"
  ) %>%
  mutate(
    Method=factor(case_when(
      Method == "n_clusters_ddip" ~ "Double Dipping",
      Method == "n_clusters_countsplit" ~ "Countsplit",
      Method == "n_clusters_stabcell" ~ "StabCell",
      TRUE ~ Method
    ), levels=c("Double Dipping", "Countsplit", "StabCell"))
  )

p1 <- ggplot(
  data=df_temp %>% filter(n_celltypes == 1),
  mapping=aes(x=n_clusters, color=Method, fill=Method)
) +
  geom_histogram(position="identity", alpha=0.2, binwidth=1) +
  geom_vline(xintercept=1, linetype="dashed") +
  # geom_density(alpha=0.2, bounds=c(0, 12)) +
  facet_nested(n_celltypes ~ lib.loc, labeller=custom_labeller) +
  labs(
    x=element_blank(),
    y=element_blank(),
    color=element_blank(),
    fill=element_blank()
  ) +
  theme(
    plot.background=element_blank(),
    panel.background=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(color="lightgrey"),
    legend.position="top",
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
  ) +
  coord_cartesian(xlim=c(1, 12)) +
  scale_fill_manual(values=my_colors[c(1,2,4)]) +
  scale_color_manual(values=my_colors[c(1,2,4)])

p2 <- ggplot(
  data=df_temp %>% filter(n_celltypes != 1),
  mapping=aes(x=n_clusters, color=Method, fill=Method)
) +
  geom_histogram(position="identity", alpha=0.2, binwidth=1) +
  geom_vline(
    xintercept=if(use_simulation == "2ClusterSimulation") 2 else 10,
    linetype="dashed"
  ) +
  # geom_density(alpha=0.2, bounds=c(0, 12)) +
  facet_nested(n_celltypes + de.prob + de.facLoc ~ lib.loc, labeller=custom_labeller, scales="free_y") +
  labs(
    x="Number of Clusters",
    y=element_blank()
  ) +
  theme(
    plot.background=element_blank(),
    panel.background=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(color="lightgrey"),
    strip.background.x=element_blank(),
    strip.text.x=element_blank(),
    legend.position="none"
  ) +
  coord_cartesian(xlim=c(1, 12)) +
  scale_x_continuous(breaks=2*0:6) +
  scale_fill_manual(values=my_colors[c(1,2,4)]) +
  scale_color_manual(values=my_colors[c(1,2,4)])

p <- (p1 / p2) + plot_layout(heights=c(1, 4))
plot(p)

df_temp %>%
  group_by(n_celltypes, Method) %>%
  summarize(
    min=min(n_clusters),
    mean=mean(n_clusters),
    max=max(n_clusters),
    .groups="drop"
  )

ggsave(
  plot=p,
  filename="cluster_histogram_numofclusters.jpg",
  path=paste0("plots/", use_simulation),
  width=19,
  height=19,
  units="cm",
  dpi=350
)
ggsave(
  plot=p,
  filename="cluster_histogram_numofclusters.tif",
  path=paste0("plots/", use_simulation),
  width=17,
  height=17,
  units="cm",
  dpi=350
)
}

# plot adjusted rand index
n_celltypes <- if(use_simulation == "2Cluster_Simulation") 2 else 10
{
p_ari <- ggplot(
  df_clust %>%
    filter(n_celltypes == !!n_celltypes) %>%
    select(-mean_cell_stability) %>%
    reshape2::melt(
      id.vars=c("complementary_pairs", "n_celltypes", "n_iter", "lib.loc", "de.prob", "de.facLoc", "iter"),
      variable.name="Method",
      value.name="ARI"
    ) %>%
    filter(stringr::str_detect(Method, "^ari")) %>%
    mutate(
      Method=factor(case_when(
        Method == "ari_ddip" ~ "Double Dipping",
        Method == "ari_countsplit" ~ "Countsplit",
        Method == "ari_stabcell" ~ "StabCell",
        TRUE ~ Method
      ), levels=c("Double Dipping", "Countsplit", "StabCell")),
      lib.loc=factor(
        case_when(lib.loc == 9.3 ~ "low", lib.loc == 10.3 ~ "medium", lib.loc == 11 ~ "large", TRUE ~ lib.loc),
        levels=c("low", "medium", "large")
      )
    ),
  aes(x=lib.loc, y=ARI, color=Method)
) +
  geom_boxplot() +
  facet_grid(de.prob ~ de.facLoc, labeller=custom_labeller) +
  theme(
    legend.position="bottom",
    plot.background=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(color="lightgrey"),
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
  ) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) +
  coord_cartesian(ylim=c(0, 1)) +
  scale_color_manual(values=my_colors[c(1, 2, 4)]) +
  labs(
    title=paste0("Cluster Quality (", n_celltypes, " Cell Type(s))"),
    x="Library Size",
    y="Adjusted Rand Index",
    color=element_blank()
  )
plot(p_ari)

ggsave(
  plot=p_ari + labs(title=element_blank()) + theme(legend.position="top"),
  filename="ARI_boxplot.eps",
  path=paste0("plots/", use_simulation),
  width=21,
  height=11,
  units="cm",
  dpi=350
)
ggsave(
  plot=p_ari + labs(title=element_blank()) + theme(legend.position="top"),
  filename="ARI_boxplot.tif",
  path=paste0("plots/", use_simulation),
  width=21,
  height=11,
  units="cm",
  dpi=350
)
}


# plot cell stability scores
{
p1 <- ggplot(
  data=df_clust %>% filter(n_celltypes == 1),
  mapping=aes(x=lib.loc, y=mean_cell_stability)
) +
  # geom_violin() +
  geom_boxplot() +
  facet_nested(. ~ n_celltypes, labeller=custom_labeller) +
  labs(
    x=element_blank(),
    y="Mean Cell Stability Score"
  ) +
  theme(
    plot.background=element_blank(),
    panel.background=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(color="lightgrey"),
    legend.position="none"
  ) +
  coord_cartesian(ylim=c(0, 1))

p2 <- ggplot(
  data=df_clust %>% filter(n_celltypes != 1),
  mapping=aes(x=lib.loc, y=mean_cell_stability)
) +
  # geom_violin() +
  geom_boxplot() +
  facet_nested(. ~ n_celltypes + de.prob + de.facLoc, labeller=custom_labeller) +
  labs(
    x=element_blank(),
    y=element_blank()
  ) +
  theme(
    plot.background=element_blank(),
    panel.background=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(color="lightgrey"),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank(),
    legend.position="none"
  ) +
  coord_cartesian(ylim=c(0, 1))

x_label <- ggplot() +
  theme_void() +
  annotate("text", x=0.5, y=0.5, label="Library Size", size=5)

p <- ((p1 + p2) + plot_layout(widths=c(1, 4))) / x_label + plot_layout(heights=c(12, 1))
plot(p)

ggsave(
  plot=p,
  filename="cluster_boxplot_cellstability.eps",
  path=paste0("plots/", use_simulation),
  width=22,
  height=11,
  units="cm",
  dpi=350
)
ggsave(
  plot=p,
  filename="cluster_boxplot_cellstability.tif",
  path=paste0("plots/", use_simulation),
  width=22,
  height=11,
  units="cm",
  dpi=350
)
}

{
p_cell_stab <- ggplot(
  bind_rows(
    df_clust %>% filter(n_celltypes %in% c(2, 10)),
    df_clust %>% filter(n_celltypes == 1),

    # copy 1 cell simulation to other settings
    df_clust %>% filter(n_celltypes == 1) %>% mutate(de.prob=factor(0.05)),
    df_clust %>% filter(n_celltypes == 1) %>% mutate(de.facLoc=factor(0.4)),
    df_clust %>% filter(n_celltypes == 1) %>% mutate(de.prob=factor(0.05), de.facLoc=factor(0.4)),
  )%>%
    mutate(lib.loc=factor(lib.loc)),
  aes(x=lib.loc, y=mean_cell_stability, color=as.factor(n_celltypes))
) +
  geom_boxplot() +
  facet_grid(de.prob ~ de.facLoc, labeller=custom_labeller) +
  theme(
    legend.position="top",
    plot.background=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(color="lightgrey"),
    panel.grid.minor=element_blank(),
    panel.background=element_blank()
  ) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) +
  coord_cartesian(ylim=c(0, 1)) +
  scale_color_manual(values=my_colors[6:7]) +
  labs(
    x="Library Size",
    y="Mean Cell Stability Score",
    color="Number of True Cell Types"
  )
p_cell_stab
}


# combination of ARI and cell stability for individual settings
for(de.prob in c(0.01, 0.05)){
  for(de.facLoc in c(0.1, 0.4)){
    p <- (

      # ari plot for multi cluster scenario
      ggplot(
        df_clust %>%
          filter(
            n_celltypes != 1,
            de.prob == !!de.prob,
            de.facLoc == !!de.facLoc
          ) %>%
          select(-c("mean_cell_stability", "n_clusters_true", "n_clusters_ddip", "n_clusters_countsplit", "n_clusters_stabcell")) %>%
          reshape2::melt(
            id.vars=c("complementary_pairs", "n_celltypes", "n_iter", "lib.loc", "de.prob", "de.facLoc", "iter"),
            variable.name="Method",
            value.name="ARI"
          ) %>%
          mutate(
            Method=factor(case_when(
              Method == "ari_ddip" ~ "Double Dipping",
              Method == "ari_countsplit" ~ "Countsplit",
              Method == "ari_stabcell" ~ "StabCell",
              TRUE ~ Method
            ), levels=c("Double Dipping", "Countsplit", "StabCell")
          )),
        aes(x=lib.loc, y=ARI, color=Method)
      ) +
        geom_boxplot() +
        theme(
          legend.position="top",
          plot.background=element_blank(),
          panel.background=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_line(color="lightgrey", linetype="dashed"),
          panel.grid.minor=element_blank(),

        ) +
        scale_color_manual(values=my_colors[c(1, 2, 4)]) +
        scale_y_continuous(breaks=c(0, 0.5, 1)) +
        coord_cartesian(ylim=c(0, 1)) +
        labs(
          x="Library Size",
          y="Adjusted Rand Index",
          color=element_blank()
        )

      +

        # cell stability scores from stabcell
        ggplot(
          data=df_clust %>%
            filter(
              (n_celltypes == 1)
              | ((de.prob == !!de.prob) & (de.facLoc == !!de.facLoc))
            ) %>%
            mutate(n_celltypes=ifelse(n_celltypes == 1, "1 Celltype", paste(n_celltypes, "Celltypes"))),
          mapping=aes(x=lib.loc, y=mean_cell_stability, color=as.factor(n_celltypes))
        ) +
        geom_boxplot() +
        theme(
          legend.position="top",
          plot.background=element_blank(),
          panel.background=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_line(color="lightgrey", linetype="dashed"),
          panel.grid.minor=element_blank(),

        ) +
        scale_color_manual(values=my_colors[6:7]) +
        scale_y_continuous(breaks=c(0, 0.5, 1)) +
        coord_cartesian(ylim=c(0, 1)) +
        labs(
          x="Library Size",
          y="Mean Cell Stability",
          color=element_blank()
        )

    ) +

      # subplot annotation
      plot_annotation(tag_levels='a', tag_prefix="", tag_suffix=")") &
      theme(plot.tag.position=c(0.025, 0.97),  # Top-left corner
            plot.tag=element_text(hjust=1, vjust=1))  # Align text in corner

    ggsave(
      plot=p,
      filename=paste0("cluster_boxplot__de.prob=", de.prob, "_de.facLoc=", de.facLoc, ".eps"),
      path=paste0("plots/", use_simulation),
      width=22,
      height=11,
      units="cm",
      dpi=350
    )
    ggsave(
      plot=p,
      filename=paste0("cluster_boxplot__de.prob=", de.prob, "_de.facLoc=", de.facLoc, ".tif"),
      path=paste0("plots/", use_simulation),
      width=22,
      height=11,
      units="cm",
      dpi=350
    )
    plot(p)
  }
}









########################################
#### investigate marker performance ####
########################################


#### single simulation setting, 2 metrics #####

use_celltypes <- if(use_simulation == "2Cluster_Simulation") 2 else 10
assumption_filter <- c("none", "unimodal", "r-concavity")
PFER_filter <- c(5, 25, 100)
de.prob <- 0.01
de.facLoc <- 0.1
lib.loc <- c("low", "medium", "large")
exceeds_one <- c(TRUE, FALSE, NA) # if selection freq thresh > 1 we can set it to 1 or omit the result

plot_metrics <- c("False Positives", "True Positives")

global <- TRUE
cluster_id <- 0
cluster_id2 <- NA # 1

{
  plot_list <- list()
  for(plot_metric in plot_metrics){

    # StabCell
    {
      p1_low <- ggplot(
        data=(
          if(global){
            df_plot_global %>% filter(id == !!cluster_id)
          } else{
            df_plot_local %>% filter(id1 == !!cluster_id, id2 == !!cluster_id2)
          }) %>%
          filter(
            Method == "StabCell",
            exceeds_one %in% !!exceeds_one,
            Metric == !!plot_metric,
            n_celltypes == !!n_celltypes,
            complementary_pairs == TRUE,
            PFER_max %in% !!PFER_filter,
            lib.loc %in% !!lib.loc,
            !!(if(n_celltypes == 1) TRUE else expr(
              (de.prob == !!de.prob)
              & (de.facLoc == !!de.facLoc)
            )),
          ),
        mapping=aes(x=lib.loc, y=value, color=assumption)
      ) +
        geom_boxplot() +
        facet_nested(
          de.prob + de.facLoc ~ PFER_max,
          labeller=custom_labeller
        ) +
        geom_hline(
          mapping=aes(yintercept=ifelse(Metric == "False Positives", as.numeric(as.character(PFER_max)), NaN)),
          color="black", linetype="dotted"
        ) +
        labs(x=element_blank(), y=element_blank(), color="Assumption") +
        theme(
          panel.background=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_line(linewidth=0.3, color="gray90"),
          panel.grid.minor=element_blank(),
          axis.line=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text=element_text(color="black"),
          axis.title=element_text(color="black"),
          plot.title=element_text(size=rel(1.2), face="bold", hjust=0),
          legend.position="top",
          legend.key=element_blank(),
          strip.background.x=element_rect(fill="gray90", color=NA),
          strip.background.y=element_blank(),
          strip.text.y=element_blank(),
        ) +
        scale_color_manual(values=my_colors[3:5])
    }

    # other approaches
    {p2_low <- ggplot(
      data=(
        if(global){
          df_plot_global %>% filter(id == !!cluster_id)
        } else{
          df_plot_local %>%filter(id1 == !!cluster_id, id2 == !!cluster_id2)
        }) %>%
        filter(
          Method != "StabCell",
          exceeds_one %in% !!exceeds_one,
          Metric == !!plot_metric,
          n_celltypes == !!n_celltypes,
          complementary_pairs == TRUE, # StabCell settings don't matter here --> pick arbitrary one
          PFER_max == 1,
          assumption == "none",
          lib.loc %in% !!lib.loc,
          !!(if(n_celltypes == 1) TRUE else expr(
            (de.prob == !!de.prob)
            & (de.facLoc == !!de.facLoc)
          )),
        ),
      mapping=aes(x=lib.loc, y=value, color=Method)
    ) +
        geom_boxplot() +
        facet_nested(
          de.prob + de.facLoc ~ "",
          labeller=labeller(.rows=custom_labeller, .cols=label_value)
        ) +
        scale_color_manual(values=my_colors) +
        labs(x=element_blank(), y=element_blank(), color=element_blank()) +
        theme(
          panel.background=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_line(linewidth=0.3, color="gray90"),
          panel.grid.minor=element_blank(),
          axis.line=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.text=element_text(color="black"),
          axis.title=element_text(color="black"),
          plot.title=element_text(size=rel(1.2), face="bold", hjust=0),
          legend.position="top",
          legend.key=element_blank(),
          strip.background=element_blank(),
          strip.text.y=element_blank(),
        )}

    if(plot_metric != plot_metrics[1]){
      p1_low <- p1_low + theme(
        legend.position="none",
        strip.background.x=element_blank(),
        strip.text.x=element_blank()
      )
      p2_low <- p2_low + theme(
        legend.position="none",
        strip.background.x=element_blank(),
        strip.text.x=element_blank()
      )
    }

    # set y limits
    if(use_simulation == "2Cluster_Simulation"){
      if(n_celltypes == 1){
        ylim_low <- ylim_high <- switch(
          plot_metric,
          "Selected"=c(0, 180),
          "False Positives"=c(0, 150),
          "True Positives"=c(0, 150),
          NULL
        )
      } else {
        ylim_low <- switch(
          plot_metric,
          "Selected"=c(0, 200),
          "False Positives"=c(0, 300),
          "True Positives"=c(0, 150),
          NULL
        )
        ylim_high <- switch(
          plot_metric,
          "Selected"=c(0, 750),
          "False Positives"=c(0, 110),
          "True Positives"=c(0, 150),
          NULL
        )
      }
    } else if(use_simulation == "NCluster_Simulation"){
      if(n_celltypes == 1){
        ylim_low <- ylim_high <- switch(
          plot_metric,
          "Selected"=c(0, 100),
          "False Positives"=c(0, 110),
          NULL
        )
      } else {
        ylim_low <- switch(
          plot_metric,
          "Selected"=c(0, 230),
          "False Positives"=c(0, 230),
          NULL
        )
        ylim_high <- switch(
          plot_metric,
          "Selected"=c(0, 350),
          "False Positives"=c(0, 320),
          NULL
        )
      }
    }
    p1_low  <- p1_low  + coord_cartesian(ylim=ylim_low)
    p2_low  <- p2_low  + coord_cartesian(ylim=ylim_low)


    x_label <- ggplot() +
      theme_void() +
      annotate("text", x=0.5, y=0.5, label="Library Size", size=5)

    y_label <- ggplot() +
      theme_void() +
      annotate("text", x=0.5, y=0.5, angle=90, label=plot_metric, size=5)

    if(n_celltypes == 1){
      p1 <- (y_label + (p1_low)) + plot_layout(widths=c(1, 25))
      p2 <- p2_low
    } else {
      p1 <- (y_label + (p1_low)) + plot_layout(widths=c(1, 25))
      p2 <- p2_low
      # p1 <- (y_label + (p1_low / p1_high)) + plot_layout(widths=c(1, 25))
      # p2 <- p2_low / p2_high
    }

    p <- (p1 | p2) + plot_layout(widths=c(length(PFER_filter), 1))
    if(plot_metric == plot_metrics[length(plot_metrics)]){
      p <- p / x_label + plot_layout(heights=c(1, 0.1))
    } else {
      p <- p + plot_layout(heights=c(1, 0.1))
    }

    plot_list[[length(plot_list) + 1]] <- p
  }
}
p <- (plot_list[[1]] / plot_list[[2]]) + plot_layout(nrow=2)
plot(p)

ggsave(
  plot=p,
  filename=paste0("marker_boxplot__de.prob=", de.prob, "_de.facLoc=", de.facLoc, ".eps"),
  path=paste0("plots/", use_simulation),
  width=22,
  height=16,
  units="cm",
  dpi=350
)
ggsave(
  plot=p,
  filename=paste0("marker_boxplot__de.prob=", de.prob, "_de.facLoc=", de.facLoc, ".tif"),
  path=paste0("plots/", use_simulation),
  width=22,
  height=16,
  units="cm",
  dpi=350
)


#### all settings for single metric #####

use_celltypes <- if(use_simulation == "2Cluster_Simulation") c(1, 2) else c(1, 10)
assumption_filter <- c("none", "unimodal", "r concavity")
PFER_filter <- c(1, 5, 10, 25, 50, 100)
de.prob <- c(0.01, 0.05)
de.facLoc <- c(0.1, 0.4)
lib.loc <- c("low", "medium", "large")
plot_metrics <- c(
  "Selected"
  , "False Positives"
  , "True Positives"
  , "False Negatives"
  , "True Negatives"
  , "Accuracy"
  , "Precision"
  , "Recall"
  , "True Negative Rate"
  , "F1 Score"
  , "Mean Cell Stability"
  , "Top 1 Agreement"
  , "Top 5 Agreement"
  , "Top 10 Agreement"
)[2:3]
exceeds_one <- c(TRUE, FALSE, NA) # if selection freq thresh > 1 we can set it to 1 or omit the result

global <- TRUE
cluster_id <- 0
cluster_id2 <- 1
{
for(n_celltypes in use_celltypes){
  for(plot_metric in plot_metrics){
    if(n_celltypes == 1 & plot_metric == "True Positives") { next }

    # StabCell
    {
    p1_low <- ggplot(
      data=(
        if(global){
          df_plot_global %>% filter(id == !!cluster_id)
        } else{
          df_plot_local %>%filter(id1 == !!cluster_id, id2 == !!cluster_id2)
        }) %>%
        filter(
          Method == "StabCell",
          de.prob == 0.01,
          exceeds_one %in% !!exceeds_one,
          Metric == !!plot_metric,
          n_celltypes == !!n_celltypes,
          complementary_pairs == TRUE,
          PFER_max %in% !!PFER_filter,
          lib.loc %in% !!lib.loc,
          !!(if(n_celltypes == 1) TRUE else expr(
            (de.prob %in% !!de.prob)
            & (de.facLoc %in% !!de.facLoc)
          )),
        ),
      mapping=aes(x=lib.loc, y=value, color=assumption)
    ) +
      geom_boxplot() +
      facet_nested(
        de.prob + de.facLoc ~ PFER_max,
        labeller=custom_labeller
      ) +
      geom_hline(
        mapping=aes(yintercept=ifelse(Metric == "False Positives", as.numeric(as.character(PFER_max)), NaN)),
        color="black", linetype="dotted"
      ) +
      labs(x=element_blank(), y=element_blank(), color="Assumption") +
      theme(
        panel.background=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(linewidth=0.3, color="gray90"),
        panel.grid.minor=element_blank(),
        axis.line=element_blank(),
        axis.ticks.x=if(n_celltypes == 1) element_line(color="black") else element_blank(),
        axis.text.x=if(n_celltypes == 1) element_text(color="black", angle=270, vjust=0.5, hjust=0) else element_blank(),
        axis.ticks.y=element_blank(),
        axis.text=element_text(color="black"),
        axis.title=element_text(color="black"),
        plot.title=element_text(size=rel(1.2), face="bold", hjust=0),
        legend.position="top",
        legend.key=element_blank(),
        strip.background.x=element_rect(fill="gray90", color=NA),
        strip.background.y=element_blank(),
        strip.text.y=element_blank(),
      ) +
      scale_color_manual(values=my_colors[3:5])
    }
    {
    if(n_celltypes == 1){
      p1_high <- ggplot() + theme_void()
    } else {
      p1_high <- ggplot(
        data=(
          if(global){
            df_plot_global %>% filter(id == !!cluster_id)
          } else{
            df_plot_local %>%filter(id1 == !!cluster_id, id2 == !!cluster_id2)
          }) %>%
          filter(
            Method == "StabCell",
            de.prob == 0.05,
            exceeds_one %in% !!exceeds_one,
            Metric == !!plot_metric,
            n_celltypes == !!n_celltypes,
            complementary_pairs == TRUE,
            PFER_max %in% !!PFER_filter,
            lib.loc %in% !!lib.loc,
            !!(if(n_celltypes == 1) TRUE else expr(
              (de.prob %in% !!de.prob)
              & (de.facLoc %in% !!de.facLoc)
            )),
          ),
        mapping=aes(x=lib.loc, y=value, color=assumption)
      ) +
        geom_boxplot() +
        facet_nested(
          de.prob + de.facLoc ~ PFER_max,
          labeller=custom_labeller
        ) +
        geom_hline(
          mapping=aes(yintercept=ifelse(Metric == "False Positives", as.numeric(as.character(PFER_max)), NaN)),
          color="black", linetype="dotted"
        ) +
        labs(x=element_blank(), y=element_blank(), color="Assumption") +
        theme(
          panel.background=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_line(linewidth=0.3, color="gray90"),
          panel.grid.minor=element_blank(),
          axis.line=element_blank(),
          axis.ticks.x=element_line(color="black"),
          axis.ticks.y=element_blank(),
          axis.text.x=element_text(color="black", angle=270, vjust=0.5, hjust=0),
          axis.text.y=element_text(color="black"),
          axis.title=element_text(color="black"),
          plot.title=element_text(size=rel(1.2), face="bold", hjust=0),
          legend.position="none",
          strip.background.x=element_blank(),
          strip.background.y=element_blank(),
          strip.text.x=element_blank(),
          strip.text.y=element_blank(),
        ) +
        scale_color_manual(values=my_colors[3:5])
    }
    }


    # other approaches
    {p2_low <- ggplot( # low de.prob to set individual ylim
      data=(
        if(global){
          df_plot_global %>% filter(id == !!cluster_id)
        } else{
          df_plot_local %>%filter(id1 == !!cluster_id, id2 == !!cluster_id2)
        }) %>%
        filter(
          Method != "StabCell",
          de.prob == 0.01,
          exceeds_one %in% !!exceeds_one,
          Metric == !!plot_metric,
          n_celltypes == !!n_celltypes,
          complementary_pairs == TRUE,
          PFER_max == 1,
          assumption == "none",
          lib.loc %in% !!lib.loc,
          !!(if(n_celltypes == 1) TRUE else expr((de.facLoc %in% !!de.facLoc))),
        ),
      mapping=aes(x=lib.loc, y=value, color=Method)
    ) +
      geom_boxplot() +
      facet_nested(
        de.prob + de.facLoc ~ '',
        labeller=labeller(.rows=custom_labeller, .cols=label_value)
      ) +
      labs(x=element_blank(), y=element_blank(), color=element_blank()) +
      theme(
        panel.background=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(linewidth=0.3, color="gray90"),
        panel.grid.minor=element_blank(),
        axis.line=element_blank(),
        axis.ticks.x=if(n_celltypes == 1) element_line(color="black") else element_blank(),
        axis.text.x=if(n_celltypes == 1) element_text(color="black", angle=270, vjust=0.5, hjust=0) else element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.text=element_text(color="black"),
        axis.title=element_text(color="black"),
        plot.title=element_text(size=rel(1.2), face="bold", hjust=0),
        legend.position="top",
        legend.key=element_blank(),
        strip.background.x=element_blank(),
        strip.background.y=if(n_celltypes == 1) element_blank() else element_rect(fill="gray90", color=NA),
        strip.text.y=if(n_celltypes == 1) element_blank() else element_text(),
      ) +
        scale_color_manual(values=my_colors)}
    {
    if(n_celltypes == 1){
      p2_high <- ggplot() + theme_void()
    } else {
      p2_high <- ggplot(  # high de.prob to set individual ylim
        data=(
          if(global){
            df_plot_global %>% filter(id == !!cluster_id)
          } else{
            df_plot_local %>%filter(id1 == !!cluster_id, id2 == !!cluster_id2)
          }) %>%
          filter(
            Method != "StabCell",
            de.prob == 0.05,
            exceeds_one %in% !!exceeds_one,
            Metric == !!plot_metric,
            n_celltypes == !!n_celltypes,
            complementary_pairs == TRUE,
            PFER_max == 1,
            assumption == "none",
            lib.loc %in% !!lib.loc,
            !!(if(n_celltypes == 1) TRUE else expr((de.facLoc %in% !!de.facLoc))),
          ),
        mapping=aes(x=lib.loc, y=value, color=Method)
      ) +
          geom_boxplot() +
          facet_nested(
            de.prob + de.facLoc ~ '',
            labeller=labeller(.rows=custom_labeller, .cols=label_value)
          ) +
          labs(x=element_blank(), y=element_blank(), color=element_blank()) +
          theme(
            panel.background=element_blank(),
            panel.grid.major.x=element_blank(),
            panel.grid.major.y=element_line(linewidth=0.3, color="gray90"),
            panel.grid.minor=element_blank(),
            axis.line=element_blank(),
            axis.ticks.x=element_line(color="black"),
            axis.ticks.y=element_blank(),
            axis.text.x=element_text(color="black", angle=270, vjust=0.5, hjust=0),
            axis.text.y=element_blank(),
            axis.title=element_text(color="black"),
            plot.title=element_text(size=rel(1.2), face="bold", hjust=0),
            legend.position="none",
            legend.key=element_blank(),
            strip.background=element_rect(fill="gray90", color=NA),
            strip.background.x=element_blank(),
            strip.text.x=element_blank(),
          ) +
        scale_color_manual(values=my_colors)
    }
    }



    # set y limits
    if(use_simulation == "2Cluster_Simulation"){
      if(n_celltypes == 1){
        ylim_low <- ylim_high <- switch(
          plot_metric,
          "Selected"=c(0, 180),
          "False Positives"=c(0, 300),
          "True Positives"=c(0, 180),
          NULL
        )
      } else {
        ylim_low <- switch(
          plot_metric,
          "Selected"=c(0, 200),
          "False Positives"=c(0, 300),
          "True Positives"=c(0, 170),
          NULL
        )
        ylim_high <- switch(
          plot_metric,
          "Selected"=c(0, 750),
          "False Positives"=c(0, 110),
          "True Positives"=c(0, 750),
          NULL
        )
      }
    } else if(use_simulation == "NCluster_Simulation"){
      if(n_celltypes == 1){
        ylim_low <- ylim_high <- switch(
          plot_metric,
          "Selected"=c(0, 100),
          "False Positives"=c(0, 110),
          "True Positives"=c(0, 10),
          NULL
        )
      } else {
        ylim_low <- switch(
          plot_metric,
          "Selected"=c(0, 230),
          "False Positives"=c(0, 230),
          "True Positives"=c(0, 60),
          NULL
        )
        ylim_high <- switch(
          plot_metric,
          "Selected"=c(0, 350),
          "False Positives"=c(0, 320),
          "True Positives"=c(0, 330),
          NULL
        )
      }
    }
    p1_low  <- p1_low  + coord_cartesian(ylim=ylim_low)
    p1_high <- p1_high + coord_cartesian(ylim=ylim_high)
    p2_low  <- p2_low  + coord_cartesian(ylim=ylim_low)
    p2_high <- p2_high + coord_cartesian(ylim=ylim_high)


    # plot size
    if(n_celltypes == 1){
      width <- 21
      height <- 17
    } else {
      width <- 21
      height <- 21
    }

    x_label <- ggplot() +
      theme_void() +
      annotate("text", x=0.5, y=0.5, label="Library Size", size=5)

    y_label <- ggplot() +
      theme_void() +
      annotate("text", x=0.5, y=0.5, angle=90, label=plot_metric, size=5)

    if(n_celltypes == 1){
      p1 <- (y_label + (p1_low)) + plot_layout(widths=c(1, 25))
      p2 <- p2_low + guides(color=guide_legend(nrow=2))
    } else {
      p1 <- (y_label + (p1_low / p1_high)) + plot_layout(widths=c(1, 25))
      p2 <- p2_low / p2_high
    }
    p <- (
      ((p1 | p2) + plot_layout(widths=c(length(PFER_filter), 1))) /
      x_label + plot_layout(heights=c(1, 0.05))
    )

    ggsave(
      plot=p,
      filename=paste0("marker_boxplot__n_celltypes=", n_celltypes, "_metric=", gsub(" ", "", plot_metric), "_", ifelse(global, "global", "local"),".eps"),
      path=paste0("plots/", use_simulation),
      width=width,
      height=height,
      units="cm",
      dpi=350
    )
    ggsave(
      plot=p,
      filename=paste0("marker_boxplot__n_celltypes=", n_celltypes, "_metric=", gsub(" ", "", plot_metric), "_", ifelse(global, "global", "local"),".tif"),
      path=paste0("plots/", use_simulation),
      width=width,
      height=height,
      units="cm",
      dpi=350
    )
    plot(p)
  }
}
}

