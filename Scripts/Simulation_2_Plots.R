setwd("/mnt/mikasa_data/nlueck/StabCell_paper")

library(ggplot2)
library(ggh4x)
library(patchwork)
library(dplyr)

source("Scripts/helpers/plot.R")


use_simulation <- c(
  "2Cluster_Simulation",
  "NCluster_Simulation"
)[1]

plot_dir <- paste0("plots/", use_simulation)
save_plots <- FALSE
save_format <- "eps"

load(paste0("results/", use_simulation, "/metrics/results_metrics_total.Rdata"))

df_markers <- bind_rows(
  df_oracle_global %>% mutate(method="Oracle"),
  df_ddip_global %>% mutate(method="Double Dipping"),
  df_countsplit_global %>% mutate(method="Countsplit"),
  df_clusterde_global %>% mutate(method="ClusterDE"),
  df_tntest_global %>% mutate(method="TN Test"),
  df_stabcell_global %>% mutate(method="StabCell")
) %>% mutate(
    method=factor(method, levels=c("StabCell", "Oracle", "Double Dipping", "Countsplit", "ClusterDE", "TN Test")),
    assumption=factor(
      ifelse(assumption == "r_concavity", "r-concavity", as.character(assumption)),
      levels=c("none", "unimodal", "r-concavity")
    )
  ) %>%
  filter(
    id == "0",
    method != "Oracle",
    p_val_thresh == 0.05 | method != "StabCell"
  )

df_clust <- df_clust %>%
  tidyr::pivot_longer(
    cols=matches("^mean_cell_stability$|^(ari|n_clusters)_(ddip|countsplit|tntest|stabcell)$"),
    names_pattern="^(mean_cell_stability|ari|n_clusters)(?:_(ddip|countsplit|tntest|stabcell))?$",
    names_to=c("Metric", "Method"),
    values_to="value"
  ) %>%
  mutate(
    lib.size=as.factor(lib.size),
    Method=factor(case_when(
      Metric == "mean_cell_stability" ~ "StabCell",
      Method == "ddip" ~ "Double Dipping\nClusterDE",
      Method == "countsplit" ~ "Countsplit",
      Method == "tntest" ~ "TN Test",
      Method == "stabcell" ~ "StabCell",
      .default=Method
    ), levels=c("StabCell", "Double Dipping\nClusterDE", "Countsplit", "TN Test"))
  )

# see which settings do not have 100 sets of iterations
# -> should be empty
print(
  df_clust %>%
    group_by(n_celltypes, lib.size, de.prob, de.facLoc, Metric, Method) %>%
    summarize(n=n(), .groups="drop") %>%
    filter(n != 100),
  n=500
)



#########################################
#### investigate cluster performance ####
#########################################

# number of clusters boxplot
ggplot(
  data=df_clust %>% filter(Metric == "n_clusters"),
  mapping=aes(x=lib.size, y=value, color=Method)
) +
  geom_boxplot() +
  facet_nested(de.facLoc + de.prob ~ n_celltypes, labeller=custom_labeller)


# number of clusters histogram
{
p1 <- ggplot(
  data=df_clust %>%
    filter(
      n_celltypes == 1,
      Metric == "n_clusters"
    ),
  mapping=aes(x=value, color=Method, fill=Method)
) +
  geom_histogram(position="dodge2", alpha=0.2, binwidth=1) +
  geom_vline(mapping=aes(xintercept=n_celltypes), linetype="dashed") +
  # geom_density(alpha=0.2, bounds=c(0, 12)) +
  facet_nested(n_celltypes ~ lib.size, labeller=custom_labeller) +
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
  scale_fill_manual(values=my_colors) +
  scale_color_manual(values=my_colors)

p2 <- ggplot(
  data=df_clust %>%
    filter(
      n_celltypes != 1,
      Metric == "n_clusters"
    ),
  mapping=aes(x=value, color=Method, fill=Method)
) +
  geom_histogram(position="dodge2", alpha=0.2, binwidth=1) +
  geom_vline(mapping=aes(xintercept=n_celltypes), linetype="dashed") +
  # geom_density(alpha=0.2, bounds=c(0, 12)) +
  facet_nested(n_celltypes + de.prob + de.facLoc ~ lib.size, labeller=custom_labeller, scales="free_y") +
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
  scale_fill_manual(values=my_colors) +
  scale_color_manual(values=my_colors)

p <- (p1 / p2) + plot_layout(heights=c(1, 4))
plot(p)

df_clust %>%
  filter(Metric == "n_clusters") %>%
  group_by(n_celltypes, Method) %>%
  summarize(
    min=min(value),
    mean=mean(value),
    max=max(value),
    .groups="drop"
  )

if(save_plots){
  ggsave(
    plot=p,
    filename=paste0("cluster_histogram_numofclusters.", save_format),
    path=plot_dir,
    width=19,
    height=19,
    units="cm",
    dpi=350,
    device=dynamic_device(save_format)
  )
}
}

# plot adjusted rand index
{
p_ari <- ggplot(
  data=df_clust %>%
    filter(
      n_celltypes > 1, # ari for one true cluster makes no sense
      Metric == "ari"
    ),
  mapping=aes(x=lib.size, y=value, color=Method)
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
  scale_color_manual(values=my_colors) +
  labs(
    x="Library Size",
    y="Adjusted Rand Index",
    color=element_blank()
  )
plot(p_ari)

if(save_plots){
  ggsave(
    plot=p_ari + labs(title=element_blank()) + theme(legend.position="top"),
    filename=paste0("ARI_boxplot.", save_format),
    path=plot_dir,
    width=21,
    height=11,
    units="cm",
    dpi=350,
    device=dynamic_device(save_format)
  )
}
}


# plot cell stability scores
{
p1 <- ggplot(
  data=df_clust %>% filter(n_celltypes == 1, Metric == "mean_cell_stability"),
  mapping=aes(x=lib.size, y=value, color=Method)
) +
  # geom_violin() +
  geom_boxplot() +
  facet_nested(. ~ n_celltypes, labeller=custom_labeller) +
  labs(
    x=element_blank(),
    y="Mean Cell Stability Score"
  ) +
  scale_color_manual(values=my_colors) +
  theme(
    plot.background=element_blank(),
    panel.background=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(color="lightgrey"),
    legend.position="none"
  ) +
  coord_cartesian(ylim=c(0, 1))

p2 <- ggplot(
  data=df_clust %>% filter(n_celltypes != 1, Metric == "mean_cell_stability"),
  mapping=aes(x=lib.size, y=value, color=Method)
) +
  # geom_violin() +
  geom_boxplot() +
  facet_nested(. ~ n_celltypes + de.prob + de.facLoc, labeller=custom_labeller) +
  labs(
    x=element_blank(),
    y=element_blank()
  ) +
  scale_color_manual(values=my_colors) +
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

if(save_plots){
  ggsave(
    plot=p,
    filename=paste0("cluster_boxplot_cellstability.", save_format),
    path=plot_dir,
    width=22,
    height=11,
    units="cm",
    dpi=350,
    device=dynamic_device(save_format)
  )
}
}


# combination of ARI and cell stability for individual settings
for(de.prob in c(0.01, 0.05)){
  for(de.facLoc in c(0.1, 0.4)){
    p <- (

      # ari plot for multi cluster scenario
      ggplot(
        df_clust %>%
          filter(
            Metric == "ari",
            n_celltypes != 1,
            de.prob == !!de.prob,
            de.facLoc == !!de.facLoc
          ),
        aes(x=lib.size, y=value, color=Method)
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
        scale_color_manual(values=my_colors) +
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
              Metric == "mean_cell_stability",
              (n_celltypes == 1) | ((de.prob == !!de.prob) & (de.facLoc == !!de.facLoc))
            ) %>%
            mutate(n_celltypes=ifelse(n_celltypes == 1, "1 Celltype", paste(n_celltypes, "Celltypes"))),
          mapping=aes(x=lib.size, y=value, color=as.factor(n_celltypes))
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
        scale_color_manual(values=my_colors) +
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
      theme(plot.tag.position=c(0.07, 0.97),  # Top-left corner
            plot.tag=element_text(hjust=1, vjust=1))  # Align text in corner

    plot(p)
  
    if(save_plots){
      ggsave(
        plot=p,
        filename=paste0("cluster_boxplot__de.prob=", de.prob, "_de.facLoc=", de.facLoc, ".", save_format),
        path=plot_dir,
        width=22,
        height=11,
        units="cm",
        dpi=350,
        device=dynamic_device(save_format)
      )
    }
  }
}




########################################
#### investigate marker performance ####
########################################

# view settings with PFER_max breach
df_markers %>%
  filter(method == "StabCell") %>%
  group_by(n_celltypes, id, assumption, PFER_max, de.prob, de.facLoc, lib.size) %>%
  summarize(
    fp_median=median(fp_markers),
    fp_mean=mean(fp_markers),
    fp_se=sd(fp_markers) / sqrt(n()),
    fp_low.ci=fp_mean - 1.96 * fp_se,
    fp_up.ci=fp_mean + 1.96 * fp_se,
    breach=case_when(fp_low.ci > as.numeric(as.character(dplyr::first(PFER_max))) ~ "X", TRUE ~ ""),
    .groups="drop"
  ) %>%
  filter(fp_mean > as.numeric(as.character(PFER_max))) %>%
  arrange(n_celltypes, PFER_max, de.prob, de.facLoc, lib.size, assumption)


# compare FP & TP in single setting
PFER_max_filter <- c("5", "25", "100")

for(n_celltypes in unique(df_markers$n_celltypes)){
  for(de.prob in unique(df_markers$de.prob)){
    for(de.facLoc in unique(df_markers$de.facLoc)){
      
      if(n_celltypes == 1 && (de.prob != 0.01 || de.facLoc != 0.1)) next
      
      (
        p <- plot_marker_metric(
          yvar="fp_markers",
          data=df_markers %>%
            filter(
              n_celltypes == !!n_celltypes,
              de.prob == if(!!n_celltypes == 1) 0.01 else !!de.prob,
              de.facLoc == if(!!n_celltypes == 1) 0.1 else !!de.facLoc,
              is.na(PFER_max) | PFER_max %in% !!PFER_max_filter
            ),
          plot_mean=FALSE
        ) +
          theme(
            axis.title.x=element_blank(),
            strip.background.y=element_blank(),
            strip.text.y=element_text( # element_blank() generates warning???
              colour=NA,
              size=0,
              margin=margin(0),
              lineheight=0
            ),
          ) 
      )
      if(n_celltypes > 1){
        p <- p /
          (
            plot_marker_metric(
              yvar="tp_markers",
              data=df_markers %>%
                filter(
                  n_celltypes == !!n_celltypes,
                  de.prob == !!de.prob,
                  de.facLoc == !!de.facLoc,
                  is.na(PFER_max) | PFER_max %in% !!PFER_max_filter
                ),
              plot_mean=FALSE
            ) +
              theme(
                legend.position="none",
                strip.background=element_blank(),
                strip.text=element_text( # element_blank() generates warning???
                  colour=NA,
                  size=0,
                  margin=margin(0),
                  lineheight=0
                ),
              )
          )
      }
      
      plot(p)
      
      if(save_plots){
        ggsave(
          plot=p,
          filename=paste0("marker_boxplot__n_celltypes=", n_celltypes, "__de.prob=", de.prob, "_de.facLoc=", de.facLoc, ".", save_format),
          path=plot_dir,
          width=22,
          height=11,
          units="cm",
          dpi=350,
          device=dynamic_device(save_format)
        )
      }
    }
  }
}




# all settings single metric
plot_metrics <- c("fp_markers", "tp_markers")

for(n_celltypes in unique(df_markers$n_celltypes)){
  for(plot_metric in plot_metrics){
    if(n_celltypes == 1 & plot_metric == "tp_markers") { next }

    p <- plot_marker_metric(
      data=df_markers %>%
        filter(
          n_celltypes == !!n_celltypes,
          # is.na(exceeds_one) | !exceeds_one
        ),
      yvar=plot_metric,
      plot_mean=FALSE
    )
    
    if(n_celltypes == 1){
      p <- p + theme(
        strip.background.y=element_blank(),
        strip.text.y=element_text( # element_blank() generates warning???
          colour=NA,
          size=0,
          margin=margin(0),
          lineheight=0
        )
      )
    }
    
    if(use_simulation == "NCluster_Simulation" & plot_metric == "fp_markers"){
      p <- p + coord_cartesian(ylim=c(0, 300))
    }
    
    plot(p)
    
    if(save_plots){
      ggsave(
        plot=p,
        filename=paste0("marker_boxplot__n_celltypes=", n_celltypes, "__metric=", plot_metric, ".", save_format),
        path=plot_dir,
        width=22,
        height=if(n_celltypes == 1) 11 else 22,
        units="cm",
        dpi=350,
        device=dynamic_device(save_format)
      )
    }
  }
}



##########################
#### stability scores ####
##########################


if(use_simulation == "2Cluster_Simulation"){
  
  for(curve in c("ROC", "PRC")){
    p <- ggplot(
      data=df_auc %>%
        filter(curvetypes == curve) %>%
        mutate(modnames=factor(modnames, levels=c("Double Dipping", "Countsplit", "ClusterDE", "TN Test", "StabCell"))),
      mapping=aes(x=modnames, y=aucs, color=modnames)
    ) +
      geom_boxplot() +
      scale_color_manual(values=my_colors) +
      facet_nested(
        de.prob + de.facLoc ~ lib.size,
        labeller=custom_labeller
      ) +
      labs(y=paste0("AU", curve)) +
      theme(
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        
        panel.background=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(color="lightgrey"),
        
        legend.position="top",
        legend.direction="horizontal",
        legend.box="vertical",
        legend.box.just="center",
        legend.title=element_blank(),
        legend.margin=margin(t=-5, b=-2)
      )
    
    plot(p)
    
    if(save_plots){
      ggsave(
        plot=p,
        filename=paste0("marker_boxplot__n_celltypes=2__metric=AU", curve, ".", save_format),
        path=plot_dir,
        width=22,
        height=22,
        units="cm",
        dpi=350,
        device=dynamic_device(save_format)
      )
    }
  }

    load(paste0("results/", use_simulation, "/results_simulation__celltypes=2__lib.size=25000__de.prob=0.01__de.facLoc=0.1__iter=5.Rdata"))
  
    seur_obj <- CreateSeuratObject(results$sim$count_data) %>%
      NormalizeData(verbose=FALSE) %>%
      FindVariableFeatures(verbose=FALSE) %>%
      ScaleData(verbose=FALSE) %>%
      RunPCA(verbose=FALSE) %>%
      RunUMAP(dims=1:50, verbose=FALSE, seed.use=1234)
    
    seur_obj$true_idents <- as.factor(results$sim$true_idents)
    seur_obj$stability_score <- results$stabcell$cell_stability_score
    seur_obj$correct_cluster <- results$stabcell$stable_clustering == as.character(results$sim$true_idents - 1)
    seur_obj$correct_cluster_verbose <- ifelse(seur_obj$correct_cluster, "Correct\nCluster Assignment", "Incorrect\nCluster Assignment")

    df_temp <- bind_rows(
      tibble(
        `Correct Cluster`=seur_obj$correct_cluster,
        `Stability Score`=seur_obj$stability_score,
        Method="StabCell"
      ),
      tibble(
        `Correct Cluster`=results$double_dipping$identities == as.character(results$sim$true_idents - 1),
        `Stability Score`=seur_obj$stability_score,
        Method="Double Dipping"
      )
    ) %>%
      mutate(`Correct Cluster`=ifelse(`Correct Cluster`, "Correct\nCluster Assignment", "Incorrect\nCluster Assignment"))
    
    p <- (
      DimPlot(seur_obj, group.by="true_idents", shape.by="correct_cluster_verbose", pt.size=1) +
        labs(color="Cell Type", shape="") +
        scale_color_manual(values=unname(my_colors[c(3, 11)])) +
        scale_shape_manual(values=c(24, 4)) +
        guides(
          colour=guide_legend(nrow=1, order=1),
          shape=guide_legend(nrow=1, order=2)
        ) +
        coord_cartesian(xlim=c(-3, 3), ylim=c(-3, 3)) +
        theme(
          plot.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          axis.line=element_blank(),
          legend.position="top",
          legend.title=element_text(margin=margin(l=20)),
          legend.box="vertical",
          legend.direction="horizontal",
          legend.box.just="left",
          panel.background=element_blank(),
          plot.background=element_blank(),
        ) |
      FeaturePlot(seur_obj, features="stability_score") +
        labs(color="Stability Score", x="UMAP 1", y="UMAP 2") +
        scale_color_gradientn(
          colors=c("red", "#DDD", "blue"),
          values=c(0, 0.5, 1),
          breaks=c(0, 0.5, 1)
        ) +
        coord_cartesian(xlim=c(-3, 3), ylim=c(-3, 3)) +
        theme(
          plot.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          axis.line=element_blank(),
          legend.position="top",
          legend.title=element_text(margin=margin(r=20, l=20)),
          legend.key.width=unit(20, "points"),
          legend.key.height=unit(4, "points"),
          panel.background=element_blank(),
          plot.background=element_blank(),
        ) |
      ggplot(
        data=df_temp %>% filter(Method == "StabCell"),
        mapping=aes(x=`Correct Cluster`, y=`Stability Score`, color=`Correct Cluster`)
      ) +
        geom_boxplot(position=position_dodge2(width=0.75)) +
        # stat_summary(
        #   fun.data=function(x) {
        #     data.frame(
        #       y=1.1,  # max(x, na.rm=TRUE),
        #       label=paste0("n=", sum(!is.na(x)))
        #     )
        #   },
        #   geom="text",
        #   vjust=-0.5,
        #   position=position_dodge2(width=0.75),
        #   show.legend=FALSE
        # ) +
        scale_color_manual(values=c(`Correct\nCluster Assignment`=unname(my_colors[9]), `Incorrect\nCluster Assignment`=unname(my_colors[1]))) +
        theme(
          panel.background=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_line(color="lightgrey"),
          legend.position="none",
          axis.title.x=element_blank(),
          # axis.text.x=element_blank(),
          # axis.ticks.length.x=unit(0, "pt")
        ) # |
      # ggplot(
      #   data=df_temp %>%
      #     filter(`Correct Cluster` == "Incorrect\nCluster Assignment") %>%
      #     group_by(Method) %>%
      #     summarize(n=n(), .groups="drop"),
      #   mapping=aes(x=Method, y=n, color=Method)
      # ) +
      #   geom_col(fill="white") +
      #   scale_color_manual(values=my_colors) +
      #   labs(y="# Incorrectly Assigned Cells") +
      #   theme(
      #     panel.background=element_blank(),
      #     panel.grid.major.x=element_blank(),
      #     panel.grid.major.y=element_line(color="lightgrey"),
      #     legend.position="none",
      #     axis.title.x=element_blank(),
      #   )
    ) +
      plot_annotation(tag_levels='a', tag_prefix="", tag_suffix=")") &
      theme(
        plot.tag.position=c(0.025, 0.985),  # Top-left corner
        plot.tag=element_text(hjust=1, vjust=1, face="plain"),   # Align text in corner
        legend.text=element_text(size=8.8),
        legend.title=element_text(size=11)
      )
    
    plot(p)
    
    df_temp %>%
      filter(`Correct Cluster` == "Incorrect\nCluster Assignment") %>%
      group_by(Method) %>%
      summarize(n=n(), .groups="drop")
    

    if(save_plots){
      ggsave(
        plot=p,
        filename=paste0("UMAP_stability_score.jpg"),
        path=plot_dir,
        width=21,
        height=9,
        units="cm",
        scale=1.4,
        dpi=350,
        device=dynamic_device("jpg")
      )
    }
    
    ##############################
    #### ROC & PrecRec Curves ####
    ##############################
    
    # investigate single scenario
    load(paste0("results/", use_simulation, "/results_simulation__celltypes=2__lib.size=25000__de.prob=0.01__de.facLoc=0.1__iter=38.Rdata"))
    
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
    
    (
    precrec::autoplot(curves, "ROC") +
      scale_color_manual(values=my_colors) +
      theme(
        plot.title=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(color="lightgrey"),
        panel.grid.minor=element_blank(),
        legend.position="top",
      ) |
    precrec::autoplot(curves, "PRC") +
      scale_color_manual(values=my_colors) +
      theme(
        plot.title=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(color="lightgrey"),
        panel.grid.minor=element_blank(),
        legend.position="top",
      )
    ) +
      plot_annotation(tag_levels='a', tag_prefix="", tag_suffix=")") +
      plot_layout(guides="collect", tag_level="keep") &
      theme(legend.position="top")
    
    precrec::auc(curves)

}

