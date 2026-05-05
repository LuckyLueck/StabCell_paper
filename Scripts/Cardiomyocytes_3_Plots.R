library(ggplot2)
library(patchwork)
library(paletteer)
library(ggh4x)
library(dplyr)
library(Matrix)
library(xtable)

save_plots <- FALSE

data_dir <- "data/cardiomyocyte_differentiation/"
results_dir <- "results/cardiomyocyte_differentiation/"
plots_dir <- "plots/cardiomyocyte_differentiation/"

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

load(paste0(results_dir, "stabcell_results__subsampling_factor=0.Rdata"))
load(paste0(results_dir, "ddip_results.Rdata"))
load(paste0(results_dir, "countsplit_results.Rdata"))


# reorder factor levels of cluster identities to match between clusterings
relabel_ddip_clusters <- function(x){
  mapping <- c(
    # "new name" = "old name"
    "0"="2",
    "1"="0",
    "2"="5",
    "3"="3",
    "4"="4",
    "5"="1",
    "6"="6",
    "7"="9",
    "8"="8",
    "9"="7"
  )
  factor(
    forcats::fct_recode(x, !!!mapping),
    levels=names(mapping)
  )
}
relabel_countsplit_clusters <- function(x){
  mapping <- c(
    # "new name" = "old name"
    "0"="2",
    "1"="0",
    "2"="5",
    "3"="3",
    "4"="4",
    "5"="1",
    "6"="6",
    "7"="9",
    "8"="8",
    "9"="7"
  )
  factor(
    forcats::fct_recode(x, !!!mapping),
    levels=names(mapping)
  )
}
cluster_celltype_assignment <- function(x){
  mapping <- c(
    # "new name" = "old name"
    "IPSC"="0",
    "MES"="1",
    "CMES"="2",
    "PROG / TRANS"="3",
    "CM"="4",
    "CF"="5",
    "UNK1"="6",
    "UNK2"="7",
    "SPUR1"="8",
    "SPUR2"="9"
  )
  factor(
    forcats::fct_recode(x, !!!mapping),
    levels=names(mapping)
  )
}

{
# double dipping preprocessing
ddip_results$idents <- ddip_results$idents %>%
  mutate(
    article_celltype=factor(article_celltype, levels=c("IPSC", "MES", "CMES", "PROG", "CM", "CF", "UNK")),
    low_res_clusters=relabel_ddip_clusters(low_res_clusters)
  )
ddip_results$markers <- ddip_results$markers %>%
  mutate(
    cluster1=relabel_ddip_clusters(cluster1),
    cluster2=relabel_ddip_clusters(cluster2)
  )

# countplit preprocessing
countsplit_results$idents <- relabel_countsplit_clusters(countsplit_results$idents)
countsplit_results$markers <- countsplit_results$markers %>%
  mutate(
    cluster1=relabel_countsplit_clusters(cluster1),
    cluster2=relabel_countsplit_clusters(cluster2)
  )

# stabcell preprocessing
stabcell_results$clustering_alignment <- stabcell_results$clustering_alignment %>% mutate(
  reference_cluster_id=relabel_ddip_clusters(reference_cluster_id)
)
colnames(stabcell_results$assignment_probabilities) <- cluster_celltype_assignment(colnames(stabcell_results$assignment_probabilities))
stabcell_results$marker_alignment <- stabcell_results$marker_alignment %>% mutate(
  reference_cluster_id=relabel_ddip_clusters(reference_cluster_id)
)
stabcell_results$subsampling_p_values <- stabcell_results$subsampling_p_values %>% mutate(
  id1=relabel_ddip_clusters(id1),
  id2=relabel_ddip_clusters(id2)
)
stabcell_results$selection_frequencies <- stabcell_results$selection_frequencies %>% mutate(
  id1=relabel_ddip_clusters(id1),
  id2=relabel_ddip_clusters(id2)
)
stabcell_results$stable_markers_global <- stabcell_results$stable_markers_global %>% mutate(
  id=relabel_ddip_clusters(id)
)
stabcell_results$stable_markers_local <- stabcell_results$stable_markers_local %>% mutate(
  id1=relabel_ddip_clusters(id1),
  id2=relabel_ddip_clusters(id2)
)
stabcell_results$clustering_mode <- relabel_ddip_clusters(stabcell_results$clustering_mode)
}


# find setting where selection_freq_thresh not greater than 1 for all cell type
# --> r_concavity & PFER_max = 100
print(
  stabcell_results$stable_markers_global %>%
    group_by(PFER_max, id, assumption, exceeds_one) %>%
    summarize(n=n(), .groups="drop") %>%
    filter(!exceeds_one),
  n=500
)


# UMAP
{
  my_theme <- theme(
    axis.line=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank()
  )

  legend_width <- 2
  coord_fixed_common <- coord_fixed(
    xlim = range(ddip_results$umap[, "umap_1"], na.rm = TRUE),
    ylim = range(ddip_results$umap[, "umap_2"], na.rm = TRUE),
    expand = FALSE
  )

  # color = differentiation day
  p_umap_diffday <- ggplot(data=ddip_results$umap, mapping=aes(x=umap_1, y=umap_2, color=ddip_results$idents$diffday)) +
    geom_point(size=0.1) +
    scale_color_paletteer_d("rcartocolor::TealRose") +
    labs(color="Differentiation Day") +
    guides(color = guide_legend(override.aes=list(size=2), ncol=2)) +  # Adjust the size of points in the legend
    theme(legend.key.width=unit(legend_width, "cm")) +
    my_theme

  # double dipping clusters
  p_umap_ddip <-
    ggplot(
      data=data.frame(
        x=ddip_results$umap[,"umap_1"],
        y=ddip_results$umap[,"umap_2"],
        z=ddip_results$idents$low_res_clusters
      ) %>%
        mutate(z=cluster_celltype_assignment(z)),
      mapping=aes(x=x, y=y, color=z)
    ) +
    geom_point(size=0.1) +
    scale_color_manual(values=my_colors[c(1, 5, 10, 2, 4, 8, 7, 6, 9, 3)]) +
    labs(color="Double Dipping", title="Double Dipping") +
    theme(legend.position="bottom") +
    my_theme +
    guides(color=guide_legend(override.aes=list(size=2), nrow=1)) # Adjust the size of points in the legend +

  # elorbany clustering
  p_umap_elorbany <- ggplot(
    data=data.frame(
      x=ddip_results$umap[,"umap_1"],
      y=ddip_results$umap[,"umap_2"],
      z=ddip_results$idents$article_celltype
    ),
    mapping=aes(x=x, y=y, color=z)
  ) +
    geom_point(size=0.1) +
    scale_color_manual(values=my_colors[c(1, 5, 10, 2, 4, 8, 7, 6, 9, 3)]) +
    labs(color="Reference", title="Reference") +
    guides(color = guide_legend(override.aes=list(size=2), ncol=2)) +  # Adjust the size of points in the legend +
    theme(legend.position="none") +
    my_theme

  # stable clustering
  p_umap_stable <-
    ggplot(
      data=data.frame(
        x=ddip_results$umap[,"umap_1"],
        y=ddip_results$umap[,"umap_2"],
        z=stabcell_results$clustering_mode
      ) %>%
        mutate(z=cluster_celltype_assignment(z)),
      mapping=aes(x=x, y=y, color=z)
    ) +
    geom_point(size=0.1) +
    scale_color_manual(values=my_colors[c(1, 5, 10, 2, 4, 8, 7, 6, 9, 3)]) +
    labs(color="StabCell", title="StabCell") +
    guides(color=guide_legend(override.aes=list(size=2), ncol=2)) +  # Adjust the size of points in the legend +
    theme(legend.position="none") +
    my_theme

  # countsplit clustering
  p_umap_countsplit <-
    ggplot(
      data=data.frame(
        x=ddip_results$umap[,"umap_1"],
        y=ddip_results$umap[,"umap_2"],
        z=countsplit_results$idents
      ) %>%
        mutate(z=cluster_celltype_assignment(z)),
      mapping=aes(x=x, y=y, color=z)
    ) +
    geom_point(size=0.1) +
    scale_color_manual(values=my_colors[c(1, 5, 10, 2, 4, 8, 7, 6, 9, 3)]) +
    labs(color="Countsplit", title="Countsplit") +
    guides(color=guide_legend(override.aes=list(size=2), ncol=2)) +  # Adjust the size of points in the legend +
    theme(legend.position="none") +
    my_theme
  p_umap_countsplit

  # color = stability score
  p_umap_stabscore <-
    ggplot(
      data=data.frame(
        x=ddip_results$umap[,"umap_1"],
        y=ddip_results$umap[,"umap_2"],
        z=stabcell_results$cell_stability_score
      ),
      mapping=aes(x=x, y=y, color=z)
    ) +
    geom_point(size=0.1) +
    scale_colour_gradient2(low=my_colors[4], mid=my_colors[2], high=my_colors[1], limits=c(0, 1), midpoint=0.5) +
    # scale_colour_paletteer_c("ggthemes::Blue-Green Sequential") +
    labs(color="Stability Score") +
    guides(color=guide_legend(override.aes=list(size=2), ncol=2)) +  # Adjust the size of points in the legend
    my_theme


  p_umap_clusterings <- (
    (p_umap_elorbany | p_umap_ddip) /
    (p_umap_countsplit | p_umap_stable)
  ) / guide_area() +
    plot_layout(heights=c(10, 10, 1), guides='collect') +
    plot_annotation(tag_levels='a', tag_prefix="", tag_suffix=")") &
      labs(color=element_blank()) &
      theme(
        plot.tag.position=c(0.0, 1.0),  # Top-left corner
        plot.tag=element_text(size=8, hjust=1, vjust=1),   # Align text in corner
        plot.title=element_text(size=8, hjust=0.5),
        legend.text=element_text(size=8, margin=margin(l=0., r=0.3)),
        legend.box.just = "center",
        legend.justification = "center",
        legend.margin=margin(l=0, r=20),
      )
  p_umap_clusterings
}
if(save_plots){
  ggsave(
    plot=p_umap_clusterings,
    file="cardiomyocyte_umap.jpg", path=plots_dir,
    width=17, height=12, units="cm", dpi=350
  )
  ggsave(
    plot=p_umap_clusterings,
    file="cardiomyocyte_umap.tif", path=plots_dir,
    width=17, height=12, units="cm", dpi=350
  )


  ggsave(
    plot=p_umap_diffday + theme(
      legend.text=element_text(size=8, margin=margin(l=0., r=0.3)),
      legend.title=element_text(size=10, margin=margin(l=0., r=0.3)),
      legend.margin=margin(l=0, r=20),
    ),
    file="cardiomyocyte_diff_day.jpg", path=plots_dir,
    width=15, height=5.5, units="cm", dpi=350
  )
  ggsave(
    plot=p_umap_diffday + theme(
      legend.text=element_text(size=8, margin=margin(l=0., r=0.3)),
      legend.title=element_text(size=10, margin=margin(l=0., r=0.3)),
      legend.margin=margin(l=0, r=20),
    ),
    file="cardiomyocyte_diff_day.tif", path=plots_dir,
    width=15, height=5.5, units="cm", dpi=350
  )
}

# UMAP clusters of subsample
plot_list <- list()
for(subsample_id in 1:nrow(stabcell_results$subsampling_identities)){

  cell_idx <- is.finite(stabcell_results$subsampling_identities[subsample_id, ])

  p <- ggplot(
    data=data.frame(
      x=ddip_results$umap[cell_idx, "umap_1"],
      y=ddip_results$umap[cell_idx, "umap_2"],
      z=stabcell_results$subsampling_identities[subsample_id, cell_idx]
    ) %>%
      mutate(z=factor(z)),
    mapping=aes(x=x, y=y, color=z)
  ) +
  geom_point(size=0.1) +
  scale_color_manual(values=my_colors[c(1, 5, 10, 2, 4, 8, 7, 6, 9, 3)]) +
  labs(color=element_blank(), title=paste("Subsample", subsample_id)) +
  theme(legend.position="none") +
  my_theme +
  guides(color=guide_legend(override.aes=list(size=5), nrow=1))
  plot_list[[length(plot_list) + 1]] <- p
}
p_umap_subsamples <- wrap_plots(plot_list, ncol=10)
p_umap_subsamples
if(save_plots){
  ggsave(
    plot=p_umap_subsamples & theme(plot.title=element_text(size=8)),
    filename="subsample_clusterings.jpg",
    path=plots_dir,
    width=22, height=22,
    units="cm", dpi=350,
    limitsize=FALSE
  )
  ggsave(
    plot=p_umap_subsamples & theme(plot.title=element_text(size=8)),
    filename="subsample_clusterings.tif",
    path=plots_dir,
    width=22, height=22,
    units="cm", dpi=350,
    limitsize=FALSE
  )
}

# Composition by Diffday
{
  p_comp_elorbany <- ggplot(
    ddip_results$idents,
    aes(x=diffday, fill=article_celltype)
  ) +
    geom_bar(position="fill") +
    scale_y_continuous(labels=scales::percent_format()) +
    guides(fill=guide_legend(ncol=2)) +
    labs(y=element_blank(), x=element_blank(), fill="Elorbany", title="Tissue Composition") +
    scale_fill_manual(values=my_colors[c(1, 5, 10, 2, 4, 8, 7, 6, 9, 3)]) +
    geom_text(
      data=ddip_results$idents %>% count(diffday),
      aes(x=diffday, y=1.05, label=paste0("n=",n)),
      inherit.aes = FALSE
    )

  p_comp_ddip <- ggplot(
    ddip_results$idents %>%
      mutate(celltype=cluster_celltype_assignment(low_res_clusters)),
    aes(x=diffday, fill=celltype)
  ) +
    geom_bar(position="fill") +
    scale_y_continuous(labels=scales::percent_format()) +
    guides(fill=guide_legend(ncol=2)) +
    labs(y=element_blank(), x=element_blank(), fill="DoubleDip") +
    scale_fill_manual(values=my_colors[c(1, 5, 10, 2, 4, 8, 7, 6, 9, 3)])

  p_comp_countsplit <- ggplot(
    mapping=aes(
      x=ddip_results$idents$diffday,
      fill=cluster_celltype_assignment(countsplit_results$idents)
    )
  ) +
    geom_bar(position="fill") +
    scale_y_continuous(labels=scales::percent_format()) +
    guides(fill=guide_legend(ncol=2)) +
    labs(y=element_blank(), x=element_blank(), fill="Countsplit") +
    scale_fill_manual(values=my_colors[c(1, 5, 10, 2, 4, 8, 7, 6, 9, 3)])

  p_comp_stabcell <- ggplot(
    mapping=aes(
      x=ddip_results$idents$diffday[rnd_idx],
      fill=cluster_celltype_assignment(stabcell_results$clustering_mode)
    )
  ) +
    geom_bar(position="fill") +
    scale_y_continuous(labels=scales::percent_format()) +
    guides(fill=guide_legend(ncol=2)) +
    labs(y=element_blank(), x="Diffday", fill="StabCell") +
    scale_fill_manual(values=my_colors[c(1, 5, 10, 2, 4, 8, 7, 6, 9, 3)])

  p_composition <- p_comp_elorbany / p_comp_ddip / p_comp_countsplit / p_comp_stabcell
  p_composition
}


# CoClustering Frequencies of Article vs DDip/Countsplit/StabCell
{
  # article vs double dipping
  tab1 <- table(ddip_results$idents[["article_celltype"]], ddip_results$idents[["low_res_clusters"]])
  p1 <- ggplot(
    tab1 %>%
      as.data.frame() %>%
      rename(Article=Var1, DoubleDipping=Var2),
    aes(x=Article, y=DoubleDipping)
  ) +
    geom_tile(aes(fill=Freq)) +
    geom_text(aes(label=Freq)) +
    scale_fill_gradient(low="white", high="red") +
    labs(title="Co-clustering: Article vs DoubleDipping") +
    theme(legend.position="none")

  # article vs countsplit
  tab2 <- table(ddip_results$idents[["article_celltype"]], countsplit_results$idents)
  p2 <- ggplot(
    tab2 %>%
      as.data.frame() %>%
      rename(Article=Var1, Countsplit=Var2),
    aes(x=Article, y=Countsplit)
  ) +
    geom_tile(aes(fill=Freq)) +
    geom_text(aes(label=Freq)) +
    scale_fill_gradient(low="white", high="red") +
    labs(title="Co-clustering: Article vs Countsplit") +
    theme(legend.position="none")

  # article vs StabCell
  tab3 <- table(ddip_results$idents[["article_celltype"]], stabcell_results$clustering_mode)
  p3 <- ggplot(
    tab3 %>%
      as.data.frame() %>%
      rename(Article=Var1, StabCell=Var2),
    aes(x=Article, y=StabCell)
  ) +
    geom_tile(aes(fill=Freq)) +
    geom_text(aes(label=Freq)) +
    scale_fill_gradient(low="white", high="red") +
    labs(title="Co-clustering: Article vs StabCell") +
    theme(legend.position="none")

  p1 | p2 | p3
}

# Number of Clusters in Subsamples
hist(
  table(stabcell_results$clustering_alignment$subsample),
  main="Number of Clusters", xlab=""
)

# Cell Stability Scores
{
p_stabscore <- ggplot(
  data=data.frame(
    cluster=stabcell_results$clustering_mode,
    cell_stability_score=stabcell_results$cell_stability_score
  ) %>%
    mutate(
      cluster=recode(
        cluster,
        "0" = "IPSC",
        "1" = "MES",
        "2" = "CMES",
        "3" = "PROG / TRANS",
        "4" = "CM",
        "5" = "CF",
        "6" = "UNK1",
        "7" = "UNK2",
        "8" = "SPUR1",
        "9" = "SPUR2",
      )
    ),
  mapping=aes(fill=cluster, y=cell_stability_score, x=cluster)
) +
  geom_violin() +
  scale_fill_manual(values=my_colors[c(1, 5, 10, 2, 4, 8, 7, 6, 9, 3)]) +
  labs(y="Cell Stability Score", x="StabCell Cluster") +
  coord_cartesian(ylim=c(0, 1)) +
  theme(
    legend.position="none",
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(color="lightgrey"),
    panel.grid.minor=element_blank(),
    plot.background=element_blank()
  )
  plot(p_stabscore)
}
if(save_plots){
  ggsave(
    plot=p_stabscore,
    file="cell_stability_scores.eps",
    path=plots_dir,
    width=17, height=12, units="cm", dpi=350
  )
  ggsave(
    plot=p_stabscore,
    file="cell_stability_scores.tif",
    path=plots_dir,
    width=17, height=12, units="cm", dpi=350
  )
}

# adjusted rand index (Elorbany vs others)
{
temp_list <- list(
  "elorbany"=ddip_results$idents$article_celltype,
  "ddip"=ddip_results$idents$low_res_clusters,
  "countsplit"=countsplit_results$idents,
  "stabcell"=stabcell_results$clustering_mode
)

ari <- matrix(
  nrow=length(temp_list), ncol=length(temp_list),
  dimnames=list(names(temp_list), names(temp_list))
)
for(i in 1:length(temp_list)){
  for(j in i:length(temp_list)){
    ari[i, j] <- ari[j, i] <- mclust::adjustedRandIndex(
      temp_list[[i]],
      temp_list[[j]]
    )
  }
}

ari
}



######################
#### marker plots ####
######################

# number of markers
only_pos_lfc <- FALSE
{
# calculate n_markers for each approach
df_n_markers_subsample <- stabcell_results$subsampling_p_values %>%
  rename(id=id1) %>%
  filter(
    is.na(id2),
    p_val <= 0.05,
    if(only_pos_lfc) lfc > 0 else TRUE
  ) %>%
  group_by(id, subsample) %>%
  summarize(n_markers=n_distinct(marker), .groups="drop")

df_n_markers_ddip <- ddip_results$markers %>%
  rename(id=cluster1) %>%
  filter(
    is.na(cluster2),
    p_val_adj <= 0.05,
    if(only_pos_lfc) avg_log2FC > 0 else TRUE
  ) %>%
  group_by(id) %>%
  summarize(n_markers_ddip=n_distinct(marker), .groups="drop")

df_n_markers_countsplit <- countsplit_results$markers %>%
  rename(id=cluster1) %>%
  filter(
    is.na(cluster2),
    p_val_adj <= 0.05,
    if(only_pos_lfc) avg_log2FC > 0 else TRUE
  ) %>%
  group_by(id) %>%
  summarize(n_markers_countsplit=n_distinct(marker), .groups="drop")

df_n_markers_stabcell <- stabcell_results$stable_markers_global %>%
  group_by(id, p_val_thresh, PFER_max, assumption, exceeds_one) %>%
  filter(
    !is.na(marker),
    if(only_pos_lfc) lfc_mean > 0 else TRUE
  ) %>%
  summarize(n_markers_stabcell=n_distinct(marker), .groups="drop") %>%
  ungroup() %>%
  arrange(p_val_thresh, assumption, PFER_max, id)

# join dataframes
df_n_markers_comparison <- df_n_markers_stabcell %>%
  full_join(df_n_markers_ddip, by="id") %>%
  full_join(df_n_markers_countsplit, by="id") %>%
  full_join(
    df_n_markers_subsample %>%
      group_by(id) %>%
      summarize(n_markers_subsample_mean=mean(n_markers), .groups="drop"),
    by="id"
  ) %>% mutate(
    id=cluster_celltype_assignment(id)
  )

p <- ggplot() +

  # subsample (violin)
  geom_violin(
    data=df_n_markers_subsample,
    mapping=aes(x=id, y=n_markers, fill=id)
  ) +

  # ddip comparison (circle)
  geom_point(
    data=df_n_markers_ddip,
    mapping=aes(x=id, y=n_markers_ddip, fill=id),
    shape=16, size=2.5
  ) +

  # countsplit comparison (triangle)
  geom_point(
    data=df_n_markers_countsplit,
    mapping=aes(x=id, y=n_markers_countsplit, fill=id),
    shape=17, size=2.5
  ) +

  # stabcell comparison (diamond)
  geom_point(
    data=df_n_markers_stabcell %>%
      filter(p_val_thresh == 0.05, PFER_max == 100, assumption == "r_concavity"),
    mapping=aes(x=id, y=n_markers_stabcell, fill=id),
    shape=18, size=3
  ) +

  # styling
  scale_fill_manual(values=my_colors[c(1, 5, 10, 2, 4, 8, 7, 6, 9, 3)]) +
  theme(
    legend.position="none",
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(color="lightgrey"),
    panel.grid.minor=element_blank(),
    plot.background=element_blank()
  ) +
  labs(x="Cluster", y="Number of Selected Markers")

plot(p)

temp_tab <- df_n_markers_comparison %>%
  filter(
    is.na(p_val_thresh) | p_val_thresh == 0.05,
    is.na(PFER_max)     | PFER_max == 100,
    is.na(assumption)   | assumption == "r_concavity"
  ) %>%
  mutate(assumption=factor(assumption, levels=c("none", "unimodal", "r_concavity"))) %>%
  arrange(id, assumption) # %>%
  # tidyr::pivot_wider(
  #   id_cols=c("id", "n_markers_ddip", "n_markers_countsplit"),
  #   names_from="assumption",
  #   values_from=c("n_markers_stabcell", "n_markers_subsample_mean")
  # ) %>%
  # select(-`NA`)
  # select(id, assumption, n_markers_stabcell, n_markers_subsample_mean, n_markers_ddip, n_markers_countsplit)

print(temp_tab)
}
print(
  temp_tab %>%
    mutate(
      n_markers_subsample_mean=as.integer(round(n_markers_subsample_mean)),
      n_markers_subsample_mean=sapply(n_markers_subsample_mean, FUN=\(x){if(is.na(x)) "" else paste0("(", x, ")")})
    ) %>%
    rename(
      Cluster=id,
      `StabCell`=n_markers_stabcell,
      `(Subsample Mean)`=n_markers_subsample_mean,
      `Double Dipping`=n_markers_ddip,
      `Countsplit`=n_markers_countsplit
    ) %>%
    select(Cluster, StabCell, `(Subsample Mean)`, `Double Dipping`, Countsplit) %>%
    xtable,
  include.rownames=FALSE
)

# selection frequencies
{
  df_sel_freq <- stabcell_results$selection_frequencies %>%
    filter(is.na(id2), p_val_thresh == 0.05) %>%
    mutate(id1=cluster_celltype_assignment(id1))

  p_sel_freq <- ggplot(
    data=df_sel_freq,
    mapping=aes(x=id1, y=selection_freq, fill=id1, color=id1)
  ) +
    geom_point(size=0.1, alpha=0.1, position=position_jitter()) +
    scale_fill_manual(values=my_colors[c(1, 5, 10, 2, 4, 8, 7, 6, 9, 3)]) +
    scale_color_manual(values=my_colors[c(1, 5, 10, 2, 4, 8, 7, 6, 9, 3)]) +
    theme(
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major.x=element_blank(),
      panel.grid.major.y=element_line(color="lightgrey"),
      panel.grid.minor=element_blank(),
      plot.background=element_blank()
    ) +
    labs(x="Cluster", y="Selection Frequency")
  plot(p_sel_freq)

  lim <- c(0.1, 0.9)
  df_sel_freq %>%
    group_by(id1) %>%
    summarize(
      n_low=sum( (selection_freq >= 0)        & (selection_freq <  !!lim[1])) / n() * 100,
      n_mid=sum( (selection_freq >= !!lim[1]) & (selection_freq <= !!lim[2])) / n() * 100,
      n_high=sum((selection_freq >  !!lim[2]) & (selection_freq <= 1)       ) / n() * 100,
      .groups="drop"
    ) %>%
    arrange(as.character(id1)) %>%
    rename_with(\(x) sapply(x, \(x) switch(x,
                                           id1="Cluster",
                                           n_low=paste0("[0%, ", lim[1]*100, "%)"),
                                           n_mid=paste0("[", lim[1]*100,"%, ", lim[2]*100, "%]"),
                                           n_high=paste0("(", lim[2]*100,"%, 100%]"),
                                           x
                                           )))

}
if(save_plots){
  ggsave(
    plot=p_sel_freq,
    file="selection_frequencies.jpg", path=plots_dir,
    width=17, height=8, units="cm", dpi=350
  )
  ggsave(
    plot=p_sel_freq,
    file="selection_frequencies.tif", path=plots_dir,
    width=17, height=8, units="cm", dpi=350
  )
}


### cell type identification ###
marker_list <- list(
  "IPSC"=c("L1TD1"),
  "MES"=c("MIXL1"),
  "CMES"=c("CRABP2", "HAS2", "SERPINE2"),
  "PROG"=c("ISL1", "NKX2-5", "GATA4", "TBX5"),
  "CF"=c("COL3A1"),
  "CM"=c("MYL7")
)

df_markers <- data.frame(
  cell_type=factor(
    rep(names(marker_list), lengths(marker_list)),
    levels=names(marker_list)
  ),
  marker=purrr::flatten_chr(marker_list)
)





# calculate marker ranking
{
stabcell_marker_ranking <- stabcell_results$stable_markers_global %>%
  rename(lfc=lfc_mean) %>%
  filter(
    PFER_max == 100,
    assumption == "r_concavity",
    lfc > 0
  ) %>%
  group_by(id) %>%
  arrange(desc(selection_freq), desc(lfc)) %>%
  mutate(rank=row_number()) %>%
  ungroup() %>%
  select(id, marker, rank, selection_freq, lfc)

ddip_marker_ranking <- ddip_results$markers %>%
  rename(id=cluster1, lfc=avg_log2FC) %>%
  filter(
    is.na(cluster2),
    p_val_adj <= 0.05,
    lfc > 0
  ) %>%
  group_by(id) %>%
  arrange(p_val_adj, desc((lfc))) %>%
  mutate(rank=row_number()) %>%
  ungroup() %>%
  select(id, marker, rank, p_val_adj, lfc)

countsplit_marker_ranking <- countsplit_results$markers %>%
  rename(id=cluster1, lfc=avg_log2FC) %>%
  filter(
    is.na(cluster2),
    p_val_adj <= 0.05,
    lfc > 0
  ) %>%
  group_by(id) %>%
  arrange(p_val_adj, desc((lfc))) %>%
  mutate(rank=row_number()) %>%
  ungroup() %>%
  select(id, marker, rank, p_val_adj, lfc)
}

# top n markers per cluster and approach
{
top_n <- 20
df_top_markers <- stabcell_marker_ranking %>%
  filter(rank <= !!top_n) %>%
  full_join(
    ddip_marker_ranking %>% filter(rank <= !!top_n),
    by=c("id", "rank"),
    suffix=c("_stabcell", "_ddip")
  ) %>%
  full_join(
    countsplit_marker_ranking %>%
      filter(rank <= !!top_n) %>%
      rename(marker_countsplit=marker),
    by=c("id", "rank")
  ) %>%
  select(id, rank, marker_stabcell, marker_ddip, marker_countsplit) %>%
  arrange(id, rank)
print(df_top_markers)
if(save_plots){
  print(df_top_markers %>%
          mutate(id=cluster_celltype_assignment(id)) %>%
          xtable,
        include.rownames=FALSE)
}
}

# DotPlot
MarkerDotPlot <- function(stable_markers_global, markers_biol){
  ggplot(
    data=stable_markers_global %>%
      filter(
        marker %in% markers_biol$marker
      ) %>%
      left_join(markers_biol, by="marker") %>%
      mutate(col_header="Target Cell Type"),
    mapping=aes(x=marker, y=id, color=lfc_mean, size=abs(lfc_mean))
  ) +
    geom_point() +
    scale_color_gradient2() +
    facet_nested(~ col_header + cell_type, scales="free_x") +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
    labs(x="Gene", y="Numeric Cluster Label", size="abs. LFC", color="LFC")
}

{
p <- MarkerDotPlot(
  stable_markers_global=stabcell_results$stable_markers_global %>%
                          filter(PFER_max == 100, assumption == "r_concavity", p_val_thresh == 0.05),
  markers_biol=df_markers
)
plot(p)
}
if(save_plots){
  ggsave(
    plot=p,
    file="biological_markers.eps", path=plots_dir,
    width=17, height=10, units="cm", dpi=350
  )
  ggsave(
    plot=p,
    file="biological_markers.tif", path=plots_dir,
    width=17, height=10, units="cm", dpi=350
  )
}

celltypes_long <-c(
  "IPSC"="IPSC",
  "MES"="mesoderm",
  "CMES"="cardiac mesoderm",
  "PROG"="progenitor",
  "CM"="cardiomyocyte",
  "CF"="cardiac fibroblast",
  "UNK1"="unknown",
  "UNK2"="unknown"
)
df_celltype_assignment <- data.frame(
  id=0:7,
  cell_type=factor(celltypes_long, levels=unique(celltypes_long)),
  abbrev=factor(names(celltypes_long), levels=unique(names(celltypes_long)))
)

# compare stabcell & ddip ranking of biological markers
{

df_assignment_markers <- df_celltype_assignment %>%
    full_join(
      df_markers,
      by=c("abbrev"="cell_type")
    ) %>%
    mutate(id=as.character(id)) %>%
    left_join(
      stabcell_marker_ranking %>%
        rename(rank_stabcell=rank, selection_freq_stabcell=selection_freq, lfc_stabcell=lfc),
      by=c("id", "marker")
    ) %>%
    left_join(
      ddip_marker_ranking %>%
        rename(rank_ddip=rank, p_val_adj_ddip=p_val_adj, lfc_ddip=lfc),
      by=c("id", "marker")
    ) %>%
    left_join(
      countsplit_marker_ranking %>%
        rename(rank_countsplit=rank, p_val_adj_countsplit=p_val_adj, lfc_countsplit=lfc),
      by=c("id", "marker")
    )
      #full_join(, ddip_marker_ranking, countsplit_marker_ranking, by=c("id", "marker"), suffix=c("_stabcell", "_ddip", "_countsplit")),

df_assignment_markers
}

# LaTeX table
print(
  df_assignment_markers %>%
    mutate(cell_type=paste0(cell_type, " (", abbrev, ")")) %>%
    select(cell_type, marker, rank_stabcell, rank_ddip, rank_countsplit) %>%
    rename(
      `Cell Type`=cell_type,
      Marker=marker,
      `StabCell rank`=rank_stabcell,
      `Double Dipping rank`=rank_ddip,
      `Countsplit rank`=rank_countsplit,
    ) %>%
    xtable,
  include.rownames=FALSE
)


use_celltype <- "MES"
ggplot(
  data=stabcell_marker_ranking %>%
    rename(rank_stabcell=rank) %>%
    mutate(cell_type=cluster_celltype_assignment(id)) %>%
    filter(cell_type == !!use_celltype) %>%
    select(marker, rank_stabcell) %>%
    inner_join(
      ddip_marker_ranking %>%
        rename(rank_ddip=rank) %>%
        mutate(cell_type=cluster_celltype_assignment(id)) %>%
        filter(cell_type == !!use_celltype) %>%
        select(marker, rank_ddip),
      by="marker"
    ),
  aes(x=rank_ddip, y=rank_stabcell)
) +
  geom_point(size=0.3) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="red", linewidth=1) +
  labs(x="Double Dipping Rank", y="StabCell Rank")

