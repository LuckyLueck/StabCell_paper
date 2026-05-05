library(ggplot2)
library(dplyr)
library(tidyr)

results_dir <- "results/runtime_study"
plots_dir <- "plots/runtime_study/"

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

df_list <- list()
for(n_celltypes in c(5, 10)){
  setting_dir <- paste0(results_dir, "__n_celltypes=", n_celltypes, "/")
  files <- list.files(setting_dir, pattern="results_runtime__.*\\.Rdata$")

  n_cells <- c()
  iter <- c()
  t_clust_ddip <- c()
  t_marker_ddip <- c()
  t_total_ddip <- c()
  t_split_countsplit <- c()
  t_clust_countsplit <- c()
  t_marker_countsplit <- c()
  t_total_countsplit <- c()
  t_refclust_stabcell <- c()
  t_subsamples_stabcell <- c()
  t_total_stabcell <- c()

  for(file in files){
    load(paste0(setting_dir, file))

    n_cells <- c(n_cells, runtimes$n_cells)
    iter <- c(iter, runtimes$iter)
    t_clust_ddip <- c(t_clust_ddip, as.numeric(runtimes$ddip$duration_clust, unit="secs"))
    t_marker_ddip <- c(t_marker_ddip, as.numeric(runtimes$ddip$duration_marker, unit="secs"))
    t_total_ddip <- c(t_total_ddip, as.numeric(runtimes$ddip$duration_total, unit="secs"))
    t_split_countsplit <- c(t_split_countsplit, as.numeric(runtimes$countsplit$duration_split, unit="secs"))
    t_clust_countsplit <- c(t_clust_countsplit, as.numeric(runtimes$countsplit$duration_clust, unit="secs"))
    t_marker_countsplit <- c(t_marker_countsplit, as.numeric(runtimes$countsplit$duration_marker, unit="secs"))
    t_total_countsplit <- c(t_total_countsplit, as.numeric(runtimes$countsplit$duration_total, unit="secs"))
    t_refclust_stabcell <- c(t_refclust_stabcell, as.numeric(runtimes$stabcell$duration_refclust, unit="secs"))
    t_subsamples_stabcell <- c(t_subsamples_stabcell, as.numeric(runtimes$stabcell$duration_subsamples, unit="secs"))
    t_total_stabcell <- c(t_total_stabcell, as.numeric(runtimes$stabcell$duration_total, unit="secs"))
  }

  df_list[[length(df_list) + 1]] <- data.frame(
    n_celltypes=n_celltypes,
    n_cells=n_cells,
    iter=iter,
    t_clust_ddip=t_clust_ddip,
    t_marker_ddip=t_marker_ddip,
    t_total_ddip=t_total_ddip,
    t_split_countsplit=t_split_countsplit,
    t_clust_countsplit=t_clust_countsplit,
    t_marker_countsplit=t_marker_countsplit,
    t_total_countsplit=t_total_countsplit,
    t_refclust_stabcell=t_refclust_stabcell,
    t_subsamples_stabcell=t_subsamples_stabcell,
    t_total_stabcell=t_total_stabcell
  ) %>%
    pivot_longer(-c("n_celltypes", "n_cells", "iter"), values_to="Runtime [s]") %>%
    mutate(
      Method=factor(case_when(
        stringr::str_ends(name, "ddip") ~ "Double Dipping",
        stringr::str_ends(name, "countsplit") ~ "Countsplit",
        stringr::str_ends(name, "stabcell") ~ "StabCell",
        TRUE ~ NA
      ), levels=c("Double Dipping", "Countsplit", "StabCell")),
      Step=case_when(
        stringr::str_starts(name, "t_clust") ~ "Cluster",
        stringr::str_starts(name, "t_marker") ~ "Marker",
        stringr::str_starts(name, "t_total") ~ "Total",
        stringr::str_starts(name, "t_split") ~ "Split",
        stringr::str_starts(name, "t_refclust") ~ "Reference Cluster",
        stringr::str_starts(name, "t_subsamples") ~ "Subsamples",
        TRUE ~ NA
      )
    ) %>%
    rename(`Cell Types`=n_celltypes) %>%
    select(Method, Step, `Cell Types`, n_cells, iter, `Runtime [s]`)
}

df <- bind_rows(df_list) %>%
  mutate(
    `Runtime [min]`=`Runtime [s]` / 60,
    `Cell Types`=factor(
      paste0(`Cell Types`, " Cell Types"),
      levels=c("5 Cell Types", "10 Cell Types")
    )
  )

p <- ggplot(
  data=df %>%
    filter(
      Step == "Total",
      # `Cell Types` == "5 Cell Types",
      # n_cells < 50000
    ),
  mapping=aes(
    x=n_cells,
    y=`Runtime [min]`,
    color=Method,
    linetype=`Cell Types`
  )
) +
  stat_summary(fun=mean, geom="line") +
  stat_summary(fun=mean, geom="point") +
  stat_summary(
    fun.data=mean_se,
    geom="errorbar",
    mapping=aes(group=interaction(Method, `Cell Types`)),
    width=0.2,
    linetype="solid"
  ) +
  #facet_grid(cols=vars(`Cell Types`)) +
  theme(
    legend.position="top",
    legend.box="vertical",
    legend.spacing.y=unit(0, "mm"),
    panel.background=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(color="lightgrey"),
    panel.grid.minor=element_blank(),
    axis.ticks.y=element_blank()
  ) +
  guides(color=guide_legend(order=1), linetype=guide_legend(order=2)) +
  labs(x="Cells", color=element_blank(), linetype=element_blank()) +
  scale_color_manual(values=my_colors[c(1, 2, 4)])

plot(p)

ggsave(
  plot=p,
  filename="RuntimeStudy.eps",
  path=plots_dir,
  width=17, height=7, units="cm", dpi=350
)
ggsave(
  plot=p,
  filename="RuntimeStudy.tif",
  path=plots_dir,
  width=17, height=7, units="cm", dpi=350
)

print(
  df %>%
    filter(Step == "Total") %>%
    group_by(Method, `Cell Types`, n_cells) %>%
    summarize(t_mean=mean_se(`Runtime [min]`), .groups="drop")
  , n=50
)


# runtime over iters (should be without a trend)
ggplot(
  data=df %>% filter(Step == "Total"),
  mapping=aes(x=iter, y=`Runtime [min]`, color=`Cell Types`)
) +
  geom_line() +
  theme(legend.position="top") +
  labs(color=element_blank()) +
  facet_wrap( . ~ Method + n_cells, scales="free_y", nrow= 3)
  # facet_grid(rows=vars(n_cells), cols=vars(Method), scales="free_y")
