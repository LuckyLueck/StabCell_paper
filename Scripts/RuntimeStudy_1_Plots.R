library(ggplot2)
library(dplyr)
library(tidyr)

results_dir <- "results/runtime_study/"
plots_dir <- "plots/runtime_study/"

my_colors <- c(
  "Double Dipping"="#e33e50",
  "Countsplit"="#F2B705",
  "StabCell"="#2E7D32",
  "ClusterDE"="#8c7f9d",
  "TN Test"="#563d2c"
)


files <- list.files(results_dir, pattern="results_runtime__.*\\.Rdata$")
methods <- c("ddip", "countsplit", "clusterde", "tntest", "stabcell")
df_list <- list()
for(file in files){
  load(paste0(results_dir, file))
  for(method in methods){
    df_list[[length(df_list) + 1]] <- data.frame(
        runtimes[[method]]
      ) %>%
      mutate(method=!!method, n_cells=runtimes$n_cells, iter=runtimes$iter)
  }
}
df <- bind_rows(df_list) %>%
  tidyr::pivot_longer(
    cols=matches("^duration_[a-z]+$"),
    names_prefix="duration_",
    names_to="step",
    values_to="runtime"
  ) %>%
  filter(!is.na(runtime)) %>%
  mutate(
    `Runtime [min]`=as.numeric(runtime, unit="mins"),
    `Runtime [sec]`=as.numeric(runtime, unit="secs"),
    method=factor(case_when(
      method == "ddip" ~ "Double Dipping",
      method == "countsplit" ~ "Countsplit",
      method == "clusterde" ~ "ClusterDE",
      method == "tntest" ~ "TN Test",
      method == "stabcell" ~ "StabCell"
    ), levels=c("Double Dipping", "Countsplit", "ClusterDE", "TN Test", "StabCell"))
  )




p <- ggplot(
  data=df %>% filter(step == "total"),
  mapping=aes(
    x=n_cells,
    y=`Runtime [min]`,
    color=method
  )
) +
  stat_summary(fun=mean, geom="line") +
  stat_summary(fun=mean, geom="point") +
  stat_summary(
    fun.data=\(x) mean_se(x, mult=1.96),
    geom="errorbar",
    mapping=aes(group=method),
    width=0.2,
    linetype="solid"
  ) +
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
  scale_y_continuous(breaks=scales::breaks_pretty()) +
  scale_color_manual(values=my_colors)

plot(p)

p_zoom <- p / p + coord_cartesian(ylim=c(0,10)) +
  plot_layout(
    axes="collect_x",
    axis_titles="collect",
    guides="collect"
  ) &
  theme(legend.position="top")

plot(p_zoom)

ggsave(
  plot=p_zoom,
  filename="RuntimeStudy.eps",
  path=plots_dir,
  width=17, height=12, units="cm", dpi=350
)


print(
  df %>%
    filter(step == "total") %>%
    group_by(method, n_cells) %>%
    summarize(t_mean=mean_se(`Runtime [min]`), .groups="drop")
  , n=50
)


# runtime over iters (should be without a trend)
ggplot(
  data=df %>% filter(step == "total"),
  mapping=aes(x=iter, y=`Runtime [min]`)
) +
  geom_line() +
  theme(legend.position="top") +
  labs(color=element_blank()) +
  facet_wrap( . ~ method + n_cells, scales="free_y", nrow=length(methods))
  
