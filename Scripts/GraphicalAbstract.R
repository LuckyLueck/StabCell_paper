library(uwot) # for umap
library(dplyr)
library(ggplot2)
library(rsvd) # for digits data

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

data('digits')
set.seed(42)
idx <- sample(1:nrow(digits), size=500, replace=FALSE)
data = digits[idx, grep("pixel", colnames(digits))]
labels = digits[idx, "label"]

# Run UMAP
umap_result <- umap(data, n_neighbors=15, min_dist=1.5, metric="euclidean")

# Create a data frame for plotting
umap_df <- data.frame(
  UMAP_1=umap_result[, 1],
  UMAP_2=umap_result[, 2],
  Label=as.factor(labels)
)

factor_replace <- function(f, to_replace=c(), new_value=to_replace[1]){ # merge clusters
  if(length(to_replace) == 0) return(f)
  f_char <- as.character(f)
  f_char[f_char %in% to_replace] <- new_value
  f_new <- factor(f_char, levels=if(new_value %in% to_replace) levels(f) else c(levels(f), new_value))
  return(f_new)
}

replace_random_fraction <- function(x, x_new, fraction = 0.3) { # make fuzzy clusters
  # Calculate how many values to replace
  n_replace <- ceiling(length(x) * fraction)
  
  # Generate random indices to replace
  indices <- sample(seq_along(x), n_replace)
  
  # Replace those values
  if(length(x_new) > 1){
    x[indices] <- sample(x_new, n_replace, TRUE)
  } else {
    x[indices] <- x_new
  }
  return(x)
}

permute_levels <- function(x) {
  if (!is.factor(x)) {
    warning("Input must be a factor.")
    return(x)
  }
  levels(x) <- sample(levels(x))
  return(x)
  
  # original_levels <- levels(x)
  # permuted_levels <- sample(original_levels)
  # 
  # # Create a named mapping
  # level_map <- setNames(permuted_levels, original_levels)
  # 
  # # Apply the mapping to the values
  # new_values <- level_map[as.character(x)]
  # 
  # factor(new_values, levels = original_levels)
}

palette <- my_colors[c(1, 6, 2, 8)]
names(palette) <- levels(umap_df$Label)

# Plot using ggplot2
plot_umap <- function(
    subsample=FALSE,
    colored=FALSE,
    complement=FALSE,
    permute_colors=FALSE,
    merge_labels=c(),
    fuzzy=0.,
    seed=42
){
  set.seed(seed)
  n <- nrow(umap_df)
  
  if(subsample){
    idx <- sample(1:n, size=floor(n / 2)) * (-1)^complement
  } else {
    idx <- 1:n
  }
  
  if(fuzzy){
    zero_stable <- umap_df %>% filter(Label == 0, UMAP_1 <= -6.5)
    zero_fuzzy <- umap_df %>% 
      filter(Label == 0, UMAP_1 > -6.5) %>%
      mutate(Label=replace_random_fraction(Label, c(2, 3), fuzzy))
   
    one_fuzzy2 <- umap_df %>% 
      filter(Label == 1, UMAP_2 <= 0) %>%
      mutate(Label=replace_random_fraction(Label, 2, fuzzy))
    one_fuzzy3 <- umap_df %>% 
      filter(Label == 1, UMAP_2 > 0) %>%
      mutate(Label=replace_random_fraction(Label, 3, fuzzy))
   
    two_stable <- umap_df %>% 
      filter(Label == 2, UMAP_1 > -2.5, UMAP_1 < 0, UMAP_2 < -2.5)
    two_fuzzy0 <- umap_df %>% 
      filter(Label == 2, UMAP_1 <= -2.5) %>%
      mutate(Label=replace_random_fraction(Label, 0, fuzzy))
    two_fuzzy1 <- umap_df %>% 
      filter(Label == 2, (UMAP_1 >= 0 & UMAP_2 < -2.5) | UMAP_1 > 2.5) %>%
      mutate(Label=replace_random_fraction(Label, 1, fuzzy))
    two_fuzzy3 <- umap_df %>% 
      filter(Label == 2, UMAP_1 > -2.5, UMAP_1 <= 2.5, UMAP_2 >= -2.5) %>%
      mutate(Label=replace_random_fraction(Label, 3, fuzzy))
    
    three_stable <- umap_df %>% 
      filter(Label == 3, UMAP_2 > 2.5)
    three_fuzzy2 <- umap_df %>% 
      filter(Label == 3, UMAP_2 <= 2.5) %>%
      mutate(Label=replace_random_fraction(Label, 2, fuzzy))
    
    plot_df <- bind_rows(list(
      zero_stable,
      zero_fuzzy,
      one_fuzzy2,
      one_fuzzy3,
      two_stable,
      two_fuzzy0,
      two_fuzzy1,
      two_fuzzy3,
      three_stable,
      three_fuzzy2
    ))
  } else {
    plot_df <- umap_df
  }

  ggplot(
    data=plot_df[idx, ] %>% 
      mutate(
        Label=factor_replace(Label, merge_labels),
        Label=if(permute_colors) permute_levels(Label) else Label
      ),
    mapping=aes(x=UMAP_1, y=UMAP_2)
  ) +
    geom_point(
      mapping=if(colored) aes(color=Label) else NULL,
      size=1,
      alpha=0.8
    ) +
    theme_minimal() +
    labs(x=element_blank(), y=element_blank()) +
    scale_color_manual(values=palette) +
    theme(
      legend.position="none",
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text=element_blank()
    )
}


save_plots <- TRUE # run on windows, linux has some svg render problems
plot_dir <- "plots/GraphicalAbstract"
width <- 3
height <- 1.5

# full data
p <- plot_umap(subsample=FALSE, colored=FALSE); plot(p); if(save_plots) ggsave(file="data_full.svg", path=plot_dir, plot=p, width=width, height=height)
p <- plot_umap(subsample=FALSE, colored=TRUE); plot(p); if(save_plots) ggsave(file="cluster_full.svg", path=plot_dir, plot=p, width=width, height=height)
p <- plot_umap(subsample=FALSE, colored=TRUE, fuzzy=0.2); plot(p); if(save_plots) ggsave(file="cluster_full_fuzzy.svg", path=plot_dir, plot=p, width=width, height=height)

# subsamples
p <- plot_umap(subsample=TRUE, seed=42, colored=FALSE, complement=FALSE);                           plot(p); if(save_plots) ggsave(file="data_sub1.svg", path=plot_dir, plot=p, width=width, height=height)
p <- plot_umap(subsample=TRUE, seed=42, colored=FALSE, complement=TRUE);                            plot(p); if(save_plots) ggsave(file="data_sub2.svg", path=plot_dir, plot=p, width=width, height=height)
p <- plot_umap(subsample=TRUE, seed=43, colored=FALSE, complement=FALSE);                           plot(p); if(save_plots) ggsave(file="data_sub3.svg", path=plot_dir, plot=p, width=width, height=height)
p <- plot_umap(subsample=TRUE, seed=42, colored=TRUE,  complement=FALSE);                           plot(p); if(save_plots) ggsave(file="cluster_sub1.svg", path=plot_dir, plot=p, width=width, height=height)
p <- plot_umap(subsample=TRUE, seed=42, colored=TRUE,  complement=TRUE,  merge_labels=c("3", "2")); plot(p); if(save_plots) ggsave(file="cluster_sub2.svg", path=plot_dir, plot=p, width=width, height=height)
p <- plot_umap(subsample=TRUE, seed=43, colored=TRUE,  complement=FALSE, merge_labels=c("1", "2")); plot(p); if(save_plots) ggsave(file="cluster_sub3.svg", path=plot_dir, plot=p, width=width, height=height)
p <- plot_umap(subsample=TRUE, seed=42, colored=TRUE,  complement=FALSE, permute_colors=TRUE);                           plot(p); if(save_plots) ggsave(file="cluster_sub1_perm.svg", path=plot_dir, plot=p, width=width, height=height)
p <- plot_umap(subsample=TRUE, seed=42, colored=TRUE,  complement=TRUE,  permute_colors=TRUE, merge_labels=c("3", "2")); plot(p); if(save_plots) ggsave(file="cluster_sub2_perm.svg", path=plot_dir, plot=p, width=width, height=height)
p <- plot_umap(subsample=TRUE, seed=43, colored=TRUE,  complement=FALSE, permute_colors=TRUE, merge_labels=c("1", "2")); plot(p); if(save_plots) ggsave(file="cluster_sub3_perm.svg", path=plot_dir, plot=p, width=width, height=height)





# gene selection histogram
set.seed(123)
n_bars <- 8
df_hist <- data.frame(
  var=as.factor(1:n_bars),
  freq=runif(n=n_bars)
)
selection_thresh <- 0.8
p <- ggplot(
  data=df_hist %>% mutate(selected=as.factor(freq > !!selection_thresh)),
  mapping=aes(x=var, y=freq, fill=selected)
) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=selection_thresh, color="black", linetype="dashed") +
  scale_fill_manual(values=my_colors[c(1, 6)]) +
  labs(
    x="Gene",
    y="Selection Frequency"
  ) +
  theme_minimal() +
  theme(
    legend.position="none",
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text=element_blank()
  )
p
ggsave(file="hist.svg", path=plot_dir, plot=p, width=2, height=2)


# stability boxplot
df_box <- data.frame(
  Label=factor(rep(1:4, each=5)),
  y=c(
    0.7, 0.8, 0.9, 1.0, 1.0,
    0.2, 0.5, 0.7, 0.9, 1.0,
    0.1, 0.3, 0.6, 0.8, 0.9,
    0.5, 0.6, 0.8, 1.0, 1.0
  )
)
p <- ggplot(
  data=df_box,
  mapping=aes(x=Label, y=y, color=Label)
) +
  geom_boxplot() +
  scale_color_manual(values=my_colors[c(1, 6, 2, 8)]) +
  theme_minimal() +
  theme(
    legend.position="none",
    plot.background=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text.x=element_blank()
  ) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) +
  coord_cartesian(ylim=c(0, 1)) +
  labs(x=element_blank(), y="Cell Stability Score")
p
ggsave(file="box.svg", path=plot_dir, plot=p, width=2, height=2)
