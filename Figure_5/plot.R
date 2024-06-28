require(tidyverse)

make_plot <- function(family_name, res.dist, xmax) {
  require(dendextend)
  require(ggdendro)
  hc <- hclust(res.dist, method = "average")
  dend0 <- stats::as.dendrogram(hc)
  ddata <- dendro_data(hc, type = "rectangle")

  clust_ref <- dendextend::cutree(dend0, h = xmax)
  xmin <- 0
  all_rect <- c() # xmin, xmax, ymin, ymax

  for (a in unique(clust_ref)) {
    b <- clust_ref[clust_ref == a]
    in_plot <- subset(ddata$labels, ddata$labels$label %in% names(b))
    ymin <- min(in_plot$x) - 0.5
    ymax <- max(in_plot$x) + 0.5
    to_add <- c(xmin, xmax, ymin, ymax)
    all_rect <- rbind(all_rect, to_add)
  }

  rect_plot <- data.frame(xmin = all_rect[, 1],
                          xmax = all_rect[, 2],
                          ymin = all_rect[, 3],
                          ymax = all_rect[, 4],
                          categ = 1:length(all_rect[, 1]))

  require(MetBrewer)
  require(ggnewscale)

  focal_palette <- "Derain"

  require(ggplot2)

  tip_info <- read_tsv("../datasets/families.txt")
  to_remove <- which(tip_info$Statistic == "number_of_lineages")
  to_remove <- c(to_remove, which(tip_info$Statistic == "crown_age"))
  to_remove <- c(to_remove, which(tip_info$Statistic == "tree_height"))


  df_tip_info <- data.frame(tip_info)

  for_plot2 <- as_tibble(ddata$labels)

  fp1 <- cbind(df_tip_info$Statistic, df_tip_info$Balance_Fischer)
  colnames(fp1) <- c("label", "balance")
  for_plot3 <- left_join(for_plot2, fp1, copy = TRUE)


  fp2 <- cbind(df_tip_info$Statistic, df_tip_info$Information)
  colnames(fp2) <- c("label", "info")
  for_plot4 <- left_join(for_plot2, fp2, copy = TRUE)

  for_plot3$y <- -0.6
  for_plot4$y <- -0.8

  p6 <-
    ggplot() +
    geom_rect(data = rect_plot,
              aes(xmin = ymin, xmax = ymax, ymin = xmin, ymax = xmax,
                  fill = as.factor(categ)), alpha = 1) +
    MetBrewer::scale_fill_met_d(name = focal_palette, direction = 1) +
    new_scale_fill() +
    geom_segment(data = segment(ddata),
                 aes(x = x, y = y, xend = xend, yend = yend)
    ) +
    geom_text(data = label(ddata),
              aes(x = x, y = y, label = label, hjust = 0),
              size = 3
    ) +
    geom_tile(data = for_plot3,
              aes(x = x, y = y, fill = balance, width = 1, height = 0.2)) +
    MetBrewer::scale_fill_met_d(name = "Hiroshige", direction = 1) +
    labs(fill = "Balance statistic?") +

    new_scale_fill() +
    geom_tile(data = for_plot4,
              aes(x = x, y = y, fill = info, width = 1, height = 0.2)) +

    MetBrewer::scale_fill_met_d(name = "Morgenstern", direction = -1) +

    labs(fill = "Information used") +
    coord_flip() +
    scale_y_reverse() +
    theme_classic() +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ylab("1 - abs(correlation)") +
    xlab("") +
    ggtitle("Simulated data")


  return(list("plot" = p6,
              "clustering" = clust_ref))
}

get_plot <- function(family_name, res.dist, ref_clustering) {

  hc <- hclust(res.dist, method = "average")
  dend0 <- stats::as.dendrogram(hc)
  ddata <- dendro_data(hc, type = "rectangle")

  for_plot2 <- as_tibble(ddata$labels)

  fp1 <- cbind(names(ref_clustering), ref_clustering)
  colnames(fp1) <- c("label", "cluster")
  for_plot3 <- left_join(for_plot2, fp1, copy = TRUE)

  for_plot3$y <- -0.6

  p6 <-
    ggplot() +
    geom_segment(data = segment(ddata),
                 aes(x = x, y = y, xend = xend, yend = yend)
    ) +
    geom_text(data = label(ddata),
              aes(x = x, y = y, label = label, hjust = 0),
              size = 3
    ) +
    new_scale_fill() +
    geom_tile(data = for_plot3,
              aes(x = x, y = y, fill = cluster, width = 1, height = 0.2)) +
    MetBrewer::scale_fill_met_d(name = "Derain", direction = 1) +
    labs(fill = "Information used") +
    coord_flip() +
    scale_y_reverse() +
    theme_classic() +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(legend.position = "none") +
    ylab("1 - abs(correlation)") +
    xlab("") +
    ggtitle(family_name)
  return(p6)
}

cor1 <- read.table("overall_cor.txt")
cor1 <- as.matrix(cor1)

cor1 <- as.data.frame(cor1)
cor1 <- as_tibble(cor1)
cor1 <- cor1 %>%
  mutate_at(1:67, as.numeric)

cor.dist <- stats::as.dist(1 - abs(as.matrix(cor1)))

ref_plot <- make_plot("Simulated data", cor.dist, 0.8)

all_plots <- list()
cnt <- 1
for (x in c("BD", "DDD", "PBD","SSE")) {
  file_name <- paste0(x, "_cor.txt")
  cor2 <- read.table(file_name)
  local_names <- names(cor2)
  rownames(cor2) <- local_names
  colnames(cor2) <- local_names
  cor.dist <- stats::as.dist(1 - abs(as.matrix(cor2)))
  local_plot <- get_plot(x, cor.dist, ref_plot$clustering)
  all_plots[[cnt]] <- local_plot
  cnt <- cnt + 1
}

library(gridExtra)
pdf("Figure_5_70.pdf", width = 25, height = 12)
grid.arrange(ref_plot$plot, arrangeGrob(all_plots[[1]], all_plots[[2]],
                                              all_plots[[3]], all_plots[[4]]), ncol = 2)
dev.off()
