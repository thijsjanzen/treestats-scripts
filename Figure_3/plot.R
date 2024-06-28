require(dendextend)
require(ggdendro)


require(rphylopic)

phylopics <- c()
phylopics <- rbind(phylopics, c("Birds", "Erithacus rubecula", 0.87)) # Robin
phylopics <- rbind(phylopics, c("Mammals", "Mus musculus", 0.8))
phylopics <- rbind(phylopics, c("Ferns", "Marattia", 0.63))
phylopics <- rbind(phylopics, c("Vascular Plants", "Platanus occidentalis", 0.7))
phylopics <- rbind(phylopics, c("Cartaliginous Fish", "Isurus oxyrinchus", 0.4))
phylopics <- rbind(phylopics, c("Ray finned Fish", "Sebastes alutus", 0.7))
phylopics <- rbind(phylopics, c("Amphibians", "Micrixalus elegans", 0.7))

make_plot <- function(family_name, res.dist, xmax) {

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

  require(ggplot2)

  tip_info <- read_tsv("../datasets/families.txt")

  to_remove <- which(tip_info$Statistic == "number_of_lineages")

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

  colors_used <- ggpubr::get_palette(palette = "Paired",
                                     k = max(rect_plot$categ))

  lvls <- sort(unique(as.numeric(rect_plot$categ)))
  rect_plot$categ <- factor(rect_plot$categ, levels = lvls)


  p6 <-
    ggplot() +
    geom_rect(data = rect_plot,
              aes(xmin = ymin, xmax = ymax, ymin = xmin, ymax = xmax,
                  fill = categ), alpha = 1) +
    scale_fill_manual(values = colors_used) +
    #MetBrewer::scale_fill_met_d(name = "Derain", direction = 1) +
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

    labs(fill = "Balance") +

    new_scale_fill() +
    geom_tile(data = for_plot4,
              aes(x = x, y = y, fill = info, width = 1, height = 0.2)) +
    MetBrewer::scale_fill_met_d(name = "Morgenstern", direction = -1) +

    labs(fill = "Information used") +

    coord_flip() +
    scale_y_reverse() +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylab("1 - abs(correlation)") +
    xlab("") +
    ggtitle(family_name)

  p6

  return(list("plot" = p6,
              "clustering" = clust_ref))
}

get_plot <- function(family_name_plot, res.dist, ref_clustering,
                     phylopic_name) {

  colors_used <- ggpubr::get_palette(palette = "Paired",
                                     k = max(ref_clustering))

  hc <- hclust(res.dist, method = "average")
  dend0 <- stats::as.dendrogram(hc)
  ddata <- dendro_data(hc, type = "rectangle")

  for_plot2 <- as_tibble(ddata$labels)

  fp1 <- cbind(names(ref_clustering), ref_clustering)
  colnames(fp1) <- c("label", "cluster")
  for_plot3 <- left_join(for_plot2, fp1, copy = TRUE)

  for_plot3$y <- -0.6

  lvls <- sort(unique(as.numeric(for_plot3$cluster)))
  for_plot3$cluster <- factor(for_plot3$cluster, levels = lvls)

  phylopic_species <- phylopics[phylopics[,1] == phylopic_name, 2]
  phylopic_pos <- as.numeric(phylopics[phylopics[,1] == phylopic_name, 3])

  p6 <-
    ggplot() +
    geom_segment(data = segment(ddata),
                 aes(x = x, y = y, xend = xend, yend = yend)
    ) +
    geom_text(data = label(ddata),
              aes(x = x, y = y, label = label, hjust = 0),
              size = 3
    ) +
    geom_phylopic(aes(x = max(ddata$segments), y = phylopic_pos, name = phylopic_species), verbose = TRUE) +
    new_scale_fill() +
    geom_tile(data = for_plot3,
              aes(x = x, y = y, fill = cluster, width = 1, height = 0.2)) +
    #MetBrewer::scale_fill_met_d(name = "Derain", direction = 1) +
    scale_fill_manual(values = colors_used) +
    labs(fill = "Information used") +
    coord_flip() +
    scale_y_reverse() +
    theme_classic() +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(legend.position = "none") +
    ylab("1 - abs(correlation)") +
    xlab("") +
    ggtitle(family_name_plot)
  return(p6)
}

# correlation across entire dataset:
cor1 <- read.table("master_cor_phy.txt")
cor1 <- as.matrix(cor1)
to_remove <- which(colnames(cor1) == "number_of_lineages")

cor1 <- cor1[-to_remove, -to_remove]
cor1 <- as.data.frame(cor1)
cor1 <- as_tibble(cor1)
cor1 <- cor1 %>%
  mutate_at(1:69, as.numeric)

cor.dist <- stats::as.dist(1 - abs(as.matrix(cor1)))

ref_plot <- make_plot("Full data set (221 trees)", cor.dist, 0.8)
ref_plot$plot

all_plots <- list()
cnt <- 1

plot_names <- c("Ferns (7 trees)", "Flowering Plants (9 trees)",
                "Birds (122 trees)", "Mammals (24 trees)", "Cartilaginous Fish (5 trees)", "Ray-finned Fish (42 trees)", "Amphibia (12 trees)")

data_names <- c("Ferns", "Vascular Plants", "Birds", "Mammals", "Cartaliginous Fish", "Ray finned Fish", "Amphibians")

require(ggdendro)

for (i in 1:7) {
  # load precalculated correlations:
  file_name <- paste0("cor_emp_phy_", data_names[i], ".txt")
  cor2 <- read.table(file_name)
  local_names <- names(cor2)
  to_remove <- which(local_names == "number_of_lineages")
  local_names <- local_names[-to_remove]
  cor2 <- cor2[-to_remove, -to_remove]

  cor.dist <- stats::as.dist(1 - abs(as.matrix(cor2)))
  names(cor.dist) <- local_names
  local_plot <- get_plot(plot_names[i], cor.dist, ref_plot$clustering,
                         data_names[i])
  all_plots[[cnt]] <- local_plot
  cnt <- cnt + 1
}

library(gridExtra)
pdf("Figure_3.pdf", width = 35, height = 12)
px <- grid.arrange(ref_plot$plot,
                   arrangeGrob(all_plots[[3]], all_plots[[5]], all_plots[[1]], all_plots[[7]],
                               all_plots[[4]], all_plots[[6]], all_plots[[2]], nrow = 2), ncol = 2)
dev.off()




