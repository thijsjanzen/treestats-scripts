require(geiger)
require(phylobase)
require(tidyverse)
require(ggfortify)

require(factoextra)

# read master tree from datelife analysis (see elsewhere)
pruned_tree <- ape::read.tree("../datasets/new_data/backbone/taxon_backbone.tree")


# tree collection files
tree_collection_files <- c("../datasets/new_data/fracced/amphibia_fracced.rds",
                           "../datasets/new_data/fracced/birds_fracced.rds",
                           "../datasets/new_data/fracced/ferns_fracced.rds",
                           "../datasets/new_data/fracced/mammals_fracced.rds",
                           "../datasets/new_data/fracced/ray_finned_fish_fracced.rds",
                           "../datasets/new_data/fracced/sharks_fracced.rds",
                           "../datasets/new_data/fracced/vascular_plants_fracced.rds")

found_stats <- c()

for (i in 1:length(tree_collection_files)) {
  tree_collection <- readRDS(tree_collection_files[i])

  families <- names(tree_collection)
  cat(i, "\n")
  for (j in 1:length(tree_collection)) {

    focal_tree <- tree_collection[[j]]

    if (length(focal_tree$tip.label) >= 10) {

      all_stats <- treestats::calc_all_stats(focal_tree, FALSE)
      to_add <- unlist(all_stats)
      to_add <- c(families[j], to_add)

      found_stats <- rbind(found_stats, to_add)
    }
  }
}

colnames(found_stats) <- c("Family", names(all_stats))
found_stats <- as_tibble(found_stats)

found_stats <- found_stats[!is.na(found_stats$sackin), ]

# make a second tree with the family names:
taxa_tree <- pruned_tree

taxa_names <- c("Amphibians", "Birds", "Ferns", "Mammals", "Ray finned Fish", "Cartaliginous Fish", "Flowering Plants")

found_stats$Taxa <- NA

for (i in 1:length(tree_collection_files)) {
  tree_collection <- readRDS(tree_collection_files[i])

  families <- names(tree_collection)
  indices <- which(taxa_tree$tip.label %in% families)
  taxa_tree$tip.label[indices] <- taxa_names[i]

  indices2 <- which(found_stats$Family %in% families)
  found_stats$Taxa[indices2] <- taxa_names[i]
}

# make sure everything is numeric:
found_stats2 <- found_stats %>%
  mutate_at(2:71, as.numeric)

# select focal data:
df <- data.matrix(found_stats2[, 2:71])
rownames(df) <- found_stats2$Family
# now, we correct everything for size
for (i in 1:ncol(df)) {
  y <- df[, i]
  x <- df[, 11] # number of lineages
  A <- lm(y~x)
  df[, i] <- A$residuals
}

df <- df[, -11] # number of lineages



df <- scale(df)


# phylo corrected pca:
pca.res <- phytools::phyl.pca(tree = pruned_tree,
                              Y = df)

diag(pca.res$Eval)[which(diag(pca.res$Eval) < 0)] <- 0

vz <- phytools::as.princomp(pca.res)

color_used <- ggpubr::get_palette(palette = "RdBu", k = 6)


found_stats2$Taxa_2 <- factor(found_stats2$Taxa, levels = c("Birds", "Mammals", "Ferns", "Flowering Plants", "Cartaliginous Fish", "Ray finned Fish", "Amphibians"))

manual_colors2 <- c( "#CC500A", #birds
                     "#9F55D6", # mammals
                     "#D4E995", # ferns
                     "#145F16", # plants 2
                     "#6CA8E9", # fish 1
                     "#3A7DB1", # fish 2
                     "grey")    # Amphibians

p2 <- autoplot(vz,
               data = found_stats2,
               colour = 'Taxa_2',
               loadings = FALSE,
               loadings.label = FALSE,
               frame = TRUE,
               frame.type = 'norm') +
  theme_minimal() +
  scale_color_manual(values = manual_colors2) +
  scale_fill_manual(values = manual_colors2) +
  theme(legend.position = "top") +
  labs(colour = "") +
  labs(fill = "")
p2

ggsave(filename = "fig_2_c.pdf",
       width = 6, height = 6)


ncp <- 10
eig <- get_eigenvalue(vz)
eig <- eig[1:min(ncp, nrow(eig)), , drop = FALSE]

eig <- eig[, 2]
text_labels <- paste0(round(eig, 1), "%")
bar_width <- 0.7
barfill <- "#32648EFF"
barcolor <- barfill
df.eig <- data.frame(dim = factor(1:length(eig)), eig = eig)
p <- ggplot(df.eig, aes(dim, eig, group = 1)) +
  geom_bar(stat = "identity", fill = viridis::mako(n = 10, begin = 0.4, end = 0.7),
           width = bar_width) +
  #geom_line() +
  geom_text(label = text_labels, vjust = -0.4,
            hjust = 0.45, size = 3) +
  labs(x = "Dimensions",
       y = "Percentage of explained variances") +
  theme_classic()
p

ggsave(p, filename = "fig_2_b.pdf",
       width = 5 , height = 5)


var <- get_pca_var(vz)

for_plot <- var$cos2
to_plot <- c()
for (i in 1:6) {
  vx <- for_plot[, i]
  vy <- sort(vx, decreasing = TRUE)
  vz <- vy
  to_add <- cbind(names(vz), vz, i)
  to_plot <- rbind(to_plot, to_add)
}

colnames(to_plot) <- c("stat", "cos2", "pca_axis")
to_plot <- as_tibble(to_plot)
to_plot$cos2 <- as.numeric(to_plot$cos2)

require(tidytext)

colors_used <- ggpubr::get_palette(palette = "Paired", k = max(6))

p2 <- to_plot %>%
  mutate("pca_lab" = paste0("PCA axis ", pca_axis)) %>%
  mutate(pca_lab = as.factor(pca_lab),
         stat = reorder_within(stat, cos2, pca_lab)) %>%
  ggplot(aes(stat, cos2, fill = pca_lab)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~pca_lab, scales = "free") +
  coord_flip() +
  scale_x_reordered() +
  scale_fill_manual(values = colors_used) +
  scale_y_continuous(expand = c(0,0)) +
  theme_minimal() +
  labs(y = "Contribution to PCA axis (cos2)",
       x = NULL)

ggsave(p2, filename = "figure_S1.pdf", width = 12, height = 10)


# get all pairwise correlations corrected:

breakz <- seq(0, 1, length.out = 99)

res.cor3 <- read.table("../Figure_3/master_cor_phy.txt")
to_remove <- which(colnames(res.cor3) == "number_of_lineages")
res.cor3 <- res.cor3[-to_remove, -to_remove]

diag(res.cor3) <- NA
res.cor4 <- abs(res.cor3)

color_used <- ggpubr::get_palette(palette = "GnBu", k = 99)
hm1 <- pheatmap::pheatmap(mat = res.cor4,
                          color = color_used,
                          breaks = breakz,
                          treeheight_col = 0,
                          treeheight_row = 0,
                        #  clustering_distance_rows = "correlation",
                       #   clustering_distance_cols = "correlation",
                          clustering_method = "average",
                          fontsize_col = 8,
                          fontsize_row = 8)

pdf("fig_2_a.pdf",
    width = 12, height = 10)
hm1
dev.off()
