

pruned_tree <-  readRDS("../datasets/pruned_tree.rds")

tree_collection_files <- c("../datasets/MammalTrees.rds",
                           "../datasets/BirdTrees.rds",
                           "../datasets/AmphibiaTrees.rds",
                           "../datasets/SquamateTrees.rds")

found_stats <- c()
for (i in 1:length(tree_collection_files)) {
  tree_collection <- readRDS(tree_collection_files[i])

  families <- names(tree_collection)

  for (j in 1:length(tree_collection)) {
    focal_tree <- tree_collection[[j]]$tree

    all_stats <- treestats::calc_all_stats(focal_tree, FALSE)
    to_add <- unlist(all_stats)
    to_add <- c(families[j], to_add)

    found_stats <- rbind(found_stats, to_add)
  }
}


require(tidyverse)

colnames(found_stats) <- c("Family", names(all_stats))
found_stats <- as_tibble(found_stats)

taxa_tree <- pruned_tree

taxa_names <- c("Mammals","Birds", "Amphibians", "Squamates")

found_stats$Taxa <- NA

for (i in 1:length(tree_collection_files)) {
  tree_collection <- readRDS(tree_collection_files[i])

  families <- names(tree_collection)
  indices <- which(taxa_tree$tip.label %in% families)
  taxa_tree$tip.label[indices] <- taxa_names[i]

  indices2 <- which(found_stats$Family %in% families)
  found_stats$Taxa[indices2] <- taxa_names[i]

}

found_stats2 <- found_stats %>%
  mutate_at(2:55, as.numeric)

df <- as.data.frame(found_stats2[, 2:55])

res.cor <- stats::cor(df, method = "pearson")

res.cor2 <- res.cor

get_cor <- function(local_stats, local_tree) {
  res.cor <- stats::cor(as.data.frame(local_stats), method = "pearson")

  res.cor2 <- res.cor

  pb <- txtProgressBar(max = length(res.cor[1, ]), style = 3)

  for (i in 1:length(res.cor[1, ])) {
    setTxtProgressBar(pb, i)
    for (j in 1:length(res.cor[1, ])) {
      stat1 <- colnames(res.cor)[i]
      stat2 <- colnames(res.cor)[j]

      if (stat1 == "number_of_lineages") {
        if (stat2 != "number_of_lineages") {
          x <- unlist(as.vector(local_stats[stat1]))
          y <- unlist(as.vector(local_stats[stat2]))
          z <- unlist(as.vector(local_stats["number_of_lineages"]))

          a1 <- nlme::gls(y~z, correlation = ape::corBrownian(1, local_tree))
          res.cor2[i, j] <- a1$coefficients[[2]]
        }
      } else if (stat2 != "number_of_lineages") {
        if (i != j) {
          x <- unlist(as.vector(local_stats[stat1]))
          y <- unlist(as.vector(local_stats[stat2]))
          z <- unlist(as.vector(local_stats["number_of_lineages"]))

          a1 <- nlme::gls(y~z, correlation = ape::corBrownian(1, local_tree))
          a2 <- nlme::gls(x~z, correlation = ape::corBrownian(1, local_tree))

          found_cor <- cor(a1$residuals, a2$residuals)

          res.cor2[i, j] <- found_cor
        }
      }
    }
  }
  return(res.cor2)
}

master_cor <- get_cor(df, pruned_tree)

write.table(master_cor, "master_cor.txt", quote = FALSE)

# now we go over the taxa
for (x in unique(found_stats$Taxa)) {
    cat(x, "\n")
   stat_subset <- subset(found_stats, found_stats$Taxa == x)
   focal_families <- as.vector(stat_subset$Family)

   focal_data <- stat_subset %>%
     select(-c("Family", "Taxa"))

   focal_data <- focal_data %>%
     mutate_at(1:54, as.numeric)

   master_tree <- pruned_tree
   to_remove <- master_tree$tip.label[!(master_tree$tip.label %in% focal_families)]
   focal_tree <- ape::drop.tip(master_tree, to_remove)

   local_cor <- get_cor(focal_data, focal_tree)
   write.table(local_cor, paste0("cor_emp_", x,".txt"), quote = FALSE)
}
