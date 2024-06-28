pruned_tree <- ape::read.tree("../datasets/new_data/backbone/taxon_backbone.tree")


# tree collection files
tree_collection_files <- c("../datasets/new_data/fracced/amphibia_fracced.rds",
                           "../datasets/new_data/fracced/birds_fracced.rds",
                           "../datasets/new_data/fracced/ferns_fracced.rds",
                           "../datasets/new_data/fracced/mammals_fracced.rds",
                           "../datasets/new_data/fracced/ray_finned_fish_fracced.rds",
                           "../datasets/new_data/fracced/sharks_fracced.rds",
                           "../datasets/new_data/fracced/vascular_plants_fracced.rds")

taxa_names <- c("Amphibians", "Birds", "Ferns", "Mammals", "Ray finned Fish", "Cartaliginous Fish", "Vascular Plants")


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


require(tidyverse)

num_stats <- length(names(all_stats))

colnames(found_stats) <- c("Family", names(all_stats))
found_stats <- as_tibble(found_stats)

taxa_tree <- pruned_tree


found_stats$Taxa <- NA

for (i in 1:length(tree_collection_files)) {
  tree_collection <- readRDS(tree_collection_files[i])

  families <- names(tree_collection)
  indices <- which(taxa_tree$tip.label %in% families)
  taxa_tree$tip.label[indices] <- taxa_names[i]

  indices2 <- which(found_stats$Family %in% families)
  found_stats$Taxa[indices2] <- taxa_names[i]
}

found_stats <- found_stats %>%
  filter(!is.na(sackin))

found_stats2 <- found_stats %>%
  mutate_at(2:(num_stats + 1), as.numeric)

df <- as.data.frame(found_stats2[, 2:(num_stats + 1)])

res.cor <- stats::cor(df, method = "pearson")

res.cor2 <- res.cor

get_cor <- function(local_stats, local_tree) {

  local_stats <- local_stats[!is.na(local_stats$gamma), ]
  local_stats <- local_stats[!is.na(local_stats$sackin), ]

  to_keep <- local_stats$Family
  to_drop <- local_tree$tip.label[!(local_tree$tip.label %in% local_stats$Family)]

  if (length(to_drop)) {
    cat("dropping: " ,to_drop, "\n")
    local_tree <- ape::drop.tip(local_tree, to_drop)
  }

  testit::assert(all.equal(local_stats$Family, local_tree$tip.label))

  sp <- local_stats$Family
  bm <- ape::corBrownian(value = 1, phy = local_tree, form =~ sp)

  local_stats2 <- local_stats %>%
    mutate_at(2:(num_stats + 1), as.numeric)

  local_stats2 <- as.data.frame(local_stats2[, 2:(num_stats + 1)])

  res.cor <- stats::cor(as.data.frame(local_stats2), method = "pearson")

  res.cor2 <- res.cor

  pb <- txtProgressBar(max = length(res.cor[1, ]), style = 3)

  for (i in 1:length(res.cor[1, ])) {
    setTxtProgressBar(pb, i)
    for (j in 1:length(res.cor[1, ])) {
      stat1 <- colnames(res.cor)[i]
      stat2 <- colnames(res.cor)[j]

      if (stat1 != stat2) {
        if (stat1 != "number_of_lineages" && stat2 != "number_of_lineages") {
          x <- unlist(as.vector(local_stats2[stat1]))
          y <- unlist(as.vector(local_stats2[stat2]))
          z <- unlist(as.vector(local_stats2["number_of_lineages"]))

          a1 <- nlme::gls(y~z, correlation = bm)
          a2 <- nlme::gls(x~z, correlation = bm)

          found_cor <- cor(a1$residuals, a2$residuals)


          res.cor2[i, j] <- found_cor
        }
      }
    }
  }
  return(res.cor2)
}

master_cor <- get_cor(found_stats2, pruned_tree)

write.table(master_cor, "master_cor_phy.txt", quote = FALSE)

# now we go over the taxa
for (x in unique(found_stats$Taxa)) {
    cat(x, "\n")
   stat_subset <- subset(found_stats, found_stats$Taxa == x)
   focal_families <- as.vector(stat_subset$Family)

   master_tree <- pruned_tree
   to_remove <- master_tree$tip.label[!(master_tree$tip.label %in% focal_families)]
   focal_tree <- ape::drop.tip(master_tree, to_remove)

   local_cor <- get_cor(stat_subset, focal_tree)
   write.table(local_cor, paste0("cor_emp_phy_", x,".txt"), quote = FALSE)
}
