
require(datelife)

tree_collection_files <- c("MammalTrees.rds",
                           "BirdTrees.rds",
                           "AmphibiaTrees.rds",
                           "SquamateTrees.rds")

# we collect all tip labels
all_names <- c()

for (i in 1:length(tree_collection_files)) {
  tree_collection <- readRDS(tree_collection_files[i])

  for (j in 1:length(tree_collection)) {
    focal_tree <- tree_collection[[j]]$tree
    all_names <- c(all_names, focal_tree$tip.label)
  }
}

# and we perform a datelife search
# because there are so many species, this may take very long (~1 day).
# please note that not all species are on datelife, but we don't know
# a priori which ones will drop out. Later on, we will check that for
# each phylogeny in the dataset we have at least one representative.
res <- datelife::datelife_search(input = all_names,
                                 summary_format = "phylo_biggest")

master_tree <- res
for (i in 1:length(tree_collection_files)) {
  tree_collection <- readRDS(tree_collection_files[i])
  families <- names(tree_collection)

  for (j in 1:length(tree_collection)) {

    focal_tree <- tree_collection[[j]]$tree

    num_found <- 0
    for (k in 1:length(focal_tree$tip.label)) {
      focal_sp <- focal_tree$tip.label[k]
      index <- which(master_tree$tip.label == focal_sp)
      if (length(index)) {
        num_found <- num_found + 1
        master_tree$tip.label[index] <- families[j]
      }
    }
  }
}

all_families <- unique(master_tree$tip.label)

pruned_tree <- master_tree
# now we need to prune down the master tree
for (focal_fam in all_families) {
  indices <- which(pruned_tree$tip.label == focal_fam)
  indices <- indices[-1]
  pruned_tree <- ape::drop.tip(pruned_tree, indices)
}

saveRDS(pruned_tree, "pruned_tree.rds")
