begin_tree <- "/source_trees/"
begin_fam <- "/family_lists/"

taxa <- c("birds", "mammals", "amphibia", "ferns", "vascular_plants", "sharks", "ray_finned_fish")

for (focal_taxa in taxa) {
  super_tree <- ape::read.tree(paste0(begin_tree, focal_taxa, ".tree"))

  tree_collection <- readRDS(paste0("/fracced/", focal_taxa, "_fracced.rds"))

  tips_info <- read_tsv(file = paste0(begin_fam, focal_taxa, "_lookuplist.txt"),
                        show_col_types = FALSE)

  for (focal_fam in unique(tips_info$family)) {
    focal_tips <- subset(tips_info, tips_info$family == focal_fam)$species
    tips_to_remove <- focal_tips
    tip_to_keep <- tips_to_remove[1]
    tips_to_remove <- tips_to_remove[-1]

    super_tree <- ape::drop.tip(super_tree, tips_to_remove)
    super_tree$tip.label[super_tree$tip.label == tip_to_keep] <- focal_fam
  }

  # focal families:
  selected_families <- names(tree_collection)

  families_to_drop <- super_tree$tip.label[!(super_tree$tip.label %in% selected_families)]
  super_tree <- ape::drop.tip(super_tree, families_to_drop)
  super_tree$node.label <- NULL
  cat(focal_taxa, length(super_tree$tip.label), treestats::crown_age(super_tree), "\n")
  ape::write.tree(super_tree, file = paste0("/family_level_trees/", focal_taxa, "_fam.tree"))
}
