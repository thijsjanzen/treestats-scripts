setwd("/Users/thijsjanzen/Library/CloudStorage/Dropbox/projects/treestats/Scripts/datasets/new_data")

begin_tree <- "/Users/thijsjanzen/Library/CloudStorage/Dropbox/projects/treestats/Scripts/datasets/new_data/source_trees/"
begin_fam <- "/Users/thijsjanzen/Library/CloudStorage/Dropbox/projects/treestats/Scripts/datasets/new_data/family_lists/"

begin <- "/Users/thijsjanzen/Library/CloudStorage/Dropbox/projects/treestats/Scripts/datasets/new_data/"

taxa <- c("birds", "mammals", "amphibia", "ferns", "vascular_plants", "sharks", "ray_finned_fish")

for (focal_taxa in taxa) {
  super_tree <- ape::read.tree(paste0(begin_tree, focal_taxa, ".tree"))
  tips_info <- read_tsv(file = paste0(begin_fam, focal_taxa, "_lookuplist.txt"))

  phylos <- list()
  all_names <- c()
  all_fams <- unique(tips_info$family)

  found <- c()
  for (f in all_fams) {
    tips <- subset(tips_info, tips_info$family == f)
    num_tips <- length(tips$species)
    if (num_tips >= 10) {
      vv <- tryCatch(expr = {taxizedb::downstream(f, db = "ncbi", downto = "species")},
                     error = function(e) {return(NA) }) #nolint
      num_ncbi <- 0
      if (!sum(is.na(vv))) {
        vv2 <- subset(vv[[1]], vv[[1]]$rank == "species")

        vv2$is_species <- 0
        for (k in 1:length(vv2$childtaxa_name)) {
          ax <- strsplit(vv2$childtaxa_name[k], split = " ")
          if (length(ax[[1]]) == 2) {
            vv2$is_species[k] <- 1
          }
        }
        vv3 <- subset(vv2, vv2$is_species == 1)

        num_ncbi <- sum(vv2$is_species)
      }

      fraction <- num_tips / num_ncbi

      if (fraction >= 0.8) {
        to_remove <- !(super_tree$tip.label %in% tips$species)
        pruned_tree <- ape::drop.tip(super_tree, super_tree$tip.label[to_remove])

        phylos[[length(phylos) + 1]] <- pruned_tree
        all_names <- c(all_names, f)
        cat(f, fraction, num_tips, length(pruned_tree$tip.label), "\n")
      }
      found <- rbind(found, c(f, fraction, num_tips, num_ncbi))
    }

  }

  names(phylos) <- all_names

  saveRDS(phylos, file = paste0(begin, "fracced/",focal_taxa, "_fracced.rds"))
  length(phylos)
  write.table(found, file = paste0(begin, "fractions/", focal_taxa, "_fracs.txt"))
}


