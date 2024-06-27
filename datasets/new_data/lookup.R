require(taxizedb)
require(tibble)
require(readr)
require(ape)

setwd("/Users/thijsjanzen/Library/CloudStorage/Dropbox/projects/treestats/Scripts/datasets/new_data")

begin_tree <- "/Users/thijsjanzen/Library/CloudStorage/Dropbox/projects/treestats/Scripts/datasets/new_data/source_trees/"
begin_fam <- "/Users/thijsjanzen/Library/CloudStorage/Dropbox/projects/treestats/Scripts/datasets/new_data/family_lists/"

taxa <- c("birds", "mammals", "amphibia", "ferns", "vascular_plants", "sharks", "ray_finned_fish")

for (focal_taxa in taxa) {
  super_tree <- ape::read.tree(paste0(begin_tree, focal_taxa, ".tree"))

  all_fams <- c()

  for (i in 1:length(super_tree$tip.label)) {
    vv <- tryCatch(expr = {
      taxizedb::classification(super_tree$tip.label[i], db = "ncbi")
    },
    error = function(e) {return(NA)})

    found_fam <- NA

    if (!is.na(vv)) {
      found_fam <- vv[[1]]$name[vv[[1]]$rank == "family"]
    } else {
      vv <- tryCatch(expr = {taxize::classification(super_tree$tip.label[i],
                                                    db = "gbif", rows = 1)},
                     error = function(e) {return(NA) })
      if (!is.na(vv)) {
        ax <- vv[[1]]$name[vv[[1]]$rank == "family"]
        if (length(ax)) {
          found_fam <- vv[[1]]$name[vv[[1]]$rank == "family"]
        }
      }
    }
    to_add <- c(super_tree$tip.label[i], found_fam)
    all_fams <- rbind(all_fams, to_add)
    cat(focal_taxa, to_add, "\n")
  }

  colnames(all_fams) <- c("species", "family")
  all_fams <- tibble::as_tibble(all_fams)
  readr::write_tsv(all_fams, file = paste0(begin_fam, focal_taxa, "_lookuplist.txt"))
}



