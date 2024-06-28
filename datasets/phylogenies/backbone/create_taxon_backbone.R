dd <- read.table("mrca_backbone.txt", header = TRUE)
tip_names <- colnames(dd)[2:8]
mat_dd <- dd[1:7, 2:8]
rownames(mat_dd) <- colnames(mat_dd)
mat_dd <- 2 * mat_dd
mat_dd <- as.dist(mat_dd)

xx <- phangorn::upgma(mat_dd)
plot(xx)
ape::cophenetic.phylo(xx)
ape::branching.times(xx)

ape::write.tree(xx, file = "taxon_backbone.tree")

