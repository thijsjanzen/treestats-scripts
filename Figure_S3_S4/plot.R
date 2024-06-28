library(treestats)
library(ggplot2)
library(tidyverse)
require(picante)
require(apTreeshape)
require(castor)
require(treebalance)
require(ape)

# RUtreebalance.R was downloaded from:
# https://github.com/robjohnnoble/RUtreebalance
source("/Users/thijsjanzen/Library/CloudStorage/Dropbox/projects/treestats/RUtreebalance-main/RUtreebalance.R")

calc_max_ladder <- function(phy) {
  ladders <- phyloTop::ladderSizes(phy)$ladderSizes
  return(max(ladders))
}

calc_colless_quad <- function(phy) {
  return(treebalance::collessI(phy, method = "quadratic"))
}

calc_colless_corr <- function(phy) {
  return(treebalance::collessI(phy, method = "corrected"))
}

calc_mean_br <- function(phy) {
  mbr <- mean(phy$edge.length)
  return(mbr)
}

calc_phylodiv <- function(phy) {
  return(sum(phy$edge.length))
}

calc_psv <- function(phy) {
  n <- length(phy$tip.label)
  sample_mat <- matrix(data = 1, nrow = n, ncol = n)
  colnames(sample_mat) <- phy$tip.label

  return(picante::psv(sample_mat, phy)[1, 1])
}

calc_mntd <- function(focal_tree) {
  n <- length(focal_tree$tip.label)
  sample_mat <- matrix(data = 1, nrow = n, ncol = n)
  colnames(sample_mat) <- focal_tree$tip.label

  picante::mntd(sample_mat, cophenetic(focal_tree),
                abundance.weighted = FALSE)[[1]]
}

calc_mpd <- function(focal_tree) {
  n <- length(focal_tree$tip.label)
  sample_mat <- matrix(data = 1, nrow = n, ncol = n)
  colnames(sample_mat) <- focal_tree$tip.label

  picante::mpd(sample_mat, cophenetic(focal_tree),
               abundance.weighted = FALSE)[[1]]
}

calc_var_mpd <- function(phy) {
  a2 <- cophenetic(phy)
  a2 <- a2[lower.tri(a2)]
  return(var(a2, na.rm = TRUE))
}

nLTT_base_R <- function(phy) {  # nolint
  empty_tree <- ape::read.tree(text = "(1:4,2:4):0;")
  return(nLTT::nltt_diff(phy, empty_tree))
}

calc_beta <- function(phy) {
  return(apTreeshape::maxlik.betasplit(phy))
}

calc_colless <- function(phy) {
  return(castor::tree_imbalance(phy, "Colless"))
}

calc_sackin <- function(phy)  {
  # return(castor::tree_imbalance(phy, "Sackin")) Castor is incorrect
  return(apTreeshape::sackin(apTreeshape::as.treeshape(phy)))
}

calc_blum <- function(phy) {
  # return(castor::tree_imbalance(phy, "Blum"))
  return(treebalance::sShapeI(phy, logbase = exp(1)))
}

calc_j_one <- function(phy) {
  return(J1_index(phy))
}

ref_tree_height <- function(phy) {
  th <- max(ape::branching.times(phy))
  if (!is.null(phy$root.edge)) th <- th + phy$root.edge
  return(th)
}

calc_i <- function(phy) {
  return(treebalance::IbasedI(phy, method = "mean", correction = "prime",
                              logs = FALSE))
}

calc_minmax_adjacency <- function(phy) {
  df <- as.data.frame(cbind(phy$edge,
                            weight = phy$edge.length))
  g <- igraph::graph_from_data_frame(df, directed = FALSE)

  adj_mat <- igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
  ref <- eigen(adj_mat)$values
  ref <- round(ref, digits = 10)
  return(ref)
}

calc_minmax_laplace <- function(phy) {
  df <- as.data.frame(cbind(phy$edge,
                            weight = phy$edge.length))
  g <- igraph::graph_from_data_frame(df, directed = FALSE)

  lapl_mat <- igraph::laplacian_matrix(g, normalized = FALSE, sparse = FALSE)
  ref <- eigen(lapl_mat)$values
  ref <- round(ref, digits = 10)
  return(ref)
}



tree_list <- list()
used_ntaxa <- c()
cnt <- 1
for (ntaxa in 10^seq(from = 1, to = 3, length.out = 10)) {
  sim_trees <- list()
  for (r in 1:10) {
    sim_trees[[r]] <- ape::rphylo(n = floor(ntaxa), birth = 1, death = 0)
  }

  tree_list[[cnt]] <- sim_trees
  used_ntaxa[cnt] <- floor(ntaxa)
  cnt <- cnt + 1
}


collect_time <- function(tree, func) {
  t1 <- Sys.time()
  func(tree)
  t2 <- Sys.time()
  return( difftime(t2, t1, units = "secs"))
}

compare_speed <- function(trees, r_func, cpp_func, used_ntaxa, text_label) {

  to_add <- matrix(nrow = length(trees) * 2, ncol = 4)

  for (r in 1:length(trees)) {
    t1 <- collect_time(trees[[r]], r_func)
    t2 <- collect_time(trees[[r]], cpp_func)

    add1 <- c(used_ntaxa, "R equivalent", t1, text_label)
    add2 <- c(used_ntaxa, "treestats", t2, text_label)
    index <- 1 + (r - 1) * 2
    to_add[index, ] <- add1
    to_add[index + 1, ] <- add2
  }

  cat(used_ntaxa, text_label, "\n")

  return(to_add)
}


found <- c()
for (i in 1:length(tree_list)) {

  found <- rbind(found, compare_speed(tree_list[[i]],
                                      castor::gamma_statistic,
                                      treestats::gamma_statistic, used_ntaxa[i], "gamma"))

  found <- rbind(found, compare_speed(tree_list[[i]],
                                      calc_sackin,
                                      treestats::sackin, used_ntaxa[i], "sackin"))

  found <- rbind(found, compare_speed(tree_list[[i]], calc_colless,  treestats::colless, used_ntaxa[i], "colless"))

  found <- rbind(found, compare_speed(tree_list[[i]], calc_beta,  treestats::beta_statistic, used_ntaxa[i], "beta"))


  found <- rbind(found,
                 compare_speed(tree_list[[i]], calc_blum,  treestats::blum, used_ntaxa[i], "blum"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], beautier::get_crown_age,  treestats::crown_age, used_ntaxa[i], "crown_age"))

  found <- rbind(found, compare_speed(tree_list[[i]], ref_tree_height, treestats::tree_height, used_ntaxa[i], "tree_height") )

  # pigot rho is not available previously
  #
  # phylogenetic diversity is not available previously
  found <- rbind(found, compare_speed(tree_list[[i]], nLTT_base_R,  treestats::nLTT_base, used_ntaxa[i], "nLTT"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], phyloTop::avgLadder,  treestats::avg_ladder, used_ntaxa[i], "ladder"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], calc_max_ladder,  treestats::max_ladder, used_ntaxa[i], "max_ladder"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], phyloTop::cherries,  treestats::cherries, used_ntaxa[i], "cherries"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], phyloTop::ILnumber,  treestats::ILnumber, used_ntaxa[i], "ILnumber"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], phyloTop::pitchforks,  treestats::pitchforks, used_ntaxa[i], "pitchforks"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], phyloTop::stairs,  treestats::stairs, used_ntaxa[i], "stairs1"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], RPANDA::spectR,  treestats::laplacian_spectrum, used_ntaxa[i], "laplacian_spectrum"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], calc_j_one,  treestats::j_one, used_ntaxa[i], "j_one"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treebalance::B1I,  treestats::b1, used_ntaxa[i], "B1"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treebalance::B2I,  treestats::b2, used_ntaxa[i], "B2"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treebalance::areaPerPairI,  treestats::area_per_pair, used_ntaxa[i], "area_per_pair"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treebalance::avgLeafDepI,  treestats::average_leaf_depth, used_ntaxa[i], "average_leaf_depth"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], calc_i,  treestats::mean_i, used_ntaxa[i], "mean_I"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treebalance::ewCollessI,  treestats::ew_colless, used_ntaxa[i], "ew_Colless"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treebalance::maxDelW,  treestats::max_del_width, used_ntaxa[i], "max_del_width"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treebalance::maxDepth,  treestats::max_depth, used_ntaxa[i], "max_depth"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treebalance::maxWidth,  treestats::max_width, used_ntaxa[i], "max_width"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treebalance::rogersI,  treestats::rogers, used_ntaxa[i], "rogers"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], phyloTop::stairs,  treestats::stairs2, used_ntaxa[i], "stairs2"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treebalance::totCophI,  treestats::tot_coph, used_ntaxa[i], "tot_coph"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treebalance::varLeafDepI,  treestats::var_leaf_depth, used_ntaxa[i], "var_leaf_depth"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treebalance::symNodesI,  treestats::sym_nodes, used_ntaxa[i], "sym_nodes"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], calc_mpd,  treestats::mean_pair_dist, used_ntaxa[i], "mpd"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], calc_psv,  treestats::psv, used_ntaxa[i], "psv"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], calc_var_mpd,  treestats::var_pair_dist, used_ntaxa[i], "vpd"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], calc_mntd,  treestats::mntd, used_ntaxa[i], "mntd"))

  # entropy J not available
  #
  found <- rbind(found,
                 compare_speed(tree_list[[i]], treebalance::rQuartetI,  treestats::rquartet, used_ntaxa[i], "rquartet"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treeCentrality::computeWienerIndex,  treestats::wiener, used_ntaxa[i], "wiener"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treeCentrality::computeBetweenness,  treestats::max_betweenness, used_ntaxa[i], "max_betweenness"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treeCentrality::computeCloseness,  treestats::max_closeness, used_ntaxa[i], "max_closeness"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treeCentrality::computeDiameter,  treestats::diameter, used_ntaxa[i], "diameter"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treeCentrality::computeEigenvector,  treestats::eigen_centrality, used_ntaxa[i], "eigen_centrality"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], calc_mean_br,  treestats::mean_branch_length, used_ntaxa[i], "mean_branch_length"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treebalance::avgVertDep, treestats::avg_vert_depth, used_ntaxa[i], "avg_vert_depth"))


  found <- rbind(found,
                 compare_speed(tree_list[[i]], treebalance::totPathLen, treestats::tot_path_length, used_ntaxa[i], "tot_path_length"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treebalance::totIntPathLen, treestats::tot_internal_path, used_ntaxa[i], "tot_internal_path"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], treebalance::mWovermD, treestats::mw_over_md, used_ntaxa[i], "mw_over_md"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], calc_colless_quad, treestats::colless_quad, used_ntaxa[i], "colless_quad"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], calc_colless_corr, treestats::colless_corr, used_ntaxa[i], "colless_corr"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], calc_minmax_adjacency, treestats:::minmax_adj, used_ntaxa[i], "minmax_adj"))

  found <- rbind(found,
                 compare_speed(tree_list[[i]], calc_minmax_laplace, treestats:::minmax_laplace, used_ntaxa[i], "minmax_laplace"))
}

colnames(found) <- c("ntips", "method", "time", "treestatsfunction")
found <- tibble::as_tibble(found)

found$method <- factor(found$method,
                       levels = c("treestats",
                                  "R equivalent"))

found$time <- as.numeric(found$time)
found$ntips <- as.numeric(found$ntips)

p1 <- found %>%
  ggplot(aes(x = ntips, y = time, col = method)) +
  geom_point() +
  stat_smooth() +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_viridis_d(option = "A", begin = 0.4, end = 0.7) +
  facet_wrap(~treestatsfunction, scales = "free", ncol = 6) +
  theme_classic() +
  theme(legend.position = "top") +
  ylab("Time per tree (seconds)") +
  xlab("Number of extant tips")

ggsave(p1, file = "Figure_S3.pdf", width = 18, height = 14)

p2 <- found %>%
  filter(ntips %in% c(1000)) %>%
  group_by(ntips, method, treestatsfunction) %>%
  summarise("mean_time" = mean(time)) %>%
  spread(key = method, value = mean_time) %>%
  mutate("speedup" = `R equivalent` / treestats) %>%
  ggplot(aes(x = stats::reorder(treestatsfunction, -speedup), y = speedup)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 0, lty = 2) +
  scale_y_log10() +
  xlab("Summary Statistic") +
  ylab("Fold speed increase for tree of 1000 tips") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.35))
p2

ggsave(p2, file = "Figure_S4.pdf", width = 8, height = 8)


found %>%
  filter(ntips %in% c(1000)) %>%
  filter(method == "treestats") %>%
  group_by(ntips, method, treestatsfunction) %>%
  summarise("mean_time" = mean(time)) %>%
  arrange(desc(mean_time))
