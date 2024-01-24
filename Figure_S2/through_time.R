collect_stats <- function(focal_tree, num_steps = 100) {
  ca <- treestats::crown_age(focal_tree)
  time_points <- seq(0.01, 0.99, by = 0.01) * ca
  ltab <- treestats::phylo_to_l(focal_tree)

  found <- c()
  for (t in time_points) {
    pruned_ltab <- ltab[which(ltab[, 1] >= t), ]
    pruned_ltab[, 1] <- pruned_ltab[, 1] - t

    stats <- treestats::calc_all_stats(pruned_ltab, normalize = TRUE)

    if (sum(is.na(stats))) {
    } else {
      found <- rbind(found, c(t / ca, as.vector(unlist(stats))))
    }
  }
  colnames(found) <- c("relative_time", names(stats))
  found <- as_tibble(found)
  return(found)
}

require(tidyverse)

source("simulation_functions.R")

param_grid <- expand.grid(model = c("BD", "DDD", "PBD", "SSE"),
                          extinction = c(TRUE),
                          repl = 1:100)

args = commandArgs(trailingOnly = TRUE)
sim_number = as.numeric(args[[1]])

repl <- param_grid$repl[sim_number]

set.seed(sim_number * 104 + 1e2 + 13 * 7 + repl * 13)

num_repl <- 1000

model <- param_grid$model[sim_number]
local_use_extinct <- param_grid$extinction[sim_number]

cat(model, "\n")
cat(local_use_extinct, "\n")

max_t <- 10
num_lin <- 300

file_name_stats <- paste0(model, "_", local_use_extinct, "_", repl, "_stats.txt")

start_repl <- 1
if (file.exists(file_name_stats)) {
  temp_file <- read_tsv(file_name_stats)
  replz <- as.numeric(temp_file$repl)
  start_repl <- max(replz) + 1
}

if (start_repl < num_repl) {
  for (r in start_repl:num_repl) {
    focal_tree <- simulate_tree(model, max_t, num_lin, local_use_extinct)

    focal_stats <- collect_stats(focal_tree)
    focal_stats$repl <- r

    write_tsv(focal_stats, "\n", file = file_name_stats, append = TRUE,
              col_names = (r == 1))

    if (r %% 10 == 0) cat(r, "\n")
  }
}
