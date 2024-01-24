source("simulation_functions.R")

param_grid <- expand.grid(model = c("BD", "DDD", "PBD", "SSE"),
                          extinction = c(FALSE, TRUE),
                          repl = 1:100)

args = commandArgs(trailingOnly = TRUE)
sim_number = as.numeric(args[[1]])

repl <- param_grid$repl[sim_number]

set.seed(sim_number * 104 + 1e2 + 13 * 7 + repl * 138311)

num_repl <- 10000

model <- param_grid$model[sim_number]
local_use_extinct <- param_grid$extinction[sim_number]

cat(model, "\n")
cat(local_use_extinct, "\n")

max_t <- 10
num_lin <- 300

file_name_trees <- paste0(model, "_", local_use_extinct, "_", repl, "_trees.trees")
file_name_stats <- paste0(model, "_", local_use_extinct, "_", repl, "_stats.txt")

all_trees <- list()

start_repl <- 1

t_prev <- Sys.time()

if (start_repl < num_repl) {
  for (r in start_repl:num_repl) {
    focal_tree <- simulate_tree(model, max_t, num_lin, local_use_extinct)

    focal_stats <- treestats::calc_all_stats(focal_tree, normalize = FALSE)
    to_add <- c(as.numeric(focal_stats))

    cat(to_add, "\n", file = file_name_stats, append = TRUE)

    if (r %% 100 == 0) cat(r, "\n")
  }
}
