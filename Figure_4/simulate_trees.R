source("simulation_functions.R")

param_grid <- expand.grid(model = c("BD", "DDD", "PBD", "SSE"),
                          extinction = c(FALSE, TRUE),
                          repl = 1,
                          num_lin = seq(1, 3, length.out = 20))


args = commandArgs(trailingOnly = TRUE)
sim_number = as.numeric(args[[1]])

repl <- param_grid$repl[sim_number]

used_seed <- sim_number * 103 + 1e5 + 13 * 7 + sim_number * 13289
used_seed <- used_seed + as.numeric(Sys.time())

set.seed(used_seed)

num_repl <- 10000

model <- param_grid$model[sim_number]
local_use_extinct <- param_grid$extinction[sim_number]

cat(model, "\n")
cat(local_use_extinct, "\n")

num_lin <- floor(10^param_grid$num_lin[sim_number])

div_rate <- (spec_rate + extinct_rate * local_use_extinct - extinct_rate * local_use_extinct)

max_t <- (1 / div_rate) * log(num_lin / 2)

file_name_trees <- paste0(model, "_", local_use_extinct, "_", sim_number, "_trees.trees")
file_name_stats <- paste0(model, "_", local_use_extinct, "_", sim_number, "_stats.txt")

all_trees <- list()

start_repl <- 1

is_continued <- FALSE
if (is_continued) {
  if (file.exists(file_name_stats)) {
    vx <- readr::read_delim(file_name_stats, delim = " ", col_names = FALSE)
    start_repl <- 1 + length(vx$X1)
    all_trees <- ape::read.tree(file_name_trees)
    set.seed(used_seed)
  }
}

if (start_repl < num_repl) {
  for (r in start_repl:num_repl) {
    focal_tree <- simulate_tree(model, max_t, num_lin, local_use_extinct)

    focal_stats <- treestats::calc_all_stats(focal_tree, normalize = TRUE)
    to_add <- c(as.numeric(focal_stats))

    cat(to_add, "\n", file = file_name_stats, append = TRUE)

    if (r %% 10 == 0) cat(r, "\n")
  }
}
