get_best_cor <- function(X, Y) {

  all_r <- c()

  all_r[1] <- cor(X, Y)

  xlog <- FALSE #min(X) > 0
  ylog <- FALSE # min(Y) > 0

  if (xlog) {
    all_r <-  c(all_r, cor(log(X), Y))
    if(ylog) {
      all_r <- c(all_r, cor(log(X), log(Y)))
    }
  }
  if (ylog) {
    all_r <- c(all_r, cor(X, log(Y)))
  }

  winner <- which.max(abs(all_r))

  return(all_r[winner])
}

require(tidyverse)
require(ggnewscale)

found <- readRDS("trees.rds")


for (focal_model in c("BD", "DDD", "SSE", "PBD")) {
  focal_data <- subset(found, found$model == focal_model)

  focal_data <- focal_data %>%
    select(-c("crown_age", "tree_height", "number_of_lineages",
              "model", "extinction", "combined"))
  focal_data <- as_tibble(focal_data)
  focal_data <- focal_data %>%
    mutate_at(1:51, as.numeric)

  stat_names <- colnames(focal_data)

  get_per_row <- function(stat1) {
    x <- unlist(as.vector(focal_data[stat1]))

    get_cor <- function(stat2) {
      if(stat2 != "number_of_lineages" && stat1 != stat2) {
        y <- unlist(as.vector(focal_data[stat2]))
        return(get_best_cor(x, y))
      }
      return(NA)
    }

    local_cor <- sapply(stat_names, get_cor)
    return(local_cor)
  }

  res <- pbapply::pblapply(stat_names, get_per_row)
  res2 <- do.call("rbind", res)

  file_name <- paste0(focal_model, "_cor.txt")

  write.table(res2, file = file_name, quote = FALSE)

  cat(focal_model, "\n")
}

focal_data <- found

focal_data <- focal_data %>%
  select(-c("crown_age", "tree_height", "number_of_lineages",
            "model", "extinction", "combined"))
focal_data <- as_tibble(focal_data)
focal_data <- focal_data %>%
  mutate_at(1:51, as.numeric)

stat_names <- colnames(focal_data)

get_per_row <- function(stat1) {
  x <- unlist(as.vector(focal_data[stat1]))

  get_cor <- function(stat2) {
    if(stat2 != "number_of_lineages" && stat1 != stat2) {
      y <- unlist(as.vector(focal_data[stat2]))
      return(get_best_cor(x, y))
    }
    return(NA)
  }

  local_cor <- sapply(stat_names, get_cor)
  return(local_cor)
}

res <- pbapply::pblapply(stat_names, get_per_row)
res2 <- do.call("rbind", res)

file_name <- paste0("overall_cor.txt")

write.table(res2, file = file_name, quote = FALSE)
