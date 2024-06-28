require(tidyverse)
require(ggnewscale)

found <- readRDS("/Users/thijsjanzen/trees_small.rds")

for (focal_model in c("BD", "DDD", "SSE", "PBD")) {
  focal_data <- subset(found, found$model == focal_model)

  focal_data <- focal_data %>%
    select(-c("crown_age", "tree_height", "number_of_lineages",
              "model"))
  focal_data <- as_tibble(focal_data)
  focal_data <- focal_data %>%
    mutate_at(1:67, as.numeric)

  stat_names <- colnames(focal_data)

  get_per_row <- function(stat1) {
    x <- unlist(as.vector(focal_data[stat1]))

    get_cor <- function(stat2) {
      if (stat2 != "number_of_lineages" && stat1 != stat2) {
        y <- unlist(as.vector(focal_data[stat2]))
        return(cor(x, y))
      }
      return(NA)
    }

    local_cor <- sapply(stat_names, get_cor)
    return(local_cor)
  }

  res <- pbmcapply::pbmclapply(stat_names, get_per_row, mc.cores = 6)

  res2 <- do.call("rbind", res)

  file_name <- paste0(focal_model, "_cor.txt")

  write.table(res2, file = file_name, quote = FALSE)

  cat(focal_model, "\n")
}

focal_data <- found

focal_data <- focal_data %>%
  select(-c("crown_age", "tree_height", "number_of_lineages",
            "model"))
focal_data <- as_tibble(focal_data)
focal_data <- focal_data %>%
  mutate_at(1:67, as.numeric)

stat_names <- colnames(focal_data)

get_per_row <- function(stat1) {
  x <- unlist(as.vector(focal_data[stat1]))

  get_cor <- function(stat2) {
    if(stat2 != "number_of_lineages" && stat1 != stat2) {
      y <- unlist(as.vector(focal_data[stat2]))
      return(cor(x, y))
    }
    return(NA)
  }

  local_cor <- sapply(stat_names, get_cor)
  return(local_cor)
}

res <- pbmcapply::pbmclapply(stat_names, get_per_row, mc.cores = 8)

res2 <- do.call("rbind", res)


file_name <- paste0("overall_cor.txt")

write.table(res2, file = file_name, quote = FALSE)
