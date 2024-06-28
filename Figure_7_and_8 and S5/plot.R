require(tidyverse)
found <- readRDS("../datasets/balance_data.rds")


found2 <- found %>%
  gather(key = "statistic", value = "val", -c(repl, n, method1, method2, model, extinction))

found3 <- found2 %>%
  filter(extinction == TRUE)

found4 <- found3 %>%
  group_by(n, method1, method2, model, statistic) %>%
  summarise("mean_val" = mean(val))

found4$combined_method <- NA
for (i in 1:length(found4$n)) {
  a1 <- found4$method1[i]
  a2 <- found4$method2[i]
  if (a1 == "any") {
    if (a2 == "random") found4$combined_method[i] <- "A-R"
    if (a2 == "youngest") found4$combined_method[i] <- "A-Y"
    if (a2 == "oldest") found4$combined_method[i] <- "A-O"
  }
  if (a1 == "terminal") {
    if (a2 == "random") found4$combined_method[i] <- "T-R"
    if (a2 == "youngest") found4$combined_method[i] <- "T-Y"
    if (a2 == "oldest") found4$combined_method[i] <- "T-O"
  }
}

found4$combined_method <- factor(found4$combined_method,
                                 levels = c("A-R", "A-O", "A-Y", "T-R", "T-O", "T-Y"))


colz <- c("#ffbd49", #Imbalance
          "#4da6c5", # No index
          "#234a79", # Shape
          "#ec632b") # Balance

tip_info <- read_tsv("../datasets/families.txt")

strip_colz <- c()
to_remove <- c()
for (x in unique(found4$statistic)) {
  ax <- subset(tip_info, tip_info$Statistic == x)
  if (ax$Balance_Fischer == "Balance") {
    strip_colz <- c(strip_colz, colz[1])
  }
  if (ax$Balance_Fischer == "Imbalance") {
    strip_colz <- c(strip_colz, colz[2])
  }
  if (ax$Balance_Fischer == "Shape") {
    strip_colz <- c(strip_colz, colz[3])
  }
  if (ax$Balance_Fischer == "No index") {
    # strip_colz <- c(strip_colz, colz[4])
    to_remove <- c(to_remove, x)
  }
}

to_stay <- unique(found4$statistic)
indices <- to_stay %in% to_remove
to_stay <- to_stay[indices]

found5 <- found4 %>%
  filter(combined_method != "A-Y") %>%
  filter(model == 1) %>%
  filter(statistic %in% to_stay)


res <- c()
for (focal_stat in unique(found5$statistic)) {
  focal_data <- subset(found5, found5$statistic == focal_stat &
                               found5$model == 1)

  for (m in unique(focal_data$combined_method)) {
    local_data <- subset(focal_data, focal_data$combined_method == m)
    y <- local_data$mean_val[local_data$n]
    ll <- loess(local_data$mean_val~local_data$n)
    y2 <- predict(ll)


    z <- round(y, 2)
    if (m == "A-R") z <- round(y2, 4)
    if (m == "T-R") z <- round(y2, 4)

    a1 <- min(which.min(z))
    a2 <- max(which.max(z))
    check_min <- TRUE
    check_max <- TRUE
    check_monotonic <- TRUE

    if (a1 > a2) {
      # imbalance
      if (a1 != 120) check_min <- FALSE
      if (a2 != 1) check_max <- FALSE

      check_monotonic <- all(z == cummin(z))
    } else {
      # balance
      if (a1 != 1) check_max <- FALSE
      if (a2 != 120) check_min <- FALSE

      check_monotonic <- all(z == cummax(z))
    }
    to_add <- c(focal_stat, m, check_min, check_max, check_monotonic, a2 > a1)
    cat(to_add, "\n")
    res <- rbind(res, to_add)
  }
}

colz <- c("#ffbd49", #Imbalance
          "#4da6c5", # No index
          "#234a79", # Shape
          "#ec632b") # Balance

strip_colz <- c()
for (x in unique(found5$statistic)) {
  local_res <- subset(res, res[, 1] == x)
  vv <- c()
  for (i in 3:5) {
    vv <- c(vv, sum(local_res[, i] == "TRUE"))
  }
  if (sum(vv) == 15) {
    if (local_res[1, 6]) {
      strip_colz <- c(strip_colz, colz[1])
    } else {
      strip_colz <- c(strip_colz, colz[4])
    }
  } else {
    strip_colz <- c(strip_colz, colz[2])
  }
}

strip <- strip_themed(background_x = elem_list_rect(fill = strip_colz))

to_plot <- found5 %>%
  ggplot(aes(x = n, y = mean_val, col = combined_method)) +
  geom_line(lwd = 1.3) +
  facet_wrap2(~statistic, strip = strip, scales = "free", ncol = 7) +
  theme_classic() +
  theme(strip.text = element_text(colour = 'white')) +
  scale_color_manual(values = c("#4f4dbc", "#8c8cff", "#b54440", "#ff7837", "#fbe1c2")) +
  theme(text = element_text(family = "Helvetica")) +
  ylab("Statistic value") +
  xlab("Number of steps taken towards fully imbalanced tree") +
  labs(col = "Algorithm")
to_plot

ggsave(to_plot, filename = "Figure_8.pdf", width = 20, height = 15)



found6 <- found4 %>%
  filter(combined_method != "A-Y") %>%
  filter(statistic %in% to_stay)

res <- c()
for (focal_stat in unique(found6$statistic)) {
  for (focal_model in unique(found6$model)) {
    focal_data <- subset(found6, found6$statistic == focal_stat &
                               found6$model == focal_model)

    for (m in unique(focal_data$combined_method)) {
      local_data <- subset(focal_data, focal_data$combined_method == m)
      y <- local_data$mean_val[local_data$n]
      ll <- loess(local_data$mean_val~local_data$n)
      y2 <- predict(ll)
      z <- round(y, 4)

      a1 <- min(which.min(z))
      a2 <- max(which.max(z))
      check_min <- TRUE
      check_max <- TRUE
      check_monotonic <- TRUE

      if (m == "A-R") z <- round(y2, 4)
      if (m == "T-R") z <- round(y2, 4)

      if (a1 > a2) {
        # imbalance
        if (a1 != 120) check_min <- FALSE
        if (a2 != 1) check_max <- FALSE

        check_monotonic <- all(z == cummin(z))
      } else {
        # balance
        if (a1 != 1) check_max <- FALSE
        if (a2 != 120) check_min <- FALSE

        check_monotonic <- all(z == cummax(z))
      }
      to_add <- c(focal_stat, focal_model, m, check_min, check_max, check_monotonic, a2 > a1)
      cat(to_add, "\n")
      res <- rbind(res, to_add)
    }
  }
}




strip_colz <- c()
for (x in unique(found6$statistic)) {
  local_res <- subset(res, res[, 1] == x)
  vv <- c()
  for (i in 4:6) {
    vv <- c(vv, sum(local_res[, i] == "TRUE"))
  }
  if (sum(vv) == 60) {
    if (local_res[1, 7]) {
      strip_colz <- c(strip_colz, colz[1])
    } else {
      strip_colz <- c(strip_colz, colz[4])
    }
  } else {
    strip_colz <- c(strip_colz, colz[2])
  }
}

strip <- strip_themed(background_x = elem_list_rect(fill = strip_colz))

to_plot <- found5 %>%
  ggplot(aes(x = n, y = mean_val, col = combined_method)) +
  geom_line(lwd = 1.3) +
  facet_wrap2(~statistic, strip = strip, scales = "free", ncol = 7) +
  theme_classic() +
  theme(strip.text = element_text(colour = 'white')) +
  scale_color_manual(values = c("#4f4dbc", "#8c8cff", "#b54440", "#ff7837", "#fbe1c2")) +
  theme(text = element_text(family = "Helvetica")) +
  ylab("Statistic value") +
  xlab("Number of steps taken towards fully imbalanced tree") +
  labs(col = "Algorithm")
to_plot

found6$model[found6$model == 1] <- "BD"
found6$model[found6$model == 2] <- "DDD"
found6$model[found6$model == 3] <- "PBD"
found6$model[found6$model == 4] <- "SSE"

to_plot <- found6 %>%
  ggplot(aes(x = n, y = mean_val, col = combined_method)) +
  geom_line(lwd = 1.3) +
  facet_grid(rows = vars(statistic), cols = vars(model), scales = "free") +
  theme_minimal() +
  #theme(strip.text = element_text(colour = 'white')) +
  scale_color_manual(values = c("#4f4dbc", "#8c8cff", "#b54440", "#ff7837", "#fbe1c2")) +
  theme(text = element_text(family = "Helvetica")) +
  ylab("Statistic value") +
  xlab("Number of steps taken towards fully imbalanced tree") +
  labs(col = "Algorithm")

to_plot
ggsave(to_plot, filename = "Figure_S8.pdf", width = 10, height = 35)
