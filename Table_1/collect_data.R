
all_found <- c()

for (x in c("Mammal", "Bird", "Amphibia","Squamate")) {
  file_name <- paste0("../Figure_3/cor_emp_", x, ".txt")
  cor2 <- read.table(file_name)
  local_names <- names(cor2)
  rownames(cor2) <- local_names
  colnames(cor2) <- local_names

  cor2[lower.tri(cor2)] <- NA
  cor3 <- as.data.frame(cor2)
  cor3$rowname <- rownames(cor3)
  cor4 <- cor3 %>%
    gather(key = "statistic", value = "val", -rowname)
  cor4 <- cor4 %>% filter(!is.na(val))
  all_found <- rbind(all_found, cor4)
}

for (x in c("BD", "DDD", "PBD","SSE")) {
  file_name <- paste0("../Figure_5/", x, "_cor.txt")
  cor2 <- read.table(file_name)
  local_names <- names(cor2)
  rownames(cor2) <- local_names
  colnames(cor2) <- local_names

  cor2[lower.tri(cor2)] <- NA
  cor3 <- as.data.frame(cor2)
  cor3$rowname <- rownames(cor3)
  cor4 <- cor3 %>%
    gather(key = "statistic", value = "val", -rowname)
  cor4 <- cor4 %>% filter(!is.na(val))
  all_found <- rbind(all_found, cor4)
}

mean_vals <- all_found %>%
  group_by(rowname, statistic) %>%
  summarise("mean" = mean(abs(val)))

sd_vals <-  all_found %>%
  group_by(rowname, statistic) %>%
  summarise("sd" = sd(abs(val)),
            "mean" = mean(abs(val))) %>%
  arrange(sd)

sd_vals2 <- sd_vals %>%
  filter(sd < 0.05) %>%
  filter(rowname != statistic) %>%
  filter(rowname != "number_of_lineages") %>%
  filter(statistic != "number_of_lineages") %>%
  filter(mean > 0.8)

write.table(sd_vals2, file = "small_sd_vals.txt")

# now we only collect simulation data:
# (not the most efficient way, but everything is fast anyway)
all_found <- c()
for (x in c("BD", "DDD", "PBD","SSE")) {
  file_name <- paste0("../Figure_5/", x, "_cor.txt")

  cor2 <- read.table(file_name)
  local_names <- names(cor2)
  rownames(cor2) <- local_names
  colnames(cor2) <- local_names

  cor2[lower.tri(cor2)] <- NA
  cor3 <- as.data.frame(cor2)
  cor3$rowname <- rownames(cor3)
  cor4 <- cor3 %>%
    gather(key = "statistic", value = "val", -rowname)
  cor4 <- cor4 %>% filter(!is.na(val))
  all_found <- rbind(all_found, cor4)
}

sd_vals <-  all_found %>%
  group_by(rowname, statistic) %>%
  summarise("sd" = sd(abs(val)),
            "mean" = mean(abs(val))) %>%
  arrange(sd)

sd_vals2 <- sd_vals %>%
  filter(sd < 0.05) %>%
  filter(rowname != statistic) %>%
  filter(rowname != "number_of_lineages") %>%
  filter(statistic != "number_of_lineages") %>%
  filter(mean > 0.8)

write.table(sd_vals2, file = "small_sd_vals_sim_only.txt")
