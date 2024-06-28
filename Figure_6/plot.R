get_all_cor <- function(cor1, label) {
  cor1 <- as.matrix(cor1)
  to_remove <- which(colnames(cor1) %in%
                       c("crown_age", "tree_height", "number_of_lineages"))
  if (length(to_remove)) cor1 <- cor1[-to_remove, -to_remove]

  diag(cor1) <- NA

  found <- c()
  for (i in 1:length(colnames(cor1))) {
    focal_label <- colnames(cor1)[i]
    all_corr <- abs(cor1[i, ])
    all_corr <- all_corr[!is.na(all_corr)]
    to_add <- cbind(focal_label, all_corr, label)
    found <- rbind(found, to_add)
  }

  return(found)
}

all_found <- c()

for (x in c("Mammals", "Birds", "Ferns", "Vascular Plants", "Cartaliginous Fish", "Ray finned Fish", "Amphibians")) {
  file_name <- paste0("../Figure_3/cor_emp_phy_", x, ".txt")
  cor2 <- read.table(file_name)
  local_names <- names(cor2)

  to_add <- get_all_cor(cor2, "Empirical")
  all_found <- rbind(all_found, to_add)
}


for (x in c("BD", "DDD", "PBD","SSE")) {
  file_name <- paste0("../Figure_5/", x, "_cor.txt")
  cor2 <- read.table(file_name)
  local_names <- names(cor2)
  rownames(cor2) <- local_names
  colnames(cor2) <- local_names
  to_add <- get_all_cor(cor2, "Simulation")
  all_found <- rbind(all_found, to_add)
}

colnames(all_found) <- c("statistic", "correlation", "method")
all_found <- as_tibble(all_found)
all_found$correlation <- as.numeric(all_found$correlation)


p1 <- all_found %>%
  filter(method == "Empirical") %>%
  mutate(class = fct_reorder(statistic, correlation, .fun='mean')) %>%
    ggplot( aes(x=class, y=correlation, fill=as.numeric(class))) +
      geom_boxplot() +
  scale_fill_viridis_c(direction = -1, begin = 0.4, end = 0.8, option = "A") +
      xlab("Summary Statistic") +
  ylab("Absolute Correlation") +
  theme_classic() +
      theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      xlab("") +
  ggtitle("Empirical trees")

p2 <- all_found %>%
  filter(method == "Simulation") %>%
  mutate(class = fct_reorder(statistic, correlation, .fun='median')) %>%
  ggplot( aes(x=class, y=correlation, fill=as.numeric(class))) +
  geom_boxplot() +
  scale_fill_viridis_c(direction = -1, begin = 0.4, end = 0.8, option = "A") +
  xlab("Summary Statistic") +
  ylab("Absolute Correlation") +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") +
  ggtitle("Simulated trees")

px <- egg::ggarrange(p1, p2, nrow = 1)

ggsave(px, filename = "Figure_6.pdf", width = 14, height = 8)


