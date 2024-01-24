require(tidyverse)
require(ggh4x)


found <- readRDS("../datasets/size_data.rds")


library(ggh4x)
norm_info <- read.table("../datasets/normalize.txt", header = TRUE)

found$normalized <- FALSE
found$type <- "None"
for (i in 1:length(norm_info$statistic)) {
  indices <- which(found$statistic == norm_info$statistic[i])
  found$normalized[indices] <- norm_info$normalize[i]
  found$type[indices] <- norm_info$type[i]
}

colz <- (c("#c5a6ff",
           "#d7ca3a",
           "grey"))

strip_colz <- c()
for (x in unique(found$statistic)) {
  ax <- subset(norm_info, norm_info$statistic == x)
  if (ax$type == "Yule") {
    strip_colz <- c(strip_colz, colz[1])
  }
  if (ax$type == "Tips") {
    strip_colz <- c(strip_colz, colz[2])
  }
  if (ax$type == "None") {
    strip_colz <- c(strip_colz, colz[3])
  }
}

strip <- strip_themed(background_x = elem_list_rect(fill = strip_colz))
found %>%
  ggplot(aes(x = number_of_lineages, y = mean_val, col = model)) +
  geom_point(size = 0.5) +
  geom_line(linewidth = 0.75) +
  scale_x_log10() +
  scale_color_brewer(type = "qual", palette = 3) +
  facet_wrap2(~statistic, scales = "free", strip = strip, ncol = 6) +
  theme_classic() +
  theme(legend.position = "top") +
  theme(legend.text = element_text(size = 13)) +
  xlab("Number of Extant Tips") +
  ylab("") +
  labs(col = "")

ggsave(filename = "Figure_4.pdf", width = 14, height = 18)
