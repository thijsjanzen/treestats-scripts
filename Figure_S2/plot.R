require(readr)
require(tidyverse)
vx <- readRDS("../datasets/through_time.rds")

num_categories <- length(unique(vx$relative_time_2))

vy <- vx %>%
  group_by(model, relative_time_2, statistic) %>%
  summarise("mean_val" = mean(mean_val, na.rm = TRUE))



p1 <- vy %>%
  ggplot(aes(x = relative_time_2, y = mean_val, col = model)) +
  geom_line(linewidth = 0.75) +
  #geom_point() +
  xlim(1, 0) +
  facet_wrap(~statistic, scales = "free", ncol = 6) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(type = "qual", palette = 3) +
  xlab("Time relative to the crown age") +
  ylab("") +
  labs(col = "")

p1

ggsave(p1, file = "Figure_S2.pdf", width = 14, height = 18)

