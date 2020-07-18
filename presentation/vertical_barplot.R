#!/usr/bin/env Rscript
args <- commandArgs(T)
file_path <- args[1]
result_save_path <- args[2]

library(dplyr)
library(ggplot2)
library(forcats)

raw_data <- read.table(file_path, header = F, stringsAsFactors = F, sep = "\t")
colnames(raw_data) <- c("name1", "name2", "color_pattern", "num1", "num2", "num3")
df_plot <- raw_data %>%
  mutate(color_pattern = as.factor(color_pattern)) %>%
  mutate(name = paste0(name1, name2), average = apply(raw_data[, 4:6], 1, mean), median = apply(raw_data[, 4:6], 1, median)) %>%
  mutate(log_ave = log10(average), log_num1 = log10(num1), log_num2 = log10(num2), log_num3 = log10(num3), log_median = log10(median)) %>%
  select(name, color_pattern, log_num1, log_num2, log_num3, average, median, log_ave, log_median)

color_scale <- c("#fff173", "#a7d1b8", "#cccccc")
names(color_scale) <- c("2", "1", "0")

color_scale_point <- c("#fff173","#a7d1b8", "#cccccc")
names(color_scale_point) <- c("2", "1", "0")
df_plot %>%
  mutate(name = fct_reorder(name, log_ave)) %>%
  ggplot() +
  geom_bar(aes(x = name, y = log_ave, fill= color_pattern), stat="identity", alpha=.6, width = .65, show.legend = F) +
  scale_fill_manual(values = color_scale) +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_point(aes(x = name, y = log_num1), shape = 16, color = "#489DCC", alpha = .9, size = 1.2) +
  geom_point(aes(x = name, y = log_num2), shape = 16, color = "#C6175A", alpha = .9, size = 1.2) +
  geom_point(aes(x = name, y = log_num3), shape = 16, color = "#EF7E5B", alpha = .9, size = 1.2) +
  coord_flip(ylim = c(0, log10(max(raw_data[,4:6])*1.1)) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

nrow <- nrow(df_plot)
height <- 0.525 + nrow * 0.169
ggsave(paste0(result_save_path, ".pdf"), device = NULL, width = 6, height = height, dpi = "retina")
