#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(forcats)
library(grid)
library(gridExtra)

args <- commandArgs(T)
file_path_1 <- args[1]
file_path_2 <- args[2]
result_save_path <- args[3]

raw_data_1 <-
  read.table(
    file_path_1,
    header = F,
    stringsAsFactors = F,
    sep = "\t"
  )

raw_data_2 <-
  read.table(
    file_path_2,
    header = F,
    stringsAsFactors = F,
    sep = "\t"
  )

colnames(raw_data_1) <-
  c("name1", "name2", "color_pattern", "score")
colnames(raw_data_2) <-
  c("name1", "name2", "color_pattern", "score")

df_plot1 <- raw_data_1 %>%
  mutate(color_pattern = as.factor(color_pattern)) %>%
  mutate(name = paste0(name1, name2)) %>%
  mutate(log_score = log10(score)) %>%
  select(name,
         color_pattern,
         log_score)

df_plot2 <- raw_data_2 %>%
  mutate(color_pattern = as.factor(color_pattern)) %>%
  mutate(name = paste0(name1, name2)) %>%
  mutate(log_score = log10(score)) %>%
  select(name,
         color_pattern,
         log_score)

color_scale <- c("#A2CA86", "#C8C8C8")
names(color_scale) <- c("1", "0")

color_scale_point <- c("#A2CA86", "#C8C8C8")
names(color_scale_point) <- c("1", "0")

p1 <-
  ggplot(mutate(df_plot1, name = fct_reorder(name, log_score))) +
  geom_bar(
    aes(x = name, y = log_score, fill = color_pattern),
    stat = "identity",
    alpha = .8,
    width = .6,
    show.legend = F
  ) +
  scale_fill_manual(values = color_scale) +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip(ylim = c(0, log10(max(raw_data[, 4]) * 1.5))) +
  theme_bw(base_line_size = 0.11672, base_rect_size = 0.11672) +
  theme(
    axis.ticks = element_line(colour = "#000000", size = 0.11672),
    axis.ticks.length = unit(0.4, 'mm'),
    axis.line = element_blank(),
    axis.text.y = element_text(
      face = "bold",
      colour = "#000000",
      size = 8
    ),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
  )

p2 <-
  ggplot(mutate(df_plot2, name = fct_reorder(name, log_score))) +
  geom_bar(
    aes(x = name, y = log_score, fill = color_pattern),
    stat = "identity",
    alpha = .8,
    width = .6,
    show.legend = F
  ) +
  scale_fill_manual(values = color_scale) +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip(ylim = c(0, log10(max(raw_data[, 4]) * 1.5))) +
  theme_bw(base_line_size = 0.11672, base_rect_size = 0.11672) +
  theme(
    axis.ticks = element_line(colour = "#000000", size = 0.11672),
    axis.ticks.length = unit(0.4, 'mm'),
    axis.line = element_blank(),
    axis.text.y = element_text(
      face = "bold",
      colour = "#000000",
      size = 8
    ),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
  )

rownum <- max(nrow(df_plot1), nrow(df_plot2))
height <- 0.525 + rownum * 0.169
pdf(
  file = result_save_path,
  width = 11,
  height = height,
  useDingbats = FALSE
)
grid.arrange(p1,
             p2,
             ncol = 2,
             nrow = 1)
dev.off()