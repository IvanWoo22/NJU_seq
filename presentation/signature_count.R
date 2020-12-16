#!/usr/bin/env Rscript
args <- commandArgs(T)

library(ggplot2)
library(reshape2)
library(dplyr)

data <- read.table(args[1],
                   header = F,
                   sep = "\t",
                   quote = "")
data$V3 <-
  c(
    rep(x = "A", times = 16),
    rep(x = "G", times = 16),
    rep(x = "C", times = 16),
    rep(x = "U", times = 16)
  )
data$V4 <- c(1:64)
color_scale <- c("#109648", "#255c99", "#f7b32b", "#d62839")
sumA <- sum(data[data$V3 == "A", ]$V2)
propA <- sum(data[data$V3 == "A", ]$V2) / sum(data$V2)
propA <- round(propA, digits = 4)
propA <- paste(propA * 100, "%", sep = '')
sumG <- sum(data[data$V3 == "G", ]$V2)
propG <- sum(data[data$V3 == "G", ]$V2) / sum(data$V2)
propG <- round(propG, digits = 4)
propG <- paste(propG * 100, "%", sep = '')
sumC <- sum(data[data$V3 == "C", ]$V2)
propC <- sum(data[data$V3 == "C", ]$V2) / sum(data$V2)
propC <- round(propC, digits = 4)
propC <- paste(propC * 100, "%", sep = '')
sumU <- sum(data[data$V3 == "U", ]$V2)
propU <- sum(data[data$V3 == "U", ]$V2) / sum(data$V2)
propU <- round(propU, digits = 4)
propU <- paste(propU * 100, "%", sep = '')

result_save_path <- args[2]

text_to_plot = data.frame(
  x = c(8.5, 24.5, 40.5, 56.5),
  y = c(
    max(data$V2) + 6,
    max(data$V2) + 6,
    max(data$V2) + 6,
    max(data$V2) + 6
  ),
  text = c(
    Reduce('paste0', c("Am\n", sumA, "\n", propA)),
    Reduce('paste0', c("Gm\n", sumG, "\n", propG)),
    Reduce('paste0', c("Cm\n", sumC, "\n", propC)),
    Reduce('paste0', c("Um\n", sumU, "\n", propU))
  )
)
line_to_plot1 = data.frame(x = c(0.65, 16.5), y = c(max(data$V2) + 11, max(data$V2) + 11))
line_to_plot2 = data.frame(x = c(16.5, 32.5), y = c(max(data$V2) + 11, max(data$V2) + 11))
line_to_plot3 = data.frame(x = c(32.5, 48.5), y = c(max(data$V2) + 11, max(data$V2) + 11))
line_to_plot4 = data.frame(x = c(48.5, 64.35), y = c(max(data$V2) + 11, max(data$V2) + 11))
p <- ggplot(data, aes(x = V4, y = V2))
p + geom_bar(
  stat = "identity",
  position = "dodge",
  aes(fill = V3),
  width = .7,
  show.legend = F,
  colour = 'black',
  size = 0.04669
) +
  ggplot2:::limits(data$V1, "x") +
  scale_y_continuous(expand = c(0, 1)) +
  scale_fill_manual(values = color_scale) +
  theme_bw(base_line_size = 0.11672, base_rect_size = 0.11672) +
  xlab("") +
  ylab("") +
  geom_text(
    aes(
      label = V2,
      y = V2 + 0.5,
      fontface = "plain"
    ),
    vjust = 0,
    size = 1,
    family = "ArialMT"
  ) +
  geom_line(
    data = line_to_plot1,
    aes(x = x, y = y),
    size = 0.701,
    color = '#109648',
    show.legend = F
  ) +
  geom_line(
    data = line_to_plot2,
    aes(x = x, y = y),
    size = 0.701,
    color = '#f7b32b',
    show.legend = F
  ) +
  geom_line(
    data = line_to_plot3,
    aes(x = x, y = y),
    size = 0.701,
    color = '#255c99',
    show.legend = F
  ) +
  geom_line(
    data = line_to_plot4,
    aes(x = x, y = y),
    size = 0.701,
    color = '#d62839',
    show.legend = F
  ) +
  geom_text(
    data = text_to_plot,
    aes(
      x = x,
      y = y,
      label = text,
      fontface = "bold",
      vjust = 0.5
    ),
    size = 5 / 3,
    family = "ArialMT",
    show.legend = F
  ) +
  theme(
    plot.margin = unit(c(0.4, 0.5, -0.5, -0.5), 'lines'),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      size = 4,
      face = "bold",
      colour = "#000000"
    ),
    axis.text.y = element_text(
      face = "plain",
      colour = "#000000",
      size = 5
    ),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "#000000", size = 0.11672),
    axis.ticks.length = unit(0.4, 'mm'),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(family = "ArialMT")
  )

ggsave(
  result_save_path,
  device = "pdf",
  width = 3.15,
  height = 1.1,
  dpi = "retina"
)
dev.off()