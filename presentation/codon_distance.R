#!/usr/bin/env Rscript
library(ggplot2)
library(splines)
library(gridExtra)
library(tidyquant)

args <- commandArgs(T)

start <- read.table(args[1], header = F)
stop <- read.table(args[2], header = F)
ymax = max(start$V2, stop$V2)

p1 <- ggplot() +
  geom_ma(aes(x = start$V1, y = start$V2), ma_fun = SMA, n = 3) +
  xlab("Distance to Start Codon (nt)") +
  ylab("Nmuber of Nm Sites") +
  geom_line(
    aes(x = start$V1, y = start$V2),
    size = 0.4669,
    show.legend = F,
    colour = "#A2CA86",
    alpha = 0.8
  ) +
  geom_line(aes(x = c(0, 0), y = c(0, ymax)),
            size = 0.23343,
            colour = "#751C2E",
            linetype = 2) +
  geom_text(
    aes(
      x = 0,
      y = ymax,
      label = "Start",
      vjust = 0
    ),
    family = "ArialMT",
    size = 5 / 3,
    show.legend = F,
    colour = "#751C2E"
  ) +
  scale_x_continuous(
    limits = c(-300, 500),
    breaks = seq(-300, 500, 100),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, ceiling(ymax / 5) * 5),
    breaks = seq(0, ceiling(ymax / 5) * 5, 10),
    expand = c(0, 0)
  ) +
  theme(
    axis.text.x = element_text(
      face = "plain",
      colour = "#000000",
      size = 5
    ),
    axis.text.y = element_text(
      face = "plain",
      colour = "#000000",
      size = 5
    ),
    axis.title.x = element_text(
      face = "plain",
      colour = "#000000",
      size = 7
    ),
    axis.title.y = element_text(
      face = "plain",
      colour = "#000000",
      size = 7
    ),
    axis.line = element_line(colour = "#000000", size = 0.11672),
    axis.ticks = element_line(colour = "#000000", size = 0.11672),
    axis.ticks.length = unit(0.4, 'mm'),
    plot.margin = unit(c(0.4, 0.4, 0.2, 0.2), 'lines'),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(family = "ArialMT")
  )

p2 <- ggplot() +
  geom_ma(aes(x = stop$V1, y = stop$V2), ma_fun = SMA, n = 3) +
  xlab("Distance to Stop Codon (nt)") +
  ylab("Nmuber of Nm Sites") +
  geom_line(
    aes(x = stop$V1, y = stop$V2),
    size = 0.4669,
    show.legend = F,
    colour = "#A2CA86",
    alpha = 0.8
  ) +
  geom_line(aes(x = c(0, 0), y = c(0, ymax)),
            size = 0.23343,
            colour = "#751C2E",
            linetype = 2) +
  geom_text(
    aes(
      x = 0,
      y = ymax,
      label = "Stop",
      vjust = 0
    ),
    family = "ArialMT",
    size = 5 / 3,
    show.legend = F,
    colour = "#751C2E"
  ) +
  scale_x_continuous(
    limits = c(-500, 500),
    breaks = seq(-500, 500, 100),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, ceiling(ymax / 5) * 5),
    breaks = seq(0, ceiling(ymax / 5) * 5, 10),
    expand = c(0, 0)
  ) +
  theme(
    axis.text.x = element_text(
      face = "plain",
      colour = "#000000",
      size = 5
    ),
    axis.text.y = element_text(
      face = "plain",
      colour = "#000000",
      size = 5
    ),
    axis.title.x = element_text(
      face = "plain",
      colour = "#000000",
      size = 7
    ),
    axis.title.y = element_text(
      face = "plain",
      colour = "#000000",
      size = 7
    ),
    axis.line = element_line(colour = "#000000", size = 0.11672),
    axis.ticks = element_line(colour = "#000000", size = 0.11672),
    axis.ticks.length = unit(0.4, 'mm'),
    plot.margin = unit(c(0.4, 0.4, 0.2, 0.2), 'lines'),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(family = "ArialMT")
  )

p <- grid.arrange(p1, p2, nrow = 1, widths = c(0.5, 0.5))
result_save_path <- args[3]
ggsave(
  result_save_path,
  plot = p,
  device = "pdf",
  width = 4.9,
  height = 3.2,
  dpi = "retina"
)