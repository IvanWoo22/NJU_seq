#!/usr/bin/env Rscript
library(ggplot2)
library(splines)
library(gridExtra)

args <- commandArgs(T)

start1 <- read.table(args[1], header = F)
start2 <- read.table(args[2], header = F)
start3 <- read.table(args[3], header = F)
start <- rbind(start1, start2, start3)
stop1 <- read.table(args[4], header = F)
stop2 <- read.table(args[5], header = F)
stop3 <- read.table(args[6], header = F)
stop <- rbind(stop1, stop2, stop3)
ymax = max(start$V2, stop$V2)


p1 <- ggplot() +
  xlab("Distance to Start Codon (nt)") +
  ylab("Nmuber of Nm Sites") +
  geom_line(
    aes(x = start1$V1, y = start1$V2),
    size = 0.4669,
    show.legend = F,
    colour = "#A2CA86",
    alpha = 0.8
  ) +
  geom_line(
    aes(x = start2$V1, y = start2$V2),
    size = 0.4669,
    show.legend = F,
    colour = "#7664AD",
    alpha = 0.8
  ) +
  geom_line(
    aes(x = start3$V1, y = start3$V2),
    size = 0.4669,
    show.legend = F,
    colour = "#7995BD",
    alpha = 0.8
  ) +
  geom_smooth(
    aes(x = start$V1, y = start$V2),
    method = "glm",
    formula = y ~ ns(x, 8),
    size = 0.701,
    show.legend = F,
    colour = "#FF335C",
    alpha = 0.4
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
  xlab("Distance to Stop Codon (nt)") +
  ylab("Nmuber of Nm Sites") +
  geom_line(
    aes(x = stop1$V1, y = stop1$V2),
    size = 0.4669,
    show.legend = F,
    colour = "#A2CA86",
    alpha = 0.8
  ) +
  geom_line(
    aes(x = stop2$V1, y = stop2$V2),
    size = 0.4669,
    show.legend = F,
    colour = "#7664AD",
    alpha = 0.8
  ) +
  geom_line(
    aes(x = stop3$V1, y = stop3$V2),
    size = 0.4669,
    show.legend = F,
    colour = "#7995BD",
    alpha = 0.8
  ) +
  geom_smooth(
    aes(x = stop$V1, y = stop$V2),
    method = "glm",
    formula = y ~ ns(x, 8),
    size = 0.701,
    show.legend = F,
    colour = "#FF335C",
    alpha = 0.4
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

p <- grid.arrange(p1, p2, ncol = 1)
result_save_path <- args[7]
ggsave(
  result_save_path,
  plot = p,
  device = "pdf",
  width = 1.9,
  height = 3.2,
  dpi = "retina"
)