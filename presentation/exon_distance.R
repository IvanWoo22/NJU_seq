#!/usr/bin/env Rscript
library("ggplot2")
args <- commandArgs(T)

count <- read.table(args[1], header = F)
porta <- read.table(args[2], header = F)
p <- ggplot() +
  xlab("Distance (nt)") +
  ylab("Nmuber of Nm Sites") +
  geom_bar(
    data = count,
    mapping = aes(
      x = V2,
      y = V3,
      fill = V1,
      color = V1
    ),
    stat = 'identity',
    width = 10,
    size = 0.04669,
    show.legend = F
  ) +
  scale_fill_manual(values = c("#4D7373", "#D67272")) +
  scale_color_manual(values = c("#94ABAB", "#8B5457")) +
  geom_line(aes(x = -porta$V1, y = porta$V2),
            size = 0.4669,
            colour = "#BFB057") +
  geom_line(aes(x = porta$V1, y = porta$V2),
            size = 0.4669,
            colour = "#BFB057") +
  geom_line(aes(x = c(0, 0), y = c(0, ceiling(max(
    count$V3
  ) / 5) * 5)), size = 0.4669, colour = "red") +
  scale_x_continuous(
    limits = c(-300, 300),
    breaks = seq(-500, 500, 50),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, ceiling(max(count$V3) / 5) * 5),
    breaks = seq(0, ceiling(max(count$V3) / 5) * 5, 10),
    expand = c(0, 0)
  ) +
  theme(
    axis.text.x = element_text(
      face = "plain",
      colour = "#000000",
      size = 4
    ),
    axis.line = element_line(colour = "#000000", size = 0.11672),
    axis.ticks = element_line(colour = "#000000", size = 0.11672),
    axis.ticks.length = unit(0.4, 'mm'),
    axis.text.y = element_text(
      face = "plain",
      colour = "#000000",
      size = 5
    ),
    axis.title = element_text(
      face = "plain",
      colour = "#000000",
      size = 7
    ),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'lines'),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(family = "ArialMT")
  )

result_save_path <- args[3]
ggsave(
  result_save_path,
  plot = p,
  device = "pdf",
  width = 1.9,
  height = 1.6,
  dpi = "retina"
)