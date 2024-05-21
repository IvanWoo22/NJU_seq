#!/usr/bin/env Rscript
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
args <- commandArgs(T)
dist_raw <- read.table(args[1], sep = "\t", header = F)

x1max <-
  ceiling(max(dist_raw[dist_raw$V1 == "five_utr",]$V2) / 5) * 5
x2max <- ceiling(max(dist_raw[dist_raw$V1 == "cds",]$V2) / 5) * 5
x3max <-
  ceiling(max(dist_raw[dist_raw$V1 == "three_utr",]$V2) / 5) * 5
y1max <-
  ceiling(max(dist_raw[dist_raw$V1 == "five_utr",]$V3) / 5) * 5
y2max <- ceiling(max(dist_raw[dist_raw$V1 == "cds",]$V3) / 5) * 5
y3max <-
  ceiling(max(dist_raw[dist_raw$V1 == "three_utr",]$V3) / 5) * 5
sum1 <- sum(dist_raw[dist_raw$V1 == "five_utr",]$V3)
sum2 <- sum(dist_raw[dist_raw$V1 == "cds",]$V3)
sum3 <- sum(dist_raw[dist_raw$V1 == "three_utr",]$V3)
sum <- sum(dist_raw$V3)
prop1 <- paste(round(sum1 / sum, 4) * 100, "%", sep = '')
prop2 <- paste(round(sum2 / sum, 4) * 100, "%", sep = '')
prop3 <- paste(round(sum3 / sum, 4) * 100, "%", sep = '')
ymax <- ceiling((max(dist_raw$V3) + 1) / 5) * 5

five_utr <- dist_raw[dist_raw$V1 == "five_utr",]
cds <- dist_raw[dist_raw$V1 == "cds",]
three_utr <- dist_raw[dist_raw$V1 == "three_utr",]

text_to_plot <- data.frame(
  x = c(x1max / 2, x1max + x2max / 2, x1max + x2max + x3max / 2),
  y = c(y1max,
        y2max,
        y3max),
  text = c(
    Reduce('paste0', c(sum1, "\n", prop1)),
    Reduce('paste0', c(sum2, "\n", prop2)),
    Reduce('paste0', c(sum3, "\n", prop3))
  )
)

p <- ggplot() +
  geom_line(
    aes(x = five_utr$V2,
        y = five_utr$V3),
    colour = "#996699",
    size = 0.701,
    show.legend = F
  ) +
  geom_line(
    aes(x = cds$V2 + x1max,
        y = cds$V3),
    colour = "#FF9933",
    size = 0.701,
    show.legend = F
  ) +
  geom_line(
    aes(x = three_utr$V2 + x1max + x2max,
        y = three_utr$V3),
    colour = "#996666",
    size = 0.701,
    show.legend = F
  ) +
  geom_text(
    data = text_to_plot,
    aes(
      x = x,
      y = y,
      label = text,
      fontface = "bold",
      vjust = 0
    ),
    family = "ArialMT",
    colour = c("#996699", "#FF9933", "#996666"),
    size = 5 / 3,
    show.legend = F
  ) +
  labs(x = '', y = '') +
  scale_y_continuous(
    limits = c(0, ymax),
    breaks = seq(0, ymax, 10),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    limits = c(0, x1max + x2max + x3max),
    breaks = seq(0, x1max + x2max + x3max, 5),
    expand = c(0, 0)
  ) +
  theme(
    plot.margin = unit(c(0.4, 0.4, -0.5, -0.6), 'lines'),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "#000000", size = 0.11672),
    axis.ticks.length = unit(0.4, 'mm'),
    axis.line.y = element_line(colour = "#000000", size = 0.11672),
    axis.text.y = element_text(
      face = "plain",
      colour = "#000000",
      size = 5
    ),
    text = element_text(family = "ArialMT")
  )

result_save_path <- args[2]
ggsave(
  result_save_path,
  plot = p,
  device = "pdf",
  width = 3.3,
  height = 1.5,
  dpi = "retina"
)