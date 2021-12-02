#!/usr/bin/env Rscript
require("ggplot2")
require("ggseqlogo")

data <- read.table(args[1], header = F)

p <- ggseqlogo(data, method = 'p') +
  theme(
    axis.text.x = element_text(face = "plain",
                               colour = "#000000"),
    axis.text.y = element_text(face = "plain",
                               colour = "#000000"),
    axis.title = element_text(face = "plain",
                              colour = "#000000"),
    text = element_text(family = "ArialMT")
  )

result_save_path <- args[2]

width = nchar(data[1, 1]) * .3

ggsave(
  result_save_path,
  plot = p,
  device = "pdf",
  width = width,
  height = 4.0,
  dpi = "retina"
)