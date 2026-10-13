#!/usr/bin/env Rscript
library(ggplot2)
library(ggseqlogo)

args <- commandArgs(T)
data <- read.table(args[1], header = F)

seqs <- toupper(chartr("T", "U", data[[1]]))
stopifnot(length(seqs) > 0, length(unique(nchar(seqs))) == 1, all(grepl("^[ACGUN]+$", seqs)))
p <- ggseqlogo(seqs, method = 'p') +
  theme(
    axis.text.x = element_text(face = "plain",
                               colour = "#000000"),
    axis.text.y = element_text(face = "plain",
                               colour = "#000000"),
    axis.title = element_text(face = "plain",
                              colour = "#000000"),
    text = element_text(family = "sans")
  )

result_save_path <- args[2]
width <- nchar(data[1, 1]) * .3
ggsave(
  result_save_path,
  plot = p,
  device = "pdf",
  width = width,
  height = 4.0,
  dpi = "retina"
)