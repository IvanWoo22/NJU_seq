#!/usr/bin/env Rscript
library(VennDiagram)

args <- commandArgs(T)
stopifnot(length(args) == 7)
read_set <- function(path) { if(file.info(path)$size == 0) character() else unique(read.table(path, header=FALSE, sep="\t")[[1]]) }

datalist <- lapply(args[c(2,4,6)], read_set)
names(datalist) <- c(args[1], args[3], args[5])

venn.plot <- venn.diagram(
  x = datalist,
  disable.logging = TRUE,
  filename = args[7],
  imagetype = "png",
  lwd = 3,
  col = c("#5E95CA", "#91211C", "#F29125"),
  fill = "#F4F4F4",
  alpha = 0.6,
  label.col = "#000000",
  cex = 3,
  fontfamily = "Arial",
  cat.col = c("#5E95CA", "#91211C", "#F29125"),
  cat.cex = 2,
  cat.fontfamily = "Arial",
  margin = 0.05,
  cat.dist = c(0.06, 0.06, 0.06),
  cat.pos = c(-20, 20, 180)
)
