#!/usr/bin/env Rscript
library("VennDiagram")

args<-commandArgs(T)

data1 = read.table(args[2],header = F,sep="\t")
data2 = read.table(args[4],header = F,sep="\t")
data3 = read.table(args[6],header = F,sep="\t")

datalist <- list(data1$V1, data2$V1, data3$V1)
names(datalist) <- c(args[1], args[3], args[5])

venn.plot <- venn.diagram(
  x = datalist,
  filename = args[7], imagetype = "png",
  lwd = 3,
  col = c("#489DCC", "#C6175A", "#EF7E5B"),
  fill = "#F4F4F4",
  alpha = 0.6,
  label.col = "#000000",
  cex = 3,
  fontfamily = "sans",
  fontface = "bold",
  cat.col = c("#489DCC", "#C6175A", "#EF7E5B"),
  cat.cex = 2,
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  margin = 0.05,
  cat.dist = c(0.06, 0.06, 0.06),
  cat.pos = c(-20, 20, 180)
)
