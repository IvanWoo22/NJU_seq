#!/usr/bin/env Rscript
args <- commandArgs(T)
library(ggplot2)

dat <- read.csv(args[1],header = T,sep = "\t",quote = "")

dat <- dat[dat$PValue<0.05,]
dat <- dat[order(dat$Category),]

dat$Alpha = 1

dat[dat$Category=="GOTERM_MF_DIRECT",]$Alpha <- c(-log10(dat[dat$Category=="GOTERM_MF_DIRECT",]$PValue)/(max(-log10(dat[dat$Category=="GOTERM_MF_DIRECT",]$PValue))))
dat[dat$Category=="GOTERM_CC_DIRECT",]$Alpha <- c(-log10(dat[dat$Category=="GOTERM_CC_DIRECT",]$PValue)/(max(-log10(dat[dat$Category=="GOTERM_CC_DIRECT",]$PValue))))
dat[dat$Category=="GOTERM_BP_DIRECT",]$Alpha <- c(-log10(dat[dat$Category=="GOTERM_BP_DIRECT",]$PValue)/(max(-log10(dat[dat$Category=="GOTERM_BP_DIRECT",]$PValue))))

p <- ggplot(dat, aes(x=Term, y=Count))
p + geom_bar(stat = "identity", aes(fill=Category), alpha=dat$Alpha, show.legend = F, width = .6) +
  scale_fill_manual(values = c("#46BC91", "#EF9C62", "#7FAFEA")) +
  scale_x_discrete(limits=factor(rev(dat$Term))) +
  geom_text(aes(label=Count, y=Count+3), vjust=0.5, hjust=0, size=3.3) +
  coord_flip() + theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "#F4F4F4",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

rowlen <- nrow(dat)
height <- 0.525 + rowlen * 0.169

result_save_path <- args[2]
ggsave(result_save_path, device = "pdf", width = 10, height = height, dpi = "retina")