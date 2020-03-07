#!/usr/bin/env Rscript
library("ggplot2")
library("gridExtra")
library("RColorBrewer")
args<-commandArgs(T)
dist_raw<-read.table(args[1],sep="\t",header = F)

five_utr<-dist_raw[dist_raw$V4=="five_utr",]
p1<-ggplot(data = five_utr) +
  stat_bin(bins = 50, geom = 'line', size = 0.5, aes(x=five_utr$V5)) +
  labs( x='5\'UTR', y='')+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_y_continuous(limits = c(0,55),
                     breaks= c(0,10,20,30,40,50),
                     expand = c(0, 0))+
  scale_x_continuous(limits = c(0,1),
                     breaks= c(0,0.25,0.5,0.75,1),
                     expand = c(0, 0))

cds<-dist_raw[dist_raw$V4=="cds",]
p2<-ggplot(data = cds) +
  stat_bin(bins = 50, geom = 'line', size = 0.5, aes(x=cds$V5)) +
  labs( x='CDS', y='')+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"))+
  scale_y_continuous(limits = c(0,55),
                     breaks=250,
                     expand = c(0, 0))+
  scale_x_continuous(limits = c(0,1),
                     breaks= c(0,0.25,0.5,0.75,1),
                     expand = c(0, 0))+
  theme(axis.text.y = element_text(size = 0))

three_utr<-dist_raw[dist_raw$V4=="three_utr",]
p3<-ggplot(data = three_utr) +
  stat_bin(bins = 50, geom = 'line', size = 0.5, aes(x=three_utr$V5)) +
  labs( x='3\'UTR', y='')+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"))+
  scale_y_continuous(limits = c(0,55),
                     breaks=250,
                     expand = c(0, 0))+
  scale_x_continuous(limits = c(0,1),
                     breaks= c(0,0.25,0.5,0.75,1),
                     expand = c(0, 0))+
  theme(axis.text.y = element_text(size = 0))

grid.arrange(p1,p2,p3,ncol=3)