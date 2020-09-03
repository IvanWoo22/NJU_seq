#!/usr/bin/env Rscript
args <- commandArgs(T)

library("ggplot2")
library("reshape2")
library("dplyr")

data <- read.table(args[1], header=F, sep="\t", quote="")
data$V3 <- c(rep(x = "A",times=16),rep(x = "G",times=16),rep(x = "C",times=16),rep(x = "U",times=16))
data$V4 <- c(1:64)
color_scale <- c("#E8845D","#6ABACC","#D1B84E","#91CE6A")
sumA<-sum(data[data$V3=="A",]$V2)
propA<-sum(data[data$V3=="A",]$V2)/sum(data$V2)
round(propA,digits=4)
paste(propA*100, "%", sep='')
sumG<-sum(data[data$V3=="G",]$V2)
propG<-sum(data[data$V3=="G",]$V2)/sum(data$V2)
round(propG,digits=4)
paste(propG*100, "%", sep='')
sumC<-sum(data[data$V3=="C",]$V2)
propC<-sum(data[data$V3=="C",]$V2)/sum(data$V2)
round(propC,digits=4)
paste(propC*100, "%", sep='')
sumU<-sum(data[data$V3=="U",]$V2)
propU<-sum(data[data$V3=="U",]$V2)/sum(data$V2)
round(propU,digits=4)
paste(propU*100, "%", sep='')

result_save_path <- args[2]

text_to_plot=data.frame(x=c(8.5,24.5,40.5,56.5),
                        y=c(max(data$V2)+6,max(data$V2)+6,max(data$V2)+6,max(data$V2)+6),
                        text=c(Reduce('paste0', c("Am\n",sumA,"\n",propA)),
                               Reduce('paste0', c("Gm\n",sumG,"\n",propG)),
                               Reduce('paste0', c("Cm\n",sumC,"\n",propC)),
                               Reduce('paste0', c("Um\n",sumU,"\n",propU))))
line_to_plot1=data.frame(x=c(0.8,16.2),y=c(max(data$V2)+11,max(data$V2)+11))
line_to_plot2=data.frame(x=c(16.8,32.2),y=c(max(data$V2)+11,max(data$V2)+11))
line_to_plot3=data.frame(x=c(32.8,48.2),y=c(max(data$V2)+11,max(data$V2)+11))
line_to_plot4=data.frame(x=c(48.8,64.2),y=c(max(data$V2)+11,max(data$V2)+11))
p <- ggplot(data, aes(x=V4, y=V2))
p + geom_bar(stat="identity", position="dodge", aes(fill=V3), width = .5, show.legend = F) +
  ggplot2:::limits(data$V1, "x") +
  scale_y_continuous(expand = c(0, 1)) +
  scale_fill_manual(values = color_scale) +
  theme_bw() +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 13, face = "bold")) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  geom_text(aes(label=V2, y=V2+0.5), vjust=0, size=3.3) +
  geom_line(data = line_to_plot1,aes(x=x,y=y,size=11), color = '#E8845D', show.legend = F) +
  geom_line(data = line_to_plot2,aes(x=x,y=y,size=11), color = '#6ABACC', show.legend = F) +
  geom_line(data = line_to_plot3,aes(x=x,y=y,size=11), color = '#D1B84E', show.legend = F) +
  geom_line(data = line_to_plot4,aes(x=x,y=y,size=11), color = '#91CE6A', show.legend = F) +
  geom_text(data=text_to_plot,aes(x=x,y=y,label=text,size=12,fontface = "bold",vjust = 0.5), show.legend = F) +
  theme(plot.margin = unit(c(1,1,0,0),'lines'),
        panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "#F4F4F4",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

ggsave(result_save_path, device = "pdf", width = 20, height = 10, dpi = "retina")
dev.off()