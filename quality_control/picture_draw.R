#!/usr/bin/env Rscript
library("gridExtra")
library("grid")
library("ggpubr")

args<-commandArgs(T)

mydata1<-read.table(args[1],header = T)

p1<-ggtexttable(mydata1, rows = c("A","G","C","T","Total"),
                theme = ttheme(rownames.style = rownames_style(color = "white",fill = "#630460",face = "bold",hjust =0.5,x = 0.5,linecolor = "#ffffff"),
                colnames.style = colnames_style(color = "white",fill = "#630460",face = "bold",linecolor = "#ffffff"),
                tbody.style = tbody_style(color = "black",fill = "#ffffff",linecolor = "#630460")))

p1<-annotate_figure(p1,fig.lab = "  All Base Contain",fig.lab.face = "bold")

mydata2<-read.table(args[2],header = T)

p2<-ggtexttable(mydata2,rows = c("A","G","C","T","Total"),
                theme = ttheme(rownames.style = rownames_style(color = "white",fill = "#630460",face = "bold",hjust =0.5,x = 0.5,linecolor = "#ffffff"),
                colnames.style = colnames_style(color = "white",fill = "#630460",face = "bold",linecolor = "#ffffff"),
                tbody.style = tbody_style(color = "black",fill = "#ffffff",linecolor = "#630460")))

p2<-annotate_figure(p2,fig.lab = "  5' Base Contain",fig.lab.face = "bold")

mydata3<-read.table(args[3],header = T)

p3<-ggtexttable(mydata3,rows = c("A","G","C","T","Total"),
                theme = ttheme(rownames.style = rownames_style(color = "white",fill = "#630460",face = "bold",hjust =0.5,x = 0.5,linecolor = "#ffffff"),
                colnames.style = colnames_style(color = "white",fill = "#630460",face = "bold",linecolor = "#ffffff"),
                tbody.style = tbody_style(color = "black",fill = "#ffffff",linecolor = "#630460")))

p3<-annotate_figure(p3,fig.lab = "  3' Base Contain",fig.lab.face = "bold")

mydata4<-read.table(args[4],header = F)

p4<-ggplot(mydata4,aes(mydata4$V2,mydata4$V3,group = mydata4$V1,color = mydata4$V1)) +
    theme_bw() +
    theme(legend.position=c(0.8,1),legend.justification=c(0,1),legend.key.size=unit(0,'mm'),legend.key.width=unit(2,'mm')) +
    theme(legend.background = element_rect(fill = NULL, color = "black",size = 0.1)) +
    geom_line(size = 0.5)+scale_colour_hue("Sample") +
    labs( x='', y='') +
    geom_point(size = 0) +
    scale_x_continuous(limits = c(10,50),breaks = seq(10,50,5))

p4<-annotate_figure(p4,fig.lab = "  Length Distribution",fig.lab.face = "bold")

mydata5<-read.table(args[5],sep = "\t",header = T)

p5<-ggtexttable(mydata5,rows = c("Total Sequence","Length Range"),
                theme = ttheme(rownames.style = rownames_style(color = "white",fill = "#630460",face = "bold",hjust =0.5,x = 0.5,linecolor = "#ffffff"),
                colnames.style = colnames_style(color = "white",fill = "#630460",face = "bold",linecolor = "#ffffff"),
                tbody.style = tbody_style(color = "black",fill = "#ffffff",linecolor = "#630460")))

p5<-annotate_figure(p5,fig.lab = "  Over View",fig.lab.face = "bold")


p<-ggarrange(NULL,p5,p1,p2,p3,p4,ncol = 1,nrow = 6,heights = c(0.3,1,2,2,2,4))

pdf(args[6],width = 12,height = 15)
print(p)
dev.off()