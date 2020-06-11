#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Statistics::R;

open( my $IN, "<", $ARGV[0] );

my %sig;
while (<$IN>) {
    chomp;
    my ( undef, undef, undef, $lef, $nm, $rgt, undef ) = split;
    if ( exists( $sig{ $lef . $nm . $rgt } ) ) {
        $sig{ $lef . $nm . $rgt }++;
    }
    else {
        $sig{ $lef . $nm . $rgt } = 1;
    }
}

my $data_table;
foreach my $nm (qw(A G C T)) {
    foreach my $lef (qw(A G C T)) {
        foreach my $rgt (qw(A G C T)) {
            $data_table .= "$lef$nm$rgt\t$sig{$lef.$nm.$rgt}\n";
        }
    }
}

print("$data_table");

my $rcmd = <<'EOF';
library(ggplot2)
library(reshape2)
library(dplyr)
data <- read.table(text=data_ori, header=F, sep="\t", quote="")
data$V3 <- c(rep(x = "A",times=16),rep(x = "G",times=16),rep(x = "C",times=16),rep(x = "U",times=16))
data$V4 <- c(1:64)
color_scale <- c("#E8845D","#6ABACC","#D1B84E","#91CE6A")
sumA<-sum(data[data$V3=="A",]$V2)
sumG<-sum(data[data$V3=="G",]$V2)
sumC<-sum(data[data$V3=="C",]$V2)
sumU<-sum(data[data$V3=="U",]$V2)
Reduce('paste0', c("A\n",sumA))
Reduce('paste0', c("G\n",sumG))
Reduce('paste0', c("C\n",sumC))
Reduce('paste0', c("U\n",sumU))
text_to_plot=data.frame(x=c(8.5,24.5,40.5,56.5),
                        y=c(max(data$V2)+6,max(data$V2)+6,max(data$V2)+6,max(data$V2)+6),
                        text=c(Reduce('paste0', c("A\n",sumA)),
                               Reduce('paste0', c("G\n",sumG)),
                               Reduce('paste0', c("C\n",sumC)),
                               Reduce('paste0', c("U\n",sumU))))
line_to_plot1=data.frame(x=c(0.8,16.2),y=c(max(data$V2)+11,max(data$V2)+11))
line_to_plot2=data.frame(x=c(16.8,32.2),y=c(max(data$V2)+11,max(data$V2)+11))
line_to_plot3=data.frame(x=c(32.8,48.2),y=c(max(data$V2)+11,max(data$V2)+11))
line_to_plot4=data.frame(x=c(48.8,64.2),y=c(max(data$V2)+11,max(data$V2)+11))
p <- ggplot(data, aes(x=V4, y=V2))
p + geom_bar(stat="identity", position="dodge", aes(fill=V3), width = .5, show.legend = F) +
  xlim(data$V1) +
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
EOF

my $R = Statistics::R->new;

$R->set( 'data_ori',         $data_table );
$R->set( 'result_save_path', $ARGV[1] );

$R->run($rcmd);


