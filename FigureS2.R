#SET DIRECTORY TO SOURCE FILE LOCATION

library(tidyverse)
library(stringi)
library(doParallel)
library(pomp)
library(panelPomp)
library(ggpubr)
library(scales)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

source("POMP_GroupLevel_DataPrep.R")

(parasites <- flow |>
  filter(day<=5,box!="05",mouseid!="01-02",mouseid!="02-03") |>
  ggplot()+
  geom_line(aes(x=day,y=Pd,group=mouseid,col=pABA))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(1,10^7.5))+ 
  scale_colour_manual(values=cbPalette[2:5])+
  theme_bw()+
  xlab("Time (d post-infection)")+ylab("Parasites (density per µL)")+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    plot.title=element_text(hjust=0.5,face="bold",size=15),
    panel.grid=element_blank(),
    legend.position=c(0.35,0.3),
    legend.title=element_text(size=15),
    legend.text=element_text(size=13),
    legend.background=element_blank()
  )
)

ggsave("FigureS2.jpeg",width=16,height=16,units="cm")
