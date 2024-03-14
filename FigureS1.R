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

flow_sub <- flow |>
  filter(day<=20,box!="05",mouseid!="01-02",mouseid!="02-03") |>
  select(day,Pd,Eryth,Retic,box,mouse,mouseid,pABA)

flow_sub$mouse<-factor(flow_sub$mouse,levels=c("01","02","03"),labels=c("Mouse 1","Mouse 2","Mouse 3"))

parasites <- flow_sub |>
  filter(day<=20) |>
  ggplot()+
  geom_line(aes(x=day,y=Pd,group=mouseid,col=mouse),linewidth=1.5)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(1,10^7.5))+ 
  scale_colour_manual(values=cbPalette[6:8])+
  theme_bw()+
  facet_grid(.~pABA)+
  labs(x="",y="Parasites\n(density per µL)",col="")+
  ggtitle("")+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    strip.background=element_blank(),
    strip.text=element_text(size=15),
    #axis.title.x=element_blank(),
    plot.title=element_text(hjust=0.5,face="bold",size=15),
    panel.grid=element_blank(),
    legend.position="bottom",
    legend.title=element_text(size=14),
    legend.text=element_text(size=13),
    legend.background=element_blank()
  )
legend<-get_legend(parasites)
parasites<-parasites+theme(legend.position="none")
parasites

erythrocytes <- flow_sub |>
  filter(day<=20) |>
  ggplot()+
  geom_line(aes(x=day,y=Eryth,group=mouseid,col=mouse),linewidth=1.5)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(10^5.5,10^7))+ 
  scale_colour_manual(values=cbPalette[6:8])+
  theme_bw()+
  facet_grid(.~pABA)+
  labs(x="Time (d post-infection)",y="Erythrocytes\n(density per µL)",col="Parasite nutrient (pABA)")+
  ggtitle("")+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    strip.background=element_blank(),
    strip.text=element_text(size=15,colour="white"),
    #axis.title.x=element_blank(),
    plot.title=element_text(hjust=0.5,face="bold",size=15),
    panel.grid=element_blank(),
    legend.position="none",
    legend.title=element_text(size=14),
    legend.text=element_text(size=13),
    legend.background=element_blank()
  )
erythrocytes

reticulocytes <- flow_sub |>
  filter(day<=20) |>
  ggplot()+
  geom_line(aes(x=day,y=Retic,group=mouseid,col=mouse),linewidth=1.5)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(10^5,10^7))+ 
  scale_colour_manual(values=cbPalette[6:8])+
  theme_bw()+
  facet_grid(.~pABA)+
  labs(x="",y="Reticulcoytes\n(density per µL)",col="Parasite nutrient (pABA)")+
  ggtitle("")+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    strip.background=element_blank(),
    strip.text=element_text(size=15,colour="white"),
    #axis.title.x=element_blank(),
    plot.title=element_text(hjust=0.5,face="bold",size=15),
    panel.grid=element_blank(),
    legend.position="none",
    legend.background=element_blank()
  )
reticulocytes

gt <- arrangeGrob(parasites,reticulocytes,erythrocytes,legend,nrow=4,
                  heights=c(1,1,1,0.2))

as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A","B","C","D",
                            "E","F","G","H",
                            "I","J","K","L"), 
                  size = 15,
                  x = c(0.075,0.315,0.54,0.77,
                        0.08,0.32,0.54,0.77,
                        0.09,0.32,0.54,0.77), 
                  y = c(0.97,0.97,0.97,0.97,
                        0.66,0.66,0.66,0.66,
                        0.35,0.35,0.35,0.35))


ggsave("FigureS1.jpeg",width=25,height=25,units="cm") 
