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
  filter(pABA=="Low")
flow_sub$Mouse <- factor(flow_sub$mouse,levels=c("01","02","03"),labels=c("Mouse 1","Mouse 2","Mouse 3"))


sm1_sub <- sm1 |>
  filter(mouseid=="03-01"|mouseid=="03-02"|mouseid=="03-03")
sm1_sub$Mouse <- factor(sm1_sub$mouseid,levels=c("03-01","03-02","03-03"),labels=c("Mouse 1","Mouse 2","Mouse 3"))
sm1_sub <- sm1_sub |>
  unite("key",c(rep,Mouse),sep="_",remove=FALSE)

(parasites <- flow_sub |>
    filter(day<=20) |>
    ggplot()+
    geom_line(aes(x=day,y=Pd,group=Mouse,col=Mouse))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(1,10^7.5))+ 
    scale_colour_manual(values=cbPalette[6:8])+
    theme_bw()+
    xlab("")+ylab("Parasites\n(density per µL)")+
    ggtitle("Raw data")+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=13),
      axis.title.x=element_blank(),
      plot.title=element_text(hjust=0.5,face="bold",size=15),
      panel.grid=element_blank(),
      legend.position=c(0.2,0.3),
      legend.title=element_blank(),
      legend.text=element_text(size=13),
      legend.background=element_blank()
    )
)

(reticulocytes <- flow_sub |>
    filter(day<=20) |>
    ggplot()+
    geom_line(aes(x=day,y=Retic,group=Mouse,col=Mouse))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(10^5,10^7))+  
    scale_colour_manual(values=cbPalette[6:8])+
    theme_bw()+
    xlab("")+ylab("Reticulocyte supply\n(density per µL)")+
    #ggtitle("")+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=13),
      axis.title.x=element_blank(),
      plot.title=element_text(hjust=0.5,face="bold",size=15),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )
)

(erythrocytes <- flow_sub |>
    filter(day<=20) |>
    ggplot()+
    geom_line(aes(x=day,y=Eryth,group=Mouse,col=Mouse))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(10^5,10^7))+  
    scale_colour_manual(values=cbPalette[6:8])+
    theme_bw()+
    xlab("Day post-infection")+ylab("Erythrocytes\n(density per µL)")+
    #ggtitle("")+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=13),
      plot.title=element_text(hjust=0.5,face="bold",size=15),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )
)

data_column <- ggarrange(parasites,reticulocytes,erythrocytes, nrow = 3, labels = c("A","B","C"))

(K_model <- sm1_sub |>
    filter(variable=="K",time<=20) |>
    #mutate(value=ifelse(value<=0,1,value)) |>
    ggplot()+
    geom_line(aes(x=time,y=value,group=key,col=Mouse),alpha=0.01)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(1,10^8))+ 
    scale_colour_manual(values=cbPalette[6:8])+
    scale_fill_manual(values=cbPalette[6:8])+
    ggtitle("Data-transformation model fits")+
    theme_bw()+
    xlab("")+ylab("Parasites,\nK (density per µL)")+
    theme(
      axis.title.y=element_text(size=11),
      axis.title.x=element_text(size=15),
      axis.text=element_text(size=13),
      plot.title=element_text(size=15,face="bold",hjust=0.5),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )
)

(N_model <- sm1_sub |>
    filter(variable=="N",time<=20) |>
    ggplot()+
    geom_line(aes(x=time,y=value,group=key,col=Mouse),alpha=0.01)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(10^4,10^7))+
    scale_colour_manual(values=cbPalette[6:8])+
    scale_fill_manual(values=cbPalette[6:8])+
    theme_bw()+
    xlab("")+ylab("Indiscriminate killing,\nN (density per µL)")+
    #ggtitle("")+
    theme(
      axis.title.y=element_text(size=11),
      axis.title.x=element_text(size=15),
      axis.text=element_text(size=13),
      plot.title=element_text(hjust=0.5,face="bold",size=15),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13),
      panel.border = element_rect(linewidth=4)
    )
)

(W_model <- sm1_sub |>
    filter(variable=="W",time<=20) |>
    ggplot()+
    geom_line(aes(x=time,y=value,group=key,col=Mouse),alpha=0.01)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(10^3,10^10))+ 
    scale_colour_manual(values=cbPalette[6:8])+
    scale_fill_manual(values=cbPalette[6:8])+
    theme_bw()+
    xlab("Day post-infection")+ylab("Targeted killing,\nW (density per µL)")+
    #ggtitle("")+
    theme(
      axis.title.y=element_text(size=11),
      axis.title.x=element_text(size=15),
      axis.text=element_text(size=13),
      plot.title=element_text(hjust=0.5,face="bold",size=15),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )
)

(R_model <- sm1_sub |>
    filter(variable=="R",time<=20) |>
    ggplot()+
    geom_line(aes(x=time,y=value,group=key,col=Mouse),alpha=0.01)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(10^5,10^7))+   
    scale_colour_manual(values=cbPalette[6:8])+
    scale_fill_manual(values=cbPalette[6:8])+
    theme_bw()+
    xlab("")+ylab("Reticulocyte supply,\nR (density per µL)")+
    #ggtitle("")+
    theme(
      axis.title.y=element_text(size=11),
      axis.title.x=element_text(size=15),
      axis.text=element_text(size=13),
      plot.title=element_text(hjust=0.5,face="bold",size=15),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13),
      panel.border = element_rect(linewidth=4)
    )
)

(E_model <- sm1_sub |>
    filter(variable=="E",time<=20) |>
    ggplot()+
    geom_line(aes(x=time,y=value,group=key,col=Mouse),alpha=0.01)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(10^5,10^7))+   
    scale_colour_manual(values=cbPalette[6:8])+
    scale_fill_manual(values=cbPalette[6:8])+
    theme_bw()+
    xlab("")+ylab("Erythrocytes,\nE (density per µL)")+
    #ggtitle("")+
    theme(
      axis.title.y=element_text(size=11),
      axis.title.x=element_text(size=15),
      axis.text=element_text(size=13),
      plot.title=element_text(hjust=0.5,face="bold",size=15),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13),
      panel.border = element_rect(linewidth=4)
    )
)

output_column <- ggarrange(K_model,R_model,E_model,N_model,W_model, nrow = 5, labels = c("D","E","F","G","H"))
output_column

joint_columns <- ggarrange(data_column,output_column,ncol=2)

ggsave("FigureS2.png",plot=joint_columns,
       width=25,height=30,units="cm") 
