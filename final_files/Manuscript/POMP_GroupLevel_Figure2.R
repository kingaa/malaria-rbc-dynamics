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
  filter(day<=20,box!="05",mouseid!="01-02",mouseid!="02-03") |>
  ggplot()+
  geom_line(aes(x=day,y=Pd,group=mouseid,col=pABA))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+ 
  scale_colour_manual(values=cbPalette[2:5])+
  theme_bw()+
  xlab("")+ylab("Parasites\n(density per µL)")+
    ggtitle("Raw data")+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    axis.title.x=element_blank(),
    plot.title=element_text(hjust=0.5,face="bold",size=15),
    panel.grid=element_blank(),
    legend.position=c(0.35,0.3),
    legend.title=element_text(size=15),
    legend.text=element_text(size=13),
    legend.background=element_blank()
  )
)

(reticulocytes <- flow |>
  filter(day<=20,box!="05",mouseid!="01-02",mouseid!="02-03") |>
  ggplot()+
  geom_line(aes(x=day,y=Retic,group=mouseid,col=pABA))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+ 
  scale_colour_manual(values=cbPalette[2:5])+
  theme_bw()+
  xlab("")+ylab("Reticulocyte supply\n(density per µL)")+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    axis.title.x=element_blank(),
    panel.grid=element_blank(),
    legend.position="none",
    legend.title=element_text(size=15),
    legend.text=element_text(size=13)
  )
)

(erythrocytes <- flow |>
    filter(day<=20,box!="05",mouseid!="01-02",mouseid!="02-03") |>
    ggplot()+
    geom_line(aes(x=day,y=Eryth,group=mouseid,col=pABA))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)))+ 
    scale_colour_manual(values=cbPalette[2:5])+
    theme_bw()+
    xlab("Day post-infection")+ylab("Erythrocytes\n(density per µL)")+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=13),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )
)

data_column <- ggarrange(parasites,reticulocytes,erythrocytes, nrow = 3, labels = c("A","B","C"))

(K_model <- group_traj |>
  filter(box!="05",variable=="K",time<=20) |>
  mutate(lo=if_else(lo<=1,1,lo)) |>
  ggplot()+
  geom_line(aes(x=time,y=med,col=pABA),linewidth=2)+
  geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
              labels=trans_format('log10',math_format(10^.x)))+ 
  scale_colour_manual(values=cbPalette[2:5])+
  scale_fill_manual(values=cbPalette[2:5])+
  ggtitle("Group-level model output")+
  theme_bw()+
  xlab("")+ylab("Parasites\n(density per µL)")+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    plot.title=element_text(size=15,face="bold",hjust=0.5),
    panel.grid=element_blank(),
    legend.position="none",
    legend.title=element_text(size=15),
    legend.text=element_text(size=13)
  )
)

(N_model <- group_traj |>
    filter(box!="05",variable=="N",time<=20) |>
    mutate(lo=if_else(lo<=1,1,lo)) |>
    ggplot()+
    geom_line(aes(x=time,y=med,col=pABA),linewidth=2)+
    geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)))+ 
    scale_colour_manual(values=cbPalette[2:5])+
    scale_fill_manual(values=cbPalette[2:5])+
    theme_bw()+
    xlab("")+ylab("RBC clearance\n(density per µL)")+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=13),
      axis.title.x=element_blank(),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )
)

(R_model <- group_traj |>
    filter(box!="05",variable=="R",time<=20) |>
    mutate(lo=if_else(lo<=1,1,lo)) |>
    ggplot()+
    geom_line(aes(x=time,y=med,col=pABA),linewidth=2)+
    geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)))+ 
    scale_colour_manual(values=cbPalette[2:5])+
    scale_fill_manual(values=cbPalette[2:5])+
    theme_bw()+
    xlab("")+ylab("Reticulocyte supply\n(density per µL)")+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=13),
      panel.grid=element_blank(),
      axis.title.x=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )
)

(E_model <- group_traj |>
    filter(box!="05",variable=="E",time<=20) |>
    mutate(lo=if_else(lo<=1,1,lo)) |>
    ggplot()+
    geom_line(aes(x=time,y=med,col=pABA),linewidth=2)+
    geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)))+ 
    scale_colour_manual(values=cbPalette[2:5])+
    scale_fill_manual(values=cbPalette[2:5])+
    theme_bw()+
    xlab("Day post-infection")+ylab("Erythrocytes\n(density per µL)")+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=13),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )
)

output_column <- ggarrange(K_model,R_model,E_model, nrow = 3, labels = c("D","E","F"))

joint_columns <- ggarrange(data_column,output_column,ncol=2)  

ggsave("data_output.png",plot=joint_columns,
       width=25,height=30,units="cm") 
