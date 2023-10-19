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

(reticVStime <- group_traj |>
  filter(variable=="R",box!="05",time<20) |>
  ggplot()+
  geom_line(aes(x=time,y=med,col=pABA),linewidth=2)+
  geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)))+ 
  xlab("Day post-infection")+
  ylab("Reticulocyte supply (density per microlitre)")+
  scale_colour_manual(values=cbPalette[2:5])+
  scale_fill_manual(values=cbPalette[2:5])+
  theme_bw()+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    panel.grid=element_blank(),
    legend.position=c(0.2,0.7),
    legend.title=element_text(size=15),
    legend.text=element_text(size=13)
  )
)
ggsave("reticVStime.png",plot=reticVStime,
       width=20,height=15,units="cm") 

(nVStime <- group_traj |>
    filter(variable=="N",box!="05",time<20) |>
    ggplot()+
    geom_line(aes(x=time,y=med,col=pABA),linewidth=2)+
    geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)))+ 
    xlab("Day post-infection")+
    ylab("Indiscriminate killing (density per microlitre)")+
    scale_colour_manual(values=cbPalette[2:5])+
    scale_fill_manual(values=cbPalette[2:5])+
    theme_bw()+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=13),
      panel.grid=element_blank(),
      legend.position=c(0.2,0.75),
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )
)
ggsave("nVStime.png",plot=nVStime,
       width=20,height=15,units="cm") 

(pABA0swirl <- group_traj |>
    filter(variable%in%c("E","R"),box!="05") |>
    select(-lo,-hi) |>
    pivot_wider(names_from="variable",values_from="med") |>
    mutate(lagE=lag(E)) |> 
    filter(time!=0,time<=20,box=="04") |>
    ggplot()+
    geom_path(aes(x=lagE,y=R),col=cbPalette[2],linewidth=2)+
    geom_text(aes(x=lagE,y=R,label=time))+
    xlab("Erythrocytes (t-1)")+ylab("Reticulocyte supply (t)")+
    theme_bw()+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=13),
      panel.grid=element_blank(),
      legend.position=c(0.8,0.75),
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )
)
ggsave("pABA0swirl.png",plot=pABA0swirl,
       width=17,height=15,units="cm")

(pABAswirl <- group_traj |>
  filter(variable%in%c("E","R"),box!="05") |>
  select(-lo,-hi) |>
  pivot_wider(names_from="variable",values_from="med") |>
  mutate(lagE=lag(E)) |> 
  filter(time!=0,time<=20) |>
  ggplot()+
  geom_path(aes(x=lagE,y=R,col=pABA),linewidth=2)+
  geom_text(aes(x=lagE,y=R,label=time))+
  xlab("Erythrocytes (t-1)")+ylab("Reticulocyte supply (t)")+
  scale_colour_manual(values=cbPalette[2:5])+
  theme_bw()+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    panel.grid=element_blank(),
    legend.position=c(0.8,0.75),
    legend.title=element_text(size=15),
    legend.text=element_text(size=13)
  )
)
ggsave("pABAswirl.png",plot=pABAswirl,
       width=17,height=15,units="cm")
  
