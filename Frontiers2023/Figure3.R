#SET DIRECTORY TO SOURCE FILE LOCATION

library(tidyverse)
library(stringi)
library(doParallel)
library(pomp)
library(panelPomp)
library(ggpubr)
library(scales)
library(aakmisc)
library(gridExtra)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

source("POMP_GroupLevel_DataPrep.R")

group_traj |>
  filter(box!="05",variable%in%c("E","R")) |>
  select(-lo,-hi) |>
  pivot_wider(names_from="variable",values_from="med") |>
  mutate(lagE=lag(E)) |>
  filter(time>0,time<=20) |>
  ggplot()+
  geom_path(aes(x=lagE,y=R,col=pABA),linewidth=2)+
  geom_text(aes(x=lagE,y=R,label=time),size=5)+
  scale_x_continuous(labels=aakmisc::scinot)+
  scale_y_continuous(labels=aakmisc::scinot)+
  scale_colour_manual(values=cbPalette[2:5])+
  #xlab("Erythrocyte density at time t-1 (density per µL)")+ylab("RBC supply at time t (density per µL)")+
  
  labs(x=expression(paste("Erythrocyte density at time ", italic("t"), "-1 (density per µL)")),
       y=expression(paste("RBC supply at time ", italic("t"), " (density per µL)"))
       )+
  
  
  labs(colour="Parasite nutrient (pABA)")+
  theme_bw()+
  theme(
    axis.title=element_text(size=15),
    strip.text=element_blank(),
    axis.text=element_text(size=11),
    panel.grid=element_blank(),
    legend.position=c(0.7,0.8),
    legend.title=element_text(size=15),
    legend.text=element_text(size=13),
    strip.background=element_blank(),
    plot.title=element_text(size=17,hjust=0.5)
    
  )

ggsave("Figure3.jpeg",width=16,height=16,units="cm")
