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
library(cowplot)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

source("POMP_GroupLevel_DataPrep.R")

flow |> 
  select(time=day,box,pABA,E=Eryth,R=Retic,pABA,mouseid) |>
  mutate(RBC=E+R,
         lagRBC=lag(RBC,1),
         lag="i = 1") |>
  na.omit() |>
  filter(time<=20,box!="05",mouseid!="01-02",mouseid!="02-03") -> data1

flow |> 
  select(time=day,box,pABA,E=Eryth,R=Retic,pABA,mouseid) |>
  mutate(RBC=E+R,
         lagRBC=lag(RBC,2),
         lag="i = 2") |>
  na.omit() |>
  filter(time<=20,box!="05",mouseid!="01-02",mouseid!="02-03") -> data2

flow |> 
  select(time=day,box,pABA,E=Eryth,R=Retic,pABA,mouseid) |>
  mutate(RBC=E+R,
         lagRBC=lag(RBC,3),
         lag="i = 3") |>
  na.omit() |>
  filter(time<=20,box!="05",mouseid!="01-02",mouseid!="02-03") -> data3

flow |> 
  select(time=day,box,pABA,E=Eryth,R=Retic,pABA,mouseid) |>
  mutate(RBC=E+R,
         lagRBC=lag(RBC,4),
         lag="i = 4") |>
  na.omit() |>
  filter(time<=20,box!="05",mouseid!="01-02",mouseid!="02-03") -> data4

flow |> 
  select(time=day,box,pABA,E=Eryth,R=Retic,pABA,mouseid) |>
  mutate(RBC=E+R,
         lagRBC=lag(RBC,5),
         lag="i = 5") |>
  na.omit() |>
  filter(time<=20,box!="05",mouseid!="01-02",mouseid!="02-03") -> data5

data <- data1 |> bind_rows(data2) |> bind_rows(data3) |> bind_rows(data4) |> bind_rows(data5)

labels <- data.frame(lag=c("i = 1","i = 2","i = 3","i = 4","i = 5"),
                     label=c("A","B","C","D","E"))

data |>
  ggplot()+
  geom_path(aes(x=lagRBC,y=R,col=pABA),linewidth=1)+
  geom_text(aes(x=lagRBC,y=R,label=time),size=5)+
  geom_text(data=labels,aes(x=9750000,y=4500000,label=label),size=5,fontface='bold')+
  scale_x_continuous(labels=aakmisc::scinot)+
  scale_y_continuous(labels=aakmisc::scinot)+
  scale_colour_manual(values=cbPalette[2:5])+
  labs(x=expression(paste("RBC density at time ", italic("t"), "-i (density per µL)")),
       y=expression(paste("Reticulocyte supply at time ", italic("t"), " (density per µL)"))
  )+
  facet_wrap(lag~.,nrow=2)+
  labs(colour="Parasite nutrient (pABA)")+
  theme_bw()+
  theme(
    axis.title=element_text(size=15),
    axis.title.y=element_text(size=15,colour="black"),
    axis.text=element_text(size=11),
    panel.grid=element_blank(),
    legend.position=c(0.825,0.2),
    legend.title=element_text(size=12),
    legend.text=element_text(size=11),
    strip.background=element_blank(),
    strip.text=element_text(size=13),
    plot.title=element_text(size=17,hjust=0.5)
    
  )

ggsave("FigureS7.jpeg",width=30,height=20,units="cm")
