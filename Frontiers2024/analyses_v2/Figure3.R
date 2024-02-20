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
library(egg)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

source("POMP_GroupLevel_DataPrep.R")

df1 <- group_traj |>
  filter(box!="05",variable%in%c("E","R")) |>
  select(-lo,-hi) |>
  pivot_wider(names_from="variable",values_from="med") |>
  mutate(RBC = E+R,
         lagRBC=lag(RBC,1)) |>
  filter(time>0,time<=20) |>
  mutate(lag="i = 1")

df2 <- group_traj |>
  filter(box!="05",variable%in%c("E","R")) |>
  select(-lo,-hi) |>
  pivot_wider(names_from="variable",values_from="med") |>
  mutate(RBC = E+R,
         lagRBC=lag(RBC,2)) |>
  filter(time>0,time<=20) |>
  mutate(lag="i = 2")

df3 <- group_traj |>
  filter(box!="05",variable%in%c("E","R")) |>
  select(-lo,-hi) |>
  pivot_wider(names_from="variable",values_from="med") |>
  mutate(RBC = E+R,
         lagRBC=lag(RBC,3)) |>
  filter(time>0,time<=20) |>
  mutate(lag="i = 3")

df4 <- group_traj |>
  filter(box!="05",variable%in%c("E","R")) |>
  select(-lo,-hi) |>
  pivot_wider(names_from="variable",values_from="med") |>
  mutate(RBC = E+R,
         lagRBC=lag(RBC,4)) |>
  filter(time>0,time<=20) |>
  mutate(lag="i = 4")

df <- df1 |> bind_rows(df2) |> bind_rows(df3) |> bind_rows(df4)

labels <- data.frame(lag=c("i = 1","i = 2","i = 3","i = 4"),
                 label=c("A","B","C","D"))

facet_plot <- df |>
  ggplot()+
  geom_path(aes(x=lagRBC,y=R,col=pABA),linewidth=2)+
  geom_text(aes(x=lagRBC,y=R,label=time),size=5)+
  geom_text(data=labels,aes(x=2100000,y=4200000,label=label),size=5,fontface='bold')+
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
    legend.position="none",
    legend.title=element_text(size=15),
    legend.text=element_text(size=13),
    strip.background=element_blank(),
    strip.text=element_text(size=13),
    plot.title=element_text(size=17,hjust=0.5)
    
  )
facet_plot

ggsave("Figure3.jpeg",width=20,height=20,units="cm")
