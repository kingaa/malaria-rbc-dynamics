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
         lagRBC=lag(RBC,1)) |>
  na.omit() |>
  filter(time<=20,box!="05",mouseid!="01-02",mouseid!="02-03") -> data1

swirl1 <- data1 |>
  ggplot()+
  geom_path(aes(x=lagRBC,y=R,col=pABA),linewidth=2)+
  geom_text(aes(x=lagRBC,y=R,label=time),size=5)+
  scale_x_continuous(labels=aakmisc::scinot)+
  scale_y_continuous(labels=aakmisc::scinot)+
  scale_colour_manual(values=cbPalette[2:5])+
  labs(x=expression(paste("RBC density at time ", italic("t"), "-1 (density per µL)")),
       y=expression(paste("Reticulocyte supply at time ", italic("t"), " (density per µL)"))
  )+
  labs(colour="Parasite nutrient (pABA)")+
  theme_bw()+
  theme(
    axis.title=element_text(size=15),
    axis.title.y=element_text(size=15,colour="black"),
    strip.text=element_blank(),
    axis.text=element_text(size=11),
    panel.grid=element_blank(),
    legend.position="none",
    legend.title=element_text(size=15),
    legend.text=element_text(size=13),
    strip.background=element_blank(),
    plot.title=element_text(size=17,hjust=0.5)
    
  )


flow |> 
  select(time=day,box,pABA,E=Eryth,R=Retic,pABA,mouseid) |>
  mutate(RBC=E+R,
         lagRBC=lag(RBC,2)) |>
  na.omit() |>
  filter(time<=20,box!="05",mouseid!="01-02",mouseid!="02-03") -> data2

swirl2 <- data2 |>
  ggplot()+
  geom_path(aes(x=lagRBC,y=R,col=pABA),linewidth=2)+
  geom_text(aes(x=lagRBC,y=R,label=time),size=5)+
  scale_x_continuous(labels=aakmisc::scinot)+
  scale_y_continuous(labels=aakmisc::scinot)+
  scale_colour_manual(values=cbPalette[2:5])+
  labs(x=expression(paste("RBC density at time ", italic("t"), "-2 (density per µL)")),
       y=expression(paste("Reticulocyte supply at time ", italic("t"), " (density per µL)"))
  )+
  labs(colour="Parasite nutrient (pABA)")+
  theme_bw()+
  theme(
    axis.title=element_text(size=15),
    axis.title.y=element_text(size=15,colour="white"),
    strip.text=element_blank(),
    axis.text=element_text(size=11),
    panel.grid=element_blank(),
    legend.position="none",
    legend.title=element_text(size=15),
    legend.text=element_text(size=13),
    strip.background=element_blank(),
    plot.title=element_text(size=17,hjust=0.5)
    
  )

flow |> 
  select(time=day,box,pABA,E=Eryth,R=Retic,pABA,mouseid) |>
  mutate(RBC=E+R,
         lagRBC=lag(RBC,3)) |>
  na.omit() |>
  filter(time<=20,box!="05",mouseid!="01-02",mouseid!="02-03") -> data3

swirl3 <- data3 |>
  ggplot()+
  geom_path(aes(x=lagRBC,y=R,col=pABA),linewidth=2)+
  geom_text(aes(x=lagRBC,y=R,label=time),size=5)+
  scale_x_continuous(labels=aakmisc::scinot)+
  scale_y_continuous(labels=aakmisc::scinot)+
  scale_colour_manual(values=cbPalette[2:5])+
  labs(x=expression(paste("RBC density at time ", italic("t"), "-3 (density per µL)")),
       y=expression(paste("Reticulocyte supply at time ", italic("t"), " (density per µL)"))
  )+
  labs(colour="Parasite nutrient (pABA)")+
  theme_bw()+
  theme(
    axis.title=element_text(size=15),
    axis.title.y=element_text(size=15,colour="white"),
    strip.text=element_blank(),
    axis.text=element_text(size=11),
    panel.grid=element_blank(),
    legend.position=c(0.75,0.85),
    legend.title=element_text(size=13),
    legend.text=element_text(size=11),
    strip.background=element_blank(),
    plot.title=element_text(size=17,hjust=0.5)
    
  )

gt <- arrangeGrob(swirl1,swirl2,swirl3,nrow=1)

# Add labels to the arranged plots
p <- as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0,0.33,0.66), y = c(1,1,1))
p
ggsave("Figure3_RawData.jpeg",width=40,height=14,units="cm")
