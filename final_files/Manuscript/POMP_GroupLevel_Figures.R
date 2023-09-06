#Load required packages
library(tidyverse)
library(stringi)
library(doParallel)
library(pomp)
library(panelPomp)
library(ggpubr)
library(scales)
library(gridExtra)

#Set working directory
setwd("~/Documents/GitHub/bdd/nw11_hier/final_files/Manuscript/")

#Set colour palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Source data preparation script
source("POMP_GroupLevel_DataPrep.R")

#View group_traj data frame
head(group_traj)

###################################
#### Uninfected versus pABA 0% ####
###################################
group_UninfInf<-group_traj |>
  filter(pABA%in%c("Uninfected","0%"),
         variable%in%c("E","R","K","N","W","Qpn","Qps","Qun","Qpw","Qus")) |>
  mutate(type=case_match(pABA,
                         "Uninfected"~"Uninfected",
                         "0%"~"Infected"),
         variable_name=case_match(variable,
                                  "E"~"Erythrocytes",
                                  "R"~"Reticulocytes",
                                  "K"~"Parasites",
                                  "N"~"Indiscriminate killing",
                                  "W"~"Targeted killing",
                                  "Qpn"~"Qpn",
                                  "Qps"~"Qps",
                                  "Qun"~"Qun",
                                  "Qpw"~"Qpw",
                                  "Qus"~"Qus"))

group_UninfInf |>
  filter(variable%in%c("E","R","K","N","W")) |>
  ggplot()+
  geom_line(aes(x=time,y=med,col=type))+
  geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=type),alpha=0.2)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+ 
  scale_colour_manual(values=cbPalette[c(2,1)],name=NULL)+
  scale_fill_manual(values=cbPalette[c(2,1)],name=NULL)+
  xlab("Day post-infection")+ylab("Density per microlitre")+
  theme_bw()+
  facet_grid(variable_name~.,scales="free_y")+
  theme(strip.background=element_blank(),
        panel.grid = element_blank(),
        legend.position="right",
        axis.title=element_text(size=15),
        axis.text=element_text(size=13),
        legend.text=element_text(size=14)
  )

group_UninfInf |>
  filter(!variable%in%c("E","R","K","N","W")) |>
  ggplot()+
  geom_line(aes(x=time,y=med,col=type))+
  geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=type),alpha=0.2)+
  scale_colour_manual(values=cbPalette[c(2,1)],name=NULL)+
  scale_fill_manual(values=cbPalette[c(2,1)],name=NULL)+
  xlab("Day post-infection")+ylab("Density per microlitre")+
  theme_bw()+
  facet_grid(variable_name~.,scales="free_y")+
  theme(strip.background=element_blank(),
        panel.grid = element_blank(),
        legend.position="right",
        axis.title=element_text(size=15),
        axis.text=element_text(size=13),
        legend.text=element_text(size=14)
  )

group_UninfInf |>
  filter(variable%in%c("R","Qun")) |>
  select(-lo,-hi,-variable_name) |>
  pivot_wider(names_from="variable",values_from="med") |>
  ggplot()+
  geom_path(aes(x=Qun,y=R,col=type),arrow=arrow())+
  scale_colour_manual(values=cbPalette[c(2,1)],name=NULL)+
  xlab("Fraction uRBCs killed by indiscriminate killing")+ylab("Reticulocyte supply")+
  theme_bw()+
  theme(strip.background=element_blank(),
        panel.grid = element_blank(),
        legend.position="right",
        axis.title=element_text(size=15),
        axis.text=element_text(size=13),
        legend.text=element_text(size=14)
  )

group_UninfInf |>
  filter(variable%in%c("E","Qun")) |>
  select(-lo,-hi,-variable_name) |>
  pivot_wider(names_from="variable",values_from="med") |>
  ggplot()+
  geom_path(aes(x=Qun,y=E,col=type))+
  scale_colour_manual(values=cbPalette[c(2,1)],name=NULL)+
  xlab("Fraction uRBCs killed by indiscriminate killing")+ylab("Erythrocytes")+
  theme_bw()+
  theme(strip.background=element_blank(),
        panel.grid = element_blank(),
        legend.position="right",
        axis.title=element_text(size=15),
        axis.text=element_text(size=13),
        legend.text=element_text(size=14)
  )

group_UninfInf |>
  filter(variable%in%c("E","N")) |>
  select(-lo,-hi,-variable_name) |>
  pivot_wider(names_from="variable",values_from="med") |>
  ggplot()+
  geom_path(aes(x=E,y=N,col=type),arrow=arrow())+
  scale_colour_manual(values=cbPalette[c(2,1)],name=NULL)+
  xlab("Erythrocytes")+ylab("Indiscriminate killing")+
  theme_bw()+
  theme(strip.background=element_blank(),
        panel.grid = element_blank(),
        legend.position="right",
        axis.title=element_text(size=15),
        axis.text=element_text(size=13),
        legend.text=element_text(size=14)
  )

group_UninfInf |>
  filter(variable%in%c("R","N")) |>
  select(-lo,-hi,-variable_name) |>
  pivot_wider(names_from="variable",values_from="med") |>
  ggplot()+
  geom_path(aes(x=R,y=N,col=type),arrow=arrow())+
  scale_colour_manual(values=cbPalette[c(2,1)],name=NULL)+
  xlab("Reticulocyte supply")+ylab("Indiscriminate killing")+
  theme_bw()+
  theme(strip.background=element_blank(),
        panel.grid = element_blank(),
        legend.position="right",
        axis.title=element_text(size=15),
        axis.text=element_text(size=13),
        legend.text=element_text(size=14)
  )

###################################
#### Different pABA treatments ####
###################################
group_pABA<-group_traj |>
  filter(variable%in%c("E","R","K","N","W","Qpn","Qps","Qun","Qpw","Qus")) |>
  mutate(variable_name=case_match(variable,
                                  "E"~"Erythrocytes",
                                  "R"~"Reticulocytes",
                                  "K"~"Parasites",
                                  "N"~"Indiscriminate killing",
                                  "W"~"Targeted killing",
                                  "Qpn"~"Qpn",
                                  "Qps"~"Qps",
                                  "Qun"~"Qun",
                                  "Qpw"~"Qpw",
                                  "Qus"~"Qus"))

group_pABA |>
  filter(variable%in%c("E","R")) |>
  select(-lo,-hi,-variable_name) |>
  pivot_wider(names_from="variable",values_from="med") |>
  ggplot()+
  geom_path(aes(x=E,y=R,col=pABA))+
  facet_wrap(pABA~.,scales="free")+
  scale_colour_manual(values=cbPalette,name=NULL)+
  xlab("Erythrocytes")+ylab("Reticulocytes")+
  theme_bw()+
  theme(strip.background=element_blank(),
        panel.grid = element_blank(),
        legend.position="right",
        axis.title=element_text(size=15),
        axis.text=element_text(size=13),
        legend.text=element_text(size=14)
  )

group_pABA |>
  filter(variable%in%c("E","N")) |>
  select(-lo,-hi,-variable_name) |>
  pivot_wider(names_from="variable",values_from="med") |>
  ggplot()+
  geom_path(aes(x=E,y=N,col=pABA),arrow=arrow(type="closed"))+
  facet_wrap(pABA~.,scales="free")+
  scale_colour_manual(values=cbPalette,name=NULL)+
  xlab("Erythrocytes")+ylab("Indiscriminate killing")+
  theme_bw()+
  theme(strip.background=element_blank(),
        panel.grid = element_blank(),
        legend.position="right",
        axis.title=element_text(size=15),
        axis.text=element_text(size=13),
        legend.text=element_text(size=14)
  )

group_pABA |>
  filter(variable%in%c("Qun")) |>
  ggplot()+
  geom_line(aes(x=time,y=med,col=pABA))+
  geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
  scale_colour_manual(values=cbPalette,name=NULL)+
  scale_fill_manual(values=cbPalette,name=NULL)+
  xlab("Time")+ylab("Fraction uRBCs killed")+
  theme_bw()+
  theme(strip.background=element_blank(),
        panel.grid = element_blank(),
        legend.position="right",
        axis.title=element_text(size=15),
        axis.text=element_text(size=13),
        legend.text=element_text(size=14)
  )


