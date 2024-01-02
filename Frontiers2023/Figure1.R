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

group_traj |>
  filter(pABA=="Unsupplemented",variable%in%c("E","R","N")) |>
  select(-lo,-hi) |>
  pivot_wider(names_from="variable",values_from="med") -> fig_df

Erange <- c(0,1)
expFun <- function(E){1/(1+exp(-10*(-E+0.5)))}

(RE_assump <- ggplot()+
    
    geom_line(aes(x=Erange,y=c(0.25,0.25)),linetype="solid",linewidth=2)+
    geom_text(aes(x=Erange[1],y=0.29,label="Constant (1,3,6,7,8,10,11,12)"),hjust=0)+
    
    geom_line(aes(x=Erange,y=c(1,0)),linetype="longdash",linewidth=2)+
    geom_text(aes(x=Erange[1]+0.55,y=0.5,label="Linear (4,5,10,13,15,16)"),hjust=0)+
    
    geom_line(aes(x=seq(0,1,0.01),y=sapply(seq(0,1,0.01),expFun)),linewidth=2,linetype="dotted")+
    geom_text(aes(x=Erange[1]+0.2,y=0.99,label="Hill function (2)"),hjust=0)+
    
    scale_x_continuous(labels = c("deficit", "homeostatic\nequilibrium"), breaks = c(0.05, 0.95))+
    
    ggtitle("A")+
    
    ylim(0,1)+
    xlab("Erythrocyte density")+ylab("RBC supply")+
    theme_bw()+
    theme(axis.title=element_text(size=15),
          axis.text.x=element_text(size=12),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          panel.grid=element_blank(),
          plot.title=element_text(hjust=0,face="bold"))
)

(Qun_assump1 <- ggplot()+
    geom_line(aes(x=seq(0,1,0.01),y=sapply(seq(0,1,0.01),function(x){dnorm(x,0.5,0.35)-0.4})),linetype="dotdash",linewidth=2)+
    geom_text(aes(x=0.8,y=0.75,label="Time-varying (8,9,14,17,18)"))+
    
    geom_line(aes(x=c(0,1),y=0.2),linetype="solid",linewidth=2)+
    geom_text(aes(x=0.5,y=0.25,label="Constant (2,12)"))+
    
    scale_x_continuous(labels = c("", "\n"), breaks = c(0.05, 0.95))+
    
    ggtitle("B i")+

    xlab("Day post-infection")+ylab("RBC clearance")+
    theme_bw()+
    theme(axis.title=element_text(size=15),
          axis.text.x=element_text(size=12),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          panel.grid=element_blank(),
          plot.title=element_text(hjust=0,face="bold")))
)

(Qun_assump2 <- ggplot()+
    
    geom_line(aes(x=Erange,y=c(1,0)),linetype="longdash",linewidth=2)+
    geom_text(aes(x=Erange[1]+0.55,y=0.5,label="Linear (1,3,4,6,7,10,11,13,15,16)"),hjust=0)+
    
    scale_x_continuous(labels = c("deficit", "homeostatic\nequilibrium"), breaks = c(0.05, 0.95))+
    
    ggtitle("ii")+
    
    ylim(0,1)+
    xlab("Erythrocyte density")+ylab("")+
    theme_bw()+
    theme(axis.title=element_text(size=15),
          axis.text.x=element_text(size=12),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          panel.grid=element_blank(),
          plot.title=element_text(hjust=0,face="bold"))
)

assump <- ggarrange(RE_assump,Qun_assump1,Qun_assump2, nrow = 1)
assump

ggsave("Figure1.jpeg",plot=assump,
       width=45,height=15,units="cm") 
