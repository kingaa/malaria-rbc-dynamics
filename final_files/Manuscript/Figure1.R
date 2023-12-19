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

(REassump <- ggplot()+
    geom_line(aes(x=Erange,y=c(0.25,0.25)),col=cbPalette[6],linewidth=2)+
    geom_text(aes(x=Erange[1],y=0.29,label="Constant (1,3,6,7,8,10,11,12)"),col=cbPalette[6],hjust=0)+
    geom_line(aes(x=Erange,y=c(1,0)),col=cbPalette[4],linewidth=2)+
    geom_text(aes(x=Erange[1]+0.55,y=0.5,label="Linear (4,5,10,13,15,16)"),col=cbPalette[4],hjust=0)+
    geom_line(aes(x=seq(0,1,0.01),y=sapply(seq(0,1,0.01),expFun)),linewidth=2,col=cbPalette[5])+
    geom_text(aes(x=Erange[1]+0.2,y=0.99,label="Hill function (2)"),col=cbPalette[5],hjust=0)+
    ylim(0,1)+
    xlab("Erythrocyte density")+ylab("RBC supply")+
    theme_bw()+
    theme(axis.title=element_text(size=15),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid=element_blank())
)

(Nassump <- ggplot()+
    geom_smooth(aes(x=fig_df$time,y=fig_df$N),se=FALSE,col=cbPalette[2],linewidth=2)+
    geom_text(aes(x=24,y=quantile(fig_df$N,probs=0.85),label="Time-varying (8,9,14,17,18)"),col=cbPalette[2])+
    geom_line(aes(x=c(0,30),y=quantile(fig_df$N,probs=0.25)),col=cbPalette[6],linewidth=2)+
    geom_text(aes(x=15,y=quantile(fig_df$N,probs=0.35),label="Constant (2,12)"),col=cbPalette[6])+
    
    geom_line(aes(x=c(0,30),y=c(quantile(fig_df$N,probs=0.9),100000)),col=cbPalette[4],linewidth=2)+
    geom_text(aes(x=21,y=100000,label="Linear (1,3,4,6,7,10,11,13,15,16)"),col=cbPalette[4])+
    
    #geom_line(aes(x=c(0,30),y=0),col=cbPalette[7],linewidth=2)+
    #geom_text(aes(x=15,y=100000,label="None (8,12)"),col=cbPalette[7])+
    xlab("Day post-infection")+ylab("RBC clearance")+
    theme_bw()+
    theme(axis.title=element_text(size=15),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid=element_blank())
)

assump <- ggarrange(REassump,Nassump, nrow = 1, labels = c("A","B"))
assump

ggsave("Figure1.png",plot=assump,
       width=27,height=12,units="cm") 
