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
  filter(pABA=="0%",variable%in%c("E","R","N")) |>
  select(-lo,-hi) |>
  pivot_wider(names_from="variable",values_from="med") -> fig_df

Erange <- c(0,1)
expFun <- function(E){1/(1+exp(-10*(-E+0.5)))}

(REassump <- ggplot()+
    geom_line(aes(x=Erange,y=c(0.25,0.25)),col=cbPalette[6],linewidth=2)+
    geom_text(aes(x=Erange[1],y=0.29,label="Constant supply (a)"),col=cbPalette[6],hjust=0)+
    geom_line(aes(x=Erange,y=c(1,0)),col=cbPalette[4],linewidth=2)+
    geom_text(aes(x=Erange[1]+0.65,y=0.4,label="Linear supply (b)"),col=cbPalette[4],hjust=0)+
    geom_line(aes(x=seq(0,1,0.01),y=sapply(seq(0,1,0.01),expFun)),linewidth=2,col=cbPalette[5])+
    geom_text(aes(x=Erange[1]+0.2,y=0.99,label="Hill function supply (c)"),col=cbPalette[5],hjust=0)+
    ylim(0,1)+
    xlab("Erythrocytes")+ylab("Reticulocyte supply")+
    theme_bw()+
    theme(axis.title=element_text(size=15),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid=element_blank())
)

(Nassump <- ggplot()+
    geom_smooth(aes(x=fig_df$time,y=fig_df$N),se=FALSE,col=cbPalette[2],linewidth=2)+
    geom_text(aes(x=23,y=quantile(fig_df$N,probs=0.85),label="Time-varying clearance (b)"),col=cbPalette[2])+
    geom_line(aes(x=c(0,30),y=quantile(fig_df$N,probs=0.25)),col=cbPalette[1],linewidth=2)+
    geom_text(aes(x=15,y=quantile(fig_df$N,probs=0.35),label="Constant clearace (a)"),col=cbPalette[1])+
    #geom_path(aes(x=c(0,7,7,23,23,30),y=c(
    #quantile(fig_df$N,probs=0.5),
    #quantile(fig_df$N,probs=0.5),
    #quantile(fig_df$N,probs=0.75),
    #quantile(fig_df$N,probs=0.75),
    #quantile(fig_df$N,probs=0.5),
    #quantile(fig_df$N,probs=0.5)
    #                                     )),col=cbPalette[3],linewidth=2)+
    geom_line(aes(x=c(0,30),y=0),col=cbPalette[7],linewidth=2)+
    geom_text(aes(x=15,y=100000,label="No clearance (c)"),col=cbPalette[7])+
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
       width=25,height=12,units="cm") 
