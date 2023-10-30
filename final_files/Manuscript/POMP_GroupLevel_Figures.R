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

#############################################################################
#### Data figures ####
#############################################################################

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

(nVStime0 <- group_traj |>
    filter(variable=="N",box=="04",time<20) |>
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
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )
)
ggsave("nVStime0.png",plot=nVStime0,
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

#############################################################################
#### Model relationship figures ####
#############################################################################
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
  geom_text(aes(x=Erange[1]+0.25,y=0.8,label="Linear supply (b)"),col=cbPalette[4],hjust=0,angle=-41)+
  geom_line(aes(x=seq(0,1,0.01),y=sapply(seq(0,1,0.01),expFun)),linewidth=2,col=cbPalette[5])+
  geom_text(aes(x=Erange[1]+0.2,y=0.99,label="Hill function supply (c)"),col=cbPalette[5],hjust=0)+
  ylim(0,1)+
  xlab("Erythrocytes (t-1)")+ylab("Reticulocyte supply (t)")+
  theme_bw()+
  theme(axis.title=element_text(size=15),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid=element_blank())
)
ggsave("REassump.png",plot=REassump,
       width=17,height=15,units="cm")



(Nassump <- ggplot()+
  geom_smooth(aes(x=fig_df$time,y=fig_df$N),se=FALSE,col=cbPalette[2],linewidth=2)+
  geom_text(aes(x=21,y=quantile(fig_df$N,probs=0.85),label="Time-varying killing (b)"),col=cbPalette[2])+
  geom_line(aes(x=c(0,30),y=quantile(fig_df$N,probs=0.25)),col=cbPalette[1],linewidth=2)+
  geom_text(aes(x=15,y=quantile(fig_df$N,probs=0.35),label="Constant killing (a)"),col=cbPalette[1])+
  #geom_path(aes(x=c(0,7,7,23,23,30),y=c(
    #quantile(fig_df$N,probs=0.5),
    #quantile(fig_df$N,probs=0.5),
    #quantile(fig_df$N,probs=0.75),
    #quantile(fig_df$N,probs=0.75),
    #quantile(fig_df$N,probs=0.5),
    #quantile(fig_df$N,probs=0.5)
     #                                     )),col=cbPalette[3],linewidth=2)+
  geom_line(aes(x=c(0,30),y=0),col=cbPalette[7],linewidth=2)+
  geom_text(aes(x=15,y=90000,label="No killing (c)"),col=cbPalette[7])+
  xlab("Day post-infection")+ylab("Indiscriminate killing")+
  theme_bw()+
  theme(axis.title=element_text(size=15),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid=element_blank())
)
ggsave("Nassump.png",plot=Nassump,
       width=17,height=15,units="cm")
