#SET DIRECTORY TO SOURCE FILE LOCATION

library(tidyverse)
library(stringi)
library(doParallel)
library(pomp)
library(panelPomp)
library(ggpubr)
library(scales)
library(cowplot)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

source("POMP_GroupLevel_DataPrep.R")

###########################
#### A plots ##############
###########################
parasites <- flow |>
  filter(day<=20,mouseid=="03-01") |>
  ggplot()+
  geom_line(aes(x=day,y=Pd,group=mouseid))+
  geom_point(aes(x=day,y=Pd,group=mouseid),size=3)+
  geom_text(aes(x=18,y=10e6,label="n=1"))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(1,10^7.5))+ 
  theme_bw()+
  xlab("")+ylab("Parasites\n(density per µL)")+
  ggtitle("Data (input)")+
  theme(
    axis.title.y=element_text(size=11),
    axis.title.x=element_text(size=13),
    axis.text=element_text(size=13),
    plot.title=element_text(hjust=0.5,face="bold",size=13),
    panel.grid=element_blank(),
    legend.position=c(0.35,0.3),
    legend.title=element_text(size=15),
    legend.text=element_text(size=13),
    legend.background=element_blank()
  )

reticulocytes <- flow |>
    filter(day<=20,mouseid=="03-01") |>
    ggplot()+
    geom_line(aes(x=day,y=Retic,group=mouseid))+
    geom_point(aes(x=day,y=Retic,group=mouseid),size=3)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(10^5,10^7))+  
    theme_bw()+
    xlab("")+ylab("Reticulocytes\n(density per µL)")+
    #ggtitle("")+
    theme(
      axis.title.y=element_text(size=11),
      axis.title.x=element_text(size=13),
      axis.text=element_text(size=13),
      plot.title=element_text(hjust=0.5,face="bold",size=15),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )

erythrocytes <- flow |>
    filter(day<=20,mouseid=="03-01") |>
    ggplot()+
    geom_line(aes(x=day,y=Eryth,group=mouseid))+
    geom_point(aes(x=day,y=Eryth,group=mouseid),size=3)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(10^6,10^7))+  
    theme_bw()+
    xlab("Time (d post-infection)")+ylab("Erythrocytes\n(density per µL)")+
    #ggtitle("")+
    theme(
      axis.title.y=element_text(size=11),
      axis.title.x=element_text(size=13),
      axis.text=element_text(size=13),
      plot.title=element_text(hjust=0.5,face="bold",size=15),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )

A_column <- ggarrange(parasites,reticulocytes,erythrocytes, nrow = 3, labels = c("A i","  ii","  iii"))

###########################
#### B plots ##############
###########################
K <- sm1 |>
  filter(variable=="K",time<=20,mouseid=="03-01") |>
  ggplot()+
  geom_line(aes(x=time,y=value,group=rep),alpha=0.01)+
  geom_text(aes(x=15,y=10^7.5,label="Parasite density"))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(1,10^7.5))+ 
  theme_bw()+
  xlab("")+ylab("K (density per µL)")+
  ggtitle("Smooth trajectories (output)")+
  theme(
    axis.title.y=element_text(size=11),
    axis.title.x=element_text(size=13),
    axis.text=element_text(size=13),
    plot.title=element_text(hjust=0.5,face="bold",size=13),
    panel.grid=element_blank(),
    legend.position=c(0.35,0.3),
    legend.title=element_text(size=15),
    legend.text=element_text(size=13),
    legend.background=element_blank()
  )

R <- sm1 |>
  filter(variable=="R",time<=20,mouseid=="03-01") |>
  ggplot()+
  geom_line(aes(x=time,y=value,group=rep),alpha=0.01)+
    geom_text(aes(x=15,y=10^7,label="Reticulocyte supply"))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(10^5,10^7))+  
  theme_bw()+
  xlab("")+ylab("R (density per µL)")+
  #ggtitle("")+
  theme(
    axis.title.y=element_text(size=11),
    axis.title.x=element_text(size=13),
    axis.text=element_text(size=13),
    plot.title=element_text(hjust=0.5,face="bold",size=15),
    panel.grid=element_blank(),
    legend.position="none",
    legend.title=element_text(size=15),
    legend.text=element_text(size=13),
    #panel.border = element_rect(linewidth=4),
    panel.background = element_rect(fill="rosybrown1")
  )

E <-  sm1 |>
    filter(variable=="E",time<=20,mouseid=="03-01") |>
    ggplot()+
    geom_line(aes(x=time,y=value,group=rep),alpha=0.01)+
    geom_text(aes(x=15,y=10^7,label="Erythrocyte density"))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(10^6,10^7))+  
    theme_bw()+
    xlab("Time (d post-infection)")+ylab("E (density per µL)")+
    #ggtitle("")+
    theme(
      axis.title.y=element_text(size=11),
      axis.title.x=element_text(size=13),
      axis.text=element_text(size=13),
      plot.title=element_text(hjust=0.5,face="bold",size=15),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13),
      #panel.border = element_rect(linewidth=4),
      panel.background = element_rect(fill="rosybrown1")
    )

N <- sm1 |>
    filter(variable=="N",time<=20,mouseid=="03-01") |>
    ggplot()+
    geom_point(aes(x=time,y=value,group=rep,fill="Analyzed trajectories",colour="Analyzed trajectories"),shape = 22,size=8)+
    geom_text(aes(x=12,y=10^7,label="RBC clearance response"))+
    scale_fill_manual(name='',
                       breaks=c('Analyzed trajectories'),
                       values=c('Analyzed trajectories'='rosybrown1'))+
    scale_colour_manual(name='',
                      breaks=c('Analyzed trajectories'),
                      values=c('Analyzed trajectories'='rosybrown1'))+
    geom_line(aes(x=time,y=value,group=rep),alpha=0.01)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(10^5,10^7))+  
    theme_bw()+
    xlab("Time (d post-infection)")+ylab("N (density per µL)")+
    #ggtitle("")+
    theme(
      axis.title.y=element_text(size=11),
      axis.title.x=element_text(size=13),
      axis.text=element_text(size=13),
      plot.title=element_text(hjust=0.5,face="bold",size=15),
      panel.grid=element_blank(),
      legend.title=element_text(size=15),
      legend.text=element_text(size=20),
      #panel.border = element_rect(linewidth=4),
      panel.background = element_rect(fill="rosybrown1")
    )
legend <- get_legend(N)
N <- N+theme(legend.position="none")


W <- sm1 |>
    filter(variable=="W",time<=20,mouseid=="03-01") |>
    ggplot()+
    geom_line(aes(x=time,y=value,group=rep),alpha=0.01)+
    geom_text(aes(x=12,y=10^8,label="iRBC clearance response"))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(10^5,10^8))+  
    theme_bw()+
    xlab("")+ylab("W (density per µL)")+
    ggtitle("")+
    theme(
      axis.title.y=element_text(size=11),
      axis.title.x=element_text(size=13),
      axis.text=element_text(size=13),
      plot.title=element_text(hjust=0.5,face="bold",size=13),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )

B_column <- ggarrange(K,W,R,N,E,legend, nrow = 3, ncol=2, labels = c("B i"," iv"," ii"," v"," iii",""))

top_columns <- ggarrange(A_column,B_column,ncol=2,widths=c(0.5,1))

###########################
#### C plots ##############
###########################
sm1_sub <- sm1 |>
  as_tibble() |>
  filter(box=="04",variable=="R",time<=20) |>
  group_by(mouseid) |>
  mutate(
    lik=exp(loglik-max(loglik))
  ) |>
  ungroup() |>
  pivot_wider(names_from=time,values_from=value) |>
  sample_n(size=300) |>
  pivot_longer(cols=`0`:`20`,names_to="time",values_to="R")

sm1_sub$pABA <- factor(sm1_sub$box,levels=c("05","04","03","02","01"),
                       labels=c("Uninfected","Unsupplemented","Low","Medium","High"))
sm1_sub$time <- as.integer(sm1_sub$time)

unsupplemented<-sm1_sub |>
  ggplot()+
  geom_line(aes(x=time,y=R,group=interaction(rep,mouse),col=lik,alpha=lik))+
  geom_text(aes(x=17,y=10^6.9,label="n=3"),size=5)+
  #geom_text(aes(x=1,y=10^6.9,label="i"),size=5)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(10^5,10^7))+
  scale_colour_gradient(low="orange",high="orange4",limits=c(0,1))+
  scale_alpha_continuous(guide="none")+
  ggtitle("Smooth trajectories (input)")+
  xlab("Time (d post-infection)")+ylab("Reticulocyte supply (density per µL)")+
  labs(alpha="Relative likelihood",colour="Relative likelihood")+
  theme_bw()+
  theme(
    axis.title.y=element_text(size=13),
    axis.title.x=element_text(size=13),
    axis.text=element_text(size=13),
    plot.title=element_text(hjust=0.5,face="bold",size=13),
    panel.grid=element_blank(),
    legend.background=element_blank(),
    legend.position=c(0.7,0.15),
    legend.title=element_text(size=10),
    legend.text=element_text(size=8),
    strip.background=element_blank(),
    strip.text=element_text(size=13)
  )
unsupplemented

median <- group_traj |>
  filter(time<=20,pABA=="Unsupplemented",variable=="R") |>
  ggplot()+
  geom_line(aes(x=time,y=med,col=pABA),linewidth=2)+
  geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
  geom_text(aes(x=17,y=10^6.9,label="n=3"),size=5)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(10^5,10^7))+
  scale_colour_manual(values=cbPalette[2])+
  scale_fill_manual(values=cbPalette[2])+
  theme_bw()+
  ggtitle("Median trajectory (output)")+
  xlab("Time (d post-infection)")+ylab("")+
  labs(colour="Parasite nutrient (pABA)",fill="Parasite nutrient (pABA)")+
  theme(
    axis.title.y=element_text(size=13),
    axis.title.x=element_text(size=13),
    axis.text=element_text(size=13),
    plot.title=element_text(hjust=0.5,face="bold",size=13),
    panel.grid=element_blank(),
    legend.position="none",
    legend.title=element_text(size=10),
    legend.text=element_text(size=8),
    strip.background=element_blank(),
    strip.text=element_text(size=13),
    legend.background=element_blank()
  )

bottom_row <- ggarrange(unsupplemented, median, nrow = 1, labels = c("C","D"))

###########################
#### Overall plot #########
###########################
overall_plot <- ggarrange(top_columns,bottom_row, nrow=2, heights=c(1,1))

ggsave("Figure2.jpeg",plot=overall_plot,
       width=30,height=40,units="cm")
