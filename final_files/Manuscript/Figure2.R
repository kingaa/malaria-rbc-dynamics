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
    xlab("")+ylab("RBC supply\n(density per µL)")+
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
    xlab("Day post-infection")+ylab("Erythrocytes\n(density per µL)")+
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
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(1,10^7.5))+ 
  theme_bw()+
  xlab("")+ylab("Parasites,\nK (density per µL)")+
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
    geom_text(aes(x=15,y=10^7,label="RBC supply"))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(10^5,10^7))+  
  theme_bw()+
  xlab("")+ylab("RBC supply,\nR (density per µL)")+
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
    xlab("Day post-infection")+ylab("Erythrocytes,\nE (density per µL)")+
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
    geom_text(aes(x=15,y=10^7,label="uRBC clearance"))+
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
    xlab("Day post-infection")+ylab("uRBC clearance,\nN (density per µL)")+
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
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(10^5,10^8))+  
    theme_bw()+
    xlab("")+ylab("iRBC clearance,\nW (density per µL)")+
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
  filter(mouseid!="01-02",mouseid!="02-03",time<=20) |> #remove underdosed mice
  pivot_wider(names_from=variable,values_from=value) |>
  separate_wider_delim(cols="mouseid",delim="-",names=c("box","mouse"),cols_remove=FALSE) |>
  filter(box!="05") |>
  group_by(mouseid) |>
  mutate(
    SM=exp(-M/(R+E)),
    SN=exp(-N/(R+E)),
    Qun=SM*(1-SN),
    lik=exp(loglik-max(loglik))
  ) |>
  ungroup()
sm1_sub$pABA <- factor(sm1_sub$box,levels=c("05","04","03","02","01"),
                       labels=c("Uninfected","Unsupplemented","Low","Medium","High"))

unsupplemented<-sm1_sub |>
  filter(pABA=="Unsupplemented") |>
  ggplot()+
  geom_line(aes(x=time,y=R,group=interaction(rep,mouse),col=lik,alpha=lik))+
  geom_text(aes(x=17,y=10^6.9,label="n=3"),size=5)+
  #geom_text(aes(x=1,y=10^6.9,label="i"),size=5)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(10^5,10^7))+
  scale_colour_gradient(low="lightgoldenrod1",high="orange4")+
  scale_alpha_continuous(guide="none")+
  ggtitle("Smooth trajectories (input)")+
  xlab("Day post-infection")+ylab("RBC supply (density per µL)")+
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

group <- group_traj |>
  filter(time<=20,pABA!="Uninfected",variable=="R") |>
  ggplot()+
  geom_line(aes(x=time,y=med,col=pABA),linewidth=2)+
  geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
  geom_text(aes(x=17,y=10^6.9,label="n=10"),size=5)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(10^5,10^7))+
  scale_colour_manual(values=cbPalette[2:5])+
  scale_fill_manual(values=cbPalette[2:5])+
  theme_bw()+
  ggtitle("Median trajectories (output)")+
  xlab("Day post-infection")+ylab("RBC supply (density per µL)")+
  labs(colour="Parasite nutrient (pABA)",fill="Parasite nutrient (pABA)")+
  theme(
    axis.title.y=element_text(size=13),
    axis.title.x=element_text(size=13),
    axis.text=element_text(size=13),
    plot.title=element_text(hjust=0.5,face="bold",size=13),
    panel.grid=element_blank(),
    legend.position=c(0.7,0.15),
    legend.title=element_text(size=10),
    legend.text=element_text(size=8),
    strip.background=element_blank(),
    strip.text=element_text(size=13),
    legend.background=element_blank()
  )

bottom_row <- ggarrange(unsupplemented, group, nrow = 1, labels = c("C","D"))

###########################
#### Overall plot #########
###########################
overall_plot <- ggarrange(top_columns,bottom_row, nrow=2, heights=c(1,1))

ggsave("Figure2.png",plot=overall_plot,
       width=30,height=40,units="cm")










Cii<-sm1_sub |>
    filter(pABA=="Low") |>
    ggplot()+
    geom_line(aes(x=time,y=R,group=interaction(rep,mouse),alpha=lik,col=lik))+
    geom_text(aes(x=17,y=10^6.9,label="n=3"),size=5)+
    geom_text(aes(x=1,y=10^6.9,label="ii"),size=5)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(10^5,10^7))+
    scale_colour_gradient(low="white",high=cbPalette[3])+
    theme_bw()+
    theme(
      axis.title=element_blank(),
      axis.text=element_text(size=13),
      plot.title=element_text(hjust=0.5,face="bold",size=13),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13),
      strip.background=element_blank(),
      strip.text=element_text(size=13)
    )

Ciii<-sm1_sub |>
  filter(pABA=="Medium") |>
  ggplot()+
  geom_line(aes(x=time,y=R,group=interaction(rep,mouse),alpha=lik,col=lik))+
  geom_text(aes(x=17,y=10^6.9,label="n=3"),size=5)+
  geom_text(aes(x=1,y=10^6.9,label="iii"),size=5)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(10^5,10^7))+
  scale_colour_gradient(low="white",high=cbPalette[4])+
  theme_bw()+
  theme(
    axis.title=element_blank(),
    axis.text=element_text(size=13),
    plot.title=element_text(hjust=0.5,face="bold",size=13),
    panel.grid=element_blank(),
    legend.position="none",
    legend.title=element_text(size=15),
    legend.text=element_text(size=13),
    strip.background=element_blank(),
    strip.text=element_text(size=13)
  )

Civ<-sm1_sub |>
  filter(pABA=="Low") |>
  ggplot()+
  geom_line(aes(x=time,y=R,group=interaction(rep,mouse),alpha=lik,col=lik))+
  geom_text(aes(x=17,y=10^6.9,label="n=3"),size=5)+
  geom_text(aes(x=1,y=10^6.9,label="iv"),size=5)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(10^5,10^7))+
  scale_colour_gradient(low="white",high=cbPalette[4])+
  theme_bw()+
  theme(
    axis.title=element_blank(),
    axis.text=element_text(size=13),
    plot.title=element_text(hjust=0.5,face="bold",size=13),
    panel.grid=element_blank(),
    legend.position="none",
    legend.title=element_text(size=15),
    legend.text=element_text(size=13),
    strip.background=element_blank(),
    strip.text=element_text(size=13)
  )

facet <- arrangeGrob(Ci,Cii,Ciii,Civ,ncol=2,nrow=2,
                   layout_matrix = cbind(c(1,3), c(2,4)),
                   top="Smooth trajectories (input)",
                   left="RBC supply (density per µL)",
                   bottom="Day post-infection")

group <- group_traj |>
  filter(time<=20,pABA!="Uninfected",variable=="R") |>
  ggplot()+
  geom_line(aes(x=time,y=med,col=pABA),linewidth=2)+
  geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
  geom_text(aes(x=17,y=10^6.9,label="n=10"),size=5)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(10^5,10^7))+
  scale_colour_manual(values=cbPalette[2:5])+
  scale_fill_manual(values=cbPalette[2:5])+
  theme_bw()+
  ggtitle("Median trajectories (output)")+
  xlab("Day post-infection")+ylab("RBC supply (density per µL)")+
  labs(colour="Parasite nutrient (pABA)",fill="Parasite nutrient (pABA)")+
  theme(
    axis.title.y=element_text(size=13),
    axis.title.x=element_text(size=13),
    axis.text=element_text(size=13),
    plot.title=element_text(hjust=0.5,face="bold",size=13),
    panel.grid=element_blank(),
    legend.position=c(0.7,0.15),
    legend.title=element_text(size=10),
    legend.text=element_text(size=8),
    strip.background=element_blank(),
    strip.text=element_text(size=13),
    legend.background=element_blank()
  )

bottom_row <- ggarrange(facet, group, nrow = 1, labels = c("C","D"))
bottom_row

###########################
#### Overall plot #########
###########################
overall_plot <- ggarrange(top_columns,bottom_row, nrow=2, heights=c(1,1))

ggsave("Figure2.png",plot=overall_plot,
       width=30,height=40,units="cm") 

########################Old 2C###############################################
ann_text <- as.data.frame(
  matrix(
    c(
      c("04","03","02","01"),
      c("n=3","n=3","n=2","n=2")
    ),
    nrow=4,ncol=2,byrow=FALSE
  )
)
names(ann_text) <- c("box","lab")
ann_text$pABA <- factor(ann_text$box,levels=c("05","04","03","02","01"),
                        labels=c("Uninfected","Unsupplemented","Low","Medium","High"))

num_text <- as.data.frame(
  matrix(
    c(
      c("04","03","02","01"),
      c("i","ii","iii","iv")
    ),
    nrow=4,ncol=2,byrow=FALSE
  )
)
names(num_text) <- c("box","lab")
num_text$pABA <- factor(num_text$box,levels=c("05","04","03","02","01"),
                        labels=c("Uninfected","Unsupplemented","Low","Medium","High"))

sm1_sub |>
  select(rep,time,box,mouseid,R,lik,pABA) |>
  pivot_wider(names_from=time,values_from=R) |>
  group_by(mouseid) |>
  sample_n(size=250,replace=FALSE) |>
  pivot_longer(cols=`0`:`20`,names_to="time",values_to="R")

facet <- sm1_sub |>
  select(rep,time,box,mouseid,R,lik,pABA) |>
  pivot_wider(names_from=time,values_from=R) |>
  group_by(mouseid) |>
  sample_n(size=250,replace=FALSE) |>
  pivot_longer(cols=`0`:`20`,names_to="time",values_to="R") |>
  group_by(box) |>
  arrange(lik, .by_group = TRUE) |>
  ggplot()+
  geom_line(aes(x=as.numeric(time),y=R,group=interaction(rep,mouseid),col=lik,alpha=lik))+
  geom_text(data=ann_text,aes(x=17,y=10^6.9,label=lab),size=5)+
  geom_text(data=num_text,aes(x=1,y=10^6.9,label=lab),size=5)+
  facet_wrap(pABA~.)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(10^5,10^7))+
  scale_colour_gradient(low="lightgrey",high="black")+
  theme_bw()+
  ggtitle("Smooth trajectories (input)")+
  xlab("Day post-infection")+ylab("RBC supply (density per µL)")+
  theme(
    axis.title.y=element_text(size=13),
    axis.title.x=element_text(size=13),
    axis.text=element_text(size=13),
    plot.title=element_text(hjust=0.5,face="bold",size=13),
    panel.grid=element_blank(),
    legend.position="none",
    legend.title=element_text(size=15),
    legend.text=element_text(size=13),
    strip.background=element_blank(),
    strip.text=element_text(size=13)
  )
facet

group <- group_traj |>
  filter(time<=20,pABA!="Uninfected",variable=="R") |>
  ggplot()+
  geom_line(aes(x=time,y=med,col=pABA),linewidth=2)+
  geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
  geom_text(aes(x=17,y=10^6.9,label="n=10"),size=5)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(10^5,10^7))+
  scale_colour_manual(values=cbPalette[2:5])+
  scale_fill_manual(values=cbPalette[2:5])+
  theme_bw()+
  ggtitle("Median trajectories (output)")+
  xlab("Day post-infection")+ylab("RBC supply (density per µL)")+
  labs(colour="Parasite nutrient (pABA)",fill="Parasite nutrient (pABA)")+
  theme(
    axis.title.y=element_text(size=13),
    axis.title.x=element_text(size=13),
    axis.text=element_text(size=13),
    plot.title=element_text(hjust=0.5,face="bold",size=13),
    panel.grid=element_blank(),
    legend.position=c(0.7,0.15),
    legend.title=element_text(size=10),
    legend.text=element_text(size=8),
    strip.background=element_blank(),
    strip.text=element_text(size=13),
    legend.background=element_blank()
  )

bottom_row <- ggarrange(facet, group, nrow = 1, labels = c("C","D"))
bottom_row

###########################
#### Overall plot #########
###########################
overall_plot <- ggarrange(top_columns,bottom_row, nrow=2, heights=c(1,1))
overall_plot

ggsave("Figure2.png",plot=overall_plot,
       width=30,height=40,units="cm") 

##########################################################################
#Original figure 2

(parasites <- flow |>
  filter(day<=20,box!="05",mouseid!="01-02",mouseid!="02-03") |>
  ggplot()+
  geom_line(aes(x=day,y=Pd,group=mouseid,col=pABA))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(1,10^7.5))+ 
  scale_colour_manual(values=cbPalette[2:5])+
  theme_bw()+
  xlab("")+ylab("Parasites\n(density per µL)")+
    ggtitle("Raw data")+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    axis.title.x=element_blank(),
    plot.title=element_text(hjust=0.5,face="bold",size=15),
    panel.grid=element_blank(),
    legend.position=c(0.35,0.3),
    legend.title=element_text(size=15),
    legend.text=element_text(size=13),
    legend.background=element_blank()
  )
)

(reticulocytes <- flow |>
  filter(day<=20,box!="05",mouseid!="01-02",mouseid!="02-03") |>
  ggplot()+
  geom_line(aes(x=day,y=Retic,group=mouseid,col=pABA))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(10^5,10^7))+  
  scale_colour_manual(values=cbPalette[2:5])+
  theme_bw()+
  xlab("")+ylab("Reticulocyte supply\n(density per µL)")+
    #ggtitle("")+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    axis.title.x=element_blank(),
    plot.title=element_text(hjust=0.5,face="bold",size=15),
    panel.grid=element_blank(),
    legend.position="none",
    legend.title=element_text(size=15),
    legend.text=element_text(size=13)
  )
)

(erythrocytes <- flow |>
    filter(day<=20,box!="05",mouseid!="01-02",mouseid!="02-03") |>
    ggplot()+
    geom_line(aes(x=day,y=Eryth,group=mouseid,col=pABA))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(10^5,10^7))+  
    scale_colour_manual(values=cbPalette[2:5])+
    theme_bw()+
    xlab("Day post-infection")+ylab("Erythrocytes\n(density per µL)")+
    #ggtitle("")+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=13),
      plot.title=element_text(hjust=0.5,face="bold",size=15),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )
)

data_column <- ggarrange(parasites,reticulocytes,erythrocytes, nrow = 3, labels = c("A","B","C"))

(K_model <- group_traj |>
  filter(box!="05",variable=="K",time<=20) |>
  mutate(lo=if_else(lo<=1,1,lo)) |>
  ggplot()+
  geom_line(aes(x=time,y=med,col=pABA),linewidth=2)+
  geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
              labels=trans_format('log10',math_format(10^.x)),
              limits=c(1,10^7.5))+ 
  scale_colour_manual(values=cbPalette[2:5])+
  scale_fill_manual(values=cbPalette[2:5])+
  ggtitle("pABA-level model output")+
  theme_bw()+
  xlab("")+ylab("")+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    axis.title.x=element_blank(),
    plot.title=element_text(size=15,face="bold",hjust=0.5),
    panel.grid=element_blank(),
    legend.position="none",
    legend.title=element_text(size=15),
    legend.text=element_text(size=13)
  )
)

(N_model <- group_traj |>
    filter(box!="05",variable=="N",time<=20) |>
    mutate(lo=if_else(lo<=1,1,lo)) |>
    ggplot()+
    geom_line(aes(x=time,y=med,col=pABA),linewidth=2)+
    geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)))+ 
    scale_colour_manual(values=cbPalette[2:5])+
    scale_fill_manual(values=cbPalette[2:5])+
    theme_bw()+
    xlab("")+ylab("")+
    #ggtitle("")+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=13),
      axis.title.x=element_blank(),
      plot.title=element_text(hjust=0.5,face="bold",size=15),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )
)

(R_model <- group_traj |>
    filter(box!="05",variable=="R",time<=20) |>
    mutate(lo=if_else(lo<=1,1,lo)) |>
    ggplot()+
    geom_line(aes(x=time,y=med,col=pABA),linewidth=2)+
    geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(10^5,10^7))+   
    scale_colour_manual(values=cbPalette[2:5])+
    scale_fill_manual(values=cbPalette[2:5])+
    theme_bw()+
    xlab("")+ylab("")+
    #ggtitle("")+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=13),
      plot.title=element_text(hjust=0.5,face="bold",size=15),
      panel.grid=element_blank(),
      axis.title.x=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )
)

(E_model <- group_traj |>
    filter(box!="05",variable=="E",time<=20) |>
    mutate(lo=if_else(lo<=1,1,lo)) |>
    ggplot()+
    geom_line(aes(x=time,y=med,col=pABA),linewidth=2)+
    geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)),
                  limits=c(10^5,10^7))+   
    scale_colour_manual(values=cbPalette[2:5])+
    scale_fill_manual(values=cbPalette[2:5])+
    theme_bw()+
    xlab("Day post-infection")+ylab("")+
    #ggtitle("")+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=13),
      plot.title=element_text(hjust=0.5,face="bold",size=15),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )
)

output_column <- ggarrange(K_model,R_model,E_model, nrow = 3,heights=c(1,1,1), labels = c("D","E","F"))

joint_columns <- ggarrange(data_column,output_column,ncol=2)
