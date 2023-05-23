library(tidyverse)
library(scales)

setwd("~/Documents/GitHub/bdd/nw11_hier/")

coef_df <- read_csv("est_coefs_joint.csv")
po_df <- read_csv("est_po_joint.csv") |>
  mutate(pABA=case_match(
    box,
    "01"~"0.05",
    "02"~"0.005",
    "03"~"0.0005",
    "04"~"0",
    "05"~"control"
  ))

po_df |> 
  ggplot(aes(x=time,y=value,color=mouse,linewidth=mouse=="GROUP"))+
  geom_line()+
  scale_linewidth_manual(guide="none",values=c(`TRUE`=2,`FALSE`=1))+
  scale_y_log10()+
  xlab("Days")+ylab("Density per microlitre")+
  facet_grid(name~pABA,scales="free_y")+
  theme_bw()

po_df |>
  filter(mouse=="GROUP") |>
  ggplot(aes(x=time,y=value,color=pABA))+
  geom_line(linewidth=2)+
  scale_y_log10()+
  xlab("Days")+ylab("Density per microlitre")+
  facet_grid(name~.,scales="free_y")+
  theme_bw()

po_df |>
  filter(mouse=="GROUP",name!="logE",name!="logM") |>
  ggplot(aes(x=time,y=value,color=name))+
  geom_line(linewidth=2)+
  scale_y_log10()+
  xlab("Days")+ylab("Density per microlitre")+
  facet_grid(.~pABA)+
  theme_bw()

po_df |>
  filter(mouse=="GROUP",name!="logE",name!="logM",time<=10) |>
  ggplot(aes(x=time,y=value,color=name))+
  geom_line(linewidth=2)+
  scale_y_log10()+
  xlab("Days")+ylab("Density per microlitre")+
  facet_grid(.~pABA)+
  theme_bw()

po_df |>
  filter(mouse=="GROUP") |>
  pivot_wider(names_from="name",values_from="value") |>
  ggplot(aes(x=logR,y=logW,color=pABA))+
  geom_path(linewidth=2)+
  geom_text(aes(label=time),col="black",size=5)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  xlab("Reticulocyte supply")+ylab("Targeted killing")+
  theme_bw()

po_df |>
  filter(mouse=="GROUP") |>
  pivot_wider(names_from="name",values_from="value") |>
  ggplot(aes(x=logN,y=logW,color=pABA))+
  geom_path(linewidth=2)+
  geom_text(aes(label=time),col="black",size=5)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  xlab("Indiscriminate killing")+ylab("Targeted killing")+
  theme_bw()

po_df |>
  filter(mouse=="GROUP") |>
  pivot_wider(names_from="name",values_from="value") |>
  ggplot(aes(x=logR,y=logN,color=pABA))+
  geom_path(linewidth=2)+
  geom_text(aes(label=time),col="black",size=5)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  xlab("Reticulocyte supply")+ylab("Indiscriminate killing")+
  theme_bw()



po_df |>
  filter(mouse=="GROUP",time<=6) |>
  pivot_wider(names_from="name",values_from="value") |>
  ggplot(aes(x=time,y=logR/logW,color=pABA))+
  geom_line()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  xlab("Days")+ylab("Reticulocyte supply:Targeted killing")+
  theme_bw()

po_df |>
  filter(mouse=="GROUP",time<=20) |>
  pivot_wider(names_from="name",values_from="value") |>
  ggplot(aes(x=logE,y=lag(logR),color=pABA))+
  geom_path(linewidth=2)+
  geom_text(aes(label=time),col="black",size=5)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  xlab("Erythrocytes")+ylab("Reticulocytes")+
  theme_bw()

po_df |>
  filter(mouse=="GROUP",time<=20) |>
  pivot_wider(names_from="name",values_from="value") |>
  ggplot(aes(x=logE,y=lead(logR),color=pABA))+
  geom_path(linewidth=2)+
  geom_text(aes(label=time),col="black",size=5)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  xlab("Erythrocytes")+ylab("Reticulocytes")+
  theme_bw()

po_df |>
  filter(mouse=="GROUP",time<=7) |>
  pivot_wider(names_from="name",values_from="value") |>
  ggplot(aes(x=logE,y=lead(logR),color=pABA))+
  geom_path(linewidth=2)+
  geom_text(aes(label=time),col="black",size=5)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  xlab("Erythrocytes")+ylab("Reticulocytes")+
  theme_bw()

po_df |>
  filter(mouse=="GROUP",time<=20) |>
  pivot_wider(names_from="name",values_from="value") |>
  ggplot(aes(x=logE,y=lead(logR)))+
  geom_point(aes(color=pABA))+
  xlab("Erythrocytes (t)")+ylab("Reticulocytes (t+1)")+
  theme_bw()


data_full |>
  filter(day<=20,paba!="0.05",paba!="control") |>
  ggplot()+
  geom_point(aes(x=day,y=Pd,col=mouse))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  facet_grid(mouse~paba)+
  theme_bw()

data_full |>
  filter(day<=20,paba!="0.05",paba!="control") |>
  ggplot()+
  geom_point(aes(x=day,y=Eryth,col=mouse))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  facet_grid(mouse~paba)+
  theme_bw()

data_full |>
  filter(day<=20) |>
  select(day,mouse,paba,Pd,Eryth,Retic) |>
  View()

po_df |>
  filter(mouse=="GROUP",time<=10) |>
  pivot_wider(names_from="name",values_from="value") |>
  ggplot(aes(x=lag(logE),y=logR,color=pABA))+
  geom_point()+
  facet_wrap(time~.,scales="free")+
  xlab("Erythrocytes (t)")+ylab("Reticulocytes (t+1)")+
  theme_bw()


