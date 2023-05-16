library(tidyverse)
library(scales)

setwd("~/Documents/GitHub/bdd/nw11_hier/")

coef02 <- read_csv("est_coefs_02_box02.csv") |>
  mutate(box="02")
po02 <- read_csv("est_po_02_box02.csv") |>
  mutate(box="02")

coef03 <- read_csv("est_coefs_03.csv") |>
  mutate(box="03")
po03 <- read_csv("est_po_03.csv") |>
  mutate(box="03")

coef04 <- read_csv("est_coefs_04.csv") |>
  mutate(box="04")
po04 <- read_csv("est_po_04.csv") |>
  mutate(box="04")

coef_df <- bind_rows(coef02,coef03,coef04)
po_df <- bind_rows(po02,po03,po04) |>
  mutate(pABA=case_match(
    box,
    "01"~"0.05",
    "02"~"0.005",
    "03"~"0.0005",
    "04"~"0",
    "05"~"control"
  ))

po_df |> ggplot(aes(x=time,y=value,color=mouse,linewidth=mouse=="GROUP"))+
  geom_line()+
  scale_linewidth_manual(guide="none",values=c(`TRUE`=2,`FALSE`=1))+
  scale_y_log10()+
  xlab("Days")+ylab("Density per microlitre")+
  facet_grid(name~pABA,scales="free_y")+
  theme_bw()

po_df |>
  filter(mouse=="GROUP") |>
  ggplot(aes(x=time,y=value,color=pABA))+
  geom_line()+
  scale_y_log10()+
  xlab("Days")+ylab("Density per microlitre")+
  facet_grid(name~.,scales="free_y")+
  theme_bw()

  po_df |>
  filter(mouse=="GROUP") |>
  pivot_wider(names_from="name",values_from="value") |>
  ggplot(aes(x=logR,y=logW,color=pABA))+
  geom_path()+
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
  geom_path()+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  xlab("Indiscriminate killing")+ylab("Targeted killing")+
  theme_bw()

po_df |>
  filter(mouse=="GROUP") |>
  pivot_wider(names_from="name",values_from="value") |>
  ggplot(aes(x=time,y=logR/logW,color=pABA))+
  geom_line()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+
  xlab("Days")+ylab("Reticulocyte supply:Targeted killing")+
  theme_bw()
