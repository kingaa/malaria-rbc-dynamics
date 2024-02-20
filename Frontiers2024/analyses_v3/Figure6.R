library(tidyverse)
library(mgcv)
library(pomp)
library(foreach)
library(iterators)
library(doFuture)
library(aakmisc)
library(egg)
library(ggpubr)
library(cowplot)

#### Colour-blind friendly palettes ####
cbPalette <- c("#332288","#117733","#88CCEE","#DDCC77","#CC6677","#882255")
cbpABA <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#### Data preparation ####
preds_df <- read.csv("results_regression_preds.csv")

preds_df$model <- factor(preds_df$model,levels=c("m1","m2","m3","m4","m5","m6"),
                         labels=c("Model A","Model B","Model C","Model D","Model E","Model F"))
preds_df$b <- factor(preds_df$b,levels=c(8,9,10,11,"None"),labels=c("Day 8","Day 9","Day 10","Day 11","No breakpoint"))
preds_df$lag <- factor(preds_df$lag,levels=c(1,2,3,4),labels=c("1-day","2-day","3-day","4-day"))

top_pred_df <- preds_df |>
  select(rep,mouseid,time,lagRBC,phase,pred,b,model,lag) |>
  unite("label",c(model,b,lag),sep=", ",remove=FALSE) |>
  filter(label=="Model B, Day 10, 3-day"|label=="Model B, Day 10, 2-day"|label=="Model B, Day 9, 2-day") |>
  separate_wider_delim(mouseid,delim="-",names=c("box","mouse")) |>
  mutate(pABA=case_match(box,
                         "01"~"High",
                         "02"~"Medium",
                         "03"~"Low",
                         "04"~"Unsupplemented"),
         facet_lab=case_match(label,
                              "Model B, Day 10, 2-day"~"Best model (23.1%)",
                              "Model B, Day 10, 3-day"~"Second-best model (11.7%)",
                              "Model B, Day 9, 2-day"~"Third-best model (11.5%)"))
top_pred_df$pABA <- factor(top_pred_df$pABA,levels=c("Unsupplemented","Low","Medium","High"))

top_pred_text <- data.frame(top_pred_df$label |> unique()) |> setNames("label") |>
  mutate(text=case_match(label,
                         "Model B, Day 10, 2-day"~"Model B\nDay 10 breakpoint\nLag (i) = 2 days",
                         "Model B, Day 10, 3-day"~"Model B\nDay 10 breakpoint\nLag (i) = 3 days",
                         "Model B, Day 9, 2-day"~"Model B\nDay 9 breakpoint\nLag (i) = 2 days"),
         facet_lab=case_match(label,
                              "Model B, Day 10, 2-day"~"Best model (23.1%)",
                              "Model B, Day 10, 3-day"~"Second-best model (11.7%)",
                              "Model B, Day 9, 2-day"~"Third-best model (11.5%)"),
         tag=case_match(label,
                        "Model B, Day 10, 2-day"~"A",
                        "Model B, Day 10, 3-day"~"B",
                        "Model B, Day 9, 2-day"~"C"))


top_pred_df |>
  ggplot(aes(x=lagRBC,y=pred))+
  geom_line(aes(x=lagRBC,y=pred,group=interaction(rep, mouse,phase,pABA),col=pABA),alpha=0.1)+
  geom_label(data=top_pred_text,aes(x=7500000,y=4200000,label=text))+
  geom_text(data=top_pred_text,aes(x=300000,y=4550000,label=tag),fontface="bold",size=7)+
  scale_x_continuous(labels = aakmisc::scinot,limits=c(0,10000000))+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4600000))+
  facet_grid(.~facet_lab)+
  scale_colour_manual(values=cbpABA[2:5])+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  labs(colour="Parasite nutrient\n(pABA)",x=expression(paste("RBC density at time ", italic("t"), "-i (density per µL)")),
       y=expression(paste("Reticulocyte supply at time ", italic("t"), " (density per µL)")))+
  theme_bw()+
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=12),
    axis.text=element_text(size=10),
    axis.title=element_text(size=15),
    legend.position=c(0.07,0.2),
    legend.background=element_blank()
  )
ggsave("Figure6.jpeg",width=30,height=15,units="cm")

