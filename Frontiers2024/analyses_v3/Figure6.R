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

##########################
#### 2 x 2 facet plot ####
##########################
top_pred_df <- preds_df |>
  select(rep,mouseid,time,lagRBC,phase,pred,b,model,lag) |>
  unite("label",c(model,b,lag),sep=", ",remove=FALSE) |>
  filter(label=="Model B, Day 10, 3-day"|label=="Model B, Day 10, 2-day"|label=="Model B, Day 9, 2-day"|label=="Model E, Day 10, 3-day") |>
  separate_wider_delim(mouseid,delim="-",names=c("box","mouse")) |>
  mutate(pABA=case_match(box,
                         "01"~"High",
                         "02"~"Medium",
                         "03"~"Low",
                         "04"~"Unsupplemented"),
         facet_lab=case_match(label,
                              "Model B, Day 10, 2-day"~"Best model (23.1%)",
                              "Model B, Day 10, 3-day"~"Second-best model (11.7%)",
                              "Model B, Day 9, 2-day"~"Third-best model (11.5%)",
                              "Model E, Day 10, 3-day"~"Fourth-best model (8.0%)"))
top_pred_df$pABA <- factor(top_pred_df$pABA,levels=c("Unsupplemented","Low","Medium","High"))

top_pred_text <- data.frame(top_pred_df$label |> unique()) |> setNames("label") |>
  mutate(text=case_match(label,
                         "Model B, Day 10, 2-day"~"Model B\nDay 10 breakpoint\nLag (i) = 2 days",
                         "Model B, Day 10, 3-day"~"Model B\nDay 10 breakpoint\nLag (i) = 3 days",
                         "Model B, Day 9, 2-day"~"Model B\nDay 9 breakpoint\nLag (i) = 2 days",
                         "Model E, Day 10, 3-day"~"Model E\nDay 10 breakpoint\nLag (i) = 3 days"),
         facet_lab=case_match(label,
                              "Model B, Day 10, 2-day"~"Best model (23.1%)",
                              "Model B, Day 10, 3-day"~"Second-best model (11.7%)",
                              "Model B, Day 9, 2-day"~"Third-best model (11.5%)",
                              "Model E, Day 10, 3-day"~"Fourth-best model (8.0%)"),
         tag=case_match(label,
                        "Model B, Day 10, 2-day"~"A",
                        "Model B, Day 10, 3-day"~"B",
                        "Model B, Day 9, 2-day"~"C",
                        "Model E, Day 10, 3-day"~"D"))

top_pred_df$facet_lab <- factor(top_pred_df$facet_lab,levels=c("Best model (23.1%)",
                                                               "Second-best model (11.7%)",
                                                               "Third-best model (11.5%)",
                                                               "Fourth-best model (8.0%)"))
top_pred_text$facet_lab <- factor(top_pred_text$facet_lab,levels=c("Best model (23.1%)",
                                                               "Second-best model (11.7%)",
                                                               "Third-best model (11.5%)",
                                                               "Fourth-best model (8.0%)"))

phase_lab <- data.frame(c("Phase 1","","",""),
           c("Phase 2","","",""),
           c("Best model (23.1%)",
             "Second-best model (11.7%)",
             "Third-best model (11.5%)",
             "Fourth-best model (8.0%)")) |> setNames(c("Phase1","Phase2","facet_lab"))
phase_lab$facet_lab <- factor(phase_lab$facet_lab,levels=c("Best model (23.1%)",
                                                                   "Second-best model (11.7%)",
                                                                   "Third-best model (11.5%)",
                                                                   "Fourth-best model (8.0%)"))


top_pred_df |>
  ggplot(aes(x=lagRBC,y=pred))+
  geom_line(aes(x=lagRBC,y=pred,group=interaction(rep, mouse,phase,pABA),col=pABA),alpha=0.1)+
  geom_label(data=top_pred_text,aes(x=7500000,y=4200000,label=text))+
  geom_text(data=top_pred_text,aes(x=300000,y=4550000,label=tag),fontface="bold",size=7)+
  
  geom_text(data=phase_lab,aes(x=7000000,y=3000000,label=Phase2),size=5)+
  geom_text(data=phase_lab,aes(x=4000000,y=750000,label=Phase1),size=5)+
  
  scale_x_continuous(labels = aakmisc::scinot,limits=c(0,10000000))+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4600000))+
  facet_wrap(facet_lab~.,nrow=2)+
  scale_colour_manual(values=cbpABA[2:5])+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  labs(colour="Parasite nutrient (pABA)",x=expression(paste("RBC density at time ", italic("t"), "-i (density per µL)")),
       y=expression(paste("Reticulocyte supply at time ", italic("t"), " (density per µL)")))+
  theme_bw()+
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=14),
    axis.text=element_text(size=12),
    axis.title=element_text(size=15),
    legend.position="top",
    legend.background=element_blank()
  )
ggsave("Figure6.jpeg",width=25,height=25,units="cm")

##########################
#### 3 x 3 facet plot ####
##########################
top_pred_df <- preds_df |>
  select(rep,mouseid,time,lagRBC,phase,pred,b,model,lag) |>
  unite("label",c(model,b,lag),sep=", ",remove=FALSE) |>
  filter(label=="Model B, Day 10, 3-day"|
           label=="Model B, Day 10, 2-day"|
           label=="Model B, Day 9, 2-day"|
           label=="Model E, Day 10, 3-day"|
           label=="Model E, Day 9, 3-day"|
           label=="Model B, Day 11, 4-day"|
           label=="Model E, Day 8, 4-day"|
           label=="Model E, Day 11, 4-day"|
           label=="Model B, Day 9, 3-day") |>
  separate_wider_delim(mouseid,delim="-",names=c("box","mouse")) |>
  mutate(pABA=case_match(box,
                         "01"~"High",
                         "02"~"Medium",
                         "03"~"Low",
                         "04"~"Unsupplemented"),
         facet_lab=case_match(label,
                              "Model B, Day 10, 2-day"~"Best model (23.1%)",
                              "Model B, Day 10, 3-day"~"Second-best model (11.7%)",
                              "Model B, Day 9, 2-day"~"Third-best model (11.5%)",
                              "Model E, Day 10, 3-day"~"Fourth-best model (8.0%)",
                              "Model E, Day 9, 3-day"~"Fifth-best model (7.3%)",
                              "Model B, Day 11, 4-day"~"Sixth-best model (6.3%)",
                              "Model E, Day 8, 4-day"~"Seventh-best model (6.0%)",
                              "Model E, Day 11, 4-day"~"Eighth-best model (5.8%)",
                              "Model B, Day 9, 3-day"~"Ninth-best model (4.2%)"))
top_pred_df$pABA <- factor(top_pred_df$pABA,levels=c("Unsupplemented","Low","Medium","High"))

top_pred_text <- data.frame(top_pred_df$label |> unique()) |> setNames("label") |>
  mutate(text=case_match(label,
                         "Model B, Day 10, 2-day"~"Model B\nDay 10 breakpoint\nLag (i) = 2 days",
                         "Model B, Day 10, 3-day"~"Model B\nDay 10 breakpoint\nLag (i) = 3 days",
                         "Model B, Day 9, 2-day"~"Model B\nDay 9 breakpoint\nLag (i) = 2 days",
                         "Model E, Day 10, 3-day"~"Model E\nDay 10 breakpoint\nLag (i) = 3 days",
                         "Model E, Day 9, 3-day"~"Model E\nDay 9 breakpoint\nLag (i) = 3 days",
                         "Model B, Day 11, 4-day"~"Model B\nDay 11 breakpoint\nLag (i) = 4 days",
                         "Model E, Day 8, 4-day"~"Model E\nDay 8 breakpoint\nLag (i) = 4 days",
                         "Model E, Day 11, 4-day"~"Model E\nDay 11 breakpoint\nLag (i) = 4 days",
                         "Model B, Day 9, 3-day"~"Model B\nDay 9 breakpoint\nLag (i) = 3 days"),
         facet_lab=case_match(label,
                              "Model B, Day 10, 2-day"~"Best model (23.1%)",
                              "Model B, Day 10, 3-day"~"Second-best model (11.7%)",
                              "Model B, Day 9, 2-day"~"Third-best model (11.5%)",
                              "Model E, Day 10, 3-day"~"Fourth-best model (8.0%)",
                              "Model E, Day 9, 3-day"~"Fifth-best model (7.3%)",
                              "Model B, Day 11, 4-day"~"Sixth-best model (6.3%)",
                              "Model E, Day 8, 4-day"~"Seventh-best model (6.0%)",
                              "Model E, Day 11, 4-day"~"Eighth-best model (5.8%)",
                              "Model B, Day 9, 3-day"~"Ninth-best model (4.2%)"),
         tag=case_match(label,
                        "Model B, Day 10, 2-day"~"A",
                        "Model B, Day 10, 3-day"~"B",
                        "Model B, Day 9, 2-day"~"C",
                        "Model E, Day 10, 3-day"~"D",
                        "Model E, Day 9, 3-day"~"E",
                        "Model B, Day 11, 4-day"~"F",
                        "Model E, Day 8, 4-day"~"G",
                        "Model E, Day 11, 4-day"~"H",
                        "Model B, Day 9, 3-day"~"I"))

top_pred_df$facet_lab <- factor(top_pred_df$facet_lab,levels=c("Best model (23.1%)",
                                                               "Second-best model (11.7%)",
                                                               "Third-best model (11.5%)",
                                                               "Fourth-best model (8.0%)",
                                                               "Fifth-best model (7.3%)",
                                                               "Sixth-best model (6.3%)",
                                                               "Seventh-best model (6.0%)",
                                                               "Eighth-best model (5.8%)",
                                                               "Ninth-best model (4.2%)"))
top_pred_text$facet_lab <- factor(top_pred_text$facet_lab,levels=c("Best model (23.1%)",
                                                                   "Second-best model (11.7%)",
                                                                   "Third-best model (11.5%)",
                                                                   "Fourth-best model (8.0%)",
                                                                   "Fifth-best model (7.3%)",
                                                                   "Sixth-best model (6.3%)",
                                                                   "Seventh-best model (6.0%)",
                                                                   "Eighth-best model (5.8%)",
                                                                   "Ninth-best model (4.2%)"))

phase_lab <- data.frame(c("Phase 1","","","","","","","",""),
                        c("Phase 2","","","","","","","",""),
                        c("Best model (23.1%)",
                                 "Second-best model (11.7%)",
                                 "Third-best model (11.5%)",
                                 "Fourth-best model (8.0%)",
                                 "Fifth-best model (7.3%)",
                                 "Sixth-best model (6.3%)",
                                 "Seventh-best model (6.0%)",
                                 "Eighth-best model (5.8%)",
                                 "Ninth-best model (4.2%)")) |> setNames(c("Phase1","Phase2","facet_lab"))
phase_lab$facet_lab <- factor(phase_lab$facet_lab,levels=c("Best model (23.1%)",
                                                                  "Second-best model (11.7%)",
                                                                  "Third-best model (11.5%)",
                                                                  "Fourth-best model (8.0%)",
                                                                  "Fifth-best model (7.3%)",
                                                                  "Sixth-best model (6.3%)",
                                                                  "Seventh-best model (6.0%)",
                                                                  "Eighth-best model (5.8%)",
                                                                  "Ninth-best model (4.2%)"))


top_pred_df |>
  ggplot(aes(x=lagRBC,y=pred))+
  geom_line(aes(x=lagRBC,y=pred,group=interaction(rep, mouse,phase,pABA),col=pABA),alpha=0.1)+
  geom_label(data=top_pred_text,aes(x=7500000,y=4200000,label=text),size=3)+
  geom_text(data=top_pred_text,aes(x=300000,y=4550000,label=tag),fontface="bold",size=5)+
  
  geom_text(data=phase_lab,aes(x=7000000,y=3000000,label=Phase2),size=3)+
  geom_text(data=phase_lab,aes(x=4000000,y=750000,label=Phase1),size=3)+
  
  scale_x_continuous(labels = aakmisc::scinot,limits=c(0,10000000))+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4600000))+
  facet_wrap(facet_lab~.,nrow=3)+
  scale_colour_manual(values=cbpABA[2:5])+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  labs(colour="Parasite nutrient (pABA)",x=expression(paste("RBC density at time ", italic("t"), "-i (density per µL)")),
       y=expression(paste("Reticulocyte supply at time ", italic("t"), " (density per µL)")))+
  theme_bw()+
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=10),
    axis.text=element_text(size=10),
    axis.title=element_text(size=12),
    legend.position="top",
    legend.background=element_blank()
  )

ggsave("Figure6_3x3.jpeg",width=25,height=25,units="cm")

top_pred_df |>
  filter(phase==1) |>
  ggplot(aes(x=lagRBC,y=pred))+
  geom_line(aes(x=lagRBC,y=pred,group=interaction(rep, mouse,phase,pABA),col=pABA),alpha=0.1)+
  geom_label(data=top_pred_text,aes(x=7500000,y=4200000,label=text),size=3)+
  geom_text(data=top_pred_text,aes(x=300000,y=4550000,label=tag),fontface="bold",size=5)+
  
  geom_text(data=phase_lab,aes(x=7000000,y=3000000,label=Phase2),size=3)+
  geom_text(data=phase_lab,aes(x=4000000,y=750000,label=Phase1),size=3)+
  
  scale_x_continuous(labels = aakmisc::scinot,limits=c(0,10000000))+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4600000))+
  facet_wrap(facet_lab~.,nrow=3)+
  scale_colour_manual(values=cbpABA[2:5])+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  labs(colour="Parasite nutrient (pABA)",x=expression(paste("RBC density at time ", italic("t"), "-i (density per µL)")),
       y=expression(paste("Reticulocyte supply at time ", italic("t"), " (density per µL)")))+
  theme_bw()+
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=10),
    axis.text=element_text(size=10),
    axis.title=element_text(size=12),
    legend.position="top",
    legend.background=element_blank()
  )

top_pred_df |>
  filter(phase==2) |>
  ggplot(aes(x=lagRBC,y=pred))+
  geom_line(aes(x=lagRBC,y=pred,group=interaction(rep, mouse,phase,pABA),col=pABA),alpha=0.1)+
  geom_label(data=top_pred_text,aes(x=7500000,y=4200000,label=text),size=3)+
  geom_text(data=top_pred_text,aes(x=300000,y=4550000,label=tag),fontface="bold",size=5)+
  
  geom_text(data=phase_lab,aes(x=7000000,y=3000000,label=Phase2),size=3)+
  geom_text(data=phase_lab,aes(x=4000000,y=750000,label=Phase1),size=3)+
  
  scale_x_continuous(labels = aakmisc::scinot,limits=c(0,10000000))+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4600000))+
  facet_wrap(facet_lab~.,nrow=3)+
  scale_colour_manual(values=cbpABA[2:5])+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  labs(colour="Parasite nutrient (pABA)",x=expression(paste("RBC density at time ", italic("t"), "-i (density per µL)")),
       y=expression(paste("Reticulocyte supply at time ", italic("t"), " (density per µL)")))+
  theme_bw()+
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=10),
    axis.text=element_text(size=10),
    axis.title=element_text(size=12),
    legend.position="top",
    legend.background=element_blank()
  )
