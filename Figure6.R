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
stats_df <- read.csv("results_regression_stats.csv") |>
  filter(dataset=="sub") |>
  group_by(rep) |>
  filter(AICc==min(AICc)) |>
  ungroup() |>
  mutate(bp = X01) |>
  select(-c(X01:X04))
stats_df$bp[is.na(stats_df$bp)] <- "None"

preds_df <- read.csv("results_regression_quants.csv")
preds_df$bp[is.na(preds_df$bp)] <- "None"
preds_df$phase[is.na(preds_df$phase)] <- 1

preds_df$model <- factor(preds_df$model,levels=c("m1","m2","m3","m4","m5","m6"),
                         labels=c("Model A","Model B","Model C","Model D","Model E","Model F"))
preds_df$bp <- factor(preds_df$bp,levels=c(8,9,10,11,"None"),labels=c("Day 8","Day 9","Day 10","Day 11","No breakpoint"))
preds_df$lag <- factor(preds_df$lag,levels=c(1,2,3,4,5),labels=c("1-day","2-day","3-day","4-day","5-day"))
preds_df <- preds_df |>
  unite("label",c(model,bp,lag),sep=", ",remove=FALSE)

head(stats_df)
stats_df |>
  unite("label",c(model,bp,lag),sep=", ",remove=FALSE) |>
  select(rep,label) |>
  group_by(label) |>
  count() |> 
  arrange(-n) |>
  mutate(freq=n/1000) -> tmp

cumsum(tmp$n/1000)
##########################
#### 2 x 2 facet plot ####
##########################
top_pred_df <- preds_df |>
  filter(label=="Model B, Day 10, 3-day"|label=="Model B, Day 10, 2-day"|label=="Model E, Day 10, 3-day"|label=="Model B, Day 9, 2-day") |> 
  #separate_wider_delim(mouseid,delim="-",names=c("box","mouse")) |>
  mutate(pABA=case_match(box,
                         1~"High",
                         2~"Medium",
                         3~"Low",
                         4~"Unsupplemented"),
         facet_lab=case_match(label,
                              "Model B, Day 10, 3-day"~"Best model (22.2%)",
                              "Model B, Day 10, 2-day"~"Second-best model (10.2%)",
                              "Model E, Day 10, 3-day"~"Third-best model (9.9%)",
                              "Model B, Day 9, 2-day"~"Fourth-best model (8.9%)"))
top_pred_df$pABA <- factor(top_pred_df$pABA,levels=c("Unsupplemented","Low","Medium","High"))

top_pred_text <- data.frame(top_pred_df$label |> unique()) |> setNames("label") |>
  mutate(text=case_match(label,
                         "Model B, Day 10, 3-day"~"Model B\nDay 10 breakpoint\nLag (i) = 3 days",
                         "Model B, Day 10, 2-day"~"Model B\nDay 10 breakpoint\nLag (i) = 2 days",
                         "Model E, Day 10, 3-day"~"Model E\nDay 10 breakpoint\nLag (i) = 3 days",
                         "Model B, Day 9, 2-day"~"Model B\nDay 9 breakpoint\nLag (i) = 2 days"),
         facet_lab=case_match(label,
                              "Model B, Day 10, 3-day"~"Best model (22.2%)",
                              "Model B, Day 10, 2-day"~"Second-best model (10.2%)",
                              "Model E, Day 10, 3-day"~"Third-best model (9.9%)",
                              "Model B, Day 9, 2-day"~"Fourth-best model (8.9%)"),
         tag=case_match(label,
                        "Model B, Day 10, 3-day"~"A",
                        "Model B, Day 10, 2-day"~"B",
                        "Model E, Day 10, 3-day"~"C",
                        "Model B, Day 9, 2-day"~"D"))

top_pred_df$facet_lab <- factor(top_pred_df$facet_lab,levels=c("Best model (22.2%)",
                                                               "Second-best model (10.2%)",
                                                               "Third-best model (9.9%)",
                                                               "Fourth-best model (8.9%)"))
top_pred_text$facet_lab <- factor(top_pred_text$facet_lab,levels=c("Best model (22.2%)",
                                                                   "Second-best model (10.2%)",
                                                                   "Third-best model (9.9%)",
                                                                   "Fourth-best model (8.9%)"))

phase_lab <- data.frame(c("Phase 1","","",""),
           c("Phase 2","","",""),
           c("Best model (22.2%)",
             "Second-best model (10.2%)",
             "Third-best model (9.9%)",
             "Fourth-best model (8.9%)")) |> setNames(c("Phase1","Phase2","facet_lab"))
phase_lab$facet_lab <- factor(phase_lab$facet_lab,levels=c("Best model (22.2%)",
                                                           "Second-best model (10.2%)",
                                                           "Third-best model (9.9%)",
                                                           "Fourth-best model (8.9%)"))
  

top_pred_df |>
  ggplot()+
  geom_line(aes(x=lagRBC,y=med,col=pABA,linetype=factor(phase)),linewidth=1)+
  geom_ribbon(aes(x=lagRBC,ymin=lo,ymax=hi,fill=pABA,group=interaction(phase,pABA)),alpha=0.2)+
  geom_label(data=top_pred_text,aes(x=7500000,y=4200000,label=text))+
  geom_text(data=top_pred_text,aes(x=300000,y=4550000,label=tag),fontface="bold",size=7)+
  
  #geom_text(data=phase_lab,aes(x=7000000,y=3000000,label=Phase2),size=5)+
  #geom_text(data=phase_lab,aes(x=4000000,y=750000,label=Phase1),size=5)+
  
  scale_x_continuous(labels = aakmisc::scinot,limits=c(0,10000000))+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4700000))+
  facet_wrap(facet_lab~.,nrow=2)+
  scale_colour_manual(values=cbpABA[2:5])+
  scale_fill_manual(values=cbpABA[2:5])+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  labs(linetype="Phase",colour="Parasite nutrient (pABA)",fill="Parasite nutrient (pABA)",x=expression(paste("RBC density at time ", italic("t"), "-i (density per µL)")),
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
  filter(label=="Model B, Day 10, 3-day"|
           label=="Model B, Day 10, 2-day"|
           label=="Model E, Day 10, 3-day"|
           label=="Model B, Day 9, 2-day"|
           label=="Model E, Day 11, 4-day"|
           label=="Model B, Day 11, 4-day"|
           label=="Model E, Day 8, 4-day"|
           label=="Model E, Day 9, 3-day"|
           label=="Model B, Day 9, 3-day"
         ) |> 
  mutate(pABA=case_match(box,
                         1~"High",
                         2~"Medium",
                         3~"Low",
                         4~"Unsupplemented"),
         facet_lab=case_match(label,
                              "Model B, Day 10, 3-day"~"Best model (22.2%)",
                              "Model B, Day 10, 2-day"~"Second-best model (10.2%)",
                              "Model E, Day 10, 3-day"~"Third-best model (9.9%)",
                              "Model B, Day 9, 2-day"~"Fourth-best model (8.9%)",
                              "Model E, Day 11, 4-day"~"Fifth-model (7.6%)",
                              "Model B, Day 11, 4-day"~"Sixth-best model (6.2%)",
                              "Model E, Day 8, 4-day"~"Seventh-best model (5.0%)",
                              "Model E, Day 9, 3-day"~"Eighth-best model (4.6%)",
                              "Model B, Day 9, 3-day"~"Ninth-best model (3.0%)")
         )
top_pred_df$pABA <- factor(top_pred_df$pABA,levels=c("Unsupplemented","Low","Medium","High"))

top_pred_text <- data.frame(top_pred_df$label |> unique()) |> setNames("label") |>
  mutate(text=case_match(label,
                         "Model B, Day 10, 3-day"~"Model B\nDay 10 breakpoint\nLag (i) = 3 days",
                         "Model B, Day 10, 2-day"~"Model B\nDay 10 breakpoint\nLag (i) = 2 days",
                         "Model E, Day 10, 3-day"~"Model E\nDay 10 breakpoint\nLag (i) = 3 days",
                         "Model B, Day 9, 2-day"~"Model B\nDay 9 breakpoint\nLag (i) = 2 days",
                         "Model E, Day 11, 4-day"~"Model E\nDay 11 breakpoint\nLag (i) = 4 days",
                         "Model B, Day 11, 4-day"~"Model B\nDay 11 breakpoint\nLag (i) = 4 days",
                         "Model E, Day 8, 4-day"~"Model E\nDay 8 breakpoint\nLag (i) = 4 days",
                         "Model E, Day 9, 3-day"~"Model E\nDay 9 breakpoint\nLag (i) = 3 days",
                         "Model B, Day 9, 3-day"~"Model B\nDay 9 breakpoint\nLag (i) = 3 days"),
         facet_lab=case_match(label,
                              "Model B, Day 10, 3-day"~"Best model (22.2%)",
                              "Model B, Day 10, 2-day"~"Second-best model (10.2%)",
                              "Model E, Day 10, 3-day"~"Third-best model (9.9%)",
                              "Model B, Day 9, 2-day"~"Fourth-best model (8.9%)",
                              "Model E, Day 11, 4-day"~"Fifth-model (7.6%)",
                              "Model B, Day 11, 4-day"~"Sixth-best model (6.2%)",
                              "Model E, Day 8, 4-day"~"Seventh-best model (5.0%)",
                              "Model E, Day 9, 3-day"~"Eighth-best model (4.6%)",
                              "Model B, Day 9, 3-day"~"Ninth-best model (3.0%)"),
         tag=case_match(label,
                        "Model B, Day 10, 3-day"~"A",
                        "Model B, Day 10, 2-day"~"B",
                        "Model E, Day 10, 3-day"~"C",
                        "Model B, Day 9, 2-day"~"D",
                        "Model E, Day 11, 4-day"~"E",
                        "Model B, Day 11, 4-day"~"F",
                        "Model E, Day 8, 4-day"~"G",
                        "Model E, Day 9, 3-day"~"H",
                        "Model B, Day 9, 3-day"~"I"))

top_pred_df$facet_lab <- factor(top_pred_df$facet_lab,levels=c("Best model (22.2%)",
                                                               "Second-best model (10.2%)",
                                                               "Third-best model (9.9%)",
                                                               "Fourth-best model (8.9%)",
                                                               "Fifth-model (7.6%)",
                                                               "Sixth-best model (6.2%)",
                                                               "Seventh-best model (5.0%)",
                                                               "Eighth-best model (4.6%)",
                                                               "Ninth-best model (3.0%)"))
top_pred_text$facet_lab <- factor(top_pred_text$facet_lab,levels=c("Best model (22.2%)",
                                                                   "Second-best model (10.2%)",
                                                                   "Third-best model (9.9%)",
                                                                   "Fourth-best model (8.9%)",
                                                                   "Fifth-model (7.6%)",
                                                                   "Sixth-best model (6.2%)",
                                                                   "Seventh-best model (5.0%)",
                                                                   "Eighth-best model (4.6%)",
                                                                   "Ninth-best model (3.0%)"))

phase_lab <- data.frame(c("Phase 1","","","","","","","",""),
                        c("Phase 2","","","","","","","",""),
                        c("Best model (22.2%)",
                          "Second-best model (10.2%)",
                          "Third-best model (9.9%)",
                          "Fourth-best model (8.9%)",
                          "Fifth-model (7.6%)",
                          "Sixth-best model (6.2%)",
                          "Seventh-best model (5.0%)",
                          "Eighth-best model (4.6%)",
                          "Ninth-best model (3.0%)")) |> setNames(c("Phase1","Phase2","facet_lab"))
phase_lab$facet_lab <- factor(phase_lab$facet_lab,levels=c("Best model (22.2%)",
                                                           "Second-best model (10.2%)",
                                                           "Third-best model (9.9%)",
                                                           "Fourth-best model (8.9%)",
                                                           "Fifth-model (7.6%)",
                                                           "Sixth-best model (6.2%)",
                                                           "Seventh-best model (5.0%)",
                                                           "Eighth-best model (4.6%)",
                                                           "Ninth-best model (3.0%)"))


top_pred_df |>
  ggplot()+
  geom_line(aes(x=lagRBC,y=med,col=pABA,linetype=factor(phase)),linewidth=1)+
  geom_ribbon(aes(x=lagRBC,ymin=lo,ymax=hi,fill=pABA,group=interaction(phase,pABA)),alpha=0.2)+
  geom_label(data=top_pred_text,aes(x=7500000,y=4200000,label=text))+
  geom_text(data=top_pred_text,aes(x=300000,y=4550000,label=tag),fontface="bold",size=7)+
  
  #geom_text(data=phase_lab,aes(x=7000000,y=3000000,label=Phase2),size=5)+
  #geom_text(data=phase_lab,aes(x=4000000,y=750000,label=Phase1),size=5)+
  
  scale_x_continuous(labels = aakmisc::scinot,limits=c(0,10000000))+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4700000))+
  facet_wrap(facet_lab~.,nrow=3)+
  scale_colour_manual(values=cbpABA[2:5])+
  scale_fill_manual(values=cbpABA[2:5])+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  labs(linetype="Phase",colour="Parasite nutrient (pABA)",fill="Parasite nutrient (pABA)",x=expression(paste("RBC density at time ", italic("t"), "-i (density per µL)")),
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

ggsave("FigureS4.jpeg",width=25,height=25,units="cm")

##########################
#### 4 x 4 facet plot ####
##########################
top_pred_df <- preds_df |>
  select(rep,mouseid,time,lagRBC,phase,pred,b,model,lag) |>
  unite("label",c(model,b,lag),sep=", ",remove=FALSE) |>
  filter(label=="Model B, Day 10, 3-day"|
           label=="Model B, Day 10, 2-day"|
           label=="Model E, Day 10, 3-day"|
           label=="Model B, Day 9, 2-day"|
           label=="Model B, Day 11, 4-day"|
           label=="Model E, Day 11, 4-day"|
           label=="Model E, Day 9, 3-day"|
           label=="Model E, Day 8, 4-day"|
           label=="Model B, Day 9, 3-day"|
           label=="Model C, Day 10, 2-day"|
           label=="Model B, Day 11, 3-day"|
           label=="Model C, Day 9, 2-day"|
           label=="Model E, Day 10, 2-day"|
           label=="Model F, Day 10, 3-day"|
           label=="Model F, Day 9, 3-day"|
           label=="Model B, No breakpoint, 4-day"
  ) |> 
  separate_wider_delim(mouseid,delim="-",names=c("box","mouse")) |>
  mutate(pABA=case_match(box,
                         "01"~"High",
                         "02"~"Medium",
                         "03"~"Low",
                         "04"~"Unsupplemented"),
         facet_lab=case_match(label,
                              "Model B, Day 10, 3-day"~"Best model (23.0%)",
                              "Model B, Day 10, 2-day"~"Second-best model (11.9%)",
                              "Model E, Day 10, 3-day"~"Third-best model (10.6%)",
                              "Model B, Day 9, 2-day"~"Fourth-best model (10.4%)",
                              "Model B, Day 11, 4-day"~"Fifth-best model (7.1%)",
                              "Model E, Day 11, 4-day"~"Sixth-best model (6.8%)",
                              "Model E, Day 9, 3-day"~"Seventh-best model (6.1%)",
                              "Model E, Day 8, 4-day"~"Eighth-best model (4.8%)",
                              "Model B, Day 9, 3-day"~"Ninth-best model (3.1%)",
                              "Model C, Day 10, 2-day"~"Tenth-best model (1.9%)",
                              "Model B, Day 11, 3-day"~"Eleventh-best model (1.8%)",
                              "Model C, Day 9, 2-day"~"Twelfth-best model (1.8%)",
                              "Model E, Day 10, 2-day"~"Thirteenth-best model (1.8%)",
                              "Model F, Day 10, 3-day"~"Fourteenth-best model (1.6%)",
                              "Model F, Day 9, 3-day"~"Fifteenth-best model (1.2%)",
                              "Model B, No breakpoint, 4-day"~"Sixteenth-best model (1.1%)")
  )
top_pred_df$pABA <- factor(top_pred_df$pABA,levels=c("Unsupplemented","Low","Medium","High"))

top_pred_text <- data.frame(top_pred_df$label |> unique()) |> setNames("label") |>
  mutate(text=case_match(label,
                         "Model B, Day 10, 3-day"~"Model B\nDay 10 breakpoint\nLag (i) = 3 days",
                         "Model B, Day 10, 2-day"~"Model B\nDay 10 breakpoint\nLag (i) = 2 days",
                         "Model E, Day 10, 3-day"~"Model E\nDay 10 breakpoint\nLag (i) = 3 days",
                         "Model B, Day 9, 2-day"~"Model B\nDay 9 breakpoint\nLag (i) = 2 days",
                         "Model B, Day 11, 4-day"~"Model B\nDay 11 breakpoint\nLag (i) = 4 days",
                         "Model E, Day 11, 4-day"~"Model E\nDay 11 breakpoint\nLag (i) = 4 days",
                         "Model E, Day 9, 3-day"~"Model E\nDay 9 breakpoint\nLag (i) = 3 days",
                         "Model E, Day 8, 4-day"~"Model E\nDay 8 breakpoint\nLag (i) = 4 days",
                         "Model B, Day 9, 3-day"~"Model B\nDay 9 breakpoint\nLag (i) = 3 days",
                         "Model C, Day 10, 2-day"~"Model C\nDay 10 breakpoint\nLag (i) = 2 days",
                         "Model B, Day 11, 3-day"~"Model B\nDay 11 breakpoint\nLag (i) = 3 days",
                         "Model C, Day 9, 2-day"~"Model C\nDay 9 breakpoint\nLag (i) = 2 days",
                         "Model E, Day 10, 2-day"~"Model E\nDay 10 breakpoint\nLag (i) = 2 days",
                         "Model F, Day 10, 3-day"~"Model F\nDay 10 breakpoint\nLag (i) = 3 days",
                         "Model F, Day 9, 3-day"~"Model F\nDay 9 breakpoint\nLag (i) = 3 days",
                         "Model B, No breakpoint, 4-day"~"Model B\nNo breakpoint\nLag (i) = 4 days"),
         facet_lab=case_match(label,
                              "Model B, Day 10, 3-day"~"Best model (23.0%)",
                              "Model B, Day 10, 2-day"~"Second-best model (11.9%)",
                              "Model E, Day 10, 3-day"~"Third-best model (10.6%)",
                              "Model B, Day 9, 2-day"~"Fourth-best model (10.4%)",
                              "Model B, Day 11, 4-day"~"Fifth-best model (7.1%)",
                              "Model E, Day 11, 4-day"~"Sixth-best model (6.8%)",
                              "Model E, Day 9, 3-day"~"Seventh-best model (6.1%)",
                              "Model E, Day 8, 4-day"~"Eighth-best model (4.8%)",
                              "Model B, Day 9, 3-day"~"Ninth-best model (3.1%)",
                              "Model C, Day 10, 2-day"~"Tenth-best model (1.9%)",
                              "Model B, Day 11, 3-day"~"Eleventh-best model (1.8%)",
                              "Model C, Day 9, 2-day"~"Twelfth-best model (1.8%)",
                              "Model E, Day 10, 2-day"~"Thirteenth-best model (1.8%)",
                              "Model F, Day 10, 3-day"~"Fourteenth-best model (1.6%)",
                              "Model F, Day 9, 3-day"~"Fifteenth-best model (1.2%)",
                              "Model B, No breakpoint, 4-day"~"Sixteenth-best model (1.1%)"),
         tag=case_match(label,
                        "Model B, Day 10, 3-day"~"A",
                        "Model B, Day 10, 2-day"~"B",
                        "Model E, Day 10, 3-day"~"C",
                        "Model B, Day 9, 2-day"~"D",
                        "Model B, Day 11, 4-day"~"E",
                        "Model E, Day 11, 4-day"~"F",
                        "Model E, Day 9, 3-day"~"G",
                        "Model E, Day 8, 4-day"~"H",
                        "Model B, Day 9, 3-day"~"I",
                        "Model C, Day 10, 2-day"~"J",
                        "Model B, Day 11, 3-day"~"K",
                        "Model C, Day 9, 2-day"~"L",
                        "Model E, Day 10, 2-day"~"M",
                        "Model F, Day 10, 3-day"~"N",
                        "Model F, Day 9, 3-day"~"O",
                        "Model B, No breakpoint, 4-day"~"P"))

top_pred_df$facet_lab <- factor(top_pred_df$facet_lab,levels=c("Best model (23.0%)",
                                                               "Second-best model (11.9%)",
                                                               "Third-best model (10.6%)",
                                                               "Fourth-best model (10.4%)",
                                                               "Fifth-best model (7.1%)",
                                                               "Sixth-best model (6.8%)",
                                                               "Seventh-best model (6.1%)",
                                                               "Eighth-best model (4.8%)",
                                                               "Ninth-best model (3.1%)",
                                                               "Tenth-best model (1.9%)",
                                                               "Eleventh-best model (1.8%)",
                                                               "Twelfth-best model (1.8%)",
                                                               "Thirteenth-best model (1.8%)",
                                                               "Fourteenth-best model (1.6%)",
                                                               "Fifteenth-best model (1.2%)",
                                                               "Sixteenth-best model (1.1%)"))
top_pred_text$facet_lab <- factor(top_pred_text$facet_lab,levels=c("Best model (23.0%)",
                                                                   "Second-best model (11.9%)",
                                                                   "Third-best model (10.6%)",
                                                                   "Fourth-best model (10.4%)",
                                                                   "Fifth-best model (7.1%)",
                                                                   "Sixth-best model (6.8%)",
                                                                   "Seventh-best model (6.1%)",
                                                                   "Eighth-best model (4.8%)",
                                                                   "Ninth-best model (3.1%)",
                                                                   "Tenth-best model (1.9%)",
                                                                   "Eleventh-best model (1.8%)",
                                                                   "Twelfth-best model (1.8%)",
                                                                   "Thirteenth-best model (1.8%)",
                                                                   "Fourteenth-best model (1.6%)",
                                                                   "Fifteenth-best model (1.2%)",
                                                                   "Sixteenth-best model (1.1%)"))

phase_lab <- data.frame(c("Phase 1","","","","","","","","","","","","","","",""),
                        c("Phase 2","","","","","","","","","","","","","","",""),
                        c("Best model (23.0%)",
                          "Second-best model (11.9%)",
                          "Third-best model (10.6%)",
                          "Fourth-best model (10.4%)",
                          "Fifth-best model (7.1%)",
                          "Sixth-best model (6.8%)",
                          "Seventh-best model (6.1%)",
                          "Eighth-best model (4.8%)",
                          "Ninth-best model (3.1%)",
                          "Tenth-best model (1.9%)",
                          "Eleventh-best model (1.8%)",
                          "Twelfth-best model (1.8%)",
                          "Thirteenth-best model (1.8%)",
                          "Fourteenth-best model (1.6%)",
                          "Fifteenth-best model (1.2%)",
                          "Sixteenth-best model (1.1%)")) |> setNames(c("Phase1","Phase2","facet_lab"))
phase_lab$facet_lab <- factor(phase_lab$facet_lab,levels=c("Best model (23.0%)",
                                                           "Second-best model (11.9%)",
                                                           "Third-best model (10.6%)",
                                                           "Fourth-best model (10.4%)",
                                                           "Fifth-best model (7.1%)",
                                                           "Sixth-best model (6.8%)",
                                                           "Seventh-best model (6.1%)",
                                                           "Eighth-best model (4.8%)",
                                                           "Ninth-best model (3.1%)",
                                                           "Tenth-best model (1.9%)",
                                                           "Eleventh-best model (1.8%)",
                                                           "Twelfth-best model (1.8%)",
                                                           "Thirteenth-best model (1.8%)",
                                                           "Fourteenth-best model (1.6%)",
                                                           "Fifteenth-best model (1.2%)",
                                                           "Sixteenth-best model (1.1%)"))

top_pred_df |>
  ggplot(aes(x=lagRBC,y=pred))+
  geom_line(aes(x=lagRBC,y=pred,group=interaction(rep, mouse,phase,pABA),col=pABA),alpha=0.1)+
  geom_label(data=top_pred_text,aes(x=7500000,y=4200000,label=text),size=3)+
  geom_text(data=top_pred_text,aes(x=300000,y=4550000,label=tag),fontface="bold",size=5)+
  
  geom_text(data=phase_lab,aes(x=7000000,y=3000000,label=Phase2),size=3)+
  geom_text(data=phase_lab,aes(x=4000000,y=750000,label=Phase1),size=3)+
  
  scale_x_continuous(labels = aakmisc::scinot,limits=c(0,10000000))+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4600000))+
  facet_wrap(facet_lab~.,nrow=4)+
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
ggsave("Figure6_4x4.jpeg",width=25,height=25,units="cm")