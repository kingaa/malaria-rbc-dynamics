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

preds_df$model <- factor(preds_df$model,levels=c("m1","m2","m3","m4","m5","m6","m7","m8"),
                         labels=c("Model A","Model B","Model C","Model D","Model E","Model F","Model G","Model H"))
preds_df$bp <- factor(preds_df$bp,levels=c(8,9,10,11,"None"),labels=c("Day 8","Day 9","Day 10","Day 11","No breakpoint"))
preds_df$lag <- factor(preds_df$lag,levels=c(1,2,3,4,5),labels=c("1-day","2-day","3-day","4-day","5-day"))
preds_df <- preds_df |>
  unite("label",c(model,bp,lag),sep=", ",remove=FALSE)

head(stats_df)
stats_df |>
  unite("label",c(model,bp,lagChoice),sep=", ",remove=FALSE) |>
  select(rep,label) |>
  group_by(label) |>
  count() |> 
  arrange(-n) |>
  mutate(freq=n/1000) -> tmp

cumsum(tmp$n/1000)
##########################
#### Top 4 facet plot ####
##########################
top_pred_df <- preds_df |>
  filter(label=="Model B, Day 10, 3-day"|label=="Model E, Day 10, 3-day"|label=="Model B, Day 10, 2-day"|label=="Model B, Day 9, 2-day") |> 
  #separate_wider_delim(mouseid,delim="-",names=c("box","mouse")) |>
  mutate(pABA=case_match(box,
                         1~"High",
                         2~"Medium",
                         3~"Low",
                         4~"Unsupplemented"),
         facet_lab=case_match(label,
                              "Model B, Day 10, 3-day"~"First-best model (20.6%)",
                              "Model E, Day 10, 3-day"~"Second-best model (12.2%)",
                              "Model B, Day 10, 2-day"~"Third-best model (10.4%)",
                              "Model B, Day 9, 2-day"~"Fourth-best model (10.4%)"))
top_pred_df$pABA <- factor(top_pred_df$pABA,levels=c("Unsupplemented","Low","Medium","High"))

top_pred_text <- data.frame(top_pred_df$label |> unique()) |> setNames("label") |>
  mutate(text=case_match(label,
                         "Model B, Day 10, 3-day"~"Model B\nBpt = d10\nLag (i) = 3d",
                         "Model E, Day 10, 3-day"~"Model E\nBpt = d10\nLag (i) = 3d",
                         "Model B, Day 10, 2-day"~"Model B\nBpt = d10\nLag (i) = 2d",
                         "Model B, Day 9, 2-day"~"Model B\nBpt = d9\nLag (i) = 2d"),
         facet_lab=case_match(label,
                              "Model B, Day 10, 3-day"~"First-best model (20.6%)",
                              "Model E, Day 10, 3-day"~"Second-best model (12.2%)",
                              "Model B, Day 10, 2-day"~"Third-best model (10.4%)",
                              "Model B, Day 9, 2-day"~"Fourth-best model (10.4%)"),
         tag=case_match(label,
                        "Model B, Day 10, 3-day"~"A",
                        "Model E, Day 10, 3-day"~"B",
                        "Model B, Day 10, 2-day"~"C",
                        "Model B, Day 9, 2-day"~"D"))

top_pred_df$facet_lab <- factor(top_pred_df$facet_lab,levels=c("First-best model (20.6%)",
                                                               "Second-best model (12.2%)",
                                                               "Third-best model (10.4%)",
                                                               "Fourth-best model (10.4%)"))
top_pred_text$facet_lab <- factor(top_pred_text$facet_lab,levels=c("First-best model (20.6%)",
                                                                   "Second-best model (12.2%)",
                                                                   "Third-best model (10.4%)",
                                                                   "Fourth-best model (10.4%)"))

phase_lab <- data.frame(c("Phase 1","","",""),
           c("Phase 2","","",""),
           c("First-best model (20.6%)",
             "Second-best model (12.2%)",
             "Third-best model (10.4%)",
             "Fourth-best model (10.4%)")) |> setNames(c("Phase1","Phase2","facet_lab"))
phase_lab$facet_lab <- factor(phase_lab$facet_lab,levels=c("First-best model (20.6%)",
                                                           "Second-best model (12.2%)",
                                                           "Third-best model (10.4%)",
                                                           "Fourth-best model (10.4%)"))
  

top4 <- top_pred_df |>
  ggplot()+
  geom_line(aes(x=lagRBC,y=med,col=pABA,linetype=factor(phase)),linewidth=1)+
  geom_ribbon(aes(x=lagRBC,ymin=lo,ymax=hi,fill=pABA,group=interaction(phase,pABA)),alpha=0.2)+
  geom_label(data=top_pred_text,aes(x=7700000,y=3700000,label=text))+
  #geom_text(data=top_pred_text,aes(x=300000,y=4550000,label=tag),fontface="bold",size=7)+
  
  #geom_text(data=phase_lab,aes(x=7000000,y=3000000,label=Phase2),size=5)+
  #geom_text(data=phase_lab,aes(x=4000000,y=750000,label=Phase1),size=5)+
  
  scale_x_continuous(breaks=c(2000000,4000000,6000000,8000000),labels = aakmisc::scinot,limits=c(1500000,9000000))+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4100000))+
  facet_wrap(facet_lab~.,nrow=2)+
  scale_colour_manual(values=cbpABA[2:5])+
  scale_fill_manual(values=cbpABA[2:5])+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  labs(linetype="Phase",colour="Parasite nutrient\n(pABA)",fill="Parasite nutrient\n(pABA)",x=expression(paste("RBC density at time ", italic("t"), "-i (density per µL)")),
       y=expression(paste("Reticulocyte supply at time ", italic("t"), " (density per µL)")))+
  theme_bw()+
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=14),
    axis.text=element_text(size=12),
    axis.title=element_text(size=15),
    legend.position="right",
    legend.background=element_blank()
  )
top4

ggdraw(top4)+draw_plot_label(label = c("A", "B","C","D"), size = 15,
                              x = c(0.1,0.47,0.1,0.47), y = c(0.99,0.99,0.53,0.53))

ggsave("Figure5.jpeg",width=25,height=22,units="cm")

##########################
#### Top nine facet plot ####
##########################
top_pred_df <- preds_df |>
  filter(label=="Model B, Day 10, 3-day"|
           label=="Model E, Day 10, 3-day"|
           label=="Model B, Day 10, 2-day"|
           label=="Model B, Day 9, 2-day"|
           label=="Model E, Day 11, 4-day"|
           label=="Model B, Day 11, 4-day"|
           label=="Model E, Day 9, 3-day"|
           label=="Model E, Day 8, 4-day"|
           label=="Model B, Day 9, 3-day"
         ) |> 
  mutate(pABA=case_match(box,
                         1~"High",
                         2~"Medium",
                         3~"Low",
                         4~"Unsupplemented"),
         facet_lab=case_match(label,
                              "Model B, Day 10, 3-day"~"Best model (20.6%)",
                              "Model E, Day 10, 3-day"~"Second-best model (12.2%)",
                              "Model B, Day 10, 2-day"~"Third-best model (10.4%)",
                              "Model B, Day 9, 2-day"~"Fourth-best model (10.4%)",
                              "Model E, Day 11, 4-day"~"Fifth-model (6.2%)",
                              "Model B, Day 11, 4-day"~"Sixth-best model (5.2%)",
                              "Model E, Day 9, 3-day"~"Seventh-best model (3.9%)",
                              "Model E, Day 8, 4-day"~"Eighth-best model (3.7%)",
                              "Model B, Day 9, 3-day"~"Ninth-best model (3.6%)")
         )
top_pred_df$pABA <- factor(top_pred_df$pABA,levels=c("Unsupplemented","Low","Medium","High"))

top_pred_text <- data.frame(top_pred_df$label |> unique()) |> setNames("label") |>
  mutate(text=case_match(label,
                         "Model B, Day 10, 3-day"~"Model B\nDay 10 breakpoint\nLag (i) = 3 days",
                         "Model E, Day 10, 3-day"~"Model E\nDay 10 breakpoint\nLag (i) = 3 days",
                         "Model B, Day 10, 2-day"~"Model B\nDay 10 breakpoint\nLag (i) = 2 days",
                         "Model B, Day 9, 2-day"~"Model B\nDay 9 breakpoint\nLag (i) = 2 days",
                         "Model E, Day 11, 4-day"~"Model E\nDay 11 breakpoint\nLag (i) = 4 days",
                         "Model B, Day 11, 4-day"~"Model B\nDay 11 breakpoint\nLag (i) = 4 days",
                         "Model E, Day 9, 3-day"~"Model E\nDay 9 breakpoint\nLag (i) = 3 days",
                         "Model E, Day 8, 4-day"~"Model E\nDay 8 breakpoint\nLag (i) = 4 days",
                         "Model B, Day 9, 3-day"~"Model B\nDay 9 breakpoint\nLag (i) = 3 days"),
         facet_lab=case_match(label,
                              "Model B, Day 10, 3-day"~"Best model (20.6%)",
                              "Model E, Day 10, 3-day"~"Second-best model (12.2%)",
                              "Model B, Day 10, 2-day"~"Third-best model (10.4%)",
                              "Model B, Day 9, 2-day"~"Fourth-best model (10.4%)",
                              "Model E, Day 11, 4-day"~"Fifth-model (6.2%)",
                              "Model B, Day 11, 4-day"~"Sixth-best model (5.2%)",
                              "Model E, Day 9, 3-day"~"Seventh-best model (3.9%)",
                              "Model E, Day 8, 4-day"~"Eighth-best model (3.7%)",
                              "Model B, Day 9, 3-day"~"Ninth-best model (3.6%)"),
         tag=case_match(label,
                        "Model B, Day 10, 3-day"~"A",
                        "Model E, Day 10, 3-day"~"B",
                        "Model B, Day 10, 2-day"~"C",
                        "Model B, Day 9, 2-day"~"D",
                        "Model E, Day 11, 4-day"~"E",
                        "Model B, Day 11, 4-day"~"F",
                        "Model E, Day 9, 3-day"~"G",
                        "Model E, Day 8, 4-day"~"H",
                        "Model B, Day 9, 3-day"~"I"))

top_pred_df$facet_lab <- factor(top_pred_df$facet_lab,levels=c("Best model (20.6%)",
                                                               "Second-best model (12.2%)",
                                                               "Third-best model (10.4%)",
                                                               "Fourth-best model (10.4%)",
                                                               "Fifth-model (6.2%)",
                                                               "Sixth-best model (5.2%)",
                                                               "Seventh-best model (3.9%)",
                                                               "Eighth-best model (3.7%)",
                                                               "Ninth-best model (3.6%)"))
top_pred_text$facet_lab <- factor(top_pred_text$facet_lab,levels=c("Best model (20.6%)",
                                                                   "Second-best model (12.2%)",
                                                                   "Third-best model (10.4%)",
                                                                   "Fourth-best model (10.4%)",
                                                                   "Fifth-model (6.2%)",
                                                                   "Sixth-best model (5.2%)",
                                                                   "Seventh-best model (3.9%)",
                                                                   "Eighth-best model (3.7%)",
                                                                   "Ninth-best model (3.6%)"))

phase_lab <- data.frame(c("Phase 1","","","","","","","",""),
                        c("Phase 2","","","","","","","",""),
                        c("Best model (20.6%)",
                          "Second-best model (12.2%)",
                          "Third-best model (10.4%)",
                          "Fourth-best model (10.4%)",
                          "Fifth-model (6.2%)",
                          "Sixth-best model (5.2%)",
                          "Seventh-best model (3.9%)",
                          "Eighth-best model (3.7%)",
                          "Ninth-best model (3.6%)")) |> setNames(c("Phase1","Phase2","facet_lab"))
phase_lab$facet_lab <- factor(phase_lab$facet_lab,levels=c("Best model (20.6%)",
                                                           "Second-best model (12.2%)",
                                                           "Third-best model (10.4%)",
                                                           "Fourth-best model (10.4%)",
                                                           "Fifth-model (6.2%)",
                                                           "Sixth-best model (5.2%)",
                                                           "Seventh-best model (3.9%)",
                                                           "Eighth-best model (3.7%)",
                                                           "Ninth-best model (3.6%)"))


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

ggsave("FigureS5.jpeg",width=25,height=25,units="cm")

##########################
#### 10-24 facet plot ####
##########################
top_pred_df <- preds_df |>
  filter(label=="Model C, Day 10, 2-day"|
           label=="Model E, Day 10, 2-day"|
           label=="Model C, Day 9, 2-day"|
           label=="Model B, No breakpoint, 4-day"|
           label=="Model F, Day 9, 3-day"|
           label=="Model F, Day 10, 2-day"|
           label=="Model F, Day 10, 3-day"|
           
           label=="Model B, Day 11, 3-day"|
           label=="Model F, Day 8, 4-day"|
           label=="Model A, No breakpoint, 4-day"|
           label=="Model E, Day 11, 3-day"|
           label=="Model F, Day 11, 4-day"|
           label=="Model B, Day 9, 1-day"|
           label=="Model C, Day 10, 3-day"|
           label=="Model E, Day 9, 4-day"
           
  ) |> 
  mutate(pABA=case_match(box,
                         1~"High",
                         2~"Medium",
                         3~"Low",
                         4~"Unsupplemented"),
         facet_lab=case_match(label,
                              "Model C, Day 10, 2-day"~"Tenth-best model (3.2%)",
                              "Model E, Day 10, 2-day"~"Eleventh-best model (2.9%)",
                              "Model C, Day 9, 2-day"~"Twelfth-best model (2.5%)",
                              "Model B, No breakpoint, 4-day"~"Thirteenth-best model (2.3%)",
                              "Model F, Day 9, 3-day"~"Fourteenth-best model (2.2%)",
                              "Model F, Day 10, 2-day"~"Fifteenth-best model (2.0%)",
                              "Model F, Day 10, 3-day"~"Sixteenth-best model (1.9%)",
                              
                              "Model B, Day 11, 3-day"~"Seventeenth-best model (1.4%)",
                              "Model F, Day 8, 4-day"~"Eighteenth-best model (1.1%)",
                              "Model A, No breakpoint, 4-day"~"Nineteenth-best model (0.9%)",
                              "Model E, Day 11, 3-day"~"Twentieth-best model (0.7%)",
                              "Model F, Day 11, 4-day"~"Twenty-first-best model (0.7%)",
                              "Model B, Day 9, 1-day"~"Twenty-second-best model (0.3%)",
                              "Model C, Day 10, 3-day"~"Twenty-third-best model (0.3%)",
                              "Model E, Day 9, 4-day"~"Twenty-fourth-best model (0.3%)")
  )
top_pred_df$pABA <- factor(top_pred_df$pABA,levels=c("Unsupplemented","Low","Medium","High"))

top_pred_text <- data.frame(top_pred_df$label |> unique()) |> setNames("label") |>
  mutate(text=case_match(label,
                         "Model C, Day 10, 2-day"~"Model C\nDay 10 breakpoint\nLag (i) = 2 days",
                         "Model E, Day 10, 2-day"~"Model E\nDay 10 breakpoint\nLag (i) = 2 days",
                         "Model C, Day 9, 2-day"~"Model C\nDay 9 breakpoint\nLag (i) = 2 days",
                         "Model B, No breakpoint, 4-day"~"Model B\nNo breakpoint\nLag (i) = 4 days",
                         "Model F, Day 9, 3-day"~"Model F\nDay 9 breakpoint\nLag (i) = 3 days",
                         "Model F, Day 10, 2-day"~"Model F\nDay 10 breakpoint\nLag (i) = 2 days",
                         "Model F, Day 10, 3-day"~"Model F\nDay 10 breakpoint\nLag (i) = 3 days",
                         
                         "Model B, Day 11, 3-day"~"Model B\nDay 11 breakpoint\nLag (i) = 3 days",
                         "Model F, Day 8, 4-day"~"Model F\nDay 8 breakpoint\nLag (i) = 4 days",
                         "Model A, No breakpoint, 4-day"~"Model A\nNo breakpoint\nLag (i) = 4 days",
                         "Model E, Day 11, 3-day"~"Model E\nDay 11 breakpoint\nLag (i) = 3 days",
                         "Model F, Day 11, 4-day"~"Model F\nDay 11 breakpoint\nLag (i) = 4 days",
                         "Model B, Day 9, 1-day"~"Model B\nDay 9 breakpoint\nLag (i) = 1 day",
                         "Model C, Day 10, 3-day"~"Model C\nDay 10 breakpoint\nLag (i) = 3 days",
                         "Model E, Day 9, 4-day"~"Model E\nDay 9 breakpoint\nLag (i) = 4 days"
                         
                         ),
         facet_lab=case_match(label,
                              "Model C, Day 10, 2-day"~"Tenth-best model (3.2%)",
                              "Model E, Day 10, 2-day"~"Eleventh-best model (2.9%)",
                              "Model C, Day 9, 2-day"~"Twelfth-best model (2.5%)",
                              "Model B, No breakpoint, 4-day"~"Thirteenth-best model (2.3%)",
                              "Model F, Day 9, 3-day"~"Fourteenth-best model (2.2%)",
                              "Model F, Day 10, 2-day"~"Fifteenth-best model (2.0%)",
                              "Model F, Day 10, 3-day"~"Sixteenth-best model (1.9%)",
                              
                              "Model B, Day 11, 3-day"~"Seventeenth-best model (1.4%)",
                              "Model F, Day 8, 4-day"~"Eighteenth-best model (1.1%)",
                              "Model A, No breakpoint, 4-day"~"Nineteenth-best model (0.9%)",
                              "Model E, Day 11, 3-day"~"Twentieth-best model (0.7%)",
                              "Model F, Day 11, 4-day"~"Twenty-first-best model (0.7%)",
                              "Model B, Day 9, 1-day"~"Twenty-second-best model (0.3%)",
                              "Model C, Day 10, 3-day"~"Twenty-third-best model (0.3%)",
                              "Model E, Day 9, 4-day"~"Twenty-fourth-best model (0.3%)"),
         tag=case_match(label,
                        "Model C, Day 10, 2-day"~"A",
                        "Model E, Day 10, 2-day"~"B",
                        "Model C, Day 9, 2-day"~"C",
                        "Model B, No breakpoint, 4-day"~"D",
                        "Model F, Day 9, 3-day"~"E",
                        "Model F, Day 10, 2-day"~"F",
                        "Model F, Day 10, 3-day"~"G",
                        
                        "Model B, Day 11, 3-day"~"H",
                        "Model F, Day 8, 4-day"~"I",
                        "Model A, No breakpoint, 4-day"~"J",
                        "Model E, Day 11, 3-day"~"K",
                        "Model F, Day 11, 4-day"~"L",
                        "Model B, Day 9, 1-day"~"M",
                        "Model C, Day 10, 3-day"~"N",
                        "Model E, Day 9, 4-day"~"O"))

top_pred_df$facet_lab <- factor(top_pred_df$facet_lab,levels=c("Tenth-best model (3.2%)",
                                                               "Eleventh-best model (2.9%)",
                                                               "Twelfth-best model (2.5%)",
                                                               "Thirteenth-best model (2.3%)",
                                                               "Fourteenth-best model (2.2%)",
                                                               "Fifteenth-best model (2.0%)",
                                                               "Sixteenth-best model (1.9%)",
                                                               
                                                               "Seventeenth-best model (1.4%)",
                                                               "Eighteenth-best model (1.1%)",
                                                               "Nineteenth-best model (0.9%)",
                                                               "Twentieth-best model (0.7%)",
                                                               "Twenty-first-best model (0.7%)",
                                                               "Twenty-second-best model (0.3%)",
                                                               "Twenty-third-best model (0.3%)",
                                                               "Twenty-fourth-best model (0.3%)"))

top_pred_text$facet_lab <- factor(top_pred_text$facet_lab,levels=c("Tenth-best model (3.2%)",
                                                                   "Eleventh-best model (2.9%)",
                                                                   "Twelfth-best model (2.5%)",
                                                                   "Thirteenth-best model (2.3%)",
                                                                   "Fourteenth-best model (2.2%)",
                                                                   "Fifteenth-best model (2.0%)",
                                                                   "Sixteenth-best model (1.9%)",
                                                                   
                                                                   "Seventeenth-best model (1.4%)",
                                                                   "Eighteenth-best model (1.1%)",
                                                                   "Nineteenth-best model (0.9%)",
                                                                   "Twentieth-best model (0.7%)",
                                                                   "Twenty-first-best model (0.7%)",
                                                                   "Twenty-second-best model (0.3%)",
                                                                   "Twenty-third-best model (0.3%)",
                                                                   "Twenty-fourth-best model (0.3%)"))

top_pred_CF <- top_pred_df |>
  filter(model=="Model C"|model=="Model F")

top_pred_none <- top_pred_df |>
  filter(bp=="No breakpoint")
  



ggplot()+
  
  geom_line(data=top_pred_CF,aes(x=lagRBC,y=med,linetype=factor(phase)),col="grey",linewidth=1)+
  geom_ribbon(data=top_pred_CF,aes(x=lagRBC,ymin=lo,ymax=hi,group=interaction(phase)),alpha=0.2)+
  
  geom_line(data=top_pred_none,aes(x=lagRBC,y=med,col=pABA),linewidth=1)+
  geom_ribbon(data=top_pred_none,aes(x=lagRBC,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
  
  geom_line(data=filter(top_pred_df,model!="Model C",model!="Model F",bp!="No breakpoint"),
            aes(x=lagRBC,y=med,col=pABA,linetype=factor(phase)),linewidth=1)+
  geom_ribbon(data=filter(top_pred_df,model!="Model C",model!="Model F",bp!="No breakpoint"),
              aes(x=lagRBC,ymin=lo,ymax=hi,fill=pABA,group=interaction(phase,pABA)),alpha=0.2)+
  
  geom_label(data=top_pred_text,aes(x=7600000,y=3900000,label=text),size=3)+
  geom_text(data=top_pred_text,aes(x=300000,y=4550000,label=tag),fontface="bold",size=5)+
  
  #geom_text(data=phase_lab,aes(x=7000000,y=3000000,label=Phase2),size=5)+
  #geom_text(data=phase_lab,aes(x=4000000,y=750000,label=Phase1),size=5)+
  
  scale_x_continuous(labels = aakmisc::scinot,limits=c(0,10000000))+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4700000))+
  facet_wrap(facet_lab~.,nrow=5)+
  scale_colour_manual(values=cbpABA[2:5])+
  scale_fill_manual(values=cbpABA[2:5])+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  labs(linetype="Phase",colour="Parasite nutrient (pABA)",fill="Parasite nutrient (pABA)",x=expression(paste("RBC density at time ", italic("t"), "-i (density per µL)")),
       y=expression(paste("Reticulocyte supply at time ", italic("t"), " (density per µL)")))+
  theme_bw()+
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=11),
    axis.text=element_text(size=9),
    axis.title=element_text(size=15),
    legend.position="top",
    legend.background=element_blank()
  )

ggsave("FigureS6.jpeg",width=25,height=25,units="cm")
