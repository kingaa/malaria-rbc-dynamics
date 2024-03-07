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

#### Data preparation ####
stats_df <- read.csv("results_regression_stats.csv") |>
  filter(dataset=="sub") |>
  group_by(rep) |>
  filter(AICc==min(AICc)) |>
  ungroup() |>
  mutate(bp = X01) |>
  select(-c(X01:X04))
stats_df$bp[is.na(stats_df$bp)] <- "None"

stats_df$model <- factor(stats_df$model,levels=c("m1","m2","m3","m4","m5","m6"),
                         labels=c("Model A","Model B","Model C","Model D","Model E","Model F"))
stats_df$bp <- factor(stats_df$bp,levels=c(8,9,10,11,"None"),labels=c("Day 8","Day 9","Day 10","Day 11","No breakpoint"))
stats_df$lag <- factor(stats_df$lag,levels=c(1,2,3,4,5),labels=c("1-day","2-day","3-day","4-day","5-day"))
stats_df$pABA <- factor(stats_df$model,levels=c("Model A","Model B","Model C","Model D","Model E","Model F"),
                        labels=c("Shape/slope and intercept","Intercept only","No effect",
                                 "Shape/slope and intercept","Intercept only","No effect"))

#### Data preparation ####
stats_df |> 
  group_by(model) |> 
  count() |>
  ungroup() |>
  mutate(freq=n/sum(n))

top_df <- stats_df |> 
  group_by(bp,model,lag) |> 
  count() |>
  ungroup() |>
  mutate(freq=n/sum(n)) |>
  filter(model=="Model B"|model=="Model E")
top_df$model <- factor(top_df$model,levels=c("Model B","Model E"),
                       labels=c("Model B (54.9%)","Model E (31.1%)"))
plot <- top_df |>
  ggplot(aes(x=lag,y=freq*100,fill=bp))+
  geom_bar(stat="identity")+
  facet_grid(.~model)+
  labs(x="Reticulocyte response lag",y="Percent of models selected",fill="Breakpoint")+
  geom_text(data=filter(top_df,freq>0.01),aes(label=paste0(sprintf("%1.1f", 100*freq),"%")),
            position=position_stack(vjust=0.5),col="white",fontface="bold")+
  scale_fill_manual(values=cbPalette)+
  theme_bw()+
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=14),
    axis.text=element_text(size=10),
    axis.title=element_text(size=15),
    legend.position=c(0.65,0.8),
    legend.background=element_blank(),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14)
  )

plot <- arrangeGrob(plot)

plot <- as_ggplot(plot)+draw_plot_label(label = c("A","B"), size = 15,
                               x = c(0.04,0.525), y = c(1))

ggsave("Figure5.jpeg",width=20,height=15,units="cm")
