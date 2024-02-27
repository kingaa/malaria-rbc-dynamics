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
stats_df <- read.csv("results_regression_stats.csv")|>
  group_by(rep) |>
  filter(AICc==min(AICc)) |>
  ungroup() |> 
  select(model,lag,bp=X01,loglik,rep,AICc)
stats_df$bp[is.na(stats_df$bp)] <- "None"

stats_df$model <- factor(stats_df$model,levels=c("m1","m2","m3","m4","m5","m6"),
                         labels=c("Model A","Model B","Model C","Model D","Model E","Model F"))
stats_df$bp <- factor(stats_df$bp,levels=c(8,9,10,11,"None"),labels=c("Day 8","Day 9","Day 10","Day 11","None"))
stats_df$lag <- factor(stats_df$lag,levels=c(1,2,3,4),labels=c("1-day","2-day","3-day","4-day"))

overall_freq <- stats_df |> 
  group_by(bp,model,lag) |> 
  count() |>
  ungroup() |>
  mutate(freq=n/sum(n))

model_list <- stats_df$model |> unique()
bp_list <- stats_df$bp |> unique()
combos <- expand_grid(model_list,bp_list) |> setNames(c("model","bp"))
df_x <- anti_join(combos, overall_freq, by = c("model","bp"))
df_x$label="X"

overall_freq |>
  ggplot(aes(x=lag,y=freq))+
  geom_bar(stat="identity",fill="darkblue")+
  geom_text(aes(label = freq),vjust=-0.5)+
  geom_text(data=df_x,aes(x=2,y=0.15,label=label),size=25,col="lightgrey")+
  facet_grid(bp~model)+
  ylim(0,0.3)+
  labs(x="Reticulocyte response lag",y="Frequency")+
  scale_x_discrete(breaks=c("1-day","2-day","3-day","4-day"))+
  theme_bw()+
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=12),
    axis.text=element_text(size=10),
    axis.title=element_text(size=15)
  )

ggsave("FigureS3.jpeg",width=25,height=20,units="cm")
