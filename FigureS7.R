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

#### Data preparation ####
stats_df <- read.csv("results_regression_stats.csv") |>
  filter(dataset == "sub") |>
  group_by(rep) |>
  filter(AICc==min(AICc)) |>
  ungroup () |>
  select(dataset,model, bp=X01,lag=lagChoice,loglik)
stats_df$bp[is.na(stats_df$bp)] <- "None"

stats_df$model <- factor(stats_df$model,levels=c("m1","m2","m3","m4","m5","m6","m7","m8"),
                         labels=c("Model A","Model B","Model C","Model D","Model E","Model F","Model G","Model H"))
stats_df$bp <- factor(stats_df$bp,levels=c(8,9,10,11,"None"),labels=c("Day 8","Day 9","Day 10","Day 11","None"))
stats_df$lag <- factor(stats_df$lag,levels=c(1,2,3,4,5),labels=c("1-day","2-day","3-day","4-day","5-day"))

lik_df <- stats_df |>
  unite("label",c(model,bp,lag),sep=", ",remove=FALSE)

label_df <- lik_df |>
  group_by(label) |>
  summarize(n=n()) |>
  arrange(n) |>
  select(label)

lik_df$label <- factor(lik_df$label,levels=label_df$label)

lik_mean <- lik_df |>
  group_by(label) |>
  summarise(mean=mean(loglik))

lik_df |>
  ggplot()+
  geom_boxplot(aes(x=label,y=loglik))+
  geom_jitter(aes(x=label,y=loglik),size=2,alpha=0.2)+
  labs(x="Model structure (form, breakpoint, lag)",y="Log likelihood")+
  theme_bw()+
  theme(
    axis.text=element_text(angle=90),
    axis.title=element_text(size=15)
  )

ggsave("FigureS7.jpeg",width=25,height=15,units="cm")
