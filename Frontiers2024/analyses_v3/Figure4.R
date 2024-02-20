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
  mutate(AIC=AIC_sub,
         AICc=AIC+2*p_sub*(p_sub+1)/(n_sub-p_sub-1),
         bp=X01) |>
  group_by(rep) |>
  filter(AICc==min(AICc)) |>
  ungroup() |> 
  select(model,lag,bp,loglik_sub,rep,AICc)
stats_df$bp[is.na(stats_df$bp)] <- "None"

stats_df$model <- factor(stats_df$model,levels=c("m1","m2","m3","m4","m5","m6"),
                         labels=c("Model A","Model B","Model C","Model D","Model E","Model F"))
stats_df$bp <- factor(stats_df$bp,levels=c(8,9,10,11,"None"),labels=c("Day 8","Day 9","Day 10","Day 11","No breakpoint"))
stats_df$lag <- factor(stats_df$lag,levels=c(1,2,3,4),labels=c("1-day","2-day","3-day","4-day"))

#### Lag stacked barplot ####
stats_lag <- stats_df |> select(lag) |> group_by(lag) |> count()
stats_lag_tmp <- data.frame("1-day",0) |> setNames(c("lag","n"))
stats_lag <- bind_rows(stats_lag_tmp,stats_lag)

lag_bar <- stats_lag |>
  ggplot()+
  geom_bar(aes(fill=factor(lag),y=n/10,x=1),stat="identity")+
  labs(x="",y="Percent of models selected",fill=NULL)+
  scale_fill_manual(values=cbPalette)+
  guides(fill = guide_legend(nrow = 2))+
  ggtitle("Reticulocyte response lag")+
  theme_bw()+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x=element_blank(),
    axis.title=element_text(size=15),
    axis.text=element_text(size=12),
    legend.position="bottom",
    panel.grid=element_blank(),
    plot.title=element_text(size=16,hjust=0.5)
  )
lag_bar

#### Lag stacked barplot ####
stats_model <- stats_df |> select(model) |> group_by(model) |> count()

model_bar <- stats_model |>
  ggplot()+
  geom_bar(aes(fill=factor(model),y=n/10,x=1),stat="identity")+
  labs(x="",y="",fill=NULL)+
  scale_fill_manual(values=cbPalette)+
  guides(fill = guide_legend(nrow = 2))+
  ggtitle("Model form")+
  theme_bw()+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x=element_blank(),
    axis.title=element_text(size=15),
    axis.text=element_text(size=12),
    legend.position="bottom",
    panel.grid=element_blank(),
    plot.title=element_text(size=16,hjust=0.5)
  )
model_bar

#### Breakpoint stacked barplot ####
stats_bp <- stats_df |> select(bp) |> group_by(bp) |> count()

bp_bar <- stats_bp |>
  ggplot()+
  geom_bar(aes(fill=factor(bp),y=n/10,x=1),stat="identity")+
  labs(x="",y="",fill=NULL)+
  scale_fill_manual(values=cbPalette)+
  guides(fill = guide_legend(nrow = 2))+
  ggtitle("Breakpoint")+
  theme_bw()+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x=element_blank(),
    axis.title=element_text(size=15),
    axis.text=element_text(size=12),
    legend.position="bottom",
    panel.grid=element_blank(),
    plot.title=element_text(size=16,hjust=0.5)
  )
bp_bar

#### Join summary plots together ####
gt <- arrangeGrob(lag_bar,model_bar,bp_bar,nrow=1)

# Add labels to the arranged plots
p <- as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B","C"), size = 15,
                  x = c(0, 0.33,0.66), y = c(1, 1,1))
p
ggsave("Figure4.jpeg",width=25,height=15,units="cm")
