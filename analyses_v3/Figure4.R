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

stats_df |> count(model) |> mutate(percent = 100* n/nrow(best))
stats_df |> count(lag) |> mutate(percent = 100* n/nrow(best))
stats_df |> count(bp) |> mutate(percent = 100* n/nrow(best))

#### Lag stacked barplot ####
stats_lag <- stats_df |> select(lag) |> group_by(lag) |> count()
stats_lag_tmp <- data.frame(c("1-day","5-day"),c(0,0)) |> setNames(c("lag","n"))
stats_lag <- bind_rows(stats_lag_tmp,stats_lag)

lag_bar <- stats_lag |>
  ggplot()+
  geom_bar(aes(fill=factor(lag),y=n/10,x=1),stat="identity")+
  labs(x="",y="",fill=NULL)+
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

#### Model stacked barplot ####
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
  labs(x="",y="Percent of models selected",fill=NULL)+
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

#### Top model bar plots ####
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
    legend.position=c(0.6,0.8),
    legend.background=element_blank(),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14)
  )

#### Join summary plots together ####
gt <- arrangeGrob(bp_bar,model_bar,lag_bar,nrow=1)

gt2 <- arrangeGrob(gt,plot,nrow=2)

# Add labels to the arranged plots
p <- as_ggplot(gt2) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B","C","D","E"), size = 15,
                  x = c(0, 0.33,0.66,0.04,0.525), y = c(1, 1,1,0.5,0.5))
p
ggsave("Figure4.jpeg",width=25,height=25,units="cm")
