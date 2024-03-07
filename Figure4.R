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

stats_df$form <- factor(stats_df$model,levels=c("m1","m2","m3","m4","m5","m6"),
                        labels=c("Non-linear","Non-linear","Non-linear","Linear","Linear","Linear"))
stats_df$model <- factor(stats_df$model,levels=c("m1","m2","m3","m4","m5","m6"),
                         labels=c("Model A","Model B","Model C","Model D","Model E","Model F"))
stats_df$bp <- factor(stats_df$bp,levels=c(8,9,10,11,"None"),labels=c("Day 8","Day 9","Day 10","Day 11","No breakpoint"))
stats_df$lag <- factor(stats_df$lag,levels=c(1,2,3,4,5),labels=c("1-day","2-day","3-day","4-day","5-day"))
stats_df$pABA <- factor(stats_df$model,levels=c("Model A","Model B","Model C","Model D","Model E","Model F"),
                        labels=c("Shape/slope and intercept","Intercept only","No effect",
                                 "Shape/slope and intercept","Intercept only","No effect"))


stats_df |> count(model) |> mutate(percent = 100* n/nrow(stats_df))
stats_df |> count(lag) |> mutate(percent = 100* n/nrow(stats_df))
stats_df |> count(bp) |> mutate(percent = 100* n/nrow(stats_df))
stats_df |> count(pABA) |> mutate(percent = 100* n/nrow(stats_df))

#### Lag stacked barplot ####
stats_lag <- stats_df |> select(lag) |> group_by(lag) |> count()
stats_lag_tmp <- data.frame(c("1-day","5-day"),c(0,0)) |> setNames(c("lag","n"))
stats_lag <- bind_rows(stats_lag_tmp,stats_lag)

lagPalette <- c("#ffffcc","#c2e699","#78c679","#31a354","#006837")

plot_labs_lag <- data.frame(
  c("4-day","3-day","2-day"),
  c(22.6,48,29.4)
) |> setNames(c("lag","percent"))

lag_bar <- stats_lag |>
  ggplot()+
  geom_bar(aes(fill=factor(lag),y=n/10,x=1,col=factor(lag)),linewidth=1,stat="identity")+
  geom_text(data=plot_labs_lag,
            aes(x=1,y=percent,label=lag,alpha=c(1,0,1)),size=5,position=position_stack(vjust = 0.5))+
  geom_text(data=plot_labs_lag,
            aes(x=1,y=percent,label=lag,alpha=c(0,1,0)),fontface="bold",size=5,position=position_stack(vjust = 0.5))+
  labs(x="",y="",fill=NULL)+
  scale_fill_manual(values=lagPalette)+
  scale_colour_manual(values=rep("black",5))+
  guides(fill = guide_legend(nrow = 1,override.aes = list(fill = NA)),
         color = "none",alpha="none")+
  ggtitle("Reticulocyte response lag")+
  theme_bw()+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x=element_blank(),
    axis.title=element_text(size=15),
    axis.text=element_text(size=12),
    legend.text=element_text(size=10,colour="white"),
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(color = "white"),
    legend.position="bottom",
    panel.grid=element_blank(),
    plot.title=element_text(size=16,hjust=0.5)
  )
lag_bar

#### Model stacked barplot ####
stats_model <- stats_df |> select(model) |> group_by(model) |> count()
stats_model$form <- factor(stats_model$model,
                           levels=c("Model A","Model B","Model C",
                                    "Model D","Model E","Model F"),
                           labels=c("Non-linear","Non-linear","Non-linear",
                                    "Linear","Linear","Linear"))

plot_labs_model <- data.frame(
  c("Model F","Model E","Model C","Model B"),
  c(7.4,31.1,5.8,54.9)
) |> setNames(c("model","percent"))

model_bar <- stats_model |>
  ggplot()+
  geom_bar(aes(fill=form,y=n/10,x=1,col=model),linewidth=1,stat="identity")+
  geom_text(data=plot_labs_model,
            aes(x=1,y=percent,label=model,,alpha=c(1,1,1,0)),size=5,position=position_stack(vjust = 0.55))+
  geom_text(data=plot_labs_model,
            aes(x=1,y=percent,label=model,alpha=c(0,0,0,1)),fontface="bold",size=5,position=position_stack(vjust = 0.55))+
  labs(x="",y="",fill=NULL)+
  #scale_fill_manual(values = c("Model A"="#de2d26",
                                #"Model B"="#fb6a4a",
                                #"Model D"="#08519c",
                                #"Model C"="#fcae91",
                                #"Model E"="#3182bd",
                                #"Model F"="#6baed6"),
                    #name = "Non-linear\nLinear",
                    #breaks=c("Model A","Model D","Model B","Model E","Model C","Model F")) +
  scale_fill_manual(values = c("Non-linear" = "brown1", "Linear" = "dodgerblue2"),
                    name = "", labels = c("Non-linear","Linear"))+
  scale_colour_manual(values=rep("black",6))+
  guides(fill = guide_legend(nrow = 1),colour="none",alpha="none")+
  ggtitle("Model form")+
  theme_bw()+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x=element_blank(),
    axis.title=element_text(size=15),
    axis.text=element_text(size=12),
    legend.position="bottom",
    legend.text=element_text(size=10),
    panel.grid=element_blank(),
    plot.title=element_text(size=16,hjust=0.5),
    legend.title=element_text(size=12,lineheight=1.25,hjust=1)
  )
model_bar

#### Breakpoint stacked barplot ####
stats_bp <- stats_df |> select(bp) |> group_by(bp) |> count()
stats_bp$form <- factor(stats_bp$bp,levels=c("Day 8","Day 9","Day 10",
                                             "Day 11", "No breakpoint"),
                        labels=c(rep("Biphasic",4),"Uniphasic"))

bpPalette <- c("#fed98e","#fe9929","#d95f0e","#993404","grey")

plot_labs_bp <- data.frame(
  c("Day 11","Day 10","Day 9","Day 8"),
  c(17.9,52.5,22,9)
) |> setNames(c("bp","percent"))

bp_bar <- stats_bp |>
  ggplot()+
  geom_bar(aes(fill=form,col=factor(bp),y=n/10,x=1),linewidth=1,stat="identity")+
  geom_text(data=plot_labs_bp,
            aes(x=1,y=percent,label=bp,alpha=c(1,0,1,1)),size=5,position=position_stack(vjust = 0.55))+
  geom_text(data=plot_labs_bp,
            aes(x=1,y=percent,label=bp,alpha=c(0,1,0,0)),fontface="bold",size=5,position=position_stack(vjust = 0.55))+
  labs(x="",y="Percent of models selected",fill=NULL)+
  scale_fill_manual(values=c("#d95f0e","grey"))+
  scale_colour_manual(values=rep("black",5))+
  guides(fill = guide_legend(nrow = 1),colour="none",alpha="none")+
  ggtitle("Breakpoint")+
  theme_bw()+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x=element_blank(),
    axis.title=element_text(size=15),
    axis.text=element_text(size=12),
    legend.position="bottom",
    legend.text=element_text(size=10),
    panel.grid=element_blank(),
    plot.title=element_text(size=16,hjust=0.5)
  )
bp_bar

#### Form stacked facet ####
form_df <- stats_df |> 
  group_by(form,lag) |> 
  count() |>
  ungroup() |>
  group_by(lag) |>
  mutate(freq=n/sum(n)) |>
  ungroup()

form_df$lag <- factor(form_df$lag,levels=c("2-day","3-day","4-day"),
                      labels=c("2-day lag (29.4%)","3-day lag (48.0%)","4-day lag (22.6%)"))

form_plot <- form_df |>
  ggplot()+
  geom_bar(aes(fill=form,y=100*freq,x=1,colour=form),linewidth=1,stat="identity")+
  labs(x="",y="Percent of models selected",fill="")+
  scale_fill_manual(values = c("brown1","dodgerblue2")) +
  scale_colour_manual(values=rep("black",2))+
  guides(fill = guide_legend(ncol = 2),
         colour="none")+
  facet_grid(.~lag)+
  theme_bw()+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x=element_blank(),
    axis.title=element_text(size=15),
    axis.text=element_text(size=12),
    legend.position="bottom",
    panel.grid=element_blank(),
    plot.title=element_text(size=14,hjust=0.5),
    legend.title=element_text(size=12),
    strip.background=element_blank(),
    strip.text=element_text(size=14),
    legend.text=element_text(size=10)
  )

#### Join summary plots together ####
layout <- rbind(c(1,2,3),
                c(4,4,4))

gt <- arrangeGrob(bp_bar,lag_bar,model_bar,form_plot,
                  layout_matrix = layout,
                  heights=c(1,0.8),
                  widths=c(1,1,1))

# Add labels to the arranged plots
p <- as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B","C","D"), size = 15,
                  x = c(0.03,0.36,0.7,0.03), y = c(1,1,1,0.44))
p
ggsave("Figure4.jpeg",width=32,height=25,units="cm")

#### pABA stacked barplot ####
stats_pABA <- stats_df |> select(pABA) |> group_by(pABA) |> count()

pABAPalette <- c("#cbc9e2","#9e9ac8","#756bb1")

pABA_bar <- stats_pABA |>
  ggplot()+
  geom_bar(aes(fill=factor(pABA),y=n/10,x=1),stat="identity")+
  labs(x="",y="Percent of models selected",fill=NULL)+
  scale_fill_manual(values=pABAPalette)+
  guides(fill = guide_legend(nrow = 2))+
  ggtitle("Effect of pABA")+
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
pABA_bar


