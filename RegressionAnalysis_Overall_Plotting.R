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

preds_df <- read.csv("results_regression_preds.csv")

preds_df$model <- factor(preds_df$model,levels=c("m1","m2","m3","m4","m5","m6"),
                         labels=c("Model A","Model B","Model C","Model D","Model E","Model F"))
preds_df$b <- factor(preds_df$b,levels=c(8,9,10,11,"None"),labels=c("Day 8","Day 9","Day 10","Day 11","No breakpoint"))
preds_df$lag <- factor(preds_df$lag,levels=c(1,2,3,4),labels=c("1-day","2-day","3-day","4-day"))

#######################
#### Summary plots ####
#######################
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
  ggtitle("Model")+
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

#######################
#### Overall facet ####
#######################
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

####################
#### Top models ####
####################
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
                       labels=c("Model B (59.4%)","Model E (29.5%)"))
top_df |>
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

ggsave("Figure5.jpeg",width=20,height=15,units="cm")

#########################
#### Likelihood plot ####
#########################
lik_df <- stats_df
lik_df$bp <- factor(lik_df$bp,levels=levels(lik_df$bp),
                                   labels=c("Day 8","Day 9","Day 10","Day 11","None"))
lik_df$bp[is.na(lik_df$bp)] <- "None"

lik_df <- lik_df |>
  unite("label",c(model,bp,lag),sep=", ",remove=FALSE)

label_df <- lik_df |>
  group_by(label) |>
  summarize(n=n()) |>
  arrange(n) |>
  select(label)

lik_df$label <- factor(lik_df$label,levels=label_df$label)

lik_mean <- lik_df |>
  group_by(label) |>
  summarise(mean=mean(loglik_sub))

lik_df |>
  ggplot()+
  geom_boxplot(aes(x=label,y=loglik_sub))+
  geom_jitter(aes(x=label,y=loglik_sub),size=2,alpha=0.2)+
  labs(x="Overall model, in order of increasing frequency among selected models",y="Log likelihood")+
  theme_bw()+
  theme(
    axis.text=element_text(angle=90),
    axis.title=element_text(size=15)
  )

ggsave("FigureS4.jpeg",width=25,height=15,units="cm")

#########################
#### Prediction plot ####
#########################
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
  
