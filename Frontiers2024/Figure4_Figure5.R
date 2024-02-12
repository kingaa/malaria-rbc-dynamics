#SET DIRECTORY TO SOURCE FILE LOCATION

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

cbModels <- c("#332288","#117733","#88CCEE","#DDCC77","#CC6677","#882255")
cbBreakpoint <- c("#44AA99","#AA4499","#D55E00","#999999")
cbpABA <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

data_lag1 <- readRDS("results_df_RBC_lag1.rds") |> mutate(lag=1)
data_lag1$pABA <- factor(data_lag1$box,
                              levels=c("04","03","02","01"),
                              labels=c("Unsupplemented","Low","Medium","High"))

data_lag2 <- readRDS("results_df_RBC_lag2.rds") |> mutate(lag=2)
data_lag2$pABA <- factor(data_lag2$box,
                         levels=c("04","03","02","01"),
                         labels=c("Unsupplemented","Low","Medium","High"))

data_lag3 <- readRDS("results_df_RBC_lag3.rds") |> mutate(lag=3)
data_lag3$pABA <- factor(data_lag3$box,
                         levels=c("04","03","02","01"),
                         labels=c("Unsupplemented","Low","Medium","High"))

full_data <- data_lag1 |> bind_rows(data_lag2) |> bind_rows(data_lag3)
full_data$lag_label <- factor(full_data$lag,levels=c(1,2,3),
                              labels=c("i = 1","i = 2","i = 3"))
full_data$model <- factor(full_data$model,
                          levels=c("m1","m2","m3","m4","m5","m6"),
                          labels=c("Model A","Model B","Model C","Model D","Model E","Model F"))
full_data$b <- factor(full_data$b,
                      levels=c(9,10,11,"None"),
                      labels=c("Day 9","Day 10","Day 11","None"))

data <- data_lag1 |> 
  bind_rows(data_lag2) |> 
  bind_rows(data_lag3) |>
  select(r,b,model,lag) |>
  unique()

data_model <- data |>
  group_by(lag,model) |>
  count()
data_model$model <- factor(data_model$model,
                           levels=c("m1","m2","m3","m4","m5","m6"),
                           labels=c("Model A","Model B","Model C","Model D","Model E","Model F"))
data_model$lag <- factor(data_model$lag,levels=c(1,2,3),labels=c("1-day","2-day","3-day"))

data_model |> group_by(model) |> summarize(percent=sum(n)/3000)

data_bp <- data |>
  group_by(lag,b) |>
  count()
data_bp$b <- factor(data_bp$b,
                           levels=c(9,10,11,"None"),
                           labels=c("Day 9","Day 10","Day 11","No breakpoint"))
data_bp$lag <- factor(data_bp$lag,levels=c(1,2,3),labels=c("1-day","2-day","3-day"))

#Stacked bar plot for model percentage by lag
gt1 <- data_model |>
  ggplot()+
  geom_bar(aes(fill=model,y=n/10,x=lag),stat="identity")+
  scale_fill_manual(values=cbModels)+
  labs(x="Reticulocyte response lag",y="Percent of models selected",fill=NULL)+
  theme_bw()+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=12),
    legend.position="top"
  )

#Stacked bar plot for breakpoint percentage by lag
gt2 <- data_bp |>
  ggplot()+
  geom_bar(aes(fill=b,y=n/10,x=lag),stat="identity")+
  scale_fill_manual(values=cbBreakpoint)+
  labs(x="Reticulocyte response lag",y="Percent of breakpoints selected",fill=NULL)+
  theme_bw()+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=12),
    legend.position="top"
  )

gt <- arrangeGrob(gt2,gt1,nrow=1)

# Add labels to the arranged plots
p <- as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0.5), y = c(1, 1))

ggsave("Figure4.jpeg",width=20,height=15,units="cm")

#Prediction facet plot
data_plot <- full_data |>
  unite("model_b",c(model,b),sep="_",remove=FALSE) |>
  filter((lag==1&model_b=="Model B_Day 9")|(lag==1&model_b=="Model C_Day 9")|(lag==2&model_b=="Model B_Day 9")|(lag==2&model_b=="Model B_Day 10")|(lag==3&model_b=="Model B_Day 10")|(lag==3&model_b=="Model E_Day 10")) |>
  select(lagRBC,pred,phase,model_b,b,r,model,lag,lag_label,pABA) |>
  unite("model_key",c(model_b,lag),sep="_",remove=FALSE) 
data_plot$rank <- factor(data_plot$model_key,
                         levels=c("Model B_Day 9_1",
                                  "Model C_Day 9_1",
                                  "Model B_Day 9_2",
                                  "Model B_Day 10_2",
                                  "Model B_Day 10_3",
                                  "Model E_Day 10_3"),
                         labels=c("First","Second","First","Second","First","Second"))
data_plot$model_label <- factor(data_plot$model_key,
                          levels=c("Model B_Day 9_1",
                                   "Model C_Day 9_1",
                                   "Model B_Day 9_2",
                                   "Model B_Day 10_2",
                                   "Model B_Day 10_3",
                                   "Model E_Day 10_3"),
                          labels=c("Model B, breakpoint 9 (82.9%)",
                                   "Model C, breakpoint 9 (10.7%)",
                                   "Model B, breakpoint 9 (42.7%)",
                                   "Model B, breakpoint 10 (42.6%)",
                                   "Model B, breakpoint 10 (50.6%)",
                                   "Model E, breakpoint 10 (15.2%)"))

facet_labels <- data_plot |>
  select(lag_label,rank,model_label) |> 
  unique()

phase_labels <- data.frame(
  phase1 = c("Phase 1","Phase 1","","","",""),
  phase2 = c("Phase 2","Phase 2","","","",""),
  lag_label   = c("i = 1","i = 1","i = 2","i = 2","i = 3","i = 3"),
  rank = c("First","Second","First","Second","First","Second")
)
                          

prediction_facet <- ggplot()+
  geom_line(data=filter(data_plot,model!="Model C"),aes(x=lagRBC,y=pred,group=interaction(r,phase,pABA),col=pABA),alpha=0.05)+
  geom_line(data=filter(data_plot,model=="Model C"),aes(x=lagRBC,y=pred,group=interaction(r,phase,pABA)),alpha=0.05)+
  geom_text(data=filter(phase_labels,rank=="First"),aes(x=7500000,y=3000000,label=phase2))+
  geom_text(data=filter(phase_labels,rank!="First"),aes(x=6500000,y=3000000,label=phase2))+
  geom_text(data=phase_labels,aes(x=5000000,y=1000000,label=phase1))+
  scale_x_continuous(labels = aakmisc::scinot,limits=c(0,10000000))+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4600000))+
  scale_colour_manual(values=cbpABA[2:5])+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  geom_text(data=facet_labels,aes(x=6250000,y=4500000,label=model_label))+
  labs(colour="Parasite nutrient (pABA)",x=expression(paste("RBC density at time ", italic("t"), "-i (density per µL)")),
       y=expression(paste("Reticulocyte supply at time ", italic("t"), " (density per µL)")))+
  facet_grid(lag_label~rank)+
  theme_bw()+
  theme(
    strip.text.x=element_blank(),
    strip.text.y=element_text(size=12),
    axis.title=element_text(size=15),
    strip.background=element_blank(),
    panel.grid.minor = element_blank(),
    legend.position="bottom"
  )
prediction_facet
prediction_facet_tagged <- tag_facet(prediction_facet,tag_pool=c("A","B","C","D","E","F"),open="",close="")

ggsave("Figure5.jpeg",width=18,height=20,units="cm")


  
