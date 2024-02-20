#SET DIRECTORY TO SOURCE FILE LOCATION

library(tidyverse)
library(mgcv)
library(pomp)
library(foreach)
library(iterators)
library(doFuture)
library(aakmisc)
plan(multisession) #for faster results, use multicore outside of RStudio

seed_choice <- 851657743
set.seed(seed_choice)

##Load in PNAS trajectories
sm1name <- "m5sm1_mod.rds"
sm1 <- readRDS(sm1name)

#Create sm1_mod tibble, with columns rep, mouse, mousid, box, time, E, R, lik and key
sm1 |>
  as_tibble() |>
  filter(time<=20,mouseid!="01-02",mouseid!="02-03") |> #remove underdosed mice
  pivot_wider(names_from=variable,values_from=value) |>
  group_by(mouseid) |>
  mutate(
    lik=exp(loglik-max(loglik))
  ) |>
  ungroup() |>
  filter(box!="05") |> #remove control mice
  select(rep,mouse,mouseid,box,time,R,lik) |>
  unite("key",c(rep,box,mouse),sep="_",remove=FALSE) -> sm1_mod

#Loop over reps
rep_num <- 1000

bake(file="results_df_R.rds",{
  foreach (
    r=1:rep_num,
    .combine=rbind,
    .options.future = list(seed = seed_choice)
  ) %dofuture% {
    
    mouse_id_list <- unique(sm1_mod$mouseid)
    
    joint_mouse_df <- data.frame()
    for (id in mouse_id_list){
      
      key_opts <- sm1_mod |> filter(mouseid==id) |> select(key,lik) |> unique()
      
      key_choice <- sample(key_opts$key,size=1,prob=key_opts$lik)
      
      sm1_mod_mouse <- sm1_mod |> 
        filter(key==key_choice)
      
      joint_mouse_df <- rbind(joint_mouse_df,sm1_mod_mouse)
      
    }
    
    joint_mouse_df$box <- factor(joint_mouse_df$box)
    
    #Models
    models <- list(
      m1=gam(R~s(time, by = box)+box,data=joint_mouse_df),
      m2=gam(R~s(time),data=joint_mouse_df),
      m3=gam(R~1,data=joint_mouse_df)
    )
    
    models |>
      lapply(\(m) tibble(AIC=AIC(m),BIC=BIC(m))) |>
      bind_rows(.id="model") -> AIC_df
    
    ##Extract breakpoints from best model (lowest AIC)
    AIC_df |>
      filter(AIC==min(AIC)) -> best
    
    
    chosen_model <- models[[best$model]]
    coefs <- chosen_model$coefficients
    
    newdata <- joint_mouse_df
    
    newdata$pred <- predict(chosen_model,newdata)
    
    newdata$r <- r
    newdata$model <- best$model
    
    newdata
    
  } -> results_df
  
  results_df 
  
}) -> results_df_R

results_bar <- results_df_R |> 
  select(r,model) |> 
  unique()

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

results_df_R$pABA <- factor(results_df_R$box,
                              levels=c("04","03","02","01"),
                              labels=c("Unsupplemented","Low","Medium","High"))

table(results_bar$model)/1000

results_bar |>
  ggplot(aes(x=reorder(model, model, function(x)-length(x)),y = (after_stat(count))/sum(after_stat(count))))+
  geom_bar()+
  ylab("Frequency")+
  scale_x_discrete(labels=c("Model 2","Model 1"))+
  theme_bw()+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=11),
    axis.title.x=element_blank(),
    panel.grid=element_blank(),
    legend.position="none",
    legend.title=element_text(size=15),
    legend.text=element_text(size=13),
    strip.text=element_text(size=12),
    strip.background=element_blank()
    
  )

range(results_df_R$pred)

second_model_plot <- results_df_R |>
    filter(model=="m2") |>
    ggplot()+
    geom_line(aes(x=time,y=pred,group=r),alpha=0.075)+
    xlab("")+ylab("")+
    geom_text(aes(x=17,y=4.7e6,label="RBC supply ~ time"),size=3)+
    scale_y_continuous(label=aakmisc::scinot,lim=range(results_df_R$pred))+
    theme_bw()+
    theme(
      axis.title=element_text(size=13),
      axis.text=element_text(size=11),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13),
      strip.text=element_text(size=12),
      strip.background=element_blank(),
      plot.title=element_text(size=14,hjust=0.5)
      
    )

best_model_plot <- results_df_R |>
    filter(model=="m1") |>
    ggplot()+
    geom_line(aes(x=time,y=pred,group=interaction(r,pABA),col=pABA),alpha=0.03)+
    xlab("")+ylab("RBC supply (density per μL)")+
    labs(colour="Parasite nutrient (pABA)")+
    scale_colour_manual(values=cbPalette[2:5])+
    geom_text(aes(x=15,y=4.7e6,label="RBC supply ~ time x pABA"),size=3)+
    guides(colour = guide_legend(override.aes = list(alpha = 1)))+
    scale_y_continuous(label=aakmisc::scinot,lim=range(results_df_R$pred))+
    theme_bw()+
    theme(
      axis.title=element_text(size=13),
      axis.text=element_text(size=11),
      panel.grid=element_blank(),
      legend.position=c(0.275,0.7),
      legend.background = element_blank(),
      legend.title=element_text(size=10),
      legend.text=element_text(size=8),
      strip.text=element_text(size=12),
      strip.background=element_blank(),
      plot.title=element_text(size=14,hjust=0.5,face="bold"),
      panel.border = element_rect(linewidth=4)
      
    )

best_model_facet <- results_df_R |>
    filter(model=="m1") |>
    ggplot()+
    geom_line(aes(x=time,y=pred,group=interaction(r,pABA),col=pABA),alpha=0.03)+
    xlab("Time (d post-infection)")+ylab("")+
    scale_colour_manual(values=cbPalette[2:5])+
    scale_y_continuous(label=aakmisc::scinot,lim=range(results_df_R$pred))+
    facet_grid(.~pABA)+
    theme_bw()+
    theme(
      axis.title=element_text(size=13),
      strip.text=element_blank(),
      axis.text=element_text(size=11),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=12),
      legend.text=element_text(size=11),
      strip.background=element_blank(),
      plot.title=element_text(size=17,hjust=0.5)
      
    )

library("gridExtra")
library("cowplot")
# Arrange plots using arrangeGrob
# returns a gtable (gt)
gt <- arrangeGrob(best_model_plot,                               
                  second_model_plot, best_model_facet,                              
                  ncol = 2, nrow = 3, 
                  layout_matrix = cbind(c(1,1,3), c(2,2,3)))
# Add labels to the arranged plots
p <- as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.4))

ggsave("FigureS6.jpeg",width=22,height=16,units="cm")
