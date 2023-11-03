#SET DIRECTORY TO SOURCE FILE LOCATION

library(tidyverse)
library(mgcv)
library(pomp)
library(foreach)
library(iterators)
library(doFuture)
plan(multisession) #for faster results, use multicore outside of RStudio

seed_choice <- 851657743
set.seed(seed_choice)

#### Make data frame with all of the data ####

##Load in PNAS trajectories
sm1name <- "m5sm1.rds"
sm1 <- readRDS(sm1name)

#Create sm1_mod tibble, with columns rep, mouse, mousid, box, time, E, R, lik and key
sm1 |>
  as_tibble() |>
  filter(time<=20,mouseid!="01-02",mouseid!="02-03") |> #remove underdosed mice
  pivot_wider(names_from=variable,values_from=value) |>
  separate_wider_delim(cols="mouseid",delim="-",names=c("box","mouse"),cols_remove=FALSE) |>
  group_by(mouseid) |>
  mutate(
    lik=exp(loglik-max(loglik))
  ) |>
  ungroup() |>
  filter(box!="05") |> #remove control mice
  select(rep,mouse,mouseid,box,time,N,lik) |>
  unite("key",c(rep,box,mouse),sep="_",remove=FALSE) -> sm1_mod

#Loop over reps -> need to change to parallel version
rep_num <- 1000

#results_df <- data.frame()
bake(file="results_df_N.rds",{
foreach (
  r=1:rep_num,
  .combine=rbind,
  .options.future = list(seed = seed_choice)
) %dofuture% {
  
  #for (r in 1:rep_num){
  
  #print(r)
  
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
    m1=gam(N~s(time, by = box),data=joint_mouse_df),
    m2=gam(N~s(time),data=joint_mouse_df),
    m3=gam(N~1,data=joint_mouse_df)
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
  #results_df<-rbind(results_df,df)
  
} -> results_df
  
  results_df 
  
  }) -> results_df_N #end for statement over r

results_bar <- results_df_N |> 
  select(r,model) |> 
  unique()

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

results_df_N$pABA <- factor(results_df_N$box,
                          levels=c("04","03","02","01"),
                          labels=c("0%","0.0005%","0.005%","0.05%"))

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

(best_model_plot <- results_df_N |>
    filter(model=="m2") |>
    ggplot()+
    geom_line(aes(x=time,y=pred,group=r),alpha=0.01)+
    xlab("Day post-infection")+ylab("RBC clearance (density per microliter)")+
    ggtitle("Model 2")+
    theme_bw()+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=11),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13),
      strip.text=element_text(size=12),
      strip.background=element_blank(),
      plot.title=element_text(size=17,hjust=0.5)
      
    )
)

(second_model_plot <- results_df_N |>
    filter(model=="m1") |>
    ggplot()+
    geom_line(aes(x=time,y=pred,group=interaction(r,pABA),col=pABA),alpha=0.075)+
    xlab("Day post-infection")+ylab("RBC clearance (density per microliter)")+
    scale_colour_manual(values=cbPalette[2:5])+
    ggtitle("Model 2")+
    theme_bw()+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=11),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13),
      strip.text=element_text(size=12),
      strip.background=element_blank(),
      plot.title=element_text(size=17,hjust=0.5)
      
    )
)

(second_model_facet <- results_df_N |>
    filter(model=="m1") |>
    ggplot()+
    geom_line(aes(x=time,y=pred,group=interaction(r,pABA),col=pABA),alpha=0.075)+
    xlab("Day post-infection")+ylab("RBC clearance")+
    scale_colour_manual(values=cbPalette[2:5])+
    facet_grid(.~pABA)+
    theme_bw()+
    theme(
      axis.title=element_text(size=15),
      strip.text=element_blank(),
      axis.text=element_text(size=11),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13),
      strip.background=element_blank(),
      plot.title=element_text(size=17,hjust=0.5)
      
    )
)

library("gridExtra")
library("cowplot")
# Arrange plots using arrangeGrob
# returns a gtable (gt)
gt <- arrangeGrob(best_model_plot,                               
                  second_model_plot, second_model_facet,                              
                  ncol = 2, nrow = 3, 
                  layout_matrix = cbind(c(1,1,3), c(2,2,3)))
# Add labels to the arranged plots
p <- as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.38))
p
