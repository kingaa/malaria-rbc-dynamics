#SET DIRECTORY TO SOURCE FILE LOCATION

library(tidyverse)
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

##Create sm1_mod tibble, with columns rep, mouse, mousid, box, time, E, R, lik and key
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
  select(rep,mouse,mouseid,box,time,E,R,lik) |>
  unite("key",c(rep,box,mouse),sep="_",remove=FALSE) -> sm1_mod

##Create data frame with breakpoints for each of the four pABA boxes
box_list <- unique(sm1_mod$box)
breakpoint_grid <- as_tibble(
  matrix(
    data=rep(c(NA,8,9,10,11),each=4),
    ncol=4,byrow=TRUE,
    dimnames=list(NULL,box_list)
  )
)

#Loop over reps -> need to change to parallel version
rep_num <- 1000

#results_df <- data.frame()

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
      filter(key==key_choice) |>
      mutate(lagE=lag(E)) |>
      na.omit()
    
    joint_mouse_df <- rbind(joint_mouse_df,sm1_mod_mouse)
    
  }
  
  joint_AIC_df <- data.frame()
  for (i in seq(1,nrow(breakpoint_grid),1)){
    
    bp <- breakpoint_grid[i,]
    
    if (is.na(bp[1])){
      
      joint_mouse_df -> df
      
      ##Models
      models <- list(
        m1=lm(R~poly(lagE,2,raw=F):(box-1)+(box-1),data=df),
        m2=lm(R~(poly(lagE,2,raw=F)+box-1),data=df),
        m3=lm(R~poly(lagE,2,raw=F),data=df),
        m4=lm(R~lagE:(box-1)+(box-1),data=df),
        m5=lm(R~(lagE+box-1),data=df),
        m6=lm(R~lagE,data=df)
      )
      
      models |>
        lapply(\(m) tibble(AIC=AIC(m),BIC=BIC(m),bp)) |>
        bind_rows(.id="model") -> model_df
      
    } else {
      
      joint_mouse_df |>
        mutate(
          phase=as.factor(if_else(time<=unlist(bp)[box],1,2))
        ) -> df
      
      ##Models
      models <- list(
        m1=lm(R~poly(lagE,2,raw=F):(box-1):(phase-1)+(box-1):(phase-1),data=df),
        m2=lm(R~(poly(lagE,2,raw=F)+box-1):(phase-1)+(phase-1),data=df),
        m3=lm(R~poly(lagE,2,raw=F):(phase-1)+(phase-1),data=df),
        m4=lm(R~lagE:(box-1):(phase-1)+(box-1):(phase-1),data=df),
        m5=lm(R~(lagE+box-1):(phase-1)+(phase-1),data=df),
        m6=lm(R~lagE:(phase-1)+(phase-1),data=df)
      )
      
      models |>
        lapply(\(m) tibble(AIC=AIC(m),BIC=BIC(m),bp)) |>
        bind_rows(.id="model") -> model_df
      
    }
    
    joint_AIC_df <- rbind(joint_AIC_df,model_df)
  }
  
  ##Extract breakpoints from best model (lowest AIC)
  joint_AIC_df |>
    filter(AIC==min(AIC)) -> best
  
  ##Obtain data frame with above breakpoints specified
  joint_mouse_df |>
    mutate(
      phase=as.factor(if_else(time<=as.integer(unlist(best)[box]),1,2))
    ) -> df
  
  ##Extract best model based on lowest AIC
  models <- list(
    m1=lm(R~poly(lagE,2,raw=F):(box-1):(phase-1)+(box-1):(phase-1),data=df),
    m2=lm(R~(poly(lagE,2,raw=F)+box-1):(phase-1)+(phase-1),data=df),
    m3=lm(R~poly(lagE,2,raw=F):(phase-1)+(phase-1),data=df),
    m4=lm(R~lagE:(box-1):(phase-1)+(box-1):(phase-1),data=df),
    m5=lm(R~(lagE+box-1):(phase-1)+(phase-1),data=df),
    m6=lm(R~lagE:(phase-1)+(phase-1),data=df)
  )
  
  chosen_model <- models[[best$model]]
  coefs <- chosen_model$coefficients
  
  pred <- model.matrix(chosen_model) %*% coefs |> as.data.frame() |> select(pred=V1)
  df <- cbind(df,pred)
  df$r <- r
  df$b <- best$`01`
  df$model <- best$model
  
  df
  #results_df<-rbind(results_df,df)
  
} -> results_df #end for statement over r

results_bar <- results_df |> 
  select(r,b,model) |> 
  unique() |>
  unite("model_b",c(model,b),sep=", ",remove=FALSE)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

results_bar |>
  ggplot(aes(x=reorder(b, b, function(x)-length(x)),y = (after_stat(count))/sum(after_stat(count))))+
  geom_bar()+
  xlab("Breakpoint (day post-infection)")+ylab("Frequency")+
  theme_bw()+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=11),
    panel.grid=element_blank(),
    legend.position="none",
    legend.title=element_text(size=15),
    legend.text=element_text(size=13),
    strip.text=element_text(size=12),
    strip.background=element_blank()
    
  )

results_bar$label <- factor(results_bar$model_b,
                            levels=c("m6, 9","m3, 9","m5, 9","m2, 9","m1, 10","m1, 9","m4, 9","m4, 10","m5, 10","m6, 10"),
                            labels=c("Model 6,\nbreakpoint 9",
                                     "Model 3,\nbreakpoint 9",
                                     "Model 5,\nbreakpoint 9",
                                     "Model 2,\nbreakpoint 9",
                                     "Model 1,\nbreakpoint 10",
                                     "Model 1,\nbreakpoint 9",
                                     "Model 4,\nbreakpoint 9",
                                     "Model 4,\nbreakpoint 10",
                                     "Model 5,\nbreakpoint 10",
                                     "Model 6,\nbreakpoint 10"))

xaxis_labs <- c("Model 6,\nbreakpoint 9",
                "Model 3,\nbreakpoint 9",
                "Model 5,\nbreakpoint 9",
                "Model 2,\nbreakpoint 9",
                "Model 1,\nbreakpoint 10",
                "Model 1,\nbreakpoint 9",
                "Model 4,\nbreakpoint 9",
                "Model 4,\nbreakpoint 10",
                "Model 5,\nbreakpoint 10",
                "Model 6,\nbreakpoint 10")

results_bar |>
  ggplot(aes(x=reorder(model_b, model_b, function(x)-length(x)),y = (after_stat(count))/sum(after_stat(count))))+
  geom_bar(fill=cbPalette[7])+
  xlab("Model, breakpoint")+ylab("Frequency")+
  scale_x_discrete(labels=xaxis_labs)+
  theme_bw()+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=11),
    axis.text.x=element_text(angle=45,vjust=0.7),
    axis.title.x=element_blank(),
    panel.grid=element_blank(),
    legend.position="none",
    legend.title=element_text(size=15),
    legend.text=element_text(size=13),
    strip.text=element_text(size=12),
    strip.background=element_blank()
    
  )



results_df$pABA <- factor(results_df$box,
                          levels=c("04","03","02","01"),
                          labels=c("0%","0.0005%","0.005%","0.05%"))


(best_model_plot <- results_df |>
  unite("model_b",c(model,b),sep=", ",remove=FALSE) |>
  filter(model_b=="m6, 9") |>
  ggplot()+
  geom_line(aes(x=lagE,y=pred,group=interaction(r,phase),col=phase),alpha=0.01)+
  annotate("text",x=7500000,y=3500000,label="Phase 1",col=cbPalette[6],size=5)+
    annotate("text",x=7500000,y=3250000,label="Phase 2",col=cbPalette[8],size=5)+
  xlab("Erythrocyte density (d-1)")+ylab("Reticulocyte supply (d)")+
  scale_colour_manual(values=cbPalette[c(6,8)])+
    labs(tag = "A") +
  ggtitle("Model 6, breakpoint 9 (day post-infection)")+
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

(second_model_plot <- results_df |>
    unite("model_b",c(model,b),sep=", ",remove=FALSE) |>
    filter(model_b=="m3, 9") |>
    ggplot()+
    geom_line(aes(x=lagE,y=pred,group=interaction(r,phase),col=phase),alpha=0.05)+
    xlab("Erythrocyte density (d-1)")+ylab("Reticulocyte supply (d)")+
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
    scale_colour_manual(values=cbPalette[c(6,8)])+
    ggtitle("Model 3, breakpoint 9")+
    theme_bw()+
    theme(
      axis.title=element_text(size=12),
      axis.text=element_text(size=9),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=10),
      legend.text=element_text(size=18),
      strip.text=element_text(size=12),
      strip.background=element_blank(),
      plot.title=element_text(size=14,hjust=0.5)
      
    )
)

(third_model_plot <- results_df |>
    unite("model_b",c(model,b),sep=", ",remove=FALSE) |>
    filter(model_b=="m5, 9") |>
    ggplot()+
    geom_line(aes(x=lagE,y=pred,group=interaction(r,phase,box),col=pABA),alpha=0.05)+
    xlab("Erythrocyte density (d-1)")+ylab("Reticulocyte supply (d)")+
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
    scale_colour_manual(values=cbPalette[c(2:5)])+
    guides(colour = guide_legend(override.aes = list(alpha = 1)))+
    ggtitle("Model 5, breakpoint 9")+
    theme_bw()+
    theme(
      axis.title=element_text(size=12),
      axis.text=element_text(size=9),
      panel.grid=element_blank(),
      legend.position=c(0.8,0.7),
      legend.title=element_text(size=10),
      legend.text=element_text(size=8),
      strip.text=element_text(size=12),
      strip.background=element_blank(),
      plot.title=element_text(size=14,hjust=0.5)
      
    )
)

library("gridExtra")
library("cowplot")
# Arrange plots using arrangeGrob
# returns a gtable (gt)
gt <- arrangeGrob(best_model_plot,                               
                  second_model_plot, third_model_plot,                              
                  ncol = 3, nrow = 2, 
                  layout_matrix = cbind(c(1,1), c(1,1), c(2,3)))
# Add labels to the arranged plots
p <- as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0.66, 0.66), y = c(1, 1, 0.5))
p


