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
    RBC = R+E,
    SM=exp(-M/(R+E)),
    SN=exp(-N/(R+E)),
    Qun=SM*(1-SN),
    lik=exp(loglik-max(loglik))
  ) |>
  ungroup() |>
  filter(box!="05") |> #remove control mice
  select(rep,mouse,mouseid,box,time,Qun,RBC,lik) |>
  unite("key",c(rep,box,mouse),sep="_",remove=FALSE) -> sm1_mod

#Loop over reps -> need to change to parallel version
rep_num <- 1000

bake(file="results_df_Qun.rds",{
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
      A=gam(Qun~s(time, by = box)+box,data=joint_mouse_df),
      B=gam(Qun~s(time),data=joint_mouse_df),
      C=gam(Qun~1,data=joint_mouse_df),
      D=gam(Qun~RBC,data=joint_mouse_df),
      E=gam(Qun~RBC:box+box,data=joint_mouse_df)
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
  
}) -> results_df_Qun

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

results_df_Qun$pABA <- factor(results_df_Qun$box,
                              levels=c("04","03","02","01"),
                              labels=c("Unsupplemented","Low","Medium","High"))

results_bar <- results_df_Qun |> 
  select(r,model) |> 
  unique()

results_bar_summ <- table(results_bar$model) |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column()
names(results_bar_summ) <- c("model","freq")

models_list <- c("A","B","C","D","E")
models_inc <- results_bar_summ$model |> c()
models_exc <- models_list[!models_list %in% models_inc]
results_bar_zero <- as.data.frame(matrix(
  c(
    models_exc,
    rep(0,length(models_exc))
  ),
  ncol=2,
  byrow=FALSE
))
names(results_bar_zero) <- c("model","freq")
results_bar_zero$freq <- as.integer(results_bar_zero$freq)

results_bar_plot <- bind_rows(results_bar_summ,results_bar_zero)

results_bar_plot$name <- factor(results_bar_plot$model,
                                levels=c("A","B","C","D","E"),
                                labels=c(
                                  expression(atop("Model A","time x pABA")),
                                  expression(atop("Model B","time")),
                                  expression(atop("Model C","constant")),
                                  expression(atop("Model D",RBC[t])),
                                  expression(atop("Model E",RBC[t]~x~pABA))
                                )
)

(bar_plot <- results_bar_plot |>
    ggplot()+
    geom_bar(aes(x=name,y=freq/1000),stat="identity")+
    annotate("text",x=3,y=0,label="X",size=5)+
    annotate("text",x=4,y=0,label="X",size=5)+
    annotate("text",x=5,y=0,label="X",size=5)+
    
    #annotate("text",x=2.75,y=0.5,label=expression(paste("Model A: ",
    #     Q[t]^un%~%s(time,~by==pABA)+pABA)),size=4,hjust=0,parse=T)+ 
    #annotate("text",x=2.75,y=0.475,label=expression(paste("Model B: ",
    #     Q[t]^un%~%s(time))),size=4,hjust=0,parse=T)+
    #annotate("text",x=2.75,y=0.45,label=expression(paste("Model C: ",
    #     Q[t]^un%~%1)),size=4,hjust=0,parse=T)+
    #annotate("text",x=2.75,y=0.425,label=expression(paste("Model D: ",
    #     Q[t]^un%~%E[t])),size=4,hjust=0,parse=T)+
    #annotate("text",x=2.75,y=0.4,label=expression(paste("Model E: ",
    #     Q[t]^un%~%E[t]:pABA+pABA)),size=4,hjust=0,parse=T)+
  
  theme_bw()+
    scale_x_discrete(labels = label_parse())+
    ylab("Frequency model selected")+
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
)

(medians_plot <- group_traj |>
    filter(variable=="Qun",time<=20) |>
    ggplot()+
    geom_line(aes(x=time,y=med,col=pABA),linewidth=2)+
    geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
    scale_colour_manual(values=cbPalette[2:5])+
    scale_fill_manual(values=cbPalette[2:5])+
    theme_bw()+
    facet_grid(.~pABA)+
    ylim(0,1)+
    labs(colour="Parasite nutrient\n(pABA)",fill="Parasite nutrient\n(pABA)")+
    xlab("Time (d post-infection)")+ylab("RBC clearance (probability/day)")+
    theme(
      axis.title=element_text(size=15),
      strip.text=element_blank(),
      axis.text=element_text(size=11),
      #panel.grid=element_blank(),
      legend.position=c(0.1,0.85),
      legend.background=element_blank(),
      legend.title=element_text(size=12),
      legend.text=element_text(size=10),
      strip.background=element_blank(),
      plot.title=element_text(size=17,hjust=0.5)
      
    )
)

library("gridExtra")

ggpubr::ggarrange(medians_plot, bar_plot, nrow = 1, labels = c("A","B"), widths=c(0.66,0.33))

ggsave("Figure6.jpeg",width=35,height=15,units="cm")
