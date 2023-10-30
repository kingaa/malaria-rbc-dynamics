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

#Load in PNAS trajectories
sm1name <- "m5sm1.rds"
sm1 <- readRDS(sm1name)

#Create sm1_mod tibble, with columns rep, mouse, mousid, box, time, E, R, lik and key
sm1 |>
  as_tibble() |>
  filter(time<=21,mouseid!="01-02",mouseid!="02-03") |> #remove underdosed mice
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

#Create dataframe with sampled keys from each box and frequency of each key
mouse_id_list <- unique(sm1_mod$mouseid)
sample_size <- 1000
foreach (
  id=mouse_id_list,
  .combine=rbind,
  .options.future = list(seed = seed_choice)
) %dofuture% {
  
  key_lik_df <- sm1_mod |> filter(mouseid==id) |> select(key,lik) |> unique()
  
  key_list <- sample(key_lik_df$key,size=sample_size,replace=TRUE,prob=key_lik_df$lik)
  
  key_table <- key_list |> table() |> as.data.frame(stringsAsFactors = FALSE) 
  names(key_table) <- c("key","Freq")
  
  key_table
  
} -> joint_key_table

#Checkpoint
stopifnot(sum(joint_key_table$Freq)==length(mouse_id_list)*sample_size)

#Create data frame with R, lagE, time, mouse and box
bake(file="jrdf.rds",{
  foreach (
    i=iter(joint_key_table,"row"),
    .combine=rbind
  ) %dofuture% {
    
    key_choice <- i$key
    freq <- i$Freq
    
    df <- sm1_mod |> 
      filter(key==key_choice) |>
      mutate(lagE=lag(E)) |>
      na.omit() |>
      select(R,lagE,time,mouse,box)
    
    if (freq>1) {
      replicate(n=freq,df,simplify=FALSE) |> bind_rows() -> df
    }
    df  
  } |>
    mutate(box=paste0("box",box)) }) -> joint_repped_df

#Create data frame with breakpoints for each of the four pABA boxes
breakpoint_grid <- as_tibble(matrix(data=c(rep(NA,4),rep(8,4),rep(9,4),rep(10,4),rep(11,4)),ncol=4,byrow=TRUE))
names(breakpoint_grid) <- c("box01","box02","box03","box04")

box_list <- unique(joint_repped_df$box)

#Create data frame that will store AICs and corresponding breakpoints for each model
bake(file="jAICdf_null.rds",{
  foreach (bp=iter(breakpoint_grid,"row"),
           .combine=rbind
  ) %dofuture% {
    
    if (is.na(bp[1])){
      
      joint_repped_df -> df
      
      #Models
      models <- list(
        m1=lm(R~poly(lagE,2,raw=T):(box-1)+(box-1),data=df),
        m2=lm(R~(poly(lagE,2,raw=T)+box-1),data=df),
        m3=lm(R~poly(lagE,2,raw=T),data=df),
        m4=lm(R~lagE:(box-1)+(box-1),data=df),
        m5=lm(R~(lagE+box-1),data=df),
        m6=lm(R~lagE,data=df)
      )
      
      models |>
        lapply(\(m) tibble(AIC=AIC(m),BIC=BIC(m),bp)) |>
        bind_rows(.id="model")
      
    } else {
     
      joint_repped_df |>
        mutate(
          phase=as.factor(if_else(time<=unlist(bp)[box],1,2))
        ) -> df
      
      #Models
      models <- list(
        m1=lm(R~poly(lagE,2,raw=T):(box-1):(phase-1)+(box-1):(phase-1),data=df),
        m2=lm(R~(poly(lagE,2,raw=T)+box-1):(phase-1)+(phase-1),data=df),
        m3=lm(R~poly(lagE,2,raw=T):(phase-1)+(phase-1),data=df),
        m4=lm(R~lagE:(box-1):(phase-1)+(box-1):(phase-1),data=df),
        m5=lm(R~(lagE+box-1):(phase-1)+(phase-1),data=df),
        m6=lm(R~lagE:(phase-1)+(phase-1),data=df)
      )
      
      models |>
        lapply(\(m) tibble(AIC=AIC(m),BIC=BIC(m),bp)) |>
        bind_rows(.id="model")
      
       
    }
    
  } }) -> joint_AIC_df

#Extract breakpoints from best model (lowest AIC)
bp <- joint_AIC_df[which(joint_AIC_df$AIC==min(joint_AIC_df$AIC)),4:7]

#Obtain dataframe with above breakpoints specified
joint_repped_df |>
  mutate(
    phase=as.factor(if_else(time<=unlist(bp)[box],1,2))
  ) -> df

#Extract best model based on lowest AIC
models <- list(
  m1=lm(R~poly(lagE,2,raw=T):(box-1):(phase-1)+(box-1):(phase-1),data=df),
  m2=lm(R~(poly(lagE,2,raw=T)+box-1):(phase-1)+(phase-1),data=df),
  m3=lm(R~poly(lagE,2,raw=T):(phase-1)+(phase-1),data=df),
  m4=lm(R~lagE:(box-1):(phase-1)+(box-1):(phase-1),data=df),
  m5=lm(R~(lagE+box-1):(phase-1)+(phase-1),data=df),
  m6=lm(R~lagE:(phase-1)+(phase-1),data=df)
)

chosen_model <- models[which(paste0("m",1:6)==joint_AIC_df[which(joint_AIC_df$AIC==min(joint_AIC_df$AIC)),1])]
chosen_model[[1]]

#Obtain confidence intervals for model parameters (2.5%, 97.5%)
library(jtools)
chosen_model_summ <- summ(chosen_model[[1]],confint=getOption("summ-confint", TRUE),ci.width=0.95)
chosen_model_coeftab <- chosen_model_summ$coeftable |> 
  as.data.frame() |> 
  rownames_to_column()

bake(file="jnd_null.rds",{
  
  joint <- data.frame()  
  
  for (box_choice in box_list){
    
    #Filter data for chosen box
    df_box <- df |> filter(box==box_choice)
    
    #Create table with coefficients for chosen box
    coeftab_box <- chosen_model_coeftab |> 
      filter(grepl(box_choice,rowname)) |>
      select(rowname,`2.5%`,`97.5%`)
    
    #Create lm object whose coefficients will be manipulated
    m01 <- lm(R~poly(lagE,2):(phase-1)+(phase-1),data=df_box)
    
    #Rewrite rowname column to match coefficient names from lm
    coeftab_box$rowname <- names(coef(m01))
    
    #Run 1000 replicates where the coefficients of m01 will be sampled from CIs
    foreach (r=1:1000,
             .combine=rbind,
             .options.future = list(seed = seed_choice)
    ) %dofuture% { #for loop over reps
      
      m02 <- m01 
      
      for (i in seq(length(coef(m02)))){ #for loop over coefficients
        
        coef_name <- names(coef(m02)[i])  
        
        coeftab_box_sub <- coeftab_box |> filter(rowname==coef_name)
        a <- coeftab_box_sub$`2.5%`
        b <- coeftab_box_sub$`97.5%`
        
        samp <- runif(1,a,b)
        
        m02$coefficients[i] <- samp
        
      } #end of for statement over i (i.e., sampling for all coefficients)
      
      #Split data based on phase
      df_box_phase1 <- df_box |> filter(phase==1)
      df_box_phase2 <- df_box |> filter(phase==2)
      
      #Generate new data to give to predict (requires lagE column and phase column)
      #Get full range across all the data and put mark what minimum of 90th percentile of that treatment
      
      newdata <- rbind(
        expand.grid(
        lagE=seq(wquant(df_box_phase1$lagE,probs=0.1), wquant(df_box_phase1$lagE,probs=0.9), 10000),
        phase=1
        ),
      expand.grid(
        lagE=seq(wquant(df_box_phase2$lagE,probs=0.1), wquant(df_box_phase2$lagE,probs=0.9), 10000),
        phase=2
      )
      )
      
      #newdata_phase1 <- data.frame(lagE=seq(min(df_box_phase1$lagE), max(df_box_phase1$lagE), 10000))
      #newdata_phase1$phase <- 1
      #newdata_phase2 <- data.frame(lagE=seq(min(df_box_phase2$lagE), max(df_box_phase2$lagE), 10000))
      #newdata_phase2$phase <- 2
      
      #newdata <- rbind(newdata_phase1,newdata_phase2)
      newdata$phase <- factor(newdata$phase)
      
      #Calculate predicted reticulocyte values for both phases using coefficients from m02
      predict(m02,newdata)
      
      pred1 <- sapply(filter(newdata,phase==1)$lagE,FUN=function(x){
        m02$coefficients[1]+
          m02$coefficients[3]*x+
          m02$coefficients[4]*x^2
      })
      
      pred2 <- sapply(filter(newdata,phase==2)$lagE,FUN=function(x){
        m02$coefficients[2]+
          m02$coefficients[5]*x+
          m02$coefficients[6]*x^2
      })
      
      #Add predicted reticulocyte values, rep number and box choice to newdata
      newdata$pred <- c(pred1,pred2)
      newdata$rep <- r
      newdata$box <- box_choice
      
      #Return newdata for use with rbind
      newdata
      
    } -> box_newdata
    
    #Join together box_newdata dfs for all boxes
    joint <- rbind(joint,box_newdata)
    
  } #end of for statement over boxes
  
  #Return join dataframe
  joint
  
}) -> joint_newdata #end of baking

#Add pABA as a factor
joint_newdata$pABA <- factor(joint_newdata$box,
                             levels=c("box04","box03","box02","box01"),
                             labels=c("0%","0.0005%","0.005%","0.05%"))

#Create summarised dataframe that calculates the avg, min and max predicted
#reticulocyte value for a given pABA, lagE and phase
joint_newdata |>
  group_by(pABA,lagE,phase) |>
  summarise(
    min=min(pred),
    avg=mean(pred),
    max=max(pred)
  ) -> joint_group


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

(ModelFitting_ER <- ggplot()+
    geom_line(data=filter(joint_group, phase==1),aes(x=lagE,y=avg,col=pABA),linewidth=2)+
    geom_line(data=filter(joint_group, phase==2),aes(x=lagE,y=avg,col=pABA),linewidth=2)+
    geom_ribbon(data=filter(joint_group, phase==1),aes(x=lagE,ymin=min,ymax=max,fill=pABA),alpha=0.2)+
    geom_ribbon(data=filter(joint_group, phase==2),aes(x=lagE,ymin=min,ymax=max,fill=pABA),alpha=0.2)+
    xlab("Erythrocytes (d-1)")+ylab("Predicted reticulocyte supply (d)")+
    scale_colour_manual(values=cbPalette[2:5])+
    scale_fill_manual(values=cbPalette[2:5])+
    theme_bw()+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=13),
      panel.grid=element_blank(),
      legend.position=c(0.8,0.8),
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )
)

#ggsave("ModelFitting_ER_Null.png",plot=ModelFitting_ER,
 #      width=20,height=15,units="cm") 

source("POMP_GroupLevel_DataPrep.R")

(RE_facet <- group_traj |>
  filter(box!="05",variable%in%c("E","R")) |>
  select(-lo,-hi) |>
  pivot_wider(names_from="variable",values_from="med") |>
  mutate(lagE=lag(E)) |>
  filter(time>0,time<=20) |>
  ggplot()+
  geom_path(aes(x=lagE,y=R,col=pABA),linewidth=2)+
  geom_text(aes(x=lagE,y=R,label=time),size=5)+
  scale_colour_manual(values=cbPalette[2:5])+
  facet_wrap(pABA~.)+
  xlab("Erythrocytes (d-1)")+ylab("Reticulocyte supply (d)")+
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
)

joint <- ggarrange(RE_facet,ModelFitting_ER, nrow = 1, labels = c("A","B"))
joint

ggsave("facet_modelFitting.png",plot=joint,
       width=30,height=15,units="cm") 
