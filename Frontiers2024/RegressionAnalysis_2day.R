#SET DIRECTORY TO SOURCE FILE LOCATION

#Load required packages
library(tidyverse)
library(pomp)
library(foreach)
library(iterators)
library(doFuture)
library(aakmisc)

#Set up for parallelisation 
plan(multisession) #for faster results, use multicore outside of RStudio

#Set seed (important because random sampling occurs below)
seed_choice <- 851657743
set.seed(seed_choice)

##Load in smooth distribution samples from PNAS work
sm1name <- "m5sm1_mod.rds"
sm1 <- readRDS(sm1name)

##Create sm1_mod tibble, with columns rep, mouse, mousid, box, time, E, R, lik and key
##From sm1_mod, we will be sampling individual distribution samples
sm1 |>
  as_tibble() |>
  filter(time<=20,mouseid!="01-02",mouseid!="02-03") |> #remove underdosed mice
  pivot_wider(names_from=variable,values_from=value) |>
  group_by(mouseid) |>
  mutate(
    RBC = R+E,
    lik=exp(loglik-max(loglik))
  ) |>
  ungroup() |>
  filter(box!="05") |> #remove control mice
  select(rep,mouse,mouseid,box,time,E,R,RBC,lik) |>
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

#Loop over reps
rep_num <- 1000

bake(file="results_df_RBC_lag2.rds",{
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
        filter(key==key_choice) |>
        mutate(lagRBC=lag(RBC,2)) |>
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
          m1=lm(R~poly(lagRBC,2,raw=F):(box-1)+(box-1),data=df),
          m2=lm(R~(poly(lagRBC,2,raw=F)+box-1),data=df),
          m3=lm(R~poly(lagRBC,2,raw=F),data=df),
          m4=lm(R~lagRBC:(box-1)+(box-1),data=df),
          m5=lm(R~(lagRBC+box-1),data=df),
          m6=lm(R~lagRBC,data=df)
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
          m1=lm(R~poly(lagRBC,2,raw=F):(box-1):(phase-1)+(box-1):(phase-1),data=df),
          m2=lm(R~(poly(lagRBC,2,raw=F)+box-1):(phase-1)+(phase-1),data=df),
          m3=lm(R~poly(lagRBC,2,raw=F):(phase-1)+(phase-1),data=df),
          m4=lm(R~lagRBC:(box-1):(phase-1)+(box-1):(phase-1),data=df),
          m5=lm(R~(lagRBC+box-1):(phase-1)+(phase-1),data=df),
          m6=lm(R~lagRBC:(phase-1)+(phase-1),data=df)
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
    if (is.na(best$`01`)){
      models <- list(
        m1=lm(R~poly(lagRBC,2,raw=F):(box-1)+(box-1),data=df),
        m2=lm(R~(poly(lagRBC,2,raw=F)+box-1),data=df),
        m3=lm(R~poly(lagRBC,2,raw=F),data=df),
        m4=lm(R~lagRBC:(box-1)+(box-1),data=df),
        m5=lm(R~(lagRBC+box-1),data=df),
        m6=lm(R~lagRBC,data=df)
      )
    } else {
      models <- list(
        m1=lm(R~poly(lagRBC,2,raw=F):(box-1):(phase-1)+(box-1):(phase-1),data=df),
        m2=lm(R~(poly(lagRBC,2,raw=F)+box-1):(phase-1)+(phase-1),data=df),
        m3=lm(R~poly(lagRBC,2,raw=F):(phase-1)+(phase-1),data=df),
        m4=lm(R~lagRBC:(box-1):(phase-1)+(box-1):(phase-1),data=df),
        m5=lm(R~(lagRBC+box-1):(phase-1)+(phase-1),data=df),
        m6=lm(R~lagRBC:(phase-1)+(phase-1),data=df)
      )
    }
    
    chosen_model <- models[[best$model]]
    coefs <- chosen_model$coefficients
    
    pred <- model.matrix(chosen_model) %*% coefs |> as.data.frame() |> select(pred=V1)
    df <- cbind(df,pred)
    df$r <- r
    df$b <- if (is.na(best$`01`)){"None"} else {best$`01`}
    df$model <- best$model
    
    df
    
  } -> results_df
  
  results_df
  
})