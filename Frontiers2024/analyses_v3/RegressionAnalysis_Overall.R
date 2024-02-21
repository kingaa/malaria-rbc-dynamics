#SET DIRECTORY TO SOURCE FILE LOCATION

#Load required packages
library(tidyverse)
library(pomp)
library(aakmisc)

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

lag_list <- 1:4

#Loop over reps
rep_num <- 1000

stats_df <- data.frame()
preds_df <- data.frame()

for (r in 1:rep_num){
  
  print(r)
  
  mouse_id_list <- unique(sm1_mod$mouseid)
  
  joint_mouse_df <- data.frame()
  for (id in mouse_id_list){
    
    key_opts <- sm1_mod |> filter(mouseid==id) |> select(key,lik) |> unique()
    
    key_choice <- sample(key_opts$key,size=1,prob=key_opts$lik)
    
    sm1_mod_mouse <- sm1_mod |> 
      filter(key==key_choice) #|>
      #mutate(lagRBC=lag(RBC,3)) |>
      #na.omit()
    
    joint_mouse_df <- rbind(joint_mouse_df,sm1_mod_mouse)
    
  }
  
  loglik <- function(model){
    res <- model$residuals
    n <- length(res)
    sigma <- sqrt(sum(res^2)/n)
    loglik_manual <- -(n/2)*(1+log(2*pi*sigma^2)) 
    return(loglik_manual)
  }
  
  joint_AIC_df <- data.frame()
  for (lag in lag_list){

    df <- joint_mouse_df |>
      mutate(lagRBC=lag(RBC,lag)) |>
      na.omit()
    
    df_sub <- df |> filter(time>max(lag_list)-1)

    for (i in seq(1,nrow(breakpoint_grid),1)){
      
      bp <- breakpoint_grid[i,]
      
      if (is.na(bp[1])){
        
        ##Models
        models <- list(
          m1=lm(R~poly(lagRBC,2,raw=F):(box-1)+(box-1),data=df),
          m2=lm(R~(poly(lagRBC,2,raw=F)+box-1),data=df),
          m3=lm(R~poly(lagRBC,2,raw=F),data=df),
          m4=lm(R~lagRBC:(box-1)+(box-1),data=df),
          m5=lm(R~(lagRBC+box-1),data=df),
          m6=lm(R~lagRBC,data=df)
        )
        
        models_sub <- list(
          m1=lm(R~poly(lagRBC,2,raw=F):(box-1)+(box-1),data=df_sub),
          m2=lm(R~(poly(lagRBC,2,raw=F)+box-1),data=df_sub),
          m3=lm(R~poly(lagRBC,2,raw=F),data=df_sub),
          m4=lm(R~lagRBC:(box-1)+(box-1),data=df_sub),
          m5=lm(R~(lagRBC+box-1),data=df_sub),
          m6=lm(R~lagRBC,data=df_sub)
        )
        
        models |>
          lapply(\(m) tibble(bp,
                             lag,
                             loglik_total=loglik(m),
                             coef_total=length(m$coefficients),
                             n_total=summary(m)$df[2]+summary(m)$df[3],
                             p_total=summary(m)$df[3],
                             AIC_total=AIC(m))) |>
          bind_rows(.id="model") -> model_df
        
        models_sub |>
          lapply(\(m) tibble(
                             loglik_sub=loglik(m),
                             coef_sub=length(m$coefficients),
                             n_sub=summary(m)$df[2]+summary(m)$df[3],
                             p_sub=summary(m)$df[3],
                             AIC_sub=AIC(m))) |>
          bind_rows(.id="model") -> model_sub
        
        model_df <- full_join(model_df,model_sub,by="model")
        
      } else {
        
        df |>
          mutate(
            phase=as.factor(if_else(time<=unlist(bp)[box],1,2))
          ) -> df
        
        df_sub <- df |> filter(time>max(lag_list)-1)
        
        ##Models
        models <- list(
          m1=lm(R~poly(lagRBC,2,raw=F):(box-1):(phase-1)+(box-1):(phase-1),data=df),
          m2=lm(R~(poly(lagRBC,2,raw=F)+box-1):(phase-1)+(phase-1),data=df),
          m3=lm(R~poly(lagRBC,2,raw=F):(phase-1)+(phase-1),data=df),
          m4=lm(R~lagRBC:(box-1):(phase-1)+(box-1):(phase-1),data=df),
          m5=lm(R~(lagRBC+box-1):(phase-1)+(phase-1),data=df),
          m6=lm(R~lagRBC:(phase-1)+(phase-1),data=df)
        )
        
        models_sub <- list(
          m1=lm(R~poly(lagRBC,2,raw=F):(box-1):(phase-1)+(box-1):(phase-1),data=df_sub),
          m2=lm(R~(poly(lagRBC,2,raw=F)+box-1):(phase-1)+(phase-1),data=df_sub),
          m3=lm(R~poly(lagRBC,2,raw=F):(phase-1)+(phase-1),data=df_sub),
          m4=lm(R~lagRBC:(box-1):(phase-1)+(box-1):(phase-1),data=df_sub),
          m5=lm(R~(lagRBC+box-1):(phase-1)+(phase-1),data=df_sub),
          m6=lm(R~lagRBC:(phase-1)+(phase-1),data=df_sub)
        )
        
        models |>
          lapply(\(m) tibble(bp,
                             lag,
                             loglik_total=loglik(m),
                             coef_total=length(m$coefficients),
                             n_total=summary(m)$df[2]+summary(m)$df[3],
                             p_total=summary(m)$df[3],
                             AIC_total=AIC(m))) |>
          bind_rows(.id="model") -> model_df
        
        models_sub |>
          lapply(\(m) tibble(
            loglik_sub=loglik(m),
            coef_sub=length(m$coefficients),
            n_sub=summary(m)$df[2]+summary(m)$df[3],
            p_sub=summary(m)$df[3],
            AIC_sub=AIC(m))) |>
          bind_rows(.id="model") -> model_sub
        
        model_df <- full_join(model_df,model_sub,by="model")
        
      }
      
      joint_AIC_df <- rbind(joint_AIC_df,model_df)
    } #end loop over breakpoints
  } #end loop over lags

  joint_AIC_df$rep <- r
  joint_AIC_df <- joint_AIC_df |>
    mutate(AICc = AIC_sub+2*p_sub*(p_sub+1)/(n_sub-p_sub-1))
  
  ##Extract breakpoints from best model (lowest AIC)
  joint_AIC_df |>
    filter(AICc==min(AICc)) -> best
  
  ##Obtain data frame with above breakpoints specified
  joint_mouse_df |>
    mutate(
      lagRBC=lag(RBC,best$lag),
    ) |>
    na.omit() |>
    mutate(phase=as.factor(if_else(time<=as.integer(unlist(best)[box]),1,2))) -> df
  
  
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
  df$rep <- r
  df$b <- best$`01`
  df$model <- best$model
  df$lag <- best$lag
  df$loglik_total <- best$loglik_total
  df$loglik_sub <- best$loglik_sub
  
  stats_df <- bind_rows(stats_df,joint_AIC_df)
  preds_df <- bind_rows(preds_df,df)
    
} #end for loop over iterations

write.csv(stats_df,"results_regression_stats.csv",row.names=FALSE)
write.csv(preds_df,"results_regression_preds.csv",row.names=FALSE)

#df |>
  #ggplot()+
  #geom_line(aes(x=lagRBC,y=pred,group=interaction(phase,box),col=phase))

#df |>
  #ggplot()+
  #geom_text(aes(x=lagRBC,y=R,label=time,group=interaction(phase,mouse,box),col=phase))
