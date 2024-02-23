#SET DIRECTORY TO SOURCE FILE LOCATION
#Load required packages
library(tidyverse)
library(pomp)
library(aakmisc)

lag_list <- 1:5
rep_num <- 20

##Create data frame with breakpoints for each of the four pABA boxes
breakpoint_grid <- as_tibble(
  matrix(
    data=rep(c(NA,8,9,10,11),each=4),
    ncol=4,byrow=TRUE,
    dimnames=list(NULL,unique(sm1_mod$box))
  )
)
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

expand_grid(
  lg=lag_list,
  bp=c(NA,8,9,10,11),
  sm1_mod
) |>
  group_by(lg,mouseid) |>
  mutate(
    lagRBC=lag(RBC,n=unique(lg)),
    phase=factor(if_else(time<=bp,1,2))
  ) |>
  ungroup() |>
  filter(!is.na(lagRBC)) |>
  group_by(lag=lg,bp,box,phase) |>
  summarize(
    lo=quantile(lagRBC,probs=0.05),
    hi=quantile(lagRBC,probs=0.95)
  ) |>
  ungroup() -> ranges

#Set seed (important because random sampling occurs below)
seed_choice <- 851657743
set.seed(seed_choice)

mouse_id_list <- unique(sm1_mod$mouseid)

loglik <- function (model) {
  res <- model$residuals
  n <- length(res)
  sigma <- sqrt(sum(res^2)/n)
  loglik_manual <- -(n/2)*(1+log(2*pi*sigma^2)) 
  loglik_manual
}

stats_df <- list()
preds_df <- list()


for (r in seq_len(rep_num)) {
  
  print(r)
  
  lapply(
    mouse_id_list,
    \(id) {
      key_opts <- sm1_mod |> filter(mouseid==id) |> select(key,lik) |> unique()
      key_choice <- sample(key_opts$key, size=1, prob=key_opts$lik)
      sm1_mod_mouse <- sm1_mod |> filter(key==key_choice)
    }
  ) |>
    bind_rows() -> joint_mouse_df
  
  joint_AIC_df <- list()
  
  for (lag in lag_list) {
    
    ## create data-frame with lagged RBC
    dat <- joint_mouse_df |>
      group_by(mouseid) |>
      mutate(lagRBC=lag(RBC,lag)) |>
      ungroup () |>
      na.omit()
    
    ## subset of data common to all lags
    dat_sub <- dat |> filter(time>max(lag_list)-1)
    
    for (i in seq_len(nrow(breakpoint_grid))) {
      
      bp <- breakpoint_grid[i,]
      
      if (is.na(bp[1])) {
        
        ##Models
        models <- list(
          m1=lm(R~poly(lagRBC,2,raw=F):(box-1)+(box-1),data=dat),
          m2=lm(R~(poly(lagRBC,2,raw=F)+box-1),data=dat),
          m3=lm(R~poly(lagRBC,2,raw=F),data=dat),
          m4=lm(R~lagRBC:(box-1)+(box-1),data=dat),
          m5=lm(R~(lagRBC+box-1),data=dat),
          m6=lm(R~lagRBC,data=dat)
        )
        
        models_sub <- list(
          m1=lm(R~poly(lagRBC,2,raw=F):(box-1)+(box-1),data=dat_sub),
          m2=lm(R~(poly(lagRBC,2,raw=F)+box-1),data=dat_sub),
          m3=lm(R~poly(lagRBC,2,raw=F),data=dat_sub),
          m4=lm(R~lagRBC:(box-1)+(box-1),data=dat_sub),
          m5=lm(R~(lagRBC+box-1),data=dat_sub),
          m6=lm(R~lagRBC,data=dat_sub)
        )
        
      } else {
        
        dat |>
          mutate(
            phase=factor(if_else(time<=unlist(bp)[box],1,2))
          ) -> dat
        
        dat_sub <- dat |> filter(time>max(lag_list)-1)
        
        ##Models
        models <- list(
          m1=lm(R~poly(lagRBC,2,raw=F):(box-1):(phase-1)+(box-1):(phase-1),data=dat),
          m2=lm(R~(poly(lagRBC,2,raw=F)+box-1):(phase-1)+(phase-1),data=dat),
          m3=lm(R~poly(lagRBC,2,raw=F):(phase-1)+(phase-1),data=dat),
          m4=lm(R~lagRBC:(box-1):(phase-1)+(box-1):(phase-1),data=dat),
          m5=lm(R~(lagRBC+box-1):(phase-1)+(phase-1),data=dat),
          m6=lm(R~lagRBC:(phase-1)+(phase-1),data=dat)
        )
        
        models_sub <- list(
          m1=lm(R~poly(lagRBC,2,raw=F):(box-1):(phase-1)+(box-1):(phase-1),data=dat_sub),
          m2=lm(R~(poly(lagRBC,2,raw=F)+box-1):(phase-1)+(phase-1),data=dat_sub),
          m3=lm(R~poly(lagRBC,2,raw=F):(phase-1)+(phase-1),data=dat_sub),
          m4=lm(R~lagRBC:(box-1):(phase-1)+(box-1):(phase-1),data=dat_sub),
          m5=lm(R~(lagRBC+box-1):(phase-1)+(phase-1),data=dat_sub),
          m6=lm(R~lagRBC:(phase-1)+(phase-1),data=dat_sub)
        )
        
      }
      
      bind_rows(
        full=models |>
          lapply(
            \(m) tibble(
              bp,
              lag,
              loglik=loglik(m),
              coef=length(m$coefficients)+1,
              n=nrow(dat),
              p=summary(m)$df[3]+1
            )
          ) |>
          bind_rows(.id="model"),
        sub=models_sub |>
          lapply(
            \(m) tibble(
              bp,
              lag,
              loglik=loglik(m),
              coef=length(m$coefficients)+1,
              n=nrow(dat_sub),
              p=summary(m)$df[3]+1
            )
          ) |>
          bind_rows(.id="model"),
        .id="dataset"
      ) |>
        mutate(
          AIC=-2*loglik+2*p,
          AICc=AIC+2*p*(p+1)/(n-p-1)
        ) -> model_df
      
      joint_AIC_df <- c(joint_AIC_df,list(model_df))
    } #end loop over breakpoints
  } #end loop over lags
  
  joint_AIC_df |>
    bind_rows() -> joint_AIC_df
  
  stats_df <- c(stats_df,list(joint_AIC_df))
  
  ##Extract breakpoints from best model (lowest AIC)
  joint_AIC_df |>
    filter(dataset=="sub") |>
    filter(AICc==min(AICc)) -> best
  
  stopifnot(
    "unequal breakpoints"=length(unique(unlist(best[c("01","02","03","04")])))==1
  )
  
  ##Obtain data frame with above breakpoints specified
  joint_mouse_df |>
    group_by(mouseid) |>
    mutate(lagRBC=lag(RBC,best$lag)) |>
    ungroup () |>
    na.omit() |>
    mutate(
      phase=factor(if_else(time<=as.integer(unlist(best)[box]),1,2))
    ) -> df
  
  ##Extract best model based on lowest AIC
  if (is.na(best$`01`)) {
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
  
  print(best)
  
  chosen_model <- models[[best$model]]
  coefs <- chosen_model$coefficients
  
  ranges |>
    filter(
      lag==best$lag,
      bp==best[["01"]] | 
        (is.na(best[["01"]]) & is.na(bp))
    ) |>
    group_by(lag,bp,box,phase) |>
    reframe(
      model=best$model,
      lagRBC=seq(from=lo,to=hi,length=1000)
    ) |>
    ungroup() -> grid
  
  grid |>
    bind_cols(pred=predict(chosen_model,newdata=grid)) -> grid
  
  preds_df <- c(preds_df,list(grid))
 
} #end for loop over iterations

stats_df |> bind_rows(.id="rep") -> stats_df
preds_df |> bind_rows(.id="rep") -> preds_df

preds_df |>
  ggplot(aes(x=lagRBC,y=pred,group=interaction(phase,box),color=box))+
  geom_line()+
  facet_wrap(~model+lag,labeller=label_both)

write.csv(stats_df,"results_regression_stats_corr.csv",row.names=FALSE)
write.csv(preds_df,"results_regression_preds_corr.csv",row.names=FALSE)
