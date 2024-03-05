#SET DIRECTORY TO SOURCE FILE LOCATION
#Load required packages
library(tidyverse)
library(pomp)
library(aakmisc)
library(foreach)
library(doFuture)
plan(multisession,workers=3)

lag_list <- 1:5
rep_num <- 1000

##Load in smooth distribution samples from PNAS work
sm1 <- readRDS("m5sm1_mod.rds")

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
breakpoint_grid <- as_tibble(
  matrix(
    data=rep(c(NA,8,9,10,11),each=4),
    ncol=4,byrow=TRUE,
    dimnames=list(NULL,unique(sm1_mod$box))
  )
)

expand_grid(
  lgg=lag_list,
  brp=c(NA,8,9,10,11),
  sm1_mod
) |>
  group_by(lgg,mouseid) |>
  mutate(
    lagRBC=lag(RBC,n=unique(lgg)),
    phase=factor(if_else(time<=brp,1,2))
  ) |>
  ungroup() |>
  filter(!is.na(lagRBC)) |>
  group_by(lgg,brp,box,phase) |>
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

foreach (
  r=seq_len(rep_num),
  .options.future = list(seed = TRUE)
) %dofuture% {
  
  stats_df <- list()
  preds_df <- list()
  
  lapply(
    mouse_id_list,
    \(id) {
      key_opts <- sm1_mod |> filter(mouseid==id) |> select(key,lik) |> unique()
      key_choice <- sample(key_opts$key, size=1, prob=key_opts$lik)
      sm1_mod_mouse <- sm1_mod |> filter(key==key_choice)
    }
  ) |>
    bind_rows() -> joint_mouse_df
  
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
          m1=lm(R~poly(lagRBC,2,raw=F):(box-1)+(box-1),data=dat), #box with different shape and intercept
          m2=lm(R~(poly(lagRBC,2,raw=F)+box-1),data=dat), #same shape, different intercept
          m3=lm(R~poly(lagRBC,2,raw=F),data=dat), #no effect
          m4=lm(R~lagRBC:(box-1)+(box-1),data=dat), #slope and intercept
          m5=lm(R~(lagRBC+box-1),data=dat), #same slope, different intercept
          m6=lm(R~lagRBC,data=dat) #no effect
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
              rep=r,
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
              rep=r,
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
      
      stats_df <- c(stats_df,list(model_df))
      
      ranges |>
        filter(
          lgg==lag,
          brp==bp[["01"]] |
            (is.na(bp[["01"]]) & is.na(brp))
        ) |>
        group_by(lag=lgg,bp=brp,box,phase) |>
        reframe(
          rep=r,
          lagRBC=seq(from=lo,to=hi,length=100)
        ) |>
        ungroup() -> grid
      
      models |>
        lapply(
          \(m) bind_cols(grid,pred=predict(m,newdata=grid))
        ) |>
        bind_rows(.id="model") -> grid
      
      preds_df <- c(preds_df,list(grid))
      
    } #end loop over breakpoints
  } #end loop over lags
  
  list(
    rep=r,
    stats=bind_rows(stats_df),
    preds=bind_rows(preds_df)
  )
} -> res

res |>
  lapply(\(x) getElement(x,"stats")) |>
  bind_rows() -> stats_df
res |>
  lapply(\(x) getElement(x,"preds")) |>
  bind_rows() -> preds_df

preds_df |>
  group_by(model,lag,bp,box,phase,lagRBC) |>
  reframe(p=c(0.1,0.5,0.9),q=quantile(pred,probs=p),label=c("lo","med","hi")) |>
  ungroup() |>
  pivot_wider(id_cols=c(model,lag,bp,box,phase,lagRBC),names_from=label,values_from=q) -> quants_df

stats_df |> write_csv("results_regression_stats.csv")
## preds_df |> write_csv("results_regression_preds.csv")
quants_df |> write_csv("results_regression_quants.csv")
