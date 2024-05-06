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
  res <- residuals(model)
  n <- length(res)
  sigma <- sqrt(sum(res^2)/n)
  loglik_manual <- -(n/2)*(1+log(2*pi*sigma^2)) 
  loglik_manual
}

loglik_optim <- function(par,data){
  
  R <- data$R
  lagRBC <- data$lagRBC
  
  F0 <- par[1]
  theta <- par[2]
  k <- par[3]
  
  res <- R - sigmoid(lagRBC,F0,theta,k)
  n <- length(res)
  sigma <- sqrt(sum(res^2)/n)
  loglik_manual <- -(n/2)*(1+log(2*pi*sigma^2)) 
  loglik_manual
  
}

loglik_optim_box <- function(parlist,datalist){
  
  sapply(1:length(parlist),FUN=function(i){
    
    R <- datalist[[i]]$R
    lagRBC <- datalist[[i]]$lagRBC
    F0 <- parlist[[i]][1]
    theta <- parlist[[i]][2]
    k <- parlist[[i]][3]
    
    R - sigmoid(lagRBC,F0,theta,k)
    
  }) |> unlist() -> res
  
  n <- length(res)
  sigma <- sqrt(sum(res^2)/n)
  loglik_manual <- -(n/2)*(1+log(2*pi*sigma^2)) 
  loglik_manual
  
}

sigmoid <- function (lagRBC,F0,theta,k) {
  F0/(1+exp(k*(lagRBC-theta)/1e5))
}

objfun <- function(par,data){
  
  R <- data$R
  lagRBC <- data$lagRBC
  
  F0 <- par[1]
  theta <- par[2]
  k <- par[3]
  
  err <- R - sigmoid(lagRBC,F0,theta,k)
  
  sum(err^2)
  
}

#foreach (
  #r=seq_len(rep_num),
  #.options.future = list(seed = TRUE)
#) %dofuture% {

res <- list()

for (r in seq_len(rep_num)){
  
  print(r)
  
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
  
  for (lagChoice in lag_list) {
    
    ## create data-frame with lagged RBC
    dat <- joint_mouse_df |>
      group_by(mouseid) |>
      mutate(lagRBC=lag(RBC,lagChoice)) |>
      ungroup () |>
      na.omit()
    
    ## subset of data common to all lags
    dat_sub <- dat |> filter(time>max(lag_list)-1)
    
    for (i in seq_len(nrow(breakpoint_grid))) {
      
      bp <- breakpoint_grid[i,]
      
      if (is.na(bp[1])) {
        
        #Linear regression models
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
        
        bind_rows(
          full=models |>
            lapply(
              \(m) tibble(
                bp,
                lagChoice,
                rep=r,
                loglik=loglik(m),
                coef=length(m$coefficients),
                n=nrow(dat),
                p=summary(m)$df[3]+1,
                con=0,
              )
            ) |>
            bind_rows(.id="model"),
          sub=models_sub |>
            lapply(
              \(m) tibble(
                bp,
                lagChoice,
                rep=r,
                loglik=loglik(m),
                coef=length(m$coefficients),
                n=nrow(dat_sub),
                p=summary(m)$df[3]+1,
                con=0
              )
            ) |>
            bind_rows(.id="model"),
          .id="dataset"
        ) |>
          mutate(
            AIC=-2*loglik+2*p,
            AICc=AIC+2*p*(p+1)/(n-p-1)
          ) -> summary_lm
        
        #Sigmoid function without box
        m7 <- optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat,method="L-BFGS-B",
                    control=list(maxit=1000),
                    lower=c(F0=1e6,theta=2e6,k=0),
                    upper=c(F0=6e6,theta=8e6,k=1))
        
        if(m7$convergence!=0){break}
        
        m7_sub <- optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat_sub,method="L-BFGS-B",
                        control=list(maxit=1000),
                        lower=c(F0=1e6,theta=2e6,k=0),
                        upper=c(F0=6e6,theta=8e6,k=1))
        
        if(m7_sub$convergence!=0){break}
        
        bind_rows(
          tibble(
            dataset="full",
            model="m7",
            bp,
            lagChoice,
            rep=r,
            loglik=loglik_optim(m7$par,dat),
            coef=length(m7$par),
            n=nrow(dat),
            p=coef+1,
            con=m7$convergence
          ),
          tibble(
            dataset="sub",
            model="m7",
            bp,
            lagChoice,
            rep=r,
            loglik=loglik_optim(m7_sub$par,dat_sub),
            coef=length(m7_sub$par),
            n=nrow(dat_sub),
            p=coef+1,
            con=m7_sub$convergence
          )
        ) |>
          mutate(
            AIC=-2*loglik+2*p,
            AICc=AIC+2*p*(p+1)/(n-p-1)
          ) -> summary_m7
        
        #Hill function with box (separate variances)
        dat_1 <- filter(dat,box=="01")
        dat_2 <- filter(dat,box=="02")
        dat_3 <- filter(dat,box=="03")
        dat_4 <- filter(dat,box=="04")
        
        sig_box01 <-  optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat_1,method="L-BFGS-B",
                            control=list(maxit=1000),
                            lower=c(F0=1e6,theta=2e6,k=0),
                            upper=c(F0=6e6,theta=8e6,k=1))
        
        if(sig_box01$convergence!=0){break}
        
        sig_box02 <-  optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat_2,method="L-BFGS-B",
                            control=list(maxit=1000),
                            lower=c(F0=1e6,theta=2e6,k=0),
                            upper=c(F0=6e6,theta=8e6,k=1))
        
        if(sig_box02$convergence!=0){break}
        
        sig_box03 <-  optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat_3,method="L-BFGS-B",
                            control=list(maxit=1000),
                            lower=c(F0=1e6,theta=2e6,k=0),
                            upper=c(F0=6e6,theta=8e6,k=1))
        
        if(sig_box03$convergence!=0){break}
        
        sig_box04 <-  optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat_4,method="L-BFGS-B",
                            control=list(maxit=1000),
                            lower=c(F0=1e6,theta=2e6,k=0),
                            upper=c(F0=6e6,theta=8e6,k=1))
        
        if(sig_box04$convergence!=0){break}
        
        dat_sub1 <- filter(dat_sub,box=="01")
        dat_sub2 <- filter(dat_sub,box=="02")
        dat_sub3 <- filter(dat_sub,box=="03")
        dat_sub4 <- filter(dat_sub,box=="04")
        
        sig_sub_box01 <-  optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat_sub1,method="L-BFGS-B",
                                control=list(maxit=1000),
                                lower=c(F0=1e6,theta=2e6,k=0),
                                upper=c(F0=6e6,theta=8e6,k=1))
        
        if(sig_sub_box01$convergence!=0){break}
        
        sig_sub_box02 <-  optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat_sub2,method="L-BFGS-B",
                                control=list(maxit=1000),
                                lower=c(F0=1e6,theta=2e6,k=0),
                                upper=c(F0=6e6,theta=8e6,k=1))
        
        if(sig_sub_box02$convergence!=0){break}
        
        sig_sub_box03 <-  optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat_sub3,method="L-BFGS-B",
                                control=list(maxit=1000),
                                lower=c(F0=1e6,theta=2e6,k=0),
                                upper=c(F0=6e6,theta=8e6,k=1))
        
        if(sig_sub_box03$convergence!=0){break}
        
        sig_sub_box04 <-  optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat_sub4,method="L-BFGS-B",
                                control=list(maxit=1000),
                                lower=c(F0=1e6,theta=2e6,k=0),
                                upper=c(F0=6e6,theta=8e6,k=1))
        
        if(sig_sub_box04$convergence!=0){break}
        
        bind_rows(
          tibble(
            dataset="full",
            model="m8",
            bp,
            lagChoice,
            rep=r,
            loglik=loglik_optim_box(list(sig_box01$par,sig_box02$par,sig_box03$par,sig_box04$par),
                                    list(dat_1,dat_2,dat_3,dat_4)),
            coef=length(c(sig_box01$par,sig_box02$par,sig_box03$par,sig_box04$par)),
            n=nrow(dat_1)+nrow(dat_2)+nrow(dat_3)+nrow(dat_4),
            p=coef+4,
            con=max(sig_box01$convergence,sig_box02$convergence,sig_box03$convergence,sig_box04$convergence)
          ),
          tibble(
            dataset="sub",
            model="m8",
            bp,
            lagChoice,
            rep=r,
            loglik=loglik_optim_box(list(sig_sub_box01$par,sig_sub_box02$par,sig_sub_box03$par,sig_sub_box04$par),
                                    list(dat_sub1,dat_sub2,dat_sub3,dat_sub4)),
            coef=length(c(sig_sub_box01$par,sig_sub_box02$par,sig_sub_box03$par,sig_sub_box04$par)),
            n=nrow(dat_sub1)+nrow(dat_sub2)+nrow(dat_sub3)+nrow(dat_sub4),
            p=coef+4,
            con=max(sig_sub_box01$convergence,sig_sub_box02$convergence,sig_sub_box03$convergence,sig_sub_box04$convergence)
          )
        ) |>
          mutate(
            AIC=-2*loglik+2*p,
            AICc=AIC+2*p*(p+1)/(n-p-1)
          ) -> summary_m8
        
        model_df <- rbind(summary_lm,summary_m7,summary_m8)
        
        stats_df <- c(stats_df,list(model_df))
        
        ranges |>
          filter(
            lgg==lagChoice,
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
          bind_rows(.id="model") -> grid_lm
        
        pred <- sigmoid(grid$lagRBC,m7$par[1],m7$par[2],m7$par[3]) |> as.data.frame()
        names(pred) <- "pred"
        grid_m7 <- bind_cols(grid,pred) |>
          mutate(model="m7") |>
          select(model,lag,bp,box,phase,rep,lagRBC,pred)
        
        pred <- sigmoid(filter(grid,box=="01")$lagRBC,sig_box01$par[1],sig_box01$par[2],sig_box01$par[3]) |> as.data.frame()
        names(pred) <- "pred"
        grid_m8_01 <- bind_cols(filter(grid,box=="01"),pred) |>
          mutate(model="m8") |>
          select(model,lag,bp,box,phase,rep,lagRBC,pred)
        
        pred <- sigmoid(filter(grid,box=="02")$lagRBC,sig_box02$par[1],sig_box02$par[2],sig_box02$par[3]) |> as.data.frame()
        names(pred) <- "pred"
        grid_m8_02 <- bind_cols(filter(grid,box=="02"),pred) |>
          mutate(model="m8") |>
          select(model,lag,bp,box,phase,rep,lagRBC,pred)
        
        pred <- sigmoid(filter(grid,box=="03")$lagRBC,sig_box03$par[1],sig_box03$par[2],sig_box03$par[3]) |> as.data.frame()
        names(pred) <- "pred"
        grid_m8_03 <- bind_cols(filter(grid,box=="03"),pred) |>
          mutate(model="m8") |>
          select(model,lag,bp,box,phase,rep,lagRBC,pred)
        
        pred <- sigmoid(filter(grid,box=="04")$lagRBC,sig_box04$par[1],sig_box04$par[2],sig_box04$par[3]) |> as.data.frame()
        names(pred) <- "pred"
        grid_m8_04 <- bind_cols(filter(grid,box=="04"),pred) |>
          mutate(model="m8") |>
          select(model,lag,bp,box,phase,rep,lagRBC,pred)
        
        grid_m8 <- bind_rows(grid_m8_01,grid_m8_02,grid_m8_03,grid_m8_04) |>
          mutate(model="m8") |>
          select(model,lag,bp,box,phase,rep,lagRBC,pred)
        
        grid <- bind_rows(grid_lm,grid_m7,grid_m8)
        preds_df <- c(preds_df,list(grid))
        
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
        
        bind_rows(
          full=models |>
            lapply(
              \(m) tibble(
                bp,
                lagChoice,
                rep=r,
                loglik=loglik(m),
                coef=length(m$coefficients),
                n=nrow(dat),
                p=summary(m)$df[3]+1,
                con=0
              )
            ) |>
            bind_rows(.id="model"),
          sub=models_sub |>
            lapply(
              \(m) tibble(
                bp,
                lagChoice,
                rep=r,
                loglik=loglik(m),
                coef=length(m$coefficients),
                n=nrow(dat_sub),
                p=summary(m)$df[3]+1,
                con=0
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
            lgg==lagChoice,
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
        
      }
      
    } #end loop over breakpoints
  } #end loop over lags
  
  tmp <- list(
    rep=r,
    stats=bind_rows(stats_df),
    preds=bind_rows(preds_df)
  )
  
  res <- append(res,tmp)
  
} #-> res

lapply(which(names(res)=="stats"),\(x) res[[x]]) |>
  bind_rows() -> stats_df

lapply(which(names(res)=="preds"),\(x) res[[x]]) |>
  bind_rows() -> preds_df

#res |>
  #lapply(\(x) getElement(x,"stats")) |>
  #bind_rows() -> stats_df
#res |>
  #lapply(\(x) getElement(x,"preds")) |>
  #bind_rows() -> preds_df

preds_df |>
  group_by(model,lag,bp,box,phase,lagRBC) |>
  reframe(p=c(0.1,0.5,0.9),q=quantile(pred,probs=p),label=c("lo","med","hi")) |>
  ungroup() |>
  pivot_wider(id_cols=c(model,lag,bp,box,phase,lagRBC),names_from=label,values_from=q) -> quants_df

stats_df |> write_csv("results_regression_stats.csv")
## preds_df |> write_csv("results_regression_preds.csv")
quants_df |> write_csv("results_regression_quants.csv")
