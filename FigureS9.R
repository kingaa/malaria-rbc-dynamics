#SET DIRECTORY TO SOURCE FILE LOCATION

library(tidyverse)
library(stringi)
library(doParallel)
library(pomp)
library(panelPomp)
library(ggpubr)
library(scales)
library(aakmisc)
library(gridExtra)
library(cowplot)

cbpABA <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

source("POMP_GroupLevel_DataPrep.R")

flow |> 
  filter(day<=20,box!="05",mouseid!="01-02",mouseid!="02-03") |>
  select(time=day,R=Retic,RBC,box,mouse,pABA,mouseid) -> data


box_list <- unique(data$box)
breakpoint_grid <- as_tibble(
  matrix(
    data=rep(c(NA,8,9,10,11),each=4),
    ncol=4,byrow=TRUE,
    dimnames=list(NULL,box_list)
  )
)

lag_list <- 1:5

#Create function to calculate log likelihood of model
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
  F0/(1+exp(k*(lagRBC-theta)/1e7))
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

joint_AIC_df <- data.frame()
for (lagChoice in lag_list){
  
  df <- data |>
    group_by(mouseid) |>
    mutate(lagRBC=lag(RBC,lagChoice)) |>
    na.omit()
  
  df_sub <- df |> filter(time>max(lag_list)-1)
  
  for (i in seq(1,nrow(breakpoint_grid),1)){
    
    bp <- breakpoint_grid[i,]
    
    if (is.na(bp[1])){

      #Linear regression models
      models <- list(
        m1=lm(R~poly(lagRBC,2,raw=F):(box-1)+(box-1),data=df), #box with different shape and intercept
        m2=lm(R~(poly(lagRBC,2,raw=F)+box-1),data=df), #same shape, different intercept
        m3=lm(R~poly(lagRBC,2,raw=F),data=df), #no effect
        m4=lm(R~lagRBC:(box-1)+(box-1),data=df), #slope and intercept
        m5=lm(R~(lagRBC+box-1),data=df), #same slope, different intercept
        m6=lm(R~lagRBC,data=df) #no effect
      )
      
      models_sub <- list(
        m1=lm(R~poly(lagRBC,2,raw=F):(box-1)+(box-1),data=df_sub),
        m2=lm(R~(poly(lagRBC,2,raw=F)+box-1),data=df_sub),
        m3=lm(R~poly(lagRBC,2,raw=F),data=df_sub),
        m4=lm(R~lagRBC:(box-1)+(box-1),data=df_sub),
        m5=lm(R~(lagRBC+box-1),data=df_sub),
        m6=lm(R~lagRBC,data=df_sub)
      )
      
      bind_rows(
        full=models |>
          lapply(
            \(m) tibble(
              bp,
              lagChoice,
              loglik=loglik(m),
              coef=length(m$coefficients),
              n=nrow(df),
              p=summary(m)$df[3]+1
            )
          ) |>
          bind_rows(.id="model"),
        sub=models_sub |>
          lapply(
            \(m) tibble(
              bp,
              lagChoice,
              loglik=loglik(m),
              coef=length(m$coefficients),
              n=nrow(df_sub),
              p=summary(m)$df[3]+1
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
      m7 <- optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=df,method="L-BFGS-B",
                  lower=c(F0=1e6,theta=2e6,k=0),
                  upper=c(F0=6e6,theta=8e6,k=1))
      
      if(m7$convergence!=0){break}
      
      m7_sub <- optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=df_sub,method="L-BFGS-B",
                      lower=c(F0=1e6,theta=2e6,k=0),
                      upper=c(F0=6e6,theta=8e6,k=1))
      
      if(m7_sub$convergence!=0){break}
      
      bind_rows(
        tibble(
          dataset="full",
          model="m7",
          bp,
          lagChoice,
          loglik=loglik_optim(m7$par,df),
          coef=length(m7$par),
          n=nrow(df),
          p=coef+1,
        ),
        tibble(
          dataset="sub",
          model="m7",
          bp,
          lagChoice,
          loglik=loglik_optim(m7_sub$par,df_sub),
          coef=length(m7_sub$par),
          n=nrow(df_sub),
          p=coef+1,
        )
      ) |>
        mutate(
          AIC=-2*loglik+2*p,
          AICc=AIC+2*p*(p+1)/(n-p-1)
        ) -> summary_m7
      
      #Hill function with box (separate variances)
      dat_1 <- filter(df,box=="01")
      dat_2 <- filter(df,box=="02")
      dat_3 <- filter(df,box=="03")
      dat_4 <- filter(df,box=="04")
      
      sig_box01 <-  optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat_1,method="L-BFGS-B",
                          lower=c(F0=1e6,theta=2e6,k=0),
                          upper=c(F0=6e6,theta=8e6,k=1))
      
      if(sig_box01$convergence!=0){break}
      
      sig_box02 <-  optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat_2,method="L-BFGS-B",
                          lower=c(F0=1e6,theta=2e6,k=0),
                          upper=c(F0=6e6,theta=8e6,k=1))
      
      if(sig_box02$convergence!=0){break}
      
      sig_box03 <-  optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat_3,method="L-BFGS-B",
                          lower=c(F0=1e6,theta=2e6,k=0),
                          upper=c(F0=6e6,theta=8e6,k=1))
      
      if(sig_box03$convergence!=0){break}
      
      sig_box04 <-  optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat_4,method="L-BFGS-B",
                          lower=c(F0=1e6,theta=2e6,k=0),
                          upper=c(F0=6e6,theta=8e6,k=1))
      
      if(sig_box04$convergence!=0){break}
      
      dat_sub1 <- filter(df_sub,box=="01")
      dat_sub2 <- filter(df_sub,box=="02")
      dat_sub3 <- filter(df_sub,box=="03")
      dat_sub4 <- filter(df_sub,box=="04")
      
      sig_sub_box01 <-  optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat_sub1,method="L-BFGS-B",
                              lower=c(F0=1e6,theta=2e6,k=0),
                              upper=c(F0=6e6,theta=8e6,k=1))
      
      if(sig_sub_box01$convergence!=0){break}
      
      sig_sub_box02 <-  optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat_sub2,method="L-BFGS-B",
                              lower=c(F0=1e6,theta=2e6,k=0),
                              upper=c(F0=6e6,theta=8e6,k=1))
      
      if(sig_sub_box02$convergence!=0){break}
      
      sig_sub_box03 <-  optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat_sub3,method="L-BFGS-B",
                              lower=c(F0=1e6,theta=2e6,k=0),
                              upper=c(F0=6e6,theta=8e6,k=1))
      
      if(sig_sub_box03$convergence!=0){break}
      
      sig_sub_box04 <-  optim(c(F0=4e6,theta=5.5e6,k=0.1),objfun,data=dat_sub4,method="L-BFGS-B",
                              lower=c(F0=1e6,theta=2e6,k=0),
                              upper=c(F0=6e6,theta=8e6,k=1))
      
      if(sig_sub_box04$convergence!=0){break}
      
      bind_rows(
        tibble(
          dataset="full",
          model="m8",
          bp,
          lagChoice,
          loglik=loglik_optim_box(list(sig_box01$par,sig_box02$par,sig_box03$par,sig_box04$par),
                                  list(dat_1,dat_2,dat_3,dat_4)),
          coef=length(c(sig_box01$par,sig_box02$par,sig_box03$par,sig_box04$par)),
          n=nrow(dat_1)+nrow(dat_2)+nrow(dat_3)+nrow(dat_4),
          p=coef+4,
        ),
        tibble(
          dataset="sub",
          model="m8",
          bp,
          lagChoice,
          loglik=loglik_optim_box(list(sig_sub_box01$par,sig_sub_box02$par,sig_sub_box03$par,sig_sub_box04$par),
                                  list(dat_sub1,dat_sub2,dat_sub3,dat_sub4)),
          coef=length(c(sig_sub_box01$par,sig_sub_box02$par,sig_sub_box03$par,sig_sub_box04$par)),
          n=nrow(dat_sub1)+nrow(dat_sub2)+nrow(dat_sub3)+nrow(dat_sub4),
          p=coef+4,
        )
      ) |>
        mutate(
          AIC=-2*loglik+2*p,
          AICc=AIC+2*p*(p+1)/(n-p-1)
        ) -> summary_m8
      
      model_df <- rbind(summary_lm,summary_m7,summary_m8)
      
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
      
      bind_rows(
        full=models |>
          lapply(
            \(m) tibble(
              bp,
              lagChoice,
              loglik=loglik(m),
              coef=length(m$coefficients),
              n=nrow(df),
              p=summary(m)$df[3]+1
            )
          ) |>
          bind_rows(.id="model"),
        sub=models_sub |>
          lapply(
            \(m) tibble(
              bp,
              lagChoice,
              loglik=loglik(m),
              coef=length(m$coefficients),
              n=nrow(df_sub),
              p=summary(m)$df[3]+1
            )
          ) |>
          bind_rows(.id="model"),
        .id="dataset"
      ) |>
        mutate(
          AIC=-2*loglik+2*p,
          AICc=AIC+2*p*(p+1)/(n-p-1)
        ) -> summary_lm
      
      model_df <- summary_lm
      
    }
    
    joint_AIC_df <- rbind(joint_AIC_df,model_df)
  } #end loop over breakpoints
}

##Extract breakpoints from best model (lowest AIC)
joint_AIC_df |>
  filter(AICc==min(AICc)) -> best

##Obtain data frame with above breakpoints specified
data |>
  group_by(mouseid) |>
  mutate(
    lagRBC=lag(RBC,best$lagChoice),
  ) |>
  na.omit() |>
  mutate(phase=as.factor(if_else(time<=as.integer(unlist(best)[box]),1,2))) |>
  ungroup() -> df

##Extract best model based on lowest AIC
chosen_model <- lm(R~(lagRBC+box-1):(phase-1)+(phase-1),data=df)
coefs <- chosen_model$coefficients

pred <- model.matrix(chosen_model) %*% coefs |> as.data.frame() |> select(pred=V1)
df <- cbind(df,pred)
df$b <- best$`01`
df$model <- best$model
df$lag <- best$lagChoice
df$loglik <- best$loglik

df |>
  select(pABA,lagRBC,phase,pred) |>
  ggplot()+
  geom_line(aes(x=lagRBC,y=pred,group=interaction(phase,pABA),col=pABA,linetype=phase),linewidth=1)+
  scale_x_continuous(labels = aakmisc::scinot)+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4500000))+
  #annotate("text",x=7000000,y=3200000,label="Phase 2",size=5)+
  #annotate("text",x=6000000,y=550000,label="Phase 1",size=5)+
  #geom_segment(aes(x=6000000,xend=6000000,y=1100000,yend=600000))+
  #geom_segment(aes(x=7000000,xend=5900000,y=3000000,yend=2500000))+
  scale_colour_manual(values=cbpABA[2:5])+
  geom_label(aes(x=9000000,y=4000000,label="Model E\nBreakpoint 10\nLag = 3 days"))+
  theme_bw()+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  labs(x=expression(paste("RBC density at time ", italic("t"), "-3 (density per µL)")),
       y=expression(paste("Reticulocyte supply at time ", italic("t"), " (density per µL)")),
                    colour="Parasite nutrient (pABA)",linetype="Phase")+
  theme(
    axis.title=element_text(size=15),
    axis.title.y=element_text(size=15,colour="black"),
    axis.text=element_text(size=12),
    axis.text.x=element_text(colour="black"),
    panel.grid.minor=element_blank(),
    legend.position=c(0.225,0.25),
    legend.background=element_blank(),
    legend.title=element_text(size=13),
    legend.text=element_text(size=10),
    strip.text=element_text(size=12),
    strip.background=element_blank(),
    plot.title=element_text(size=14,hjust=0.5,face="bold")#,
    #panel.border = element_rect(linewidth=4)
    
  )
ggsave("FigureS9.jpeg",width=15,height=15,units="cm")
