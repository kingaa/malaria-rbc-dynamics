#SET DIRECTORY TO SOURCE FILE LOCATION

library(tidyverse)
library(pomp)
library(stringi)
library(doParallel)
library(panelPomp)
library(ggpubr)
library(iterators)
library(doFuture)
plan(multisession) #for faster results, use multicore outside of RStudio

source("POMP_GroupLevel_DataPrep.R")

head(flow)

#Set up data
flow |>
  select(day,mouseid,box,pABA,E=Eryth,R=Retic) |>
  mutate(lagE=lag(E)) |>
  filter(box!="05",day>0,day<=21,mouseid!="01-02",mouseid!="02-03") -> data

##Create data frame with breakpoints for each of the four pABA boxes
box_list <- unique(data$box)
breakpoint_grid <- as_tibble(
  matrix(
    data=rep(c(NA,8,9,10,11),each=4),
    ncol=4,byrow=TRUE,
    dimnames=list(NULL,box_list)
  )
)

##Model fitting
bake(file="jAICdf_data.rds",{
  foreach (
    bp=iter(breakpoint_grid,"row"),
    .combine=rbind
  ) %dofuture% {
    
    if (is.na(bp[1])){
      
      data -> df
      
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
        bind_rows(.id="model")
      
    } else {
      
      data |>
        mutate(
          phase=as.factor(if_else(day<=unlist(bp)[box],1,2))
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
        bind_rows(.id="model")
      
    }
  }
}) -> joint_AIC_df

##Extract breakpoints from best model (lowest AIC)
joint_AIC_df |>
  filter(AIC==min(AIC)) -> best

##Obtain data frame with above breakpoints specified
data |>
  mutate(
    phase=as.factor(if_else(day<=as.integer(unlist(best)[box]),1,2))
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

#Split data based on phase
df_phase1 <- df |> filter(phase==1)
df_phase2 <- df |> filter(phase==2)

#Generate new data to give to predict (requires lagE column and phase column)
newdata <- rbind(
  expand.grid(
    lagE=seq(min(df_phase1$lagE), max(df_phase1$lagE), 10000),
    phase=1
  ),
  expand.grid(
    lagE=seq(min(df_phase2$lagE), max(df_phase2$lagE), 10000),
    phase=2
  )
)
newdata$phase <- factor(newdata$phase)

pred <- predict(chosen_model,newdata=newdata, interval = "confidence", level = 0.95)

newdata_pred <- bind_cols(newdata,pred)

newdata_pred |>
  ggplot()+
  geom_line(aes(x=lagE,y=fit,group=phase))+
  geom_ribbon(aes(x=lagE,ymin=lwr,ymax=upr,group=phase),alpha=0.2)+
  theme_bw()+
  xlab("Erythrocyte density (d-1)")+ylab("Reticulocyte supply (d)")+
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
