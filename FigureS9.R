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

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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
loglik <- function(model){
  res <- model$residuals
  n <- length(res)
  sigma <- sqrt(sum(res^2)/n)
  loglik_manual <- -(n/2)*(1+log(2*pi*sigma^2)) 
  return(loglik_manual)
}

joint_AIC_df <- data.frame()
for (lag in lag_list){
  
  df <- data |>
    group_by(mouseid) |>
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
}
joint_AIC_df <- joint_AIC_df |>
  mutate(AICc = AIC_sub+2*p_sub*(p_sub+1)/(n_sub-p_sub-1))

##Extract breakpoints from best model (lowest AIC)
joint_AIC_df |>
  filter(AICc==min(AICc)) -> best

##Obtain data frame with above breakpoints specified
data |>
  group_by(mouseid) |>
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
df$b <- best$`01`
df$model <- best$model
df$lag <- best$lag
df$loglik_total <- best$loglik_total
df$loglik_sub <- best$loglik_sub

df |>
  ggplot()+
  geom_line(aes(x=lagRBC,y=pred,group=interaction(mouse,phase,pABA),col=pABA))+
  scale_x_continuous(labels = aakmisc::scinot)+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4500000))+
  annotate("text",x=7000000,y=3200000,label="Phase 2",size=5)+
  annotate("text",x=6000000,y=550000,label="Phase 1",size=5)+
  geom_segment(aes(x=6000000,xend=6000000,y=1100000,yend=600000))+
  geom_segment(aes(x=7000000,xend=5900000,y=3000000,yend=2500000))+
  xlab("")+ylab("Reticulocyte supply (t)")+
  scale_colour_manual(values=cbpABA[2:5])+
  geom_label(aes(x=9000000,y=4000000,label="Model E\nBreakpoint 10\nLag = 3 days"))+
  theme_bw()+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  labs(x=expression(paste("RBC density at time ", italic("t"), "-3 (density per µL)")),
       y=expression(paste("Reticulocyte supply at time ", italic("t"), " (density per µL)"),
                    colour="Parasite nutrient (pABA)"))+
  theme(
    axis.title=element_text(size=15),
    axis.title.y=element_text(size=15,colour="black"),
    axis.text=element_text(size=12),
    axis.text.x=element_text(colour="black"),
    panel.grid.minor=element_blank(),
    legend.position=c(0.2,0.2),
    legend.background=element_blank(),
    legend.title=element_text(size=13),
    legend.text=element_text(size=10),
    strip.text=element_text(size=12),
    strip.background=element_blank(),
    plot.title=element_text(size=14,hjust=0.5,face="bold")#,
    #panel.border = element_rect(linewidth=4)
    
  )
ggsave("FigureS9.jpeg",width=15,height=15,units="cm")
