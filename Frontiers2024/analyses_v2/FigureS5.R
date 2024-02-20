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

#### 1-day lag ####
flow |> 
  select(time=day,box,pABA,E=Eryth,R=Retic,pABA) |>
  mutate(RBC=E+R,
         lagRBC=lag(RBC,1)) |>
  na.omit() |>
  filter(time<=20,box!="05",) -> data

##Create data frame with breakpoints for each of the four pABA boxes
box_list <- unique(data$box)
breakpoint_grid <- as_tibble(
  matrix(
    data=rep(c(NA,8,9,10,11),each=4),
    ncol=4,byrow=TRUE,
    dimnames=list(NULL,box_list)
  )
)

joint_AIC_df <- data.frame()
for (i in seq(1,nrow(breakpoint_grid),1)){
  print(i)
  bp <- breakpoint_grid[i,]
  
  if (is.na(bp[1])){
    
    data -> df
    
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
    
    data |>
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
data |>
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
df$b <- if (is.na(best$`01`)){"None"} else {best$`01`}
df$model <- best$model

#Plotting
lag1 <- df |>
  ggplot()+
  geom_line(aes(x=lagRBC,y=pred,group=interaction(phase,pABA),col=pABA))+
  scale_x_continuous(labels = aakmisc::scinot)+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4000000))+
  annotate("text",x=7000000,y=3200000,label="Phase 2",size=5)+
  annotate("text",x=6000000,y=1100000,label="Phase 1",size=5)+
  geom_segment(aes(x=6000000,xend=6000000,y=1000000,yend=600000))+
  geom_segment(aes(x=7000000,xend=5900000,y=3000000,yend=2500000))+
  xlab("")+ylab("Reticulocyte supply (t)")+
  scale_colour_manual(values=cbPalette[2:5])+
  geom_text(aes(x=7500000,y=3900000,label="Model B, breakpoint 9"))+
  theme_bw()+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  labs(x=expression(paste("RBC density at time ", italic("t"), "-1 (density per µL)")),
       y=expression(paste("Reticulocyte supply at time ", italic("t"), " (density per µL)"),
                    colour="Parasite nutrient (pABA)"))+
  theme(
    axis.title=element_text(size=15),
    axis.title.y=element_text(size=15,colour="black"),
    axis.text=element_text(size=12),
    axis.text.x=element_text(colour="black"),
    panel.grid.minor=element_blank(),
    legend.position=c(0.2,0.8),
    legend.background=element_blank(),
    legend.title=element_text(size=13),
    legend.text=element_text(size=10),
    strip.text=element_text(size=12),
    strip.background=element_blank(),
    plot.title=element_text(size=14,hjust=0.5,face="bold")#,
    #panel.border = element_rect(linewidth=4)
    
  )
lag1

#### 2-day lag ####
flow |> 
  select(time=day,box,pABA,E=Eryth,R=Retic,pABA) |>
  mutate(RBC=E+R,
         lagRBC=lag(RBC,2)) |>
  na.omit() |>
  filter(time<=20,box!="05") -> data

##Create data frame with breakpoints for each of the four pABA boxes
box_list <- unique(data$box)
breakpoint_grid <- as_tibble(
  matrix(
    data=rep(c(NA,8,9,10,11),each=4),
    ncol=4,byrow=TRUE,
    dimnames=list(NULL,box_list)
  )
)

joint_AIC_df <- data.frame()
for (i in seq(1,nrow(breakpoint_grid),1)){
  print(i)
  bp <- breakpoint_grid[i,]
  
  if (is.na(bp[1])){
    
    data -> df
    
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
    
    data |>
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
data |>
  mutate(
    phase=as.factor(if_else(time<=as.integer(unlist(best)[box]),1,2))
  ) -> df2


##Extract best model based on lowest AIC
if (is.na(best$`01`)){
  models <- list(
    m1=lm(R~poly(lagRBC,2,raw=F):(box-1)+(box-1),data=df2),
    m2=lm(R~(poly(lagRBC,2,raw=F)+box-1),data=df2),
    m3=lm(R~poly(lagRBC,2,raw=F),data=df2),
    m4=lm(R~lagRBC:(box-1)+(box-1),data=df2),
    m5=lm(R~(lagRBC+box-1),data=df2),
    m6=lm(R~lagRBC,data=df2)
  )
} else {
  models <- list(
    m1=lm(R~poly(lagRBC,2,raw=F):(box-1):(phase-1)+(box-1):(phase-1),data=df2),
    m2=lm(R~(poly(lagRBC,2,raw=F)+box-1):(phase-1)+(phase-1),data=df2),
    m3=lm(R~poly(lagRBC,2,raw=F):(phase-1)+(phase-1),data=df2),
    m4=lm(R~lagRBC:(box-1):(phase-1)+(box-1):(phase-1),data=df2),
    m5=lm(R~(lagRBC+box-1):(phase-1)+(phase-1),data=df2),
    m6=lm(R~lagRBC:(phase-1)+(phase-1),data=df2)
  )
}

chosen_model <- models[[best$model]]
coefs <- chosen_model$coefficients

pred <- model.matrix(chosen_model) %*% coefs |> as.data.frame() |> select(pred=V1)
df2 <- cbind(df2,pred)
df2$b <- if (is.na(best$`01`)){"None"} else {best$`01`}
df2$model <- best$model

#Plotting
lag2 <- df2 |>
  ggplot()+
  geom_line(aes(x=lagRBC,y=pred,group=interaction(phase,pABA),col=pABA))+
  scale_x_continuous(labels = aakmisc::scinot)+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4000000))+
  annotate("text",x=7000000,y=3200000,label="Phase 2",size=5)+
  annotate("text",x=6000000,y=1100000,label="Phase 1",size=5)+
  geom_segment(aes(x=6000000,xend=6000000,y=1000000,yend=800000))+
  geom_segment(aes(x=7000000,xend=5900000,y=3000000,yend=2400000))+
  scale_colour_manual(values=cbPalette[2:5])+
  geom_text(aes(x=7500000,y=3900000,label="Model B, breakpoint 9"))+
  theme_bw()+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  labs(x=expression(paste("RBC density at time ", italic("t"), "-2 (density per µL)")),
       y=expression(paste("Reticulocyte supply at time ", italic("t"), " (density per µL)"),
                    colour="Parasite nutrient (pABA)"))+
  theme(
    axis.title=element_text(size=15),
    axis.title.y=element_text(size=15,colour="white"),
    axis.text=element_text(size=12),
    axis.text.x=element_text(colour="black"),
    panel.grid.minor=element_blank(),
    legend.position="none",
    legend.background=element_blank(),
    legend.title=element_text(size=13),
    legend.text=element_text(size=10),
    strip.text=element_text(size=12),
    strip.background=element_blank(),
    plot.title=element_text(size=14,hjust=0.5,face="bold")#,
    #panel.border = element_rect(linewidth=4)
    
  )
lag2

#### 3-day lag ####
flow |> 
  select(time=day,box,pABA,E=Eryth,R=Retic,pABA) |>
  mutate(RBC=E+R,
         lagRBC=lag(RBC,3)) |>
  na.omit() |>
  filter(time<=20,box!="05") -> data

##Create data frame with breakpoints for each of the four pABA boxes
box_list <- unique(data$box)
breakpoint_grid <- as_tibble(
  matrix(
    data=rep(c(NA,8,9,10,11),each=4),
    ncol=4,byrow=TRUE,
    dimnames=list(NULL,box_list)
  )
)

joint_AIC_df <- data.frame()
for (i in seq(1,nrow(breakpoint_grid),1)){
  print(i)
  bp <- breakpoint_grid[i,]
  
  if (is.na(bp[1])){
    
    data -> df
    
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
    
    data |>
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
data |>
  mutate(
    phase=as.factor(if_else(time<=as.integer(unlist(best)[box]),1,2))
  ) -> df3


##Extract best model based on lowest AIC
if (is.na(best$`01`)){
  models <- list(
    m1=lm(R~poly(lagRBC,2,raw=F):(box-1)+(box-1),data=df3),
    m2=lm(R~(poly(lagRBC,2,raw=F)+box-1),data=df3),
    m3=lm(R~poly(lagRBC,2,raw=F),data=df3),
    m4=lm(R~lagRBC:(box-1)+(box-1),data=df3),
    m5=lm(R~(lagRBC+box-1),data=df3),
    m6=lm(R~lagRBC,data=df3)
  )
} else {
  models <- list(
    m1=lm(R~poly(lagRBC,2,raw=F):(box-1):(phase-1)+(box-1):(phase-1),data=df3),
    m2=lm(R~(poly(lagRBC,2,raw=F)+box-1):(phase-1)+(phase-1),data=df3),
    m3=lm(R~poly(lagRBC,2,raw=F):(phase-1)+(phase-1),data=df3),
    m4=lm(R~lagRBC:(box-1):(phase-1)+(box-1):(phase-1),data=df3),
    m5=lm(R~(lagRBC+box-1):(phase-1)+(phase-1),data=df3),
    m6=lm(R~lagRBC:(phase-1)+(phase-1),data=df3)
  )
}

chosen_model <- models[[best$model]]
coefs <- chosen_model$coefficients

pred <- model.matrix(chosen_model) %*% coefs |> as.data.frame() |> select(pred=V1)
df3 <- cbind(df3,pred)
df3$b <- if (is.na(best$`01`)){"None"} else {best$`01`}
df3$model <- best$model

#Plotting
lag3 <- df3 |>
  ggplot()+
  geom_line(aes(x=lagRBC,y=pred,group=interaction(phase,pABA),col=pABA))+
  scale_x_continuous(labels = aakmisc::scinot)+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4000000))+
  annotate("text",x=7000000,y=3200000,label="Phase 2",size=5)+
  annotate("text",x=6000000,y=1100000,label="Phase 1",size=5)+
  geom_segment(aes(x=6000000,xend=6000000,y=1000000,yend=650000))+
  geom_segment(aes(x=7000000,xend=5700000,y=3000000,yend=2200000))+
  scale_colour_manual(values=cbPalette[2:5])+
  geom_text(aes(x=7500000,y=3900000,label="Model E, breakpoint 9"))+
  theme_bw()+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  labs(x=expression(paste("RBC density at time ", italic("t"), "-3 (density per µL)")),
       y=expression(paste("Reticulocyte supply at time ", italic("t"), " (density per µL)"),
                    colour="Parasite nutrient (pABA)"))+
  theme(
    axis.title=element_text(size=15),
    axis.title.y=element_text(size=15,colour="white"),
    axis.text=element_text(size=12),
    axis.text.x=element_text(colour="black"),
    panel.grid.minor=element_blank(),
    legend.position="none",
    legend.background=element_blank(),
    legend.title=element_text(size=10),
    legend.text=element_text(size=8),
    strip.text=element_text(size=12),
    strip.background=element_blank(),
    plot.title=element_text(size=14,hjust=0.5,face="bold")#,
    #panel.border = element_rect(linewidth=4)
    
  )
lag3

#Plotting together
gt <- arrangeGrob(lag1,lag2,lag3,nrow=1)

# Add labels to the arranged plots
p <- as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0,0.33,0.66), y = c(1,1,1))
p
ggsave("FigureS6.jpeg",width=40,height=14,units="cm")
