library(tidyverse)
library(stringi)
library(doParallel)
library(pomp)
library(panelPomp)
library(foreach)
library(iterators)
library(doRNG)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

setwd("~/Documents/GitHub/bdd/nw11_hier/final_files/Manuscript/")

source("POMP_GroupLevel_DataPrep.R")

group_traj |> 
  filter(variable%in%c("E","R"),time<=21) |>
  select(-lo,-hi) |>
  pivot_wider(names_from="variable",values_from="med") |>
  mutate(lagE=lag(E)) |>
  filter(pABA!="Uninfected") -> group_ER
group_ER$lagE[group_ER$time==0] <- NA
group_ER |> na.omit() -> group_ER     
breakpoint <- 9

df1 <- group_ER |> filter(time<=breakpoint)
df2 <- group_ER |> filter(time>breakpoint)

lm1 <- lm(R~poly(lagE,2),data=df1)
lm2 <- lm(R~poly(lagE,2)*pABA,data=df1)
anov1 <- anova(lm1,lm2)
if(anov1$`Pr(>F)`[2]<0.05){
  pred1<-predict(lm2) |> as.data.frame()
} else {
  pred1<-predict(lm1) |> as.data.frame()
}
pred1$lagE <- df1$lagE
pred1$pABA <- df1$pABA
names(pred1) <- c("R","lagE","pABA")

lm3 <- lm(R~poly(lagE,2),data=df2)
lm4 <- lm(R~poly(lagE,2)*pABA,data=df2)
anov2 <- anova(lm3,lm4)
if(anov2$`Pr(>F)`[2]<0.05){
  pred2<-predict(lm4) |> as.data.frame()
} else {
  pred2<-predict(lm3) |> as.data.frame()
}
pred2$lagE <- df2$lagE
pred2$pABA <- df2$pABA
names(pred2) <- c("R","lagE","pABA")

pred1$lm<-1
pred2$lm<-2

ggplot()+
  geom_line(data=pred1,aes(x=lagE,y=R,col=pABA),linewidth=2)+
  geom_line(data=pred2,aes(x=lagE,y=R),linewidth=2)+
  scale_colour_manual(values=cbPalette[2:5])+
  theme_bw()+
  xlab("Erythrocytes (t-1)")+ylab("Reticulcoytes (t)")+
  theme(
    strip.background=element_blank(),
    panel.grid = element_blank(),
    legend.position = c(0.8,0.8),
    legend.title = element_text(size=15),
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    legend.text=element_text(size=14))

#############################################################################
group_traj |> 
  filter(variable%in%c("Qun"),time<=21,pABA!="Uninfected") -> group_Qun

group_Qun |>
  ggplot()+
  geom_line(aes(x=time,y=med,col=pABA),linewidth=2)+
  geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
  scale_colour_manual(values=cbPalette[2:5])+
  scale_fill_manual(values=cbPalette[2:5])+
  theme_bw()+
  xlab("Day post-infection")+ylab("Probability uRBC removed\nby indiscriminate killing")+
  theme(
    strip.background=element_blank(),
    panel.grid = element_blank(),
    legend.position=c(0.2,0.8),
    legend.title=element_text(size=15),
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    legend.text=element_text(size=14))
