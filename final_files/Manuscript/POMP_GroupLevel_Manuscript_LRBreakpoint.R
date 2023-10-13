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

#Load in PNAS trajectories
sm1name <- "m5sm1.rds"
sm1 <- readRDS(sm1name)

#Convert sm1 to tibble, remove underdosed mice, calculate relative likelihood, select for pABA 0% (box 4)
sm1 |>
  as_tibble() |>
  filter(time<=21,mouseid!="01-02",mouseid!="02-03") |> #remove underdosed mice
  pivot_wider(names_from=variable,values_from=value) |>
  separate_wider_delim(cols="mouseid",delim="-",names=c("box","mouse"),cols_remove=FALSE) |>
  group_by(mouseid) |>
  mutate(
    lik=exp(loglik-max(loglik))
  ) |>
  ungroup() |> #|> select(rep,mouse,box,lik) |> distinct() ->sm1_tmp
  filter(box=="04") |>
  select(rep,mouse,time,E,R,lik) |>
  unite("key",rep:mouse,sep="_",remove=FALSE) -> sm1_mod

#Create df of key and likelihoods to sample
sm1_mod |> filter(mouse=="01") |> select(key,lik) |> unique() -> sample_mouse01
sm1_mod |> filter(mouse=="02") |> select(key,lik) |> unique() -> sample_mouse02
sm1_mod |> filter(mouse=="03") |> select(key,lik) |> unique() -> sample_mouse03
nrow(sample_mouse01) #should be nrow = 2000

#Create list of keys based on likelihoods
sample_list01 <- sample(sample_mouse01$key,size=1000,replace=TRUE,prob=sample_mouse01$lik)
sample_list02 <- sample(sample_mouse02$key,size=1000,replace=TRUE,prob=sample_mouse02$lik)
sample_list03 <- sample(sample_mouse03$key,size=1000,replace=TRUE,prob=sample_mouse03$lik)
length(sample_list01) #should be length 1000

#Tabulate how often a key shows up in sample_list
sample_table01 <- sample_list01 |> table() |> as.data.frame() 
sample_table02 <- sample_list02 |> table() |> as.data.frame() 
sample_table03 <- sample_list03 |> table() |> as.data.frame() 
names(sample_table01) <- c("sample_list","Freq")
names(sample_table02) <- c("sample_list","Freq")
names(sample_table03) <- c("sample_list","Freq")

sum(sample_table01$Freq) #total of frequency column should be 1000

#Join together sample_tables
sample_table <- bind_rows(sample_table01,sample_table02,sample_table03)

#Create empty data frame to store output
pred<-data.frame()

#Create list of potential breakpoint days
breakpoint_list <- c(8,9,10,11,12,13,14)

for (i in 1:nrow(sample_table)){
  
  #Assign key from sample_list to s
  s<-sample_table$sample_list[i]
  
  #Filter sm1_sample for assigned key, calculate lagE and omit time = 0 (which is NA for lagE)
  df<- sm1_mod |> 
    filter(key==s) |>
    mutate(lagE=lag(E)) |>
    na.omit()
  
  #Create empty list to store log likelihoods
  #create vector of length of the breakpoint list then fill in
  ll_list<-c()
  
  for (b in breakpoint_list){ #loop over all potential breakpoints
    
    #Create dfs that split the data based on the breakpoint b    
    df1<- df |> filter(time<=b) #data frame for first regression
    df2<- df |> filter(time>b) #data frame for second regression
    
    #Run quadratic regression
    lm1<-lm(R~poly(lagE,2),data=df1)
    lm2<-lm(R~poly(lagE,2),data=df2)
    #Run linear regression for the second df
    lm3<-lm(R~lagE,data=df2)
    
    #Compare lm2 and lm3 using ANOVA
    anov<-anova(lm3,lm2)
    
    #Calculate log likelihood for each regression
    ll1<-stats::logLik(lm1)
    if(anov$`Pr(>F)`[2]<0.05){
      ll2<-stats::logLik(lm2)
    } else
    {
      ll2<-stats::logLik(lm3)
    }
    
    #Sum the log likelihoods of both regressions
    llsum<-ll1+ll2
    
    #Append summed log likelihoods to ll_list
    ll_list<-append(ll_list,llsum)
    
  } #end of loop over breakpoint list
  
  #Select breakpoint based on largest log likelihood
  breakpoint<-breakpoint_list[which(ll_list==max(ll_list))]
  
  #Create dfs that split the data based on the selected breakpoint
  df1<- df |> filter(time<=breakpoint)
  df2<- df |> filter(time>breakpoint)
  
  #Run quadratic regression
  lm1<-lm(R~poly(lagE,2),data=df1)
  lm2<-lm(R~poly(lagE,2),data=df2)
  #Run linear regression for the second df
  lm3<-lm(R~lagE,data=df2)
  
  #Compare lm2 and lm3 using ANOVA
  anov<-anova(lm3,lm2)
  
  #Calculate the predicted values for R given the two regressions
  pred1<-predict(lm1) |> as.data.frame()
  if(anov$`Pr(>F)`[2]<0.05){
    pred2<-predict(lm2) |> as.data.frame()
  } else
  {
    pred2<-predict(lm3) |> as.data.frame()
  }
  
  #Create lagE column in predicted R dfs
  pred1$lagE<-df1$lagE
  pred2$lagE<-df2$lagE
  
  names(pred1)<-c("predR","lagE")
  names(pred2)<-c("predR","lagE") 
  
  pred1$lm<-1
  pred2$lm<-2
  
  pred_tmp<-bind_rows(pred1,pred2)
  pred_tmp$mouse<-df$mouse |> unique()
  pred_tmp$key<-s
  pred_tmp$bp<-breakpoint
  pred_tmp$i<-i
  
  reps<-sample_table$Freq[which(sample_table$sample_list==s)]
  
  pred_tmp2<-data.frame()
  for (r in 1:reps){
    pred_tmp$rep <- r
    pred_tmp2<-bind_rows(pred_tmp2,pred_tmp) 
  }
  
  pred<-bind_rows(pred,pred_tmp2)
  #print(pred |> select(i,bp,rep) |> unique() |> nrow())

} #end of loop over sampled trajectories

pred<-unite(pred,"group_ID",i:rep,sep="_",remove=FALSE)

stopifnot(pred |> select(i) |> unique() |> nrow()==nrow(sample_table)) #should equal nrow(sample_table)
stopifnot(pred |> select(i,rep) |> unique() |> nrow()==3000) #should equal number of sampled keys (default 3000)

pred |> select(i,bp,rep) |> unique() |>
  ggplot()+
  geom_bar(aes(x=bp))

table(unique(select(pred,i,bp,rep))$bp)/3000

ggplot()+
  geom_line(data=filter(pred,lm==1),aes(x=lagE,y=predR,group=group_ID,col=factor(mouse)),alpha=0.01)+
  geom_line(data=filter(pred,lm==2),aes(x=lagE,y=predR,group=group_ID,col=factor(mouse)),alpha=0.01)+
  #geom_path(data=med_pABA0,aes(x=lagE,y=R))+
  #geom_text(data=med_pABA0,aes(x=lagE,y=R,label=time),size=5)+
  theme_bw()+
  xlab("Erythrocytes (t-1)")+ylab("Reticulcoytes (t)")+
  theme(
    strip.background=element_blank(),
    panel.grid = element_blank(),
    legend.position="none",
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    legend.text=element_text(size=14))

#Prepare dataframe for plotting group level R vs. lag E for pABA 0%
group_traj |>
  filter(pABA=="0%",variable%in%c("E","R")) |>
  select(-lo,-hi,-pABA,-box) |>
  pivot_wider(names_from="variable",values_from="med") |>
  mutate(lagE=lag(E)) |>
  filter(time<=20) -> med_pABA0

ggplot()+
  geom_path(data=med_pABA0,aes(x=lagE,y=R),col=cbPalette[2],linewidth=2)+
  geom_text(data=med_pABA0,aes(x=lagE,y=R,label=time),size=5)+
  theme_bw()+
  xlab("Erythrocytes (t-1)")+ylab("Reticulcoytes (t)")+
  theme(
    strip.background=element_blank(),
    panel.grid = element_blank(),
    legend.position="none",
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    legend.text=element_text(size=14))

#Prepare dataframe for plotting group level Qun vs. time for pABA 0%
group_traj |>
  filter(pABA=="0%",variable%in%c("Qun")) |>
  filter(time<=20) -> Qun_pABA0

ggplot()+
  geom_line(data=Qun_pABA0,aes(x=time,y=med),col=cbPalette[2],linewidth=2)+
  geom_ribbon(data=Qun_pABA0,aes(x=time,ymin=lo,ymax=hi),fill=cbPalette[2],alpha=0.2)+
  #geom_text(data=med_pABA0,aes(x=lagE,y=R,label=time),size=5)+
  theme_bw()+
  xlab("Day post-infection")+ylab("Probability uRBC removed\nby indiscriminate killing")+
  theme(
    strip.background=element_blank(),
    panel.grid = element_blank(),
    legend.position="none",
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    legend.text=element_text(size=14))

#Prepare dataframe for plotting group level R vs. time for pABA 0%
group_traj |>
  filter(pABA=="0%",variable%in%c("R")) |>
  filter(time<=20) -> R_pABA0

ggplot()+
  geom_line(data=R_pABA0,aes(x=time,y=med),col=cbPalette[2],linewidth=2)+
  geom_ribbon(data=R_pABA0,aes(x=time,ymin=lo,ymax=hi),fill=cbPalette[2],alpha=0.2)+
  #geom_text(data=med_pABA0,aes(x=lagE,y=R,label=time),size=5)+
  theme_bw()+
  xlab("Day post-infection")+ylab("Reticulocyte supply (density per microlitre)")+
  theme(
    strip.background=element_blank(),
    panel.grid = element_blank(),
    legend.position="none",
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    legend.text=element_text(size=14))
