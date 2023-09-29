library(tidyverse)
library(stringi)
library(doParallel)
library(pomp)
library(panelPomp)
library(foreach)
library(iterators)
library(doRNG)
setwd("~/Documents/GitHub/bdd/nw11_hier/final_files/Manuscript/")

source("POMP_GroupLevel_DataPrep.R")

#Prepare dataframe for plotting group level R vs. lag E for pABA 0%
group_traj |>
  filter(pABA=="0%",variable%in%c("E","R")) |>
  select(-lo,-hi,-pABA,-box) |>
  pivot_wider(names_from="variable",values_from="med") |>
  mutate(lagE=lag(E)) |>
  filter(time<=20) -> med_pABA0

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

sm1_mod |> select(rep,mouse,lik) |> distinct() -> tmp

tmp |>
  ggplot()+
  geom_histogram(aes(x=log(lik)))+
  facet_grid(.~mouse)

#Create df of key and likelihoods to sample
sm1_mod |> select(key,lik) |> unique() -> sample_df
nrow(sample_df) #should be nrow = 6000

#Create list of keys based on likelihoods
sample_list <- sample(sample_df$key,size=2000,replace=TRUE,prob=sample_df$lik)
length(sample_list) #should be length 2000

#Tabulate how often a key shows up in sample_list
sample_table <- sample_list |> table() |> as.data.frame() 
sum(sample_table$Freq) #total of frequency column should be 2000

#Create empty data frame to store output
pred<-data.frame()

#Create list of potential breakpoint days
breakpoint_list <- c(8,9,10,11,12,13,14)

#Rep counter
rep_count<-0

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
    
    #Calculate log likelihood for each regression
    ll1<-stats::logLik(lm1)
    ll2<-stats::logLik(lm2)
    
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
  
  #Calculate the predicted values for R given the two regressions
  pred1<-predict(lm1) |> as.data.frame()
  pred2<-predict(lm2) |> as.data.frame()
  
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

pred |> select(i) |> unique() |> nrow() #should equal nrow(sample_table)
pred |> select(i,rep) |> unique() |> nrow() #should equal number of sampled keys (default 2000)

pred |> select(i,bp,rep) |> unique() |>
  ggplot()+
  geom_bar(aes(x=bp))

table(unique(select(pred,i,bp,rep))$bp)/2000

ggplot()+
  geom_line(data=filter(pred,lm==1),aes(x=lagE,y=predR,group=group_ID,col=factor(mouse)),alpha=0.05)+
  geom_line(data=filter(pred,lm==2),aes(x=lagE,y=predR,group=group_ID,col=factor(mouse)),alpha=0.05)+
  geom_path(data=med_pABA0,aes(x=lagE,y=R))+
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

ggplot()+
  geom_path(data=pred,aes(x=lagE,y=predR,group=group_ID,col=factor(mouse)),alpha=0.005)+
  #geom_path(data=med_pABA0,aes(x=lagE,y=R))+
  #geom_text(data=med_pABA0,aes(x=lagE,y=R,label=time),size=5)+
  theme_bw()+
  xlab("Erythrocytes (t-1)")+ylab("Reticulcoytes (t)")+
  facet_wrap(mouse~.)+
  theme(
    strip.background=element_blank(),
    panel.grid = element_blank(),
    legend.position="none",
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    legend.text=element_text(size=14))

ggplot()+
  geom_line(data=filter(pred,lm==1),aes(x=lagE,y=predR,group=group_ID,col=factor(mouse)),alpha=0.05)+
  geom_line(data=filter(pred,lm==2),aes(x=lagE,y=predR,group=group_ID,col=factor(mouse)),alpha=0.05)+
  theme_bw()+
  facet_wrap(mouse~.)+
  xlab("Erythrocytes (t-1)")+ylab("Reticulcoytes (t)")
