library(tidyverse)
library(stringi)
library(doParallel)
library(pomp)
library(panelPomp)

setwd("~/Documents/GitHub/bdd/nw11_hier/final_files/Manuscript/")

#Load in PNAS trajectories
sm1name <- "m5sm1.rds"
sm1 <- readRDS(sm1name)

sm1 |>
  as_tibble() |>
  filter(time<=21,mouseid!="01-02",mouseid!="02-03") |> #remove underdosed mice
  pivot_wider(names_from=variable,values_from=value) |>
  separate_wider_delim(cols="mouseid",delim="-",names=c("box","mouse"),cols_remove=FALSE) |>
  group_by(mouseid) |>
  mutate(
    lagE=lag(E),
    lik=exp(loglik-max(loglik))
  ) |>
  ungroup() |>
  filter(box=="04") |>
  select(rep,mouse,time,E,R,lik) |>
  unite("key",rep:mouse,sep="_",remove=FALSE) -> sm1_mod

sm1_mod |> select(key,lik) |> unique() -> sample_df

sample_list <- sample(sample_df$key,size=2000,replace=TRUE,prob=sample_df$lik)

sm1_mod |> filter(key%in%sample_list) -> sm1_sample

pred<-data.frame()

breakpoint_list <- c(8,9,10,11,12,13,14)

for (i in 1:2000){
  
  s<-sample_list[i]
  
  df<- sm1_sample |> 
    filter(key==s) |>
    mutate(lagE=lag(E)) |>
    na.omit()
  
  ll_list<-c()
  
  for (b in breakpoint_list){
    
    df1<- df |> filter(time<=b)
    df2<- df |> filter(time>b)
    
    df1$lagE2<-df1$lagE^2
    df2$lagE2<-df2$lagE^2
    
    lm1<-lm(R~lagE+lagE2,data=df1)
    lm2<-lm(R~lagE+lagE2,data=df2)
    
    ll1<-stats::logLik(lm1)
    ll2<-stats::logLik(lm2)
    
    llsum<-ll1+ll2
    
    ll_list<-append(ll_list,llsum)
    
  } #end of loop over breakpoint list
  
  breakpoint<-breakpoint_list[which(ll_list==max(ll_list))]
  
  df1<- df |> filter(time<=breakpoint)
  df2<- df |> filter(time>breakpoint)
  
  df1$lagE2<-df1$lagE^2
  df2$lagE2<-df2$lagE^2
  
  lm1<-lm(R~lagE+lagE2,data=df1)
  lm2<-lm(R~lagE+lagE2,data=df2)
  
  pred1<-predict(lm1) |> as.data.frame()
  pred2<-predict(lm2) |> as.data.frame()
  
  pred1$lagE<-df1$lagE
  pred2$lagE<-df2$lagE
  
  names(pred1)<-c("predR","lagE")
  names(pred2)<-c("predR","lagE") 
  
  pred1$key<-s
  pred1$lm<-1
  pred1$bp<-breakpoint
  pred1$i<-i
  
  pred2$key<-s
  pred2$lm<-2
  pred2$bp<-breakpoint
  pred2$i<-i
  
  pred_tmp<-bind_rows(pred1,pred2)
  
  pred<-bind_rows(pred,pred_tmp)

} #end of loop over sampled trajectories


pred |> select(i,bp) |> unique() |>
  ggplot()+
  geom_histogram(aes(x=bp))

ggplot()+
  geom_line(data=filter(pred,lm==1),aes(x=lagE,y=predR,group=i),alpha=0.01)+
  geom_line(data=filter(pred,lm==2),aes(x=lagE,y=predR,group=i),alpha=0.01)+
  theme_bw()+
  xlab("Erythrocytes (t-1)")+ylab("Reticulcoytes (t)")
