library(tidyverse)

setwd("~/Documents/GitHub/bdd/nw11_hier/final_files/Manuscript/") #sorry Aaron

#### Make data frame with all of the data ####

#Load in PNAS trajectories
sm1name <- "m5sm1.rds"
sm1 <- readRDS(sm1name)

#Create sm1_mod tibble, with columns rep, mouse, mousid, box, time, E, R, lik and key
sm1 |>
  as_tibble() |>
  filter(time<=21,mouseid!="01-02",mouseid!="02-03") |> #remove underdosed mice
  pivot_wider(names_from=variable,values_from=value) |>
  separate_wider_delim(cols="mouseid",delim="-",names=c("box","mouse"),cols_remove=FALSE) |>
  group_by(mouseid) |>
  mutate(
    lik=exp(loglik-max(loglik))
  ) |>
  ungroup() |>
  filter(box!="01") |> #remove control mice
  select(rep,mouse,mouseid,box,time,E,R,lik) |>
  unite("key",c(rep,box,mouse),sep="_",remove=FALSE) -> sm1_mod

#Create dataframe with sampled keys from each box and frequency of each key
joint_key_table<-data.frame()
mouse_id_list <- unique(sm1_mod$mouseid)
sample_size <- 1000
for(id in mouse_id_list){
  
  key_lik_df <- sm1_mod |> filter(mouseid==id) |> select(key,lik) |> unique()
  
  key_list <- sample(key_lik_df$key,size=sample_size,replace=TRUE,prob=key_lik_df$lik)
  
  key_table <- key_list |> table() |> as.data.frame() 
  names(key_table) <- c("key","Freq")
  
  joint_key_table <- bind_rows(joint_key_table,key_table)
  
}
#Checkpoint
stopifnot(sum(joint_key_table$Freq)==length(mouse_id_list)*sample_size)

#Create empty data frame to store output
pred<-data.frame()

#Create list of potential breakpoint days
breakpoint_list <- c(8,9,10,11,12,13,14)

for (i in 1:nrow(joint_key_table)){
  
  #Assign key from sample_list to s
  s<-joint_key_table$key[i] |> as.character()
  
  #Filter sm1_sample for assigned key, calculate lagE and omit time = 0 (which is NA for lagE)
  df <- sm1_mod |> 
    filter(key==s) |>
    mutate(lagE=lag(E)) |>
    na.omit() 
  
  #Create empty list to store log likelihoods
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
  } else {
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
  
  reps<-joint_key_table$Freq[which(joint_key_table$key==s)]
  
  pred_tmp2<-data.frame()
  for (r in 1:reps){
    pred_tmp$rep <- r
    pred_tmp2<-bind_rows(pred_tmp2,pred_tmp) 
  }
  
  pred<-bind_rows(pred,pred_tmp2)
  
} #end of for statement over rows in joint_key_table
