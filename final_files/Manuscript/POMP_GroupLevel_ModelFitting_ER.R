#SET DIRECTORY TO SOURCE FILE LOCATION

library(tidyverse)
library(pomp)
library(foreach)
library(iterators)
library(doFuture)
plan(multisession)

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
  filter(box!="05") |> #remove control mice
  select(rep,mouse,mouseid,box,time,E,R,lik) |>
  unite("key",c(rep,box,mouse),sep="_",remove=FALSE) -> sm1_mod

#Create dataframe with sampled keys from each box and frequency of each key
mouse_id_list <- unique(sm1_mod$mouseid)
sample_size <- 1000
foreach (
  id=mouse_id_list,
  .combine=rbind,
  .options.future = list(seed = 851657743)
) %dofuture% {
  
  key_lik_df <- sm1_mod |> filter(mouseid==id) |> select(key,lik) |> unique()
  
  key_list <- sample(key_lik_df$key,size=sample_size,replace=TRUE,prob=key_lik_df$lik)
  
  key_table <- key_list |> table() |> as.data.frame(stringsAsFactors = FALSE) 
  names(key_table) <- c("key","Freq")
  
  key_table
  
} -> joint_key_table

#Checkpoint
stopifnot(sum(joint_key_table$Freq)==length(mouse_id_list)*sample_size)

### Ask Aaron about "baking" data frame below###

#Create data frame with R, lagE, time, mouse and box
bake(file="jrdf.rds",{
foreach (
  i=iter(joint_key_table,"row"),
  .combine=rbind
) %dofuture% {
  
  key_choice <- i$key
  freq <- i$Freq
  
  df <- sm1_mod |> 
    filter(key==key_choice) |>
    mutate(lagE=lag(E)) |>
    na.omit() |>
    select(R,lagE,time,mouse,box)
  
  if (freq>1) {
    replicate(n=freq,df,simplify=FALSE) |> bind_rows() -> df
  }
  df  
} |>
  mutate(box=paste0("box",box)) }) -> joint_repped_df

#Create data frame with breakpoints for each of the four pABA boxes
breakpoint_range <- 8:10
breakpoint_grid <- expand_grid(
  box01=breakpoint_range,
  box02=breakpoint_range,
  box03=breakpoint_range,
  box04=breakpoint_range
)

box_list <- unique(joint_repped_df$box)

#Create data frame that will store AICs and corresponding breakpoints for each model
bake(file="jAICdf.rds",{
foreach (bp=iter(breakpoint_grid,"row"),
         .combine=rbind
) %dofuture% {
  
  joint_repped_df |>
    mutate(
      phase=as.factor(if_else(time<=unlist(bp)[box],1,2))
    ) -> df
  
  #Models
  models <- list(
    m1=lm(R~poly(lagE,2):(box-1):(phase-1)+(box-1):(phase-1),data=df),
    m2=lm(R~(poly(lagE,2)+box-1):(phase-1)+(phase-1),data=df),
    m3=lm(R~poly(lagE,2):(phase-1)+(phase-1),data=df),
    m4=lm(R~lagE:(box-1):(phase-1)+(box-1):(phase-1),data=df),
    m5=lm(R~(lagE+box-1):(phase-1)+(phase-1),data=df),
    m6=lm(R~lagE:(phase-1)+(phase-1),data=df)
  )
  
  models |>
    lapply(\(m) tibble(AIC=AIC(m),BIC=BIC(m),bp)) |>
    bind_rows(.id="model")
  
} }) -> joint_AIC_df

bp <- joint_AIC_df[which(joint_AIC_df$AIC==min(joint_AIC_df$AIC)),4:7]

joint_repped_df |>
  mutate(
    phase=as.factor(if_else(time<=unlist(bp)[box],1,2))
  ) -> df

models <- list(
  m1=lm(R~poly(lagE,2):(box-1):(phase-1)+(box-1):(phase-1),data=df),
  m2=lm(R~(poly(lagE,2)+box-1):(phase-1)+(phase-1),data=df),
  m3=lm(R~poly(lagE,2):(phase-1)+(phase-1),data=df),
  m4=lm(R~lagE:(box-1):(phase-1)+(box-1):(phase-1),data=df),
  m5=lm(R~(lagE+box-1):(phase-1)+(phase-1),data=df),
  m6=lm(R~lagE:(phase-1)+(phase-1),data=df)
)

chosen_model <- models[which(paste0("m",1:6)==joint_AIC_df[which(joint_AIC_df$AIC==min(joint_AIC_df$AIC)),1])]
chosen_model

