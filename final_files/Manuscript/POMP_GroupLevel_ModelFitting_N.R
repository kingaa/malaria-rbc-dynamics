#SET DIRECTORY TO SOURCE FILE LOCATION

library(tidyverse)
library(mgcv)
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
  select(rep,mouse,mouseid,box,time,N,lik) |>
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

#Create data frame with N
foreach (
  i=iter(joint_key_table,"row"),
  .combine=rbind
) %dofuture% {
  
  key_choice <- i$key
  freq <- i$Freq
  
  df <- sm1_mod |> 
    filter(key==key_choice) |>
    select(N,time,mouse,box)
  
  if (freq>1) {
    replicate(n=freq,df,simplify=FALSE) |> bind_rows() -> df
  }
  df  
} |>
  mutate(box=paste0("box",box)) -> joint_repped_df

#Models
models <- list(
  m1=gam(N~s(time*box),data=joint_repped_df), #ASK AARON
  m2=gam(N~s(time),data=joint_repped_df),
  m3=gam(N~1,data=joint_repped_df)
)

models |>
  lapply(\(m) tibble(AIC=AIC(m),BIC=BIC(m))) |>
  bind_rows(.id="model")
