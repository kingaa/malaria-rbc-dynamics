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

### Ask Aaron about "baking" data frame below###

#Create data frame with R, lagE, time, mouse and box
joint_repped_df <- data.frame()
for (i in 1:nrow(joint_key_table)){
  
  key_choice <- joint_key_table$key[i] |> droplevels()
  freq <- joint_key_table$Freq[i]

  df <- sm1_mod |> 
    filter(key==key_choice) |>
    mutate(lagE=lag(E)) |>
    na.omit() |>
    select(R,lagE,time,mouse,box)
  
  if (freq>1){
    
    repped_df <- data.frame()
    for (rep in 1:freq){
      
      repped_df <- bind_rows(repped_df,df)
      
    } #end of for statement over frequency
    
  } else {
    
    repped_df <- df
    
  } #end of if else statement
  
  joint_repped_df <- bind_rows(joint_repped_df,repped_df)
    
} #end of for statement over rows in joint_key_table

#Create column "phase" to be filled in below
joint_repped_df$phase <- NA

#Create data frame with breakpoints for each of the four pABA boxes
breakpoint_range <- 8:9
breakpoint_grid <- expand_grid(breakpoint_range,
                               breakpoint_range,
                               breakpoint_range,
                               breakpoint_range)
names(breakpoint_grid) <- c("box02","box03","box04","box05")

#Create data frame that will store AICs and corresponding breakpoints for each model
joint_AIC_df <- data.frame()
for (i in 1:nrow(breakpoint_grid)){
  
  breakpoint_choices <- breakpoint_grid[i,]
  
  box_list <- unique(joint_repped_df$box)
  for(box_choice in box_list){
    
    joint_repped_subset <- filter(joint_repped_df,box==box_choice)
    breakpoint_choice_index <- grep(paste0("box",box_choice),
                                    names(breakpoint_choices))
    
    breakpoint_choice <- as.integer(breakpoint_choices[breakpoint_choice_index])
    
    joint_repped_df$phase[joint_repped_df$box==box_choice&joint_repped_df$time<=breakpoint_choice] <- 1
    joint_repped_df$phase[joint_repped_df$box==box_choice&joint_repped_df$time>breakpoint_choice] <- 2
    
  } #end for statement over box_list
  
  #Models
  m1 <- lm(R~(poly(lagE,2)*box):phase,data=joint_repped_df)
  m2 <- lm(R~(poly(lagE,2)+box):phase,data=joint_repped_df)
  m3 <- lm(R~(poly(lagE,2)):phase,data=joint_repped_df)
  m4 <- lm(R~(lagE*box):phase,data=joint_repped_df)
  m5 <- lm(R~(lagE+box):phase,data=joint_repped_df)
  m6 <- lm(R~lagE:phase,data=joint_repped_df)
  
  #ASK AARON!!!!
  AIC1<-extractAIC(m1)[2]
  AIC2<-extractAIC(m2)[2]
  AIC3<-extractAIC(m3)[2]
  AIC4<-extractAIC(m4)[2]
  AIC5<-extractAIC(m5)[2]
  AIC6<-extractAIC(m6)[2]
  
  AIC_df <- as.data.frame(matrix(nrow=6,ncol=6))
  names(AIC_df) <- c("model","AIC",paste0("bp0",2:5))
  
  AIC_df$model <- 1:6
  AIC_df$AIC <- c(AIC1,AIC2,AIC3,AIC4,AIC5,AIC6)
  AIC_df$bp02 <- rep(breakpoint_choices[1],6)
  AIC_df$bp03 <- rep(breakpoint_choices[2],6)
  AIC_df$bp04 <- rep(breakpoint_choices[3],6)
  AIC_df$bp05 <- rep(breakpoint_choices[4],6)
  
  joint_AIC_df <- bind_rows(joint_AIC_df,AIC_df)
  
  
} #end of for statement over rows in breakpoint_grid




