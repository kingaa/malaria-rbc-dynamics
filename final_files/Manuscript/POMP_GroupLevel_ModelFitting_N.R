#SET DIRECTORY TO SOURCE FILE LOCATION

library(tidyverse)
library(mgcv)
library(foreach)
library(iterators)
library(doFuture)
plan(multisession) #for faster results, use multicore outside of RStudio

seed_choice <- 851657743
set.seed(seed_choice)

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

joint_repped_df$box <- factor(joint_repped_df$box) #important that this is factor so gam works with "by = box"

#Models
models <- list(
  m1=gam(N~(box-1) + s(time, by = box),data=joint_repped_df),
  m2=gam(N~s(time),data=joint_repped_df),
  m3=gam(N~1,data=joint_repped_df)
)

models |>
  lapply(\(m) tibble(AIC=AIC(m),BIC=BIC(m))) |>
  bind_rows(.id="model") -> AIC_df

chosen_model <- models[which(paste0("m",1:3)==AIC_df[which(AIC_df$AIC==min(AIC_df$AIC)),1])]
m01 <- chosen_model[[1]]

beta <- coef(m01)

Vb <- vcov(m01, unconditional = TRUE)
se <- sqrt(diag(Vb))

CI <- sapply(names(beta),FUN=function(name){
  
  i <- which(names(beta) == name)
  x <- beta[i] + (c(-1,1) * (2 * se[i]))
  x
  
})


#foreach (r=1:1000,
 #        .combine=rbind,
  #       .options.future = list(seed = seed_choice)
#) %dofuture% { #for loop over reps
bake(file="jndN.rds",{
  
  joint <- data.frame()
  for (r in 1:1000){
    m02 <- m01
    
    for (co in names(coef(m02))){
      
      lwr <- CI[1,co]
      upr <- CI[2,co]
      
      samp <- runif(1,lwr,upr)
      
      m02$coefficients[co] <- samp
      
    }
    
    newdata <- expand.grid(
      time=0:20,
      box=c("box01","box02","box03","box04"))
    
    newdata$pred <- predict(m02,newdata)
    newdata$rep <- r
    
    #newdata
    joint<- rbind(joint,newdata)

  }
  
  joint
  
}) -> joint_newdata

joint_newdata$pABA <- factor(joint_newdata$box,
                             levels=c("box04","box03","box02","box01"),
                             labels=c("0%","0.0005%","0.005%","0.05%"))

joint_newdata |>
  group_by(pABA,time) |>
  summarise(
    min=min(pred),
    avg=mean(pred),
    max=max(pred)
  ) -> joint_group

(ModelFitting_N <- ggplot()+
  geom_line(data=joint_group,aes(x=time,y=avg,col=pABA),linewidth=2)+
  geom_ribbon(data=joint_group,aes(x=time,ymin=min,ymax=max,fill=pABA),alpha=0.2)+
  xlab("Day post-infection")+ylab("Indiscriminate killing (density per microlitre)")+
  scale_colour_manual(values=cbPalette[2:5])+
  scale_fill_manual(values=cbPalette[2:5])+
  theme_bw()+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    panel.grid=element_blank(),
    legend.position=c(0.2,0.7),
    legend.title=element_text(size=15),
    legend.text=element_text(size=13)
  )
)
ggsave("ModelFitting_N.png",plot=ModelFitting_N,
       width=20,height=15,units="cm") 
