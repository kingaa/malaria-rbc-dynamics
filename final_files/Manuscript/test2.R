library(tidyverse)
library(pomp)
library(foreach)
library(iterators)
library(doFuture)
plan(multisession) #for faster results, use multicore outside of RStudio

seed_choice <- 851657743
set.seed(seed_choice)

#### Make data frame with all of the data ####

##Load in PNAS trajectories
sm1name <- "m5sm1.rds"
sm1 <- readRDS(sm1name)

##Create sm1_mod tibble, with columns rep, mouse, mousid, box, time, E, R, lik and key
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

##Create dataframe with sampled keys from each box and frequency of each key
mouse_id_list <- unique(sm1_mod$mouseid)
sample_size <- 1000
foreach (
  id=mouse_id_list,
  .combine=rbind,
  .options.future = list(seed = seed_choice)
) %dofuture% {
  key_lik_df <- sm1_mod |> filter(mouseid==id) |> select(key,lik) |> unique()
  key_list <- sample(key_lik_df$key,size=sample_size,replace=TRUE,prob=key_lik_df$lik)
  key_table <- key_list |> table() |> as.data.frame(stringsAsFactors = FALSE)
  names(key_table) <- c("key","Freq")
  key_table
} -> joint_key_table

##Checkpoint
stopifnot(sum(joint_key_table$Freq)==length(mouse_id_list)*sample_size)

##Create data frame with R, lagE, time, mouse and box
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
    mutate(box=paste0("box",box))
}
) -> joint_repped_df

##Create data frame with breakpoints for each of the four pABA boxes
box_list <- unique(joint_repped_df$box)
breakpoint_grid <- as_tibble(
  matrix(
    data=rep(c(NA,8,9,10,11),each=4),
    ncol=4,byrow=TRUE,
    dimnames=list(NULL,box_list)
  )
)

##Model fitting
bake(file="jAICdf_null.rds",{
  foreach (
    bp=iter(breakpoint_grid,"row"),
    .combine=rbind
  ) %dofuture% {

    if (is.na(bp[1])){

      joint_repped_df -> df

      ##Models
      models <- list(
        m1=lm(R~poly(lagE,2,raw=F):(box-1)+(box-1),data=df),
        m2=lm(R~(poly(lagE,2,raw=F)+box-1),data=df),
        m3=lm(R~poly(lagE,2,raw=F),data=df),
        m4=lm(R~lagE:(box-1)+(box-1),data=df),
        m5=lm(R~(lagE+box-1),data=df),
        m6=lm(R~lagE,data=df)
      )

      models |>
        lapply(\(m) tibble(AIC=AIC(m),BIC=BIC(m),bp)) |>
        bind_rows(.id="model")

    } else {

      joint_repped_df |>
        mutate(
          phase=as.factor(if_else(time<=unlist(bp)[box],1,2))
        ) -> df

      ##Models
      models <- list(
        m1=lm(R~poly(lagE,2,raw=F):(box-1):(phase-1)+(box-1):(phase-1),data=df),
        m2=lm(R~(poly(lagE,2,raw=F)+box-1):(phase-1)+(phase-1),data=df),
        m3=lm(R~poly(lagE,2,raw=F):(phase-1)+(phase-1),data=df),
        m4=lm(R~lagE:(box-1):(phase-1)+(box-1):(phase-1),data=df),
        m5=lm(R~(lagE+box-1):(phase-1)+(phase-1),data=df),
        m6=lm(R~lagE:(phase-1)+(phase-1),data=df)
      )

      models |>
        lapply(\(m) tibble(AIC=AIC(m),BIC=BIC(m),bp)) |>
        bind_rows(.id="model")

    }
  }
}) -> joint_AIC_df

##Extract breakpoints from best model (lowest AIC)
joint_AIC_df |>
  filter(AIC==min(AIC)) -> best

##Obtain data frame with above breakpoints specified
joint_repped_df |>
  mutate(
    phase=as.factor(if_else(time<=as.integer(unlist(best)[box]),1,2))
  ) -> df

##Extract best model based on lowest AIC
models <- list(
  m1=lm(R~poly(lagE,2,raw=F):(box-1):(phase-1)+(box-1):(phase-1),data=df),
  m2=lm(R~(poly(lagE,2,raw=F)+box-1):(phase-1)+(phase-1),data=df),
  m3=lm(R~poly(lagE,2,raw=F):(phase-1)+(phase-1),data=df),
  m4=lm(R~lagE:(box-1):(phase-1)+(box-1):(phase-1),data=df),
  m5=lm(R~(lagE+box-1):(phase-1)+(phase-1),data=df),
  m6=lm(R~lagE:(phase-1)+(phase-1),data=df)
)

chosen_model <- models[[best$model]]

##Obtain confidence intervals for model parameters (2.5%, 97.5%)
library(jtools)
chosen_model_summ <- summ(
  chosen_model,
  confint=TRUE,
  ci.width=0.95
)

coeftab <- chosen_model_summ$coeftable
coefs <- runif_design(
  lower=coeftab[,"2.5%"],
  upper=coeftab[,"97.5%"],
  nseq=1000
)

model.matrix(chosen_model) %*% t(coefs) |>
  apply(1L,quantile,probs=c(lo=0.025,med=0.5,hi=0.975),names=FALSE) |>
  t() -> predq
colnames(predq) <- c("lo","med","hi")

df |>
  bind_cols(predq) -> df

df$pABA <- factor(df$box,
                  levels=c("box04","box03","box02","box01"),
                  labels=c("0%","0.0005%","0.005%","0.05%"))

(ModelFitting_ER <- df |> filter(lo>0) |>
    arrange(lagE) |>
    ggplot(aes(x=lagE,y=med,ymin=lo,ymax=hi,color=pABA,fill=pABA,
               group=interaction(pABA,phase)))+
    geom_line()+
    geom_ribbon(alpha=0.2,color=NA)+
    theme_bw()+
    xlab("Erythrocytes (d-1)")+ylab("Predicted reticulocyte supply (d)")+
    scale_colour_manual(values=cbPalette[2:5])+
    scale_fill_manual(values=cbPalette[2:5])+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=13),
      panel.grid=element_blank(),
      legend.position=c(0.8,0.8),
      legend.title=element_text(size=15),
      legend.text=element_text(size=13)
    )
)

library(stringi)
library(doParallel)
library(panelPomp)
library(ggpubr)
source("POMP_GroupLevel_DataPrep.R")

(RE_facet <- group_traj |>
    filter(box!="05",variable%in%c("E","R")) |>
    select(-lo,-hi) |>
    pivot_wider(names_from="variable",values_from="med") |>
    mutate(lagE=lag(E)) |>
    filter(time>0,time<=20) |>
    ggplot()+
    geom_path(aes(x=lagE,y=R,col=pABA),linewidth=2)+
    geom_text(aes(x=lagE,y=R,label=time),size=5)+
    scale_colour_manual(values=cbPalette[2:5])+
    facet_wrap(pABA~.)+
    xlab("Erythrocytes (d-1)")+ylab("Reticulocyte supply (d)")+
    theme_bw()+
    theme(
      axis.title=element_text(size=15),
      axis.text=element_text(size=11),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13),
      strip.text=element_text(size=12),
      strip.background=element_blank()
      
    )
)

joint <- ggarrange(RE_facet,ModelFitting_ER, nrow = 1, labels = c("A","B"))
joint

ggsave("facet_modelFitting.png",plot=joint,
       width=30,height=15,units="cm") 

