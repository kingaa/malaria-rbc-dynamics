#Load required packages
library(tidyverse)
library(stringi)
library(doParallel)
library(pomp)
library(panelPomp)
library(ggpubr)
library(scales)
library(gridExtra)
library(strucchange)
library(fUnitRoots)
library(forecast)
library(foreach)
library(iterators)
library(doRNG)
#Set working directory
setwd("~/Documents/GitHub/bdd/nw11_hier/final_files/Manuscript/")

read_csv("data.csv",
         col_types="iiinnnn"
) %>%
  mutate(
    mouseid=sprintf("%02d-%02d",box,mouse),
    box=sprintf("%02d",box),
    mouse=sprintf("%02d",mouse),
    paba=as.factor(box),
    rbc_density=rbc_density/1000
  ) %>%
  mutate(
    paba=recode(
      paba,
      "01"="0.05","02"="0.005","03"="0.0005","04"="0","05"="control"
    )
  ) %>%
  select(
    day,
    Pd=ama_density,
    RBC=rbc_density,
    Ter119=ter119_density,
    CD71=cd71_density,
    mouseid,
    paba,box,mouse
  ) %>%
  arrange(mouseid,day) %>%
  mutate(
    paba=as.character(paba),
    Ter119=ifelse(Ter119==0,NA,Ter119),
    CD71=ifelse(CD71==0,NA,CD71)
  ) %>%
  mutate(
    Eryth=(1-CD71/Ter119)*RBC,
    Retic=CD71/Ter119*RBC
  ) -> flow

flow %>%
  filter(day<=4) %>%
  lm(log(Pd)~box:day+mouseid-1,data=.) -> fit2

coef(fit2) %>% 
  bind_rows() %>%
  gather(var,val) %>%
  mutate(
    var=stri_replace_all_regex(var,"mouseid(\\d{2})-(\\d{2})","dose[$1-$2]"),
    var=stri_replace_all_regex(var,"box(\\d{2}):day","Beta[$1]"),
    val=exp(val)
  ) -> theta1

expand.grid(
  box=sprintf("%02d",1:5),
  mouse=sprintf("%02d",1:3)
) %>%
  mutate(
    mouseid=paste0(box,"-",mouse),
    betavar=sprintf("Beta[%s]",box),
    dosevar=sprintf("dose[%s]",mouseid)
  ) %>%
  left_join(theta1,by=c("betavar"="var")) %>%
  dplyr::rename(Beta=val) %>%
  left_join(theta1,by=c("dosevar"="var")) %>%
  dplyr::rename(dose=val) %>%
  select(-betavar,-dosevar) %>%
  arrange(box,mouse) %>%
  mutate(
    Beta=coalesce(Beta,0),
    dose=coalesce(dose,0)
  ) -> theta

library(doParallel)
registerDoParallel()

theta %>%
  mutate(
    sigmaPd = 2, 
    sigmaRBC = 0.1, 
    sigmaRetic = 0.3, 
    sigmaW = 1, 
    sigmaN = 0.5, 
    sigmaR = 0.5,
    E_0 = 8.0e6, 
    R_0 = 3e5, 
    W_0 = 8.8e4, 
    N_0 = 8e5
  ) -> theta

foreach (m = iter(theta,"row"),.inorder=TRUE,.combine=c) %dopar% {
  
  flow %>%
    filter(mouseid==m$mouseid) %>%
    select(day,Pd,RBC,Retic) %>%
    mutate(Retic=if_else(day %in% c(0,14),NA_real_,Retic)) %>%
    pomp(
      params=select(m,-mouseid,-box,-mouse) %>% unlist(),
      times="day",
      t0=0,
      rmeasure=Csnippet("
  Retic = rlnorm(log(1+R),sigmaRetic)-1;
  RBC = rlnorm(log(1+E+R),sigmaRBC)-1;
  Pd = rlnorm(log(1+K),sigmaPd)-1;"),
      dmeasure=Csnippet("
  double l1, l2, l3;
  l1 = (R_FINITE(Retic)) ? dlnorm(1+Retic,log(1+R),sigmaRetic,1) : 0;
  l2 = (R_FINITE(RBC)) ? dlnorm(1+RBC,log(1+E+R),sigmaRBC,1) : 0;
  l3 = (R_FINITE(Pd) && Pd>0) ? dlnorm(1+Pd,log(1+K),sigmaPd,1) : 0;
  lik = (give_log) ? l1+l2+l3 : exp(l1+l2+l3);"),
      rprocess=discrete_time(
        step.fun=Csnippet("
  double Mold = M;
  M = Beta*K*exp(-(W+N)/(R+E));
  E = (R+E)*exp(-(Mold+N)/(R+E));
  N = rlnorm(log(N),sigmaN);
  W = rlnorm(log(W),sigmaW);
  R = rlnorm(log(R),sigmaR);
  K = (R+E>0) ? (R+E)*(1-exp(-M/(R+E))): 0;
  "),
        delta.t=1
      ),
      partrans=parameter_trans(
        log=c("sigmaW","sigmaR","sigmaN",
              "sigmaPd","sigmaRBC","sigmaRetic",
              "N_0","W_0","E_0","R_0")
      ),
      rinit=Csnippet("
  E = E_0;
  R = R_0;
  N = N_0;
  W = W_0;
  M = 0;
  K = dose;"),
      statenames=c("E","R","W","N","M","K"),
      paramnames=c(
        "Beta","dose",
        "sigmaPd","sigmaRBC","sigmaRetic",
        "sigmaW","sigmaN","sigmaR",
        "E_0","R_0","W_0","N_0"
      )
    )
} %>% 
  set_names(theta$mouseid) -> pos

pf4 <- readRDS("m5pf4.rds")
pf4 %>%
  select(loglik,starts_with("sigma"),ends_with("_0")) %>%
  filter(loglik==max(loglik)) -> mle

bake(file="m5sm1_200.rds",{
  
  registerDoRNG(1510516029)
  
  foreach (i=1:200, 
           .inorder=FALSE,.combine=bind_rows,
           .packages=c("panelPomp","magrittr","reshape2","plyr","tibble","tidyr","dplyr")
  ) %dopar% 
    {
      
      mle %>% select(-loglik) %>% unlist() -> p
      panelPomp(pos,shared=p) -> m
      
      m %>%
        pfilter(Np=5e5,filter.traj=TRUE) %>%
        as("list") -> pfs
      
      pfs %>%
        plyr::llply(filter_traj) %>%
        melt() %>%
        dplyr::rename(mouseid=L1) %>%
        mutate(rep=i) %>%
        left_join(
          data.frame(mouseid=names(pfs),loglik=sapply(pfs,logLik)),
          by="mouseid"
        )
      
    } -> res
  attr(res,"ncpu") <- getDoParWorkers()
  res
  
}) -> sm1

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
  select(box,rep,mouse,loglik,lik) |> distinct() -> sm1_lik

sm1_lik |> filter(box=="04") |>
  ggplot()+
  geom_histogram(aes(x=log(lik)))+
  facet_grid(.~mouse)
