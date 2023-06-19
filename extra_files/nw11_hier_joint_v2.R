## -----------------------------------------------------------------------------
library(tidyverse)
library(stringi)
library(pomp)
options(dplyr.summarise.inform=FALSE)


setwd("~/Documents/GitHub/bdd/nw11_hier/")
## -----------------------------------------------------------------------------

read_csv(
  "data.csv",
  col_types="iiinnnn"
) |>
  mutate(
    mouseid=sprintf("%02d-%02d",box,mouse),
    box=sprintf("%02d",box),
    mouse=sprintf("%02d",mouse),
    paba=as.factor(box),
    rbc_density=rbc_density/1000,
    paba=case_match(
      paba,
      "01"~"0.05",
      "02"~"0.005",
      "03"~"0.0005",
      "04"~"0",
      "05"~"control"
    )
  ) |>
  select(
    day,
    mouseid,
    paba,
    box,
    mouse,
    Pd=ama_density,
    RBC=rbc_density,
    Ter119=ter119_density,
    CD71=cd71_density
  ) |>
  arrange(mouseid,day) |>
  mutate(
    paba=as.character(paba),
    Ter119=if_else(Ter119==0,NA_real_,Ter119),
    CD71=if_else(CD71==0,NA_real_,CD71),
    Eryth=(1-CD71/Ter119)*RBC,
    Retic=CD71/Ter119*RBC,
    mouse=case_when(
      mouse=="01"~"A",
      mouse=="02"~"B",
      mouse=="03"~"C",
      TRUE~NA_character_
    )
  ) -> dat

stopifnot(!is.na(dat$mouse))

## -----------------------------------------------------------------------------
#Pomp object for box 2 (pABA = 0.005)
dat |>
  filter(box=="02",mouse!="C") |>
  select(day,mouse,Retic,RBC,Pd) |>
  filter(day<=20) |>
  pivot_wider(values_from=c(Retic,RBC,Pd),names_from=mouse) |>
  pomp(
    times="day",t0=0,
    dmeasure=Csnippet(r"{
      double rbc;
      double K;
      double tweak = 0.01;
      lik = 0;
      // first mouse
      rbc = exp(logR_A)+exp(logE_A);
      K = rbc*(1-exp(-exp(logM_A)/rbc));
      lik += (R_FINITE(Retic_A)) ? dnorm(log(Retic_A),logR_A,sigmaRetic+tweak,1) : 0;
      lik += (R_FINITE(RBC_A)) ? dnorm(log(RBC_A),log(rbc),sigmaRBC+tweak,1) : 0;
      lik += (R_FINITE(Pd_A)) ? dnorm(log(Pd_A),log(K),sigmaPd+tweak,1) : 0;
      // second mouse
      rbc = exp(logR_B)+exp(logE_B);
      K = rbc*(1-exp(-exp(logM_B)/rbc));
      lik += (R_FINITE(Retic_B)) ? dnorm(log(Retic_B),logR_B,sigmaRetic+tweak,1) : 0;
      lik += (R_FINITE(RBC_B)) ? dnorm(log(RBC_B),log(rbc),sigmaRBC+tweak,1) : 0;
      lik += (R_FINITE(Pd_B)) ? dnorm(log(Pd_B),log(K),sigmaPd+tweak,1) : 0;
      if (!give_log) lik = exp(lik);
    }"
    ),
    dprocess=Csnippet(r"{
      double rbc, K, md, ed;
      // GAUSSIAN MARKOV RF MODEL:
      loglik = dnorm(logR_2,logR_1,sigmaR,1) +
               dnorm(logN_2,logN_1,sigmaN,1) +
               dnorm(logW_2,logW_1,sigmaW,1);
      // FIRST MOUSE (A):
      rbc = exp(logR_A_1)+exp(logE_A_1);
      K = rbc*(1-exp(-exp(logM_A_1)/rbc));
      // unlawfulness penalty
      md = logM_A_2-log(Beta*K*exp(-(exp(logW_A_1)+exp(logN_A_1))/rbc));
      ed = logE_A_2-log(rbc*exp(-(exp(logM_A_1)+exp(logN_A_1))/rbc));
      loglik -= lambdaM*md*md + lambdaE*ed*ed;
      // AR(1) individual deviation
      if (t_1 < 1) {
        loglik += dnorm(logW_A_1-logW_1,0,sigmaw/sqrt(1-alphaw),1) +
                  dnorm(logN_A_1-logN_1,0,sigman/sqrt(1-alphan),1) +
                  dnorm(logR_A_1-logR_1,0,sigmar/sqrt(1-alphar),1);
      }
      loglik += dnorm(logW_A_2-logW_2,alphaw*(logW_A_1-logW_1),sigmaw,1) +
                dnorm(logN_A_2-logN_2,alphan*(logN_A_1-logN_1),sigman,1) +
                dnorm(logR_A_2-logR_2,alphar*(logR_A_1-logR_1),sigmar,1);
      // SECOND MOUSE (B):
      rbc = exp(logR_B_1)+exp(logE_B_1);
      K = rbc*(1-exp(-exp(logM_B_1)/rbc));
      // unlawfulness penalty
      md = logM_B_2-log(Beta*K*exp(-(exp(logW_B_1)+exp(logN_B_1))/rbc));
      ed = logE_B_2-log(rbc*exp(-(exp(logM_B_1)+exp(logN_B_1))/rbc));
      loglik -= lambdaM*md*md + lambdaE*ed*ed;
      // AR(1) individual deviation
      if (t_1 < 1) {
        loglik += dnorm(logW_B_1-logW_1,0,sigmaw/sqrt(1-alphaw),1) +
                  dnorm(logN_B_1-logN_1,0,sigman/sqrt(1-alphan),1) +
                  dnorm(logR_B_1-logR_1,0,sigmar/sqrt(1-alphar),1);
      }
      loglik += dnorm(logW_B_2-logW_2,alphaw*(logW_B_1-logW_1),sigmaw,1) +
                dnorm(logN_B_2-logN_2,alphan*(logN_B_1-logN_1),sigman,1) +
                dnorm(logR_B_2-logR_2,alphar*(logR_B_1-logR_1),sigmar,1);
    }"
    ),
    partrans=parameter_trans(
      log=c(
        "sigmaRBC","sigmaPd","sigmaRetic",
        "lambdaM","lambdaE",
        "sigmaR","sigmaW","sigmaN",
        "sigmaw","sigman","sigmar"
      ),
      logit=c(
        "alphaw","alphan","alphar"
      )        
    ),
    params=c(
      Beta=6.09, # fixed for now
      lambdaM=1000,lambdaE=1000, # unlawfulness penalty: chosen to be large
      sigmaR=0.75,sigmaW=0.75,sigmaN=0.75, # fixed at these values in example 9
      ## based on example 9 estimates above
      alphan=0.5,sigman=0.6,
      alphar=0.8,sigmar=0.25,
      alphaw=0.85,sigmaw=0.3,
      sigmaRBC=0.1,sigmaRetic=0.1,sigmaPd=0.6 # chosen to be bigger than needed before
    ),
    paramnames=c(
      "Beta",
      "lambdaM","lambdaE",
      "sigmaR","sigmaW","sigmaN",
      "alphaw","alphan","alphar",
      "sigmaw","sigman","sigmar",
      "sigmaRBC","sigmaPd","sigmaRetic"
    ),
    statenames=c(
      "logM","logE","logR","logW","logN",
      outer(
        c("logE","logM","logR","logW","logN"),
        c("A","B"),
        paste,sep="_"
      )
    )
  ) -> po02

## -----------------------------------------------------------------------------
#Pomp object for box 3 (pABA = 0.0005)
dat |>
  filter(box=="03") |>
  select(day,mouse,Retic,RBC,Pd) |>
  filter(day<=20) |>
  pivot_wider(values_from=c(Retic,RBC,Pd),names_from=mouse) |>
  pomp(
    times="day",t0=0,
    dmeasure=Csnippet(r"{
      double rbc;
      double K;
      double tweak = 0.01;
      lik = 0;
      // first mouse
      rbc = exp(logR_A)+exp(logE_A);
      K = rbc*(1-exp(-exp(logM_A)/rbc));
      lik += (R_FINITE(Retic_A)) ? dnorm(log(Retic_A),logR_A,sigmaRetic+tweak,1) : 0;
      lik += (R_FINITE(RBC_A)) ? dnorm(log(RBC_A),log(rbc),sigmaRBC+tweak,1) : 0;
      lik += (R_FINITE(Pd_A)) ? dnorm(log(Pd_A),log(K),sigmaPd+tweak,1) : 0;
      // second mouse
      rbc = exp(logR_B)+exp(logE_B);
      K = rbc*(1-exp(-exp(logM_B)/rbc));
      lik += (R_FINITE(Retic_B)) ? dnorm(log(Retic_B),logR_B,sigmaRetic+tweak,1) : 0;
      lik += (R_FINITE(RBC_B)) ? dnorm(log(RBC_B),log(rbc),sigmaRBC+tweak,1) : 0;
      lik += (R_FINITE(Pd_B)) ? dnorm(log(Pd_B),log(K),sigmaPd+tweak,1) : 0;
      // third mouse
      rbc = exp(logR_C)+exp(logE_C);
      K = rbc*(1-exp(-exp(logM_C)/rbc));
      lik += (R_FINITE(Retic_C)) ? dnorm(log(Retic_C),logR_C,sigmaRetic+tweak,1) : 0;
      lik += (R_FINITE(RBC_C)) ? dnorm(log(RBC_C),log(rbc),sigmaRBC+tweak,1) : 0;
      lik += (R_FINITE(Pd_C)) ? dnorm(log(Pd_C),log(K),sigmaPd+tweak,1) : 0;
      if (!give_log) lik = exp(lik);
    }"
    ),
    dprocess=Csnippet(r"{
      double rbc, K, md, ed;
      // GAUSSIAN MARKOV RF MODEL:
      loglik = dnorm(logR_2,logR_1,sigmaR,1) +
               dnorm(logN_2,logN_1,sigmaN,1) +
               dnorm(logW_2,logW_1,sigmaW,1);
      // FIRST MOUSE (A):
      rbc = exp(logR_A_1)+exp(logE_A_1);
      K = rbc*(1-exp(-exp(logM_A_1)/rbc));
      // unlawfulness penalty
      md = logM_A_2-log(Beta*K*exp(-(exp(logW_A_1)+exp(logN_A_1))/rbc));
      ed = logE_A_2-log(rbc*exp(-(exp(logM_A_1)+exp(logN_A_1))/rbc));
      loglik -= lambdaM*md*md + lambdaE*ed*ed;
      // AR(1) individual deviation
      if (t_1 < 1) {
        loglik += dnorm(logW_A_1-logW_1,0,sigmaw/sqrt(1-alphaw),1) +
                  dnorm(logN_A_1-logN_1,0,sigman/sqrt(1-alphan),1) +
                  dnorm(logR_A_1-logR_1,0,sigmar/sqrt(1-alphar),1);
      }
      loglik += dnorm(logW_A_2-logW_2,alphaw*(logW_A_1-logW_1),sigmaw,1) +
                dnorm(logN_A_2-logN_2,alphan*(logN_A_1-logN_1),sigman,1) +
                dnorm(logR_A_2-logR_2,alphar*(logR_A_1-logR_1),sigmar,1);
      // SECOND MOUSE (B):
      rbc = exp(logR_B_1)+exp(logE_B_1);
      K = rbc*(1-exp(-exp(logM_B_1)/rbc));
      // unlawfulness penalty
      md = logM_B_2-log(Beta*K*exp(-(exp(logW_B_1)+exp(logN_B_1))/rbc));
      ed = logE_B_2-log(rbc*exp(-(exp(logM_B_1)+exp(logN_B_1))/rbc));
      loglik -= lambdaM*md*md + lambdaE*ed*ed;
      // AR(1) individual deviation
      if (t_1 < 1) {
        loglik += dnorm(logW_B_1-logW_1,0,sigmaw/sqrt(1-alphaw),1) +
                  dnorm(logN_B_1-logN_1,0,sigman/sqrt(1-alphan),1) +
                  dnorm(logR_B_1-logR_1,0,sigmar/sqrt(1-alphar),1);
      }
      loglik += dnorm(logW_B_2-logW_2,alphaw*(logW_B_1-logW_1),sigmaw,1) +
                dnorm(logN_B_2-logN_2,alphan*(logN_B_1-logN_1),sigman,1) +
                dnorm(logR_B_2-logR_2,alphar*(logR_B_1-logR_1),sigmar,1);
      // THIRD MOUSE (C):
      rbc = exp(logR_C_1)+exp(logE_C_1);
      K = rbc*(1-exp(-exp(logM_C_1)/rbc));
      // unlawfulness penalty
      md = logM_C_2-log(Beta*K*exp(-(exp(logW_C_1)+exp(logN_C_1))/rbc));
      ed = logE_C_2-log(rbc*exp(-(exp(logM_C_1)+exp(logN_C_1))/rbc));
      loglik -= lambdaM*md*md + lambdaE*ed*ed;
      // AR(1) individual deviation
      if (t_1 < 1) {
        loglik += dnorm(logW_C_1-logW_1,0,sigmaw/sqrt(1-alphaw),1) +
                  dnorm(logN_C_1-logN_1,0,sigman/sqrt(1-alphan),1) +
                  dnorm(logR_C_1-logR_1,0,sigmar/sqrt(1-alphar),1);
      }
      loglik += dnorm(logW_C_2-logW_2,alphaw*(logW_C_1-logW_1),sigmaw,1) +
                dnorm(logN_C_2-logN_2,alphan*(logN_C_1-logN_1),sigman,1) +
                dnorm(logR_C_2-logR_2,alphar*(logR_C_1-logR_1),sigmar,1);
    }"
    ),
    partrans=parameter_trans(
      log=c(
        "sigmaRBC","sigmaPd","sigmaRetic",
        "lambdaM","lambdaE",
        "sigmaR","sigmaW","sigmaN",
        "sigmaw","sigman","sigmar"
      ),
      logit=c(
        "alphaw","alphan","alphar"
      )        
    ),
    params=c(
      Beta=5.76, # fixed for now
      lambdaM=1000,lambdaE=1000, # unlawfulness penalty: chosen to be large
      sigmaR=0.75,sigmaW=0.75,sigmaN=0.75, # fixed at these values in example 9
      ## based on example 9 estimates above
      alphan=0.5,sigman=0.6,
      alphar=0.8,sigmar=0.25,
      alphaw=0.85,sigmaw=0.3,
      sigmaRBC=0.1,sigmaRetic=0.1,sigmaPd=0.6 # chosen to be bigger than needed before
    ),
    paramnames=c(
      "Beta",
      "lambdaM","lambdaE",
      "sigmaR","sigmaW","sigmaN",
      "alphaw","alphan","alphar",
      "sigmaw","sigman","sigmar",
      "sigmaRBC","sigmaPd","sigmaRetic"
    ),
    statenames=c(
      "logM","logE","logR","logW","logN",
      outer(
        c("logE","logM","logR","logW","logN"),
        c("A","B","C"),
        paste,sep="_"
      )
    )
  ) -> po03

## -----------------------------------------------------------------------------
#Pomp object for box 4 (pABA = 0)
dat |>
  filter(box=="04") |>
  select(day,mouse,Retic,RBC,Pd) |>
  filter(day<=20) |>
  pivot_wider(values_from=c(Retic,RBC,Pd),names_from=mouse) |>
  pomp(
    times="day",t0=0,
    dmeasure=Csnippet(r"{
      double rbc;
      double K;
      double tweak = 0.01;
      lik = 0;
      // first mouse
      rbc = exp(logR_A)+exp(logE_A);
      K = rbc*(1-exp(-exp(logM_A)/rbc));
      lik += (R_FINITE(Retic_A)) ? dnorm(log(Retic_A),logR_A,sigmaRetic+tweak,1) : 0;
      lik += (R_FINITE(RBC_A)) ? dnorm(log(RBC_A),log(rbc),sigmaRBC+tweak,1) : 0;
      lik += (R_FINITE(Pd_A)) ? dnorm(log(Pd_A),log(K),sigmaPd+tweak,1) : 0;
      // second mouse
      rbc = exp(logR_B)+exp(logE_B);
      K = rbc*(1-exp(-exp(logM_B)/rbc));
      lik += (R_FINITE(Retic_B)) ? dnorm(log(Retic_B),logR_B,sigmaRetic+tweak,1) : 0;
      lik += (R_FINITE(RBC_B)) ? dnorm(log(RBC_B),log(rbc),sigmaRBC+tweak,1) : 0;
      lik += (R_FINITE(Pd_B)) ? dnorm(log(Pd_B),log(K),sigmaPd+tweak,1) : 0;
      // third mouse
      rbc = exp(logR_C)+exp(logE_C);
      K = rbc*(1-exp(-exp(logM_C)/rbc));
      lik += (R_FINITE(Retic_C)) ? dnorm(log(Retic_C),logR_C,sigmaRetic+tweak,1) : 0;
      lik += (R_FINITE(RBC_C)) ? dnorm(log(RBC_C),log(rbc),sigmaRBC+tweak,1) : 0;
      lik += (R_FINITE(Pd_C)) ? dnorm(log(Pd_C),log(K),sigmaPd+tweak,1) : 0;
      if (!give_log) lik = exp(lik);
    }"
    ),
    dprocess=Csnippet(r"{
      double rbc, K, md, ed;
      // GAUSSIAN MARKOV RF MODEL:
      loglik = dnorm(logR_2,logR_1,sigmaR,1) +
               dnorm(logN_2,logN_1,sigmaN,1) +
               dnorm(logW_2,logW_1,sigmaW,1);
      // FIRST MOUSE (A):
      rbc = exp(logR_A_1)+exp(logE_A_1);
      K = rbc*(1-exp(-exp(logM_A_1)/rbc));
      // unlawfulness penalty
      md = logM_A_2-log(Beta*K*exp(-(exp(logW_A_1)+exp(logN_A_1))/rbc));
      ed = logE_A_2-log(rbc*exp(-(exp(logM_A_1)+exp(logN_A_1))/rbc));
      loglik -= lambdaM*md*md + lambdaE*ed*ed;
      // AR(1) individual deviation
      if (t_1 < 1) {
        loglik += dnorm(logW_A_1-logW_1,0,sigmaw/sqrt(1-alphaw),1) +
                  dnorm(logN_A_1-logN_1,0,sigman/sqrt(1-alphan),1) +
                  dnorm(logR_A_1-logR_1,0,sigmar/sqrt(1-alphar),1);
      }
      loglik += dnorm(logW_A_2-logW_2,alphaw*(logW_A_1-logW_1),sigmaw,1) +
                dnorm(logN_A_2-logN_2,alphan*(logN_A_1-logN_1),sigman,1) +
                dnorm(logR_A_2-logR_2,alphar*(logR_A_1-logR_1),sigmar,1);
      // SECOND MOUSE (B):
      rbc = exp(logR_B_1)+exp(logE_B_1);
      K = rbc*(1-exp(-exp(logM_B_1)/rbc));
      // unlawfulness penalty
      md = logM_B_2-log(Beta*K*exp(-(exp(logW_B_1)+exp(logN_B_1))/rbc));
      ed = logE_B_2-log(rbc*exp(-(exp(logM_B_1)+exp(logN_B_1))/rbc));
      loglik -= lambdaM*md*md + lambdaE*ed*ed;
      // AR(1) individual deviation
      if (t_1 < 1) {
        loglik += dnorm(logW_B_1-logW_1,0,sigmaw/sqrt(1-alphaw),1) +
                  dnorm(logN_B_1-logN_1,0,sigman/sqrt(1-alphan),1) +
                  dnorm(logR_B_1-logR_1,0,sigmar/sqrt(1-alphar),1);
      }
      loglik += dnorm(logW_B_2-logW_2,alphaw*(logW_B_1-logW_1),sigmaw,1) +
                dnorm(logN_B_2-logN_2,alphan*(logN_B_1-logN_1),sigman,1) +
                dnorm(logR_B_2-logR_2,alphar*(logR_B_1-logR_1),sigmar,1);
      // THIRD MOUSE (C):
      rbc = exp(logR_C_1)+exp(logE_C_1);
      K = rbc*(1-exp(-exp(logM_C_1)/rbc));
      // unlawfulness penalty
      md = logM_C_2-log(Beta*K*exp(-(exp(logW_C_1)+exp(logN_C_1))/rbc));
      ed = logE_C_2-log(rbc*exp(-(exp(logM_C_1)+exp(logN_C_1))/rbc));
      loglik -= lambdaM*md*md + lambdaE*ed*ed;
      // AR(1) individual deviation
      if (t_1 < 1) {
        loglik += dnorm(logW_C_1-logW_1,0,sigmaw/sqrt(1-alphaw),1) +
                  dnorm(logN_C_1-logN_1,0,sigman/sqrt(1-alphan),1) +
                  dnorm(logR_C_1-logR_1,0,sigmar/sqrt(1-alphar),1);
      }
      loglik += dnorm(logW_C_2-logW_2,alphaw*(logW_C_1-logW_1),sigmaw,1) +
                dnorm(logN_C_2-logN_2,alphan*(logN_C_1-logN_1),sigman,1) +
                dnorm(logR_C_2-logR_2,alphar*(logR_C_1-logR_1),sigmar,1);
    }"
    ),
    partrans=parameter_trans(
      log=c(
        "sigmaRBC","sigmaPd","sigmaRetic",
        "lambdaM","lambdaE",
        "sigmaR","sigmaW","sigmaN",
        "sigmaw","sigman","sigmar"
      ),
      logit=c(
        "alphaw","alphan","alphar"
      )        
    ),
    params=c(
      Beta=4.08, # fixed for now
      lambdaM=1000,lambdaE=1000, # unlawfulness penalty: chosen to be large
      sigmaR=0.75,sigmaW=0.75,sigmaN=0.75, # fixed at these values in example 9
      ## based on example 9 estimates above
      alphan=0.5,sigman=0.6,
      alphar=0.8,sigmar=0.25,
      alphaw=0.85,sigmaw=0.3,
      sigmaRBC=0.1,sigmaRetic=0.1,sigmaPd=0.6 # chosen to be bigger than needed before
    ),
    paramnames=c(
      "Beta",
      "lambdaM","lambdaE",
      "sigmaR","sigmaW","sigmaN",
      "alphaw","alphan","alphar",
      "sigmaw","sigman","sigmar",
      "sigmaRBC","sigmaPd","sigmaRetic"
    ),
    statenames=c(
      "logM","logE","logR","logW","logN",
      outer(
        c("logE","logM","logR","logW","logN"),
        c("A","B","C"),
        paste,sep="_"
      )
    )
  ) -> po04


## -----------------------------------------------------------------------------
#Initial trajectories for box 2 (pABA = 0.005)
po_df <- read_csv(
  "po_df_Example9.csv",
  col_types="nnnnnnnnnncccnnnnnnnnnnn"
) |>
  mutate(
    mouse=case_when(
      mouse=="01"~"A",
      mouse=="02"~"B",
      mouse=="03"~"C",
      TRUE~NA_character_
    )
  ) |>
  filter(box=="02",mouse!="C") |>
  select(day,mouse,logM=M,logE=E,logR=R,logW=W,logN=N)

po_df |> pivot_longer(-c(day,mouse)) |>
  mutate(value=exp(value)) |>
  group_by(day,name) |>
  summarise(mean=log(mean(value))) |>
  ungroup() |>
  pivot_wider(names_from=name,values_from=mean) -> po_avg

po_df |> pivot_longer(-c(day,mouse)) |>
  pivot_wider(names_from=c(name,mouse),values_from=value) -> po_ind

full_join(po_avg,po_ind,by="day") |>
  ungroup() |>
  select(-day) |>
  t() -> coefs02

## -----------------------------------------------------------------------------
#Initial trajectories for box 3 (pABA = 0.0005)
po_df <- read_csv(
  "po_df_Example9.csv",
  col_types="nnnnnnnnnncccnnnnnnnnnnn"
) |>
  mutate(
    mouse=case_when(
      mouse=="01"~"A",
      mouse=="02"~"B",
      mouse=="03"~"C",
      TRUE~NA_character_
    )
  ) |>
  filter(box=="03") |>
  select(day,mouse,logM=M,logE=E,logR=R,logW=W,logN=N)

po_df |> pivot_longer(-c(day,mouse)) |>
  mutate(value=exp(value)) |>
  group_by(day,name) |>
  summarise(mean=log(mean(value))) |>
  ungroup() |>
  pivot_wider(names_from=name,values_from=mean) -> po_avg

po_df |> pivot_longer(-c(day,mouse)) |>
  pivot_wider(names_from=c(name,mouse),values_from=value) -> po_ind

full_join(po_avg,po_ind,by="day") |>
  ungroup() |>
  select(-day) |>
  t() -> coefs03

## -----------------------------------------------------------------------------
#Initial trajectories for box 4 (pABA = 0)
po_df <- read_csv(
  "po_df_Example9.csv",
  col_types="nnnnnnnnnncccnnnnnnnnnnn"
) |>
  mutate(
    mouse=case_when(
      mouse=="01"~"A",
      mouse=="02"~"B",
      mouse=="03"~"C",
      TRUE~NA_character_
    )
  ) |>
  filter(box=="04") |>
  select(day,mouse,logM=M,logE=E,logR=R,logW=W,logN=N)

po_df |> pivot_longer(-c(day,mouse)) |>
  mutate(value=exp(value)) |>
  group_by(day,name) |>
  summarise(mean=log(mean(value))) |>
  ungroup() |>
  pivot_wider(names_from=name,values_from=mean) -> po_avg

po_df |> pivot_longer(-c(day,mouse)) |>
  pivot_wider(names_from=c(name,mouse),values_from=value) -> po_ind

full_join(po_avg,po_ind,by="day") |>
  ungroup() |>
  select(-day) |>
  t() -> coefs04


## -----------------------------------------------------------------------------

create_objfun <- function (
    object1, params1 = coef(object1), coefs1,
    object2, params2 = coef(object2), coefs2,
    object3, params3 = coef(object3), coefs3,
    est,
    control = list(reltol = 1e-3, maxit = 1e6)
) {
  
  params1 <- partrans(object1,params1,dir="toEst")
  params2 <- partrans(object2,params2,dir="toEst")
  params3 <- partrans(object3,params3,dir="toEst")
  
  idx <- match(est,names(params1)) #vector of positions
  
  ofun <- function (par) {
    
    params1[idx] <- par
    params2[idx] <- par
    params3[idx] <- par
    
    coef(object1,transform=TRUE) <<- params1
    coef(object2,transform=TRUE) <<- params2
    coef(object3,transform=TRUE) <<- params3
    
    fit1 <- optim(
      fn = function (x) {
        
        dim(x)<-dim(coefs1)
        dimnames(x) <- dimnames(coefs1)
        
        object1@states <<- x
        
        
        object1 |> dprocess(log=TRUE) |> sum() -> ss
        object1 |> dmeasure(log=TRUE) |> sum() -> ll
        
        -(ss+ll)
        
      },
      par=coefs1,
      method="BFGS",
      control=control
    )
    
    fit2 <- optim(
      fn = function (x) {
        
        dim(x)<-dim(coefs2)
        dimnames(x) <- dimnames(coefs2)
        
        object2@states <<- x
        
        object2 |> dprocess(log=TRUE) |> sum() -> ss
        object2 |> dmeasure(log=TRUE) |> sum() -> ll
        
        -(ss+ll)
        
      },
      par=coefs2,
      method="BFGS",
      control=control
    )
    
    fit3 <- optim(
      fn = function (x) {
        
        dim(x)<-dim(coefs3)
        dimnames(x) <- dimnames(coefs3)
        
        object3@states <<- x
        
        object3 |> dprocess(log=TRUE) |> sum() -> ss
        object3 |> dmeasure(log=TRUE) |> sum() -> ll
        
        -(ss+ll)
        
      },
      par=coefs3,
      method="BFGS",
      control=control
    )
    
    coefs1 <<- fit1$par
    coefs2 <<- fit2$par
    coefs3 <<- fit3$par
    
    fit1$value+fit2$value+fit3$value
  }
  environment(ofun) <- list2env(
    list(object1=object1,params1=params1,coefs1=coefs1,
         object2=object2,params2=params2,coefs2=coefs2,
         object3=object3,params3=params3,coefs3=coefs3,
         idx=idx,control=control),
    parent=parent.frame(2)
  )
  ofun
}

## -----------------------------------------------------------------------------

stew(
  file="nw11_hier_joint_v2.rda",
  info=TRUE,
  {
    create_objfun(
      object1=po02,
      object2=po03,
      object3=po04,
      coefs1=coefs02,
      coefs2=coefs03,
      coefs3=coefs04,
      est=c("alphaw","alphan","alphar",
            "sigmaw","sigman","sigmar",
            "sigmaRBC","sigmaPd","sigmaRetic"),
    ) -> ofun
    try({
      optim(
        fn=ofun,
        par=coef(po02,c("alphaw","alphan","alphar",
                        "sigmaw","sigman","sigmar",
                        "sigmaRBC","sigmaPd","sigmaRetic"),transform=TRUE),
        control=list(maxit=500,reltol=1e-4,trace=1)
      ) -> fit
      ofun(fit$par)
    })
  })

evalq(coef(object1),envir=environment(ofun)) |> melt() -> est.coefs02
est.coefs02$box = "02"
evalq(coef(object2),envir=environment(ofun)) |> melt() -> est.coefs03
est.coefs03$box = "03"
evalq(coef(object3),envir=environment(ofun)) |> melt() -> est.coefs04
est.coefs04$box = "04"

est.coefs<-bind_rows(est.coefs02,est.coefs03,est.coefs04)

write.csv(est.coefs,"est_coefs_joint.csv",row.names=FALSE)

evalq(coefs1,envir=environment(ofun)) |>
  melt() |>
  as_tibble() |>
  separate(Var1,into=c("name","mouse")) |>
  mutate(
    time=time(po02)[Var2],
    mouse=if_else(is.na(mouse),"GROUP",mouse),
    value=exp(value),
    box="02"
  ) |>
  select(-Var2) -> est.po02

evalq(coefs2,envir=environment(ofun)) |>
  melt() |>
  as_tibble() |>
  separate(Var1,into=c("name","mouse")) |>
  mutate(
    time=time(po02)[Var2],
    mouse=if_else(is.na(mouse),"GROUP",mouse),
    value=exp(value),
    box="03"
  ) |>
  select(-Var2) -> est.po03

evalq(coefs3,envir=environment(ofun)) |>
  melt() |>
  as_tibble() |>
  separate(Var1,into=c("name","mouse")) |>
  mutate(
    time=time(po02)[Var2],
    mouse=if_else(is.na(mouse),"GROUP",mouse),
    value=exp(value),
    box="04"
  ) |>
  select(-Var2) -> est.po04

est.po<-bind_rows(est.po02,est.po03,est.po04)
write.csv(est.po,"est_po_joint.csv",row.names=FALSE)
