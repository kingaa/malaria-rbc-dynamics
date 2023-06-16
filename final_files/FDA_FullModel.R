## -----------------------------------------------------------------------------
library(tidyverse)
library(mgcv)
library(stringi)
library(pomp)
options(dplyr.summarise.inform=FALSE)


#setwd("~/Documents/GitHub/bdd/nw11_hier/final_files/")

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
#Pomp object for box 1 (pABA = 0.05)
dat |>
  filter(box=="01",mouse!="B") |>
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
      if (t < 9) {
        rbc = exp(logR_A)+exp(logE_A);
        K = rbc*(1-exp(-exp(logM_A)/rbc));
        lik += (R_FINITE(Retic_A)) ? dnorm(log(Retic_A),logR_A,sigmaRetic+tweak,1) : 0;
        lik += (R_FINITE(RBC_A)) ? dnorm(log(RBC_A),log(rbc),sigmaRBC+tweak,1) : 0;
        lik += (R_FINITE(Pd_A)) ? dnorm(log(Pd_A),log(K),sigmaPd+tweak,1) : 0;
      }
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
      double tweak2 = 0.0001;
      // GAUSSIAN MARKOV RF MODEL:
      loglik = dnorm(logR_2,logR_1,sigmaR,1) +
               dnorm(logN_2,logN_1,sigmaN,1) +
               dnorm(logW_2,logW_1,sigmaW,1);
      // FIRST MOUSE (A):
      if (t_1 < 8) {
        rbc = exp(logR_A_1)+exp(logE_A_1);
        K = rbc*(1-exp(-exp(logM_A_1)/rbc));
        // unlawfulness penalty
        md = logM_A_2-log(Beta*K*exp(-(exp(logW_A_1)+exp(logN_A_1))/rbc));
        ed = logE_A_2-log(rbc*exp(-(exp(logM_A_1)+exp(logN_A_1))/rbc));
        loglik -= lambdaM*md*md + lambdaE*ed*ed;
        // AR(1) individual deviation
        if (t_1 < 1) {
          loglik += dnorm(logW_A_1-logW_1,0,sigmaw/sqrt(1-(alphaw-tweak2)),1) +
                    dnorm(logN_A_1-logN_1,0,sigman/sqrt(1-(alphan-tweak2)),1) +
                    dnorm(logR_A_1-logR_1,0,sigmar/sqrt(1-(alphar-tweak2)),1);
        }
        loglik += dnorm(logW_A_2-logW_2,(alphaw-tweak2)*(logW_A_1-logW_1),sigmaw,1) +
                  dnorm(logN_A_2-logN_2,(alphan-tweak2)*(logN_A_1-logN_1),sigman,1) +
                  dnorm(logR_A_2-logR_2,(alphar-tweak2)*(logR_A_1-logR_1),sigmar,1);
      }
      // THIRD MOUSE (C):
      rbc = exp(logR_C_1)+exp(logE_C_1);
      K = rbc*(1-exp(-exp(logM_C_1)/rbc));
      // unlawfulness penalty
      md = logM_C_2-log(Beta*K*exp(-(exp(logW_C_1)+exp(logN_C_1))/rbc));
      ed = logE_C_2-log(rbc*exp(-(exp(logM_C_1)+exp(logN_C_1))/rbc));
      loglik -= lambdaM*md*md + lambdaE*ed*ed;
      // AR(1) individual deviation
      if (t_1 < 1) {
        loglik += dnorm(logW_C_1-logW_1,0,sigmaw/sqrt(1-(alphaw-tweak2)),1) +
                  dnorm(logN_C_1-logN_1,0,sigman/sqrt(1-(alphaw-tweak2)),1) +
                  dnorm(logR_C_1-logR_1,0,sigmar/sqrt(1-(alphaw-tweak2)),1);
      }
      loglik += dnorm(logW_C_2-logW_2,(alphaw-tweak2)*(logW_C_1-logW_1),sigmaw,1) +
                dnorm(logN_C_2-logN_2,(alphan-tweak2)*(logN_C_1-logN_1),sigman,1) +
                dnorm(logR_C_2-logR_2,(alphar-tweak2)*(logR_C_1-logR_1),sigmar,1);
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
      Beta=6.51, # fixed for now
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
        c("A","C"),
        paste,sep="_"
      )
    )
  ) -> po01

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
      double tweak2 = 0.0001;
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
        loglik += dnorm(logW_A_1-logW_1,0,sigmaw/sqrt(1-(alphaw-tweak2)),1) +
                  dnorm(logN_A_1-logN_1,0,sigman/sqrt(1-(alphan-tweak2)),1) +
                  dnorm(logR_A_1-logR_1,0,sigmar/sqrt(1-(alphar-tweak2)),1);
      }
      loglik += dnorm(logW_A_2-logW_2,(alphaw-tweak2)*(logW_A_1-logW_1),sigmaw,1) +
                dnorm(logN_A_2-logN_2,(alphan-tweak2)*(logN_A_1-logN_1),sigman,1) +
                dnorm(logR_A_2-logR_2,(alphar-tweak2)*(logR_A_1-logR_1),sigmar,1);
      // SECOND MOUSE (B):
      rbc = exp(logR_B_1)+exp(logE_B_1);
      K = rbc*(1-exp(-exp(logM_B_1)/rbc));
      // unlawfulness penalty
      md = logM_B_2-log(Beta*K*exp(-(exp(logW_B_1)+exp(logN_B_1))/rbc));
      ed = logE_B_2-log(rbc*exp(-(exp(logM_B_1)+exp(logN_B_1))/rbc));
      loglik -= lambdaM*md*md + lambdaE*ed*ed;
      // AR(1) individual deviation
      if (t_1 < 1) {
        loglik += dnorm(logW_B_1-logW_1,0,sigmaw/sqrt(1-(alphaw-tweak2)),1) +
                  dnorm(logN_B_1-logN_1,0,sigman/sqrt(1-(alphan-tweak2)),1) +
                  dnorm(logR_B_1-logR_1,0,sigmar/sqrt(1-(alphar-tweak2)),1);
      }
      loglik += dnorm(logW_B_2-logW_2,(alphaw-tweak2)*(logW_B_1-logW_1),sigmaw,1) +
                dnorm(logN_B_2-logN_2,(alphan-tweak2)*(logN_B_1-logN_1),sigman,1) +
                dnorm(logR_B_2-logR_2,(alphar-tweak2)*(logR_B_1-logR_1),sigmar,1);
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
      sigmaRBC=0.1,sigmaRetic=0.1,sigmaPd=0.2 # chosen to be bigger than needed before
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
      if (t < 19) {
        rbc = exp(logR_B)+exp(logE_B);
        K = rbc*(1-exp(-exp(logM_B)/rbc));
        lik += (R_FINITE(Retic_B)) ? dnorm(log(Retic_B),logR_B,sigmaRetic+tweak,1) : 0;
        lik += (R_FINITE(RBC_B)) ? dnorm(log(RBC_B),log(rbc),sigmaRBC+tweak,1) : 0;
        lik += (R_FINITE(Pd_B)) ? dnorm(log(Pd_B),log(K),sigmaPd+tweak,1) : 0;
      }
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
      double tweak2 = 0.0001;
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
        loglik += dnorm(logW_A_1-logW_1,0,sigmaw/sqrt(1-(alphaw-tweak2)),1) +
                  dnorm(logN_A_1-logN_1,0,sigman/sqrt(1-(alphan-tweak2)),1) +
                  dnorm(logR_A_1-logR_1,0,sigmar/sqrt(1-(alphar-tweak2)),1);
      }
      loglik += dnorm(logW_A_2-logW_2,(alphaw-tweak2)*(logW_A_1-logW_1),sigmaw,1) +
                dnorm(logN_A_2-logN_2,(alphan-tweak2)*(logN_A_1-logN_1),sigman,1) +
                dnorm(logR_A_2-logR_2,(alphar-tweak2)*(logR_A_1-logR_1),sigmar,1);
      // SECOND MOUSE (B):
      if (t_1 < 18) {
        rbc = exp(logR_B_1)+exp(logE_B_1);
        K = rbc*(1-exp(-exp(logM_B_1)/rbc));
        // unlawfulness penalty
        md = logM_B_2-log(Beta*K*exp(-(exp(logW_B_1)+exp(logN_B_1))/rbc));
        ed = logE_B_2-log(rbc*exp(-(exp(logM_B_1)+exp(logN_B_1))/rbc));
        loglik -= lambdaM*md*md + lambdaE*ed*ed;
        // AR(1) individual deviation
        if (t_1 < 1) {
          loglik += dnorm(logW_B_1-logW_1,0,sigmaw/sqrt(1-(alphaw-tweak2)),1) +
                    dnorm(logN_B_1-logN_1,0,sigman/sqrt(1-(alphan-tweak2)),1) +
                    dnorm(logR_B_1-logR_1,0,sigmar/sqrt(1-(alphar-tweak2)),1);
        }
        loglik += dnorm(logW_B_2-logW_2,(alphaw-tweak2)*(logW_B_1-logW_1),sigmaw,1) +
                  dnorm(logN_B_2-logN_2,(alphan-tweak2)*(logN_B_1-logN_1),sigman,1) +
                  dnorm(logR_B_2-logR_2,(alphar-tweak2)*(logR_B_1-logR_1),sigmar,1);
      }
      // THIRD MOUSE (C):
      rbc = exp(logR_C_1)+exp(logE_C_1);
      K = rbc*(1-exp(-exp(logM_C_1)/rbc));
      // unlawfulness penalty
      md = logM_C_2-log(Beta*K*exp(-(exp(logW_C_1)+exp(logN_C_1))/rbc));
      ed = logE_C_2-log(rbc*exp(-(exp(logM_C_1)+exp(logN_C_1))/rbc));
      loglik -= lambdaM*md*md + lambdaE*ed*ed;
      // AR(1) individual deviation
      if (t_1 < 1) {
        loglik += dnorm(logW_C_1-logW_1,0,sigmaw/sqrt(1-(alphaw-tweak2)),1) +
                  dnorm(logN_C_1-logN_1,0,sigman/sqrt(1-(alphan-tweak2)),1) +
                  dnorm(logR_C_1-logR_1,0,sigmar/sqrt(1-(alphar-tweak2)),1);
      }
      loglik += dnorm(logW_C_2-logW_2,(alphaw-tweak2)*(logW_C_1-logW_1),sigmaw,1) +
                dnorm(logN_C_2-logN_2,(alphan-tweak2)*(logN_C_1-logN_1),sigman,1) +
                dnorm(logR_C_2-logR_2,(alphar-tweak2)*(logR_C_1-logR_1),sigmar,1);
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
      if (t < 17) {
        rbc = exp(logR_A)+exp(logE_A);
        K = rbc*(1-exp(-exp(logM_A)/rbc));
        lik += (R_FINITE(Retic_A)) ? dnorm(log(Retic_A),logR_A,sigmaRetic+tweak,1) : 0;
        lik += (R_FINITE(RBC_A)) ? dnorm(log(RBC_A),log(rbc),sigmaRBC+tweak,1) : 0;
        lik += (R_FINITE(Pd_A)) ? dnorm(log(Pd_A),log(K),sigmaPd+tweak,1) : 0;
      }
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
      double tweak2 = 0.0001;
      // GAUSSIAN MARKOV RF MODEL:
      loglik = dnorm(logR_2,logR_1,sigmaR,1) +
               dnorm(logN_2,logN_1,sigmaN,1) +
               dnorm(logW_2,logW_1,sigmaW,1);
      // FIRST MOUSE (A):
      if (t_1 < 16) {
        rbc = exp(logR_A_1)+exp(logE_A_1);
        K = rbc*(1-exp(-exp(logM_A_1)/rbc));
        // unlawfulness penalty
        md = logM_A_2-log(Beta*K*exp(-(exp(logW_A_1)+exp(logN_A_1))/rbc));
        ed = logE_A_2-log(rbc*exp(-(exp(logM_A_1)+exp(logN_A_1))/rbc));
        loglik -= lambdaM*md*md + lambdaE*ed*ed;
        // AR(1) individual deviation
        if (t_1 < 1) {
          loglik += dnorm(logW_A_1-logW_1,0,sigmaw/sqrt(1-(alphaw-tweak2)),1) +
                    dnorm(logN_A_1-logN_1,0,sigman/sqrt(1-(alphan-tweak2)),1) +
                    dnorm(logR_A_1-logR_1,0,sigmar/sqrt(1-(alphar-tweak2)),1);
        }
        loglik += dnorm(logW_A_2-logW_2,(alphaw-tweak2)*(logW_A_1-logW_1),sigmaw,1) +
                  dnorm(logN_A_2-logN_2,(alphan-tweak2)*(logN_A_1-logN_1),sigman,1) +
                  dnorm(logR_A_2-logR_2,(alphar-tweak2)*(logR_A_1-logR_1),sigmar,1);
      }
      // SECOND MOUSE (B):
      rbc = exp(logR_B_1)+exp(logE_B_1);
      K = rbc*(1-exp(-exp(logM_B_1)/rbc));
      // unlawfulness penalty
      md = logM_B_2-log(Beta*K*exp(-(exp(logW_B_1)+exp(logN_B_1))/rbc));
      ed = logE_B_2-log(rbc*exp(-(exp(logM_B_1)+exp(logN_B_1))/rbc));
      loglik -= lambdaM*md*md + lambdaE*ed*ed;
      // AR(1) individual deviation
      if (t_1 < 1) {
        loglik += dnorm(logW_B_1-logW_1,0,sigmaw/sqrt(1-(alphaw-tweak2)),1) +
                  dnorm(logN_B_1-logN_1,0,sigman/sqrt(1-(alphan-tweak2)),1) +
                  dnorm(logR_B_1-logR_1,0,sigmar/sqrt(1-(alphar-tweak2)),1);
      }
      loglik += dnorm(logW_B_2-logW_2,(alphaw-tweak2)*(logW_B_1-logW_1),sigmaw,1) +
                dnorm(logN_B_2-logN_2,(alphan-tweak2)*(logN_B_1-logN_1),sigman,1) +
                dnorm(logR_B_2-logR_2,(alphar-tweak2)*(logR_B_1-logR_1),sigmar,1);
      // THIRD MOUSE (C):
      rbc = exp(logR_C_1)+exp(logE_C_1);
      K = rbc*(1-exp(-exp(logM_C_1)/rbc));
      // unlawfulness penalty
      md = logM_C_2-log(Beta*K*exp(-(exp(logW_C_1)+exp(logN_C_1))/rbc));
      ed = logE_C_2-log(rbc*exp(-(exp(logM_C_1)+exp(logN_C_1))/rbc));
      loglik -= lambdaM*md*md + lambdaE*ed*ed;
      // AR(1) individual deviation
      if (t_1 < 1) {
        loglik += dnorm(logW_C_1-logW_1,0,sigmaw/sqrt(1-(alphaw-tweak2)),1) +
                  dnorm(logN_C_1-logN_1,0,sigman/sqrt(1-(alphan-tweak2)),1) +
                  dnorm(logR_C_1-logR_1,0,sigmar/sqrt(1-(alphar-tweak2)),1);
      }
      loglik += dnorm(logW_C_2-logW_2,(alphaw-tweak2)*(logW_C_1-logW_1),sigmaw,1) +
                dnorm(logN_C_2-logN_2,(alphan-tweak2)*(logN_C_1-logN_1),sigman,1) +
                dnorm(logR_C_2-logR_2,(alphar-tweak2)*(logR_C_1-logR_1),sigmar,1);
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
#Pomp object for box 5 (control)
dat |>
  filter(box=="05") |>
  select(day,mouse,Retic,RBC) |>
  filter(day<=20) |>
  pivot_wider(values_from=c(Retic,RBC),names_from=mouse) |>
  pomp(
    times="day",t0=0,
    dmeasure=Csnippet(r"{
      double rbc;
      double tweak = 0.01;
      lik = 0;
      // first mouse
      rbc = exp(logR_A)+exp(logE_A);
      lik += (R_FINITE(Retic_A)) ? dnorm(log(Retic_A),logR_A,sigmaRetic+tweak,1) : 0;
      lik += (R_FINITE(RBC_A)) ? dnorm(log(RBC_A),log(rbc),sigmaRBC+tweak,1) : 0;
      // second mouse
      rbc = exp(logR_B)+exp(logE_B);
      lik += (R_FINITE(Retic_B)) ? dnorm(log(Retic_B),logR_B,sigmaRetic+tweak,1) : 0;
      lik += (R_FINITE(RBC_B)) ? dnorm(log(RBC_B),log(rbc),sigmaRBC+tweak,1) : 0;
      // third mouse
      rbc = exp(logR_C)+exp(logE_C);
      lik += (R_FINITE(Retic_C)) ? dnorm(log(Retic_C),logR_C,sigmaRetic+tweak,1) : 0;
      lik += (R_FINITE(RBC_C)) ? dnorm(log(RBC_C),log(rbc),sigmaRBC+tweak,1) : 0;
      if (!give_log) lik = exp(lik);
    }"
    ),
    dprocess=Csnippet(r"{
      double rbc, ed;
      double tweak2 = 0.0001;
      // GAUSSIAN MARKOV RF MODEL:
      loglik = dnorm(logR_2,logR_1,sigmaR,1) +
               dnorm(logN_2,logN_1,sigmaN,1);
      // FIRST MOUSE (A):
      rbc = exp(logR_A_1)+exp(logE_A_1);
      // unlawfulness penalty
      ed = logE_A_2-log(rbc*exp(-(exp(logN_A_1))/rbc));
      loglik -= lambdaE*ed*ed;
      // AR(1) individual deviation
      if (t_1 < 1) {
        loglik += dnorm(logN_A_1-logN_1,0,sigman/sqrt(1-(alphan-tweak2)),1) +
                  dnorm(logR_A_1-logR_1,0,sigmar/sqrt(1-(alphar-tweak2)),1);
      }
      loglik += dnorm(logN_A_2-logN_2,(alphan-tweak2)*(logN_A_1-logN_1),sigman,1) +
                dnorm(logR_A_2-logR_2,(alphar-tweak2)*(logR_A_1-logR_1),sigmar,1);
      // SECOND MOUSE (B):
      rbc = exp(logR_B_1)+exp(logE_B_1);
      // unlawfulness penalty
      ed = logE_B_2-log(rbc*exp(-(exp(logN_B_1))/rbc));
      loglik -= lambdaE*ed*ed;
      // AR(1) individual deviation
      if (t_1 < 1) {
        loglik += dnorm(logN_B_1-logN_1,0,sigman/sqrt(1-(alphan-tweak2)),1) +
                  dnorm(logR_B_1-logR_1,0,sigmar/sqrt(1-(alphar-tweak2)),1);
      }
      loglik += dnorm(logN_B_2-logN_2,(alphan-tweak2)*(logN_B_1-logN_1),sigman,1) +
                dnorm(logR_B_2-logR_2,(alphar-tweak2)*(logR_B_1-logR_1),sigmar,1);
      // THIRD MOUSE (C):
      rbc = exp(logR_C_1)+exp(logE_C_1);
      // unlawfulness penalty
      ed = logE_C_2-log(rbc*exp(-(exp(logN_C_1))/rbc));
      loglik -= lambdaE*ed*ed;
      // AR(1) individual deviation
      if (t_1 < 1) {
        loglik += dnorm(logN_C_1-logN_1,0,sigman/sqrt(1-(alphan-tweak2)),1) +
                  dnorm(logR_C_1-logR_1,0,sigmar/sqrt(1-(alphar-tweak2)),1);
      }
      loglik += dnorm(logN_C_2-logN_2,(alphan-tweak2)*(logN_C_1-logN_1),sigman,1) +
                dnorm(logR_C_2-logR_2,(alphar-tweak2)*(logR_C_1-logR_1),sigmar,1);
    }"
    ),
    partrans=parameter_trans(
      log=c(
        "sigmaRBC","sigmaRetic",
        "lambdaE",
        "sigmaR","sigmar",
        "sigmaN","sigman"
      ),
      logit=c(
        "alphar","alphan"
      )        
    ),
    params=c(
      lambdaE=1000, # unlawfulness penalty: chosen to be large
      sigmaR=0.75,sigmaN=0.75, # fixed at these values in example 9
      ## based on example 9 estimates above
      alphan=0.5,sigman=0.6,
      alphar=0.8,sigmar=0.25,
      sigmaRBC=0.1,sigmaRetic=0.1 # chosen to be bigger than needed before
    ),
    paramnames=c(
      "lambdaE",
      "sigmaR","sigmaN",
      "alphan","alphar",
      "sigman","sigmar",
      "sigmaRBC","sigmaRetic"
    ),
    statenames=c(
      "logE","logR","logN",
      outer(
        c("logE","logR","logN"),
        c("A","B","C"),
        paste,sep="_"
      )
    )
  ) -> po05

## -----------------------------------------------------------------------------
#Initial trajectories for box 1 (pABA = 0.05)
po_C <- dat |>
  filter(box=="01",mouse=="C") |>
  select(
    day,
    mouse,
    R=Retic,
    E=Eryth,
    Pd=Pd
  ) |>
  mutate(
    K=coalesce(pmin(Pd,R+E-1),50), # the choice of 50 is ad hoc
    M=-(R+E)*log(1-K/(R+E)),
    N=-(R+E)*log(lead(E)/(R+E))-M,
    W=-(R+E)*log(lead(M)/(6.51*K))-N,
    N=if_else(N>0,N,NA_real_),
    W=if_else(W>0,W,NA_real_)
  ) |>
  select(day,mouse,R,E,K,M,N,W) |>
  pivot_longer(-c(day,mouse)) |>
  group_by(name,mouse) |>
  mutate(
    value=gam(log(value)~s(day)) |>
      predict(newdata=pick(day)) |>
      exp() |>
      as.numeric()
  ) |>
  ungroup() |>
  pivot_wider() |>
  filter(day<=20)

po_A <- dat|>
  filter(box=="01",mouse=="A") |>
  select(
    day,
    mouse,
    R=Retic,
    E=Eryth,
    Pd=Pd
  ) |>
  mutate(
    K=coalesce(pmin(Pd,R+E-1),50), # the choice of 50 is ad hoc
    M=-(R+E)*log(1-K/(R+E)),
    N=-(R+E)*log(lead(E)/(R+E))-M,
    W=-(R+E)*log(lead(M)/(6.51*K))-N,
    N=if_else(N>0,N,NA_real_),
    W=if_else(W>0,W,NA_real_)
  ) |>
  select(day,mouse,R,E,K,M,N,W)

#Substitute in N and W values from mouse C to replace NAs from mouse A
po_A$N[is.na(po_A$N)] <- po_C$N[po_C$day%in%po_A$day[is.na(po_A$N)]]
po_A$W[is.na(po_A$W)] <- po_C$W[po_C$day%in%po_A$day[is.na(po_A$W)]]

po_df <- bind_rows(po_A,po_C) |>
  mutate(
    logM=log(M),
    logE=log(E),
    logR=log(R),
    logW=log(W),
    logN=log(N)
  ) |>
  select(day,mouse,logM,logE,logR,logW,logN)

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
  t() -> coefs01

## -----------------------------------------------------------------------------
#Initial trajectories for box 2 (pABA = 0.005)
po_df <- dat |>
  filter(box=="02",mouse!="C") |>
  select(
    day,
    mouse,
    R=Retic,
    E=Eryth,
    Pd=Pd
  ) |>
  mutate(
    K=coalesce(pmin(Pd,R+E-1),50), # the choice of 50 is ad hoc
    M=-(R+E)*log(1-K/(R+E)),
    N=-(R+E)*log(lead(E)/(R+E))-M,
    W=-(R+E)*log(lead(M)/(5.76*K))-N,
    N=if_else(N>0,N,NA_real_),
    W=if_else(W>0,W,NA_real_)
  ) |>
  select(day,mouse,R,E,M,N,W) |>
  pivot_longer(-c(day,mouse)) |>
  group_by(name,mouse) |>
  mutate(
    value=gam(log(value)~s(day)) |>
      predict(newdata=pick(day)) |>
      as.numeric()
  ) |>
  ungroup() |>
  pivot_wider() |>
  filter(day<=20)
names(po_df) <- c("day","mouse","logR","logE","logM","logN","logW")

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
po_df <- dat |>
  filter(box=="03") |>
  select(
    day,
    mouse,
    R=Retic,
    E=Eryth,
    Pd=Pd
  ) |>
  mutate(
    K=coalesce(pmin(Pd,R+E-1),50), # the choice of 50 is ad hoc
    M=-(R+E)*log(1-K/(R+E)),
    N=-(R+E)*log(lead(E)/(R+E))-M,
    W=-(R+E)*log(lead(M)/(5.76*K))-N,
    N=if_else(N>0,N,NA_real_),
    W=if_else(W>0,W,NA_real_)
  ) |>
  select(day,mouse,R,E,M,N,W) |>
  pivot_longer(-c(day,mouse)) |>
  group_by(name,mouse) |>
  mutate(
    value=gam(log(value)~s(day)) |>
      predict(newdata=pick(day)) |>
      as.numeric()
  ) |>
  ungroup() |>
  pivot_wider() |>
  filter(day<=20)
names(po_df) <- c("day","mouse","logR","logE","logM","logN","logW")

po_df <- slice(po_df,-which(po_df$mouse=="B"&po_df$day%in%c(19,20)))

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
po_df <- dat |>
  filter(box=="04") |>
  select(
    day,
    mouse,
    R=Retic,
    E=Eryth,
    Pd=Pd
  ) |>
  mutate(
    K=coalesce(pmin(Pd,R+E-1),50), # the choice of 50 is ad hoc
    M=-(R+E)*log(1-K/(R+E)),
    N=-(R+E)*log(lead(E)/(R+E))-M,
    W=-(R+E)*log(lead(M)/(5.76*K))-N,
    N=if_else(N>0,N,NA_real_),
    W=if_else(W>0,W,NA_real_)
  ) |>
  select(day,mouse,R,E,M,N,W) |>
  pivot_longer(-c(day,mouse)) |>
  group_by(name,mouse) |>
  mutate(
    value=gam(log(value)~s(day)) |>
      predict(newdata=pick(day)) |>
      as.numeric()
  ) |>
  ungroup() |>
  pivot_wider() |>
  filter(day<=20)
names(po_df) <- c("day","mouse","logR","logE","logM","logN","logW")

po_df <- slice(po_df,-which(po_df$mouse=="A"&po_df$day%in%c(17:20)))

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
#Initial trajectories for box 5 (control)
po_df <- dat |>
  filter(box=="05") |>
  select(day,mouse,Retic,Eryth) |>
  pivot_longer(-c(day,mouse)) |>
  group_by(name,mouse) |>
  mutate(
    value=gam(log(value)~s(day)) |>
      predict(newdata=pick(day)) |>
      exp() |>
      as.numeric()
  ) |>
  ungroup() |>
  pivot_wider() |>
  filter(day<=20) |>
  mutate(logE=log(Eryth),logR=log(Retic)) |>
  select(day,mouse,logE,logR) |>
  mutate(logN=logR)

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
  t() -> coefs05

## -----------------------------------------------------------------------------

create_objfun <- function (
    object1, params1 = coef(object1), coefs1,
    object2, params2 = coef(object2), coefs2,
    object3, params3 = coef(object3), coefs3,
    object4, params4 = coef(object4), coefs4,
    object_control, params_control = coef(object_control), coefs_control,
    est,
    est_control,
    control = list(reltol = 1e-8, maxit = 1e6)
) {
  
  params1 <- partrans(object1,params1,dir="toEst")
  params2 <- partrans(object2,params2,dir="toEst")
  params3 <- partrans(object3,params3,dir="toEst")
  params4 <- partrans(object4,params4,dir="toEst")
  params_control <- partrans(object_control,params_control,dir="toEst")
  
  idx <- match(est,names(params1)) #vector of positions
  
  idx_control <- match(est_control,names(params_control))
  
  ofun <- function (par) {
    
    params1[idx] <- par
    params2[idx] <- par
    params3[idx] <- par
    params4[idx] <- par
    params_control[idx_control] <- par[match(est_control,names(par))]
    
    pomp::coef(object1,transform=TRUE) <<- params1
    pomp::coef(object2,transform=TRUE) <<- params2
    pomp::coef(object3,transform=TRUE) <<- params3
    pomp::coef(object4,transform=TRUE) <<- params4
    pomp::coef(object_control,transform=TRUE) <<- params_control
    
    #Box 01 (pABA = 0.05)
    coefs1val <- c(coefs1)
    coefs1idx <- which(!is.na(coefs1val))
    
    fit1 <- optim(
      fn = function (x,coefs_idx) {
        
        y<-array(data=NA,dim=length(coefs1))
        y[coefs_idx] <- x
        dim(y) <- dim(coefs1)
        dimnames(y) <- dimnames(coefs1)
        
        object1@states <<- y
        
        object1 |> dprocess(log=TRUE) |> sum() -> ss
        object1 |> dmeasure(log=TRUE) |> sum() -> ll
        
        -(ss+ll)
        
      },
      par=c(coefs1)[coefs1idx],
      coefs_idx=coefs1idx,
      method="BFGS",
      control=control
    )
    
    sol01 <- array(data=NA,dim=length(coefs1))
    sol01[coefs1idx] <- fit1$par
    dim(sol01) <- dim(coefs1)
    dimnames(sol01) <- dimnames(coefs1)
    
    #Box 02 (pABA = 0.005)
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
    
    sol02 <- fit2$par
    dimnames(sol02) <- dimnames(coefs2)
    
    #Box 03 (pABA = 0.0005)
    coefs3val <- c(coefs3)
    coefs3idx <- which(!is.na(coefs3val))
    
    fit3 <- optim(
      fn = function (x,coefs_idx) {
        
        y<-array(data=NA,dim=length(coefs3))
        y[coefs_idx] <- x
        dim(y) <- dim(coefs3)
        dimnames(y) <- dimnames(coefs3)
        
        object3@states <<- y
        
        object3 |> dprocess(log=TRUE) |> sum() -> ss
        object3 |> dmeasure(log=TRUE) |> sum() -> ll
        
        -(ss+ll)
        
      },
      par=c(coefs3)[coefs3idx],
      coefs_idx=coefs3idx,
      method="BFGS",
      control=control
    )
    
    sol03 <- array(data=NA,dim=length(coefs3))
    sol03[coefs3idx] <- fit3$par
    dim(sol03) <- dim(coefs3)
    dimnames(sol03) <- dimnames(coefs3)
    
    #Box 04 (pABA = 0)
    coefs4val <- c(coefs4)
    coefs4idx <- which(!is.na(coefs4val))
    
    fit4 <- optim(
      fn = function (x,coefs_idx) {
        
        y<-array(data=NA,dim=length(coefs4))
        y[coefs_idx] <- x
        dim(y) <- dim(coefs4)
        dimnames(y) <- dimnames(coefs4)
        
        object4@states <<- y
        
        object4 |> dprocess(log=TRUE) |> sum() -> ss
        object4 |> dmeasure(log=TRUE) |> sum() -> ll
        
        -(ss+ll)
        
      },
      par=c(coefs4)[coefs4idx],
      coefs_idx=coefs4idx,
      method="BFGS",
      control=control
    )
    
    sol04 <- array(data=NA,dim=length(coefs4))
    sol04[coefs4idx] <- fit4$par
    dim(sol04) <- dim(coefs4)
    dimnames(sol04) <- dimnames(coefs4)
    
    #Box 05 (control)
    fit_control <- optim(
      fn = function (x) {
        
        dim(x)<-dim(coefs_control)
        dimnames(x) <- dimnames(coefs_control)
        
        object_control@states <<- x
        
        object_control |> dprocess(log=TRUE) |> sum() -> ss
        object_control |> dmeasure(log=TRUE) |> sum() -> ll
        
        -(ss+ll)
        
      },
      par=coefs_control,
      method="BFGS",
      control=control
    )
    
    sol_control <- fit_control$par
    dimnames(sol_control) <- dimnames(coefs_control)
    
    coefs1 <<- sol01
    coefs2 <<- sol02
    coefs3 <<- sol03
    coefs4 <<- sol04
    coefs_control <<- sol_control
    
    fit1$value+fit2$value+fit3$value+fit4$value+fit_control$value
  }
  environment(ofun) <- list2env(
    list(object1=object1,params1=params1,coefs1=coefs1,
         object2=object2,params2=params2,coefs2=coefs2,
         object3=object3,params3=params3,coefs3=coefs3,
         object4=object4,params4=params4,coefs4=coefs4,
         object_control=object_control,params_control=params_control,coefs_control=coefs_control,
         idx=idx,idx_control=idx_control,
         control=control,est_control=est_control),
    parent=parent.frame(2)
  )
  ofun
}

## -----------------------------------------------------------------------------

stew(
  file="FDA_FullModel.rda",
  info=TRUE,
  {
    create_objfun(
      object1=po01,
      object2=po02,
      object3=po03,
      object4=po04,
      object_control=po05,
      coefs1=coefs01,
      coefs2=coefs02,
      coefs3=coefs03,
      coefs4=coefs04,
      coefs_control=coefs05,
      est=c("alphaw","alphan","alphar",
            "sigmaw","sigman","sigmar",
            "sigmaRBC","sigmaRetic","sigmaPd"),
      est_control=c("alphar","alphan",
                    "sigmar","sigman",
                    "sigmaRBC","sigmaRetic"),
    ) -> ofun
    try({
      optim(
        fn=ofun,
        par=coef(po02,c("alphaw","alphan","alphar",
                        "sigmaw","sigman","sigmar",
                        "sigmaRBC","sigmaRetic","sigmaPd"),transform=TRUE),
        control=list(maxit=3000,reltol=1e-4,trace=1)
      ) -> fit
      ofun(fit$par)
    })
  })

evalq(coef(object1),envir=environment(ofun)) |> melt() -> est.coefs01
est.coefs01$box = "01"
evalq(coef(object2),envir=environment(ofun)) |> melt() -> est.coefs02
est.coefs02$box = "02"
evalq(coef(object3),envir=environment(ofun)) |> melt() -> est.coefs03
est.coefs03$box = "03"
evalq(coef(object4),envir=environment(ofun)) |> melt() -> est.coefs04
est.coefs04$box = "04"
evalq(coef(object_control),envir=environment(ofun)) |> melt() -> est.coefs05
est.coefs05$box = "05"

est.coefs<-bind_rows(est.coefs01,est.coefs02,est.coefs03,est.coefs04,est.coefs05)

write.csv(est.coefs,"parameter_results_FullModel.csv",row.names=FALSE)

evalq(coefs1,envir=environment(ofun)) |>
  melt() |>
  as_tibble() |>
  separate(Var1,into=c("name","mouse")) |>
  mutate(
    time=time(po02)[Var2],
    mouse=if_else(is.na(mouse),"Group",mouse),
    value=exp(value),
    box="01"
  ) |>
  select(-Var2) -> est.po01

evalq(coefs2,envir=environment(ofun)) |>
  melt() |>
  as_tibble() |>
  separate(Var1,into=c("name","mouse")) |>
  mutate(
    time=time(po02)[Var2],
    mouse=if_else(is.na(mouse),"Group",mouse),
    value=exp(value),
    box="02"
  ) |>
  select(-Var2) -> est.po02

evalq(coefs3,envir=environment(ofun)) |>
  melt() |>
  as_tibble() |>
  separate(Var1,into=c("name","mouse")) |>
  mutate(
    time=time(po02)[Var2],
    mouse=if_else(is.na(mouse),"Group",mouse),
    value=exp(value),
    box="03"
  ) |>
  select(-Var2) -> est.po03

evalq(coefs4,envir=environment(ofun)) |>
  melt() |>
  as_tibble() |>
  separate(Var1,into=c("name","mouse")) |>
  mutate(
    time=time(po02)[Var2],
    mouse=if_else(is.na(mouse),"Group",mouse),
    value=exp(value),
    box="04"
  ) |>
  select(-Var2) -> est.po04

evalq(coefs_control,envir=environment(ofun)) |>
  melt() |>
  as_tibble() |>
  separate(Var1,into=c("name","mouse")) |>
  mutate(
    time=time(po02)[Var2],
    mouse=if_else(is.na(mouse),"Group",mouse),
    value=exp(value),
    box="05"
  ) |>
  select(-Var2) -> est.po05

est.po<-bind_rows(est.po01,est.po02,est.po03,est.po04,est.po05)
write.csv(est.po,"trajectory_results_FullModel.csv",row.names=FALSE)
