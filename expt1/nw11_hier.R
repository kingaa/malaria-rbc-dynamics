## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## -----------------------------------------------------------------------------
library(tidyverse)
library(stringi)
library(pomp)


## -----------------------------------------------------------------------------
read_csv(
  "../data.csv",
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
read_csv("../po_df_Example9.csv") |>
  mutate(
    mouse=case_when(
      mouse=="01"~"A",
      mouse=="02"~"B",
      mouse=="03"~"C",
      TRUE~NA_character_
    )
  ) |>
  filter(box=="02") |>
  select(day,mouse,logR=R,logW=W,logN=N) |>
  pivot_longer(-c(day,mouse)) |>
  group_by(day,name) |>
  mutate(value=value-mean(value)) -> resids

resids |>
  ggplot(aes(x=day,y=value,group=mouse,color=mouse))+
  geom_line()+
  facet_wrap(~name,ncol=1)+
  theme_bw()

resids |>
  group_by(mouse,name) |>
  reframe(
    var=c("ar1","sigma"),
    value={
      arima(
        x=value,
        order=c(1,0,0),
        SSinit="Rossignol2011",
        include.mean=FALSE
      ) -> fit
      c(fit$coef,sqrt(fit$sigma2))
    }
  ) |>
  ungroup() |>
  group_by(name,var) |>
  summarize(value=mean(value)) |>
  pivot_wider(names_from=var)


## -----------------------------------------------------------------------------
read_csv("../coef_df_Example9.csv") |>
  mutate(
    mouse=case_when(
      mouse=="01"~"A",
      mouse=="02"~"B",
      mouse=="03"~"C",
      TRUE~NA_character_
    )
  ) |>
  filter(box=="02") |>
  select(mouse,lambdaM,lambdaE,sigmaR,sigmaW,sigmaN,sigmaRBC,sigmaRetic,sigmaPd)


## -----------------------------------------------------------------------------
dat |>
  filter(box=="02") |>
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
        c("A","B","C"),
        paste,sep="_"
      )
    )
  ) -> po


## -----------------------------------------------------------------------------
po_df <- read_csv(
  "../po_df_Example9.csv",
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
  filter(box=="02") |>
  select(day,mouse,logM=M,logE=E,logR=R,logW=W,logN=N)

po_df |> pivot_longer(-c(day,mouse)) |>
  mutate(value=exp(value)) |>
  group_by(day,name) |>
  summarise(mean=log(mean(value))) |>
  pivot_wider(names_from=name,values_from=mean) -> po_avg

po_df |> pivot_longer(-c(day,mouse)) |>
  pivot_wider(names_from=c(name,mouse),values_from=value) -> po_ind

full_join(po_avg,po_ind,by="day") |>
  ungroup() |>
  select(-day) |>
  t() -> coefs



## -----------------------------------------------------------------------------
create_objfun <- function (
    object, params = coef(object), est, coefs,
    control = list(reltol = 1e-8, maxit = 1e6)
) {
  params <- partrans(object,params,dir="toEst")
  idx <- match(est,names(params))
  ofun <- function (par) {
    params[idx] <- par
    coef(object,transform=TRUE) <<- params
    fit <- optim(
      fn = function (x) {
        dim(x) <- dim(coefs)
        dimnames(x) <- dimnames(coefs)
        object@states <<- x
        object |> dprocess(log=TRUE) |> sum() -> ss
        object |> dmeasure(log=TRUE) |> sum() -> ll
        -(ss+ll)
      },
      par=coefs,
      method="BFGS",
      control=control
    )
    coefs <<- fit$par
    fit$value
  }
  environment(ofun) <- list2env(
    list(object=object,params=params,coefs=coefs,
         idx=idx,control=control),
    parent=parent.frame(2)
  )
  ofun
}


## -----------------------------------------------------------------------------
po |>
  create_objfun(
    est=c("alphaw","alphan","alphar",
      "sigmaw","sigman","sigmar",
      "sigmaRBC","sigmaPd","sigmaRetic"),
    coefs=coefs,
    control=list(reltol=1e-4,maxit=1e6)
  ) -> ofun

tic <- Sys.time()
try(
  optim(
    fn=ofun,
    par=coef(po,c("alphaw","alphan","alphar",
      "sigmaw","sigman","sigmar",
      "sigmaRBC","sigmaPd","sigmaRetic"),transform=TRUE),
    control=list(maxit=200,reltol=1e-3,trace=0)
  )
) -> fit
toc <- Sys.time()
fit
toc-tic

save.image(file="nw11_hier.rda")

evalq(coef(object),envir=environment(ofun)) |> melt()

convergence <- fit$convergence
min.fit <- fit$value


## -----------------------------------------------------------------------------
## evalq(object,envir=environment(ofun)) -> po
## 
## coef(po) |>
##   as.data.frame() |>
##   t() |>
##   as_tibble() |>
##   mutate(
##     box="02",
##     paba="0.005",
##     lambdaM=1000,
##     lambdaE=1000,
##     sigmaR=0.75,
##     sigmaW=0.75,
##     sigmaN=0.75,
##     fit=min.fit,
##     con=convergence
##   ) -> coef_df
## 
## po |>
##   as_tibble() |>
##   mutate(
##     K=log((exp(R)+exp(E))*(1-exp(-exp(M)/(exp(R)+exp(E))))),
##     box="02",
##     paba="0.005",
##     lambdaM=1000,
##     lambdaE=1000,
##     sigmaR=0.75,
##     sigmaW=0.75,
##     sigmaN=0.75,
##     fit=min.fit,
##     con=convergence
##   ) -> po_df

