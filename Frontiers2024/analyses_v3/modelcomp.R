library(tidyverse)

read_csv("results_regression_stats.csv") |>
  rename(b1=`01`,b2=`02`,b3=`03`,b4=`04`) -> dat

stopifnot(
  "AICs not accurate"=dat |>
    filter(
      round(AIC_total+2*loglik_total-2*p_total-2,8)!=0 |
        round(AIC_sub+2*loglik_sub-2*p_sub-2,8)!=0 |
        p_total != coef_total |
        p_sub != coef_sub
    ) |>
    nrow()==0,
  "n_sub incorrect"=dat |> select(n_sub) |> distinct() |> nrow()==1,
  "n_total incorrect"=dat |> reframe(n_total+lag) |> distinct() |> nrow()==1,
  "breakpoints inconsistent"=dat |>
    select(model,lag,rep,b1,b2,b3,b4) |>
    filter(b1 != b2 | b2 != b3 | b3 != b4) |> nrow()==0,
  "inconsistent parameter count (1)"=dat |> filter(p_total!=p_sub) |> nrow()==0,
  "inconsistent parameter count (2)"=dat |> filter(coef_total!=p_total) |> nrow()==0,
  "inconsistent parameter count (3)"=dat |> filter(coef_sub!=p_sub) |> nrow()==0
)

dat |>
  select(-b2,-b3,-b4) |>
  rename(bp=b1) |>
  mutate(
    bp=coalesce(as.character(bp),"NA"),
    bp=ordered(bp,levels=c("NA",8,9,10,11)),
    lag=ordered(lag),
    ##    AIC=-2*loglik_total*n_sub/n_total+2*p_total,
    AIC=AIC_sub,
    AICc=AIC+2*p_sub*(p_sub+1)/(n_sub-p_sub-1)
  ) -> dat

stopifnot(
  "wrong number of models"=dat |> select(model,bp,lag) |> distinct() |> nrow()==120,
  "wrong number of replicates"=dat |> select(rep) |> distinct() |> nrow()==1000,
  "wrong reps per model"=dat |> count(model,bp,lag) |> filter(n!=1000) |> nrow()==0,
  "wrong models per rep"=dat |> count(rep) |> filter(n!=120) |> nrow()==0
)

dat |>
  select(model,bp,lag,rep,AICc) |>
  ##  select(model,bp,lag,rep,AIC) |>
  group_by(rep) |>
  filter(AICc==min(AICc)) |>
  ##  filter(AIC==min(AIC)) |>
  ungroup() |>
  count(model,bp,lag) |>
  mutate(freq=n/sum(n)) -> freqs

freqs |>
  ggplot(aes(x=lag,y=freq))+
  geom_bar(stat="identity",fill="darkblue")+
  facet_grid(bp~model,labeller=label_both)+
  theme_bw()

freqs |>
  group_by(model,lag) |>
  summarize(freq=sum(freq)) |>
  ungroup() |>
  ggplot(aes(x=model,y=freq,fill=lag))+
  geom_bar(stat="identity",position="stack")+
  theme_bw()

freqs |>
  group_by(model,bp) |>
  summarize(freq=sum(freq)) |>
  ungroup() |>
  ggplot(aes(x=model,y=freq,fill=bp))+
  geom_bar(stat="identity",position="stack")+
  theme_bw()

freqs |>
  filter(
    model %in% c("m2","m5"),
    bp!="NA"
  ) |>
  ggplot(aes(x=lag,y=freq,fill=bp))+
  geom_bar(stat="identity")+
  facet_grid(.~model,labeller=label_both)+
  theme_bw()
