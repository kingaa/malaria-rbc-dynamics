#Read in PNAS data
read_csv("data.csv",
         col_types="iiinnnn"
) |>
  mutate(
    mouseid=sprintf("%02d-%02d",box,mouse),
    box=sprintf("%02d",box),
    mouse=sprintf("%02d",mouse),
    paba=as.factor(box),
    rbc_density=rbc_density/1000
  ) |>
  mutate(
    paba=recode(
      paba,
      "01"="0.05","02"="0.005","03"="0.0005","04"="0","05"="control"
    )
  ) |>
  select(
    day,
    Pd=ama_density,
    RBC=rbc_density,
    Ter119=ter119_density,
    CD71=cd71_density,
    mouseid,
    paba,box,mouse
  ) |>
  arrange(mouseid,day) |>
  mutate(
    paba=as.character(paba),
    Ter119=ifelse(Ter119==0,NA,Ter119),
    CD71=ifelse(CD71==0,NA,CD71)
  ) |>
  mutate(
    Eryth=(1-CD71/Ter119)*RBC,
    Retic=CD71/Ter119*RBC
  ) -> flow

flow$pABA <- factor(flow$box,levels=c("05","04","03","02","01"),
                    labels=c("Uninfected","Unsupplemented","Low","Medium","High"))


#Read in PNAS trajectories
sm1name <- "m5sm1.rds"
sm1 <- readRDS(sm1name)

#Create dataframe with weighted trajectories (grouped by mouseid first, then by box)
sm1 |>
  as_tibble() |>
  filter(mouseid!="01-02",mouseid!="02-03") |> #remove underdosed mice
  pivot_wider(names_from=variable,values_from=value) |>
  #left_join(bdf,by=c("mouseid")) |>
  separate_wider_delim(cols="mouseid",delim="-",names=c("box","mouse"),cols_remove=FALSE) |>
  group_by(mouseid) |>
  mutate(
    SM=exp(-M/(R+E)),
    SN=exp(-N/(R+E)),
    Qun=SM*(1-SN),
    lik=exp(loglik-max(loglik))
  ) |>
  ungroup() |>
  select(-loglik) |>
  gather(variable,value,-rep,-time,-mouse,-box,-lik,-mouseid) |>
  filter(is.finite(value)) |>
  group_by(box,time,variable) |>
  dplyr::reframe(
    value=pomp::wquant(x=value,probs=c(0.05,0.5,0.95),weights=lik),
    name=c("lo","med","hi")
  ) |>
  ungroup() |>
  pivot_wider() -> group_traj

group_traj$pABA <- factor(group_traj$box,levels=c("05","04","03","02","01"),
                          labels=c("Uninfected","Unsupplemented","Low","Medium","High"))

#Remove estimates for W for control mice
group_traj <- group_traj |> dplyr::slice(-which(group_traj$box=="05"&group_traj$variable=="W"))
