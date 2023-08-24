#Load required packages
library(tidyverse)
library(stringi)
library(doParallel)
library(pomp)
library(panelPomp)
library(ggpubr)
library(scales)
library(gridExtra)

#Set working directory
setwd("~/Documents/GitHub/bdd/nw11_hier/final_files/")

#Set colour palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Read in PNAS data
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

#Obtain estimates for beta and dose
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
  rename(Beta=val) %>%
  left_join(theta1,by=c("dosevar"="var")) %>%
  rename(dose=val) %>%
  select(-betavar,-dosevar) %>%
  arrange(box,mouse) %>%
  mutate(
    Beta=coalesce(Beta,0),
    dose=coalesce(dose,0)
  ) -> theta

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

#Create object pos with pomp object for each mouse  
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

#Obtain MLE for sigmas, initial values, betas and dose
pf4name <- "m5pf4.rds"
pf4 <- readRDS(pf4name)

pf4 %>%
  select(loglik,starts_with("sigma"),ends_with("_0")) %>%
  filter(loglik==max(loglik)) -> mle

mle %>% unlist() -> p
pos %>% panelPomp(shared=p) %>% coef() %>% 
  rbind() %>% as_tibble() -> theta
theta %>%
  select(starts_with("Beta"),starts_with("dose")) %>%
  gather(variable, value) %>%
  tidyr::extract(variable, into=c("variable","mouseid"), 
                 regex="([[:alnum:]]+)\\[(.+?)\\]") %>%
  spread(variable,value) -> bdf

#Load in PNAS trajectories
sm1name <- "m5sm1.rds"
sm1 <- readRDS(sm1name)

#Create dataframe with weighted trajectories (grouped by mouseid first, then by box)
sm1 |>
  as_tibble() |>
  filter(mouseid!="01-02",mouseid!="02-03") |> #remove underdosed mice
  pivot_wider(names_from=variable,values_from=value) |>
  left_join(bdf,by=c("mouseid")) |>
  separate_wider_delim(cols="mouseid",delim="-",names=c("box","mouse"),cols_remove=FALSE) |>
  group_by(mouseid) |>
  mutate(
    SM=exp(-M/(R+E)),
    SN=exp(-N/(R+E)),
    SW=exp(-W/(R+E)),
    Qps=(1-SM)*SW*SN,
    Qpn=N/(N+W)*(1-SM)*(1-SW*SN),
    Qpw=W/(N+W)*(1-SM)*(1-SW*SN),
    Qun=SM*(1-SN),
    Qus=SM*SN,
    lambda_r=Beta*R/M*Qps,
    lambda_e=Beta*E/M*Qps,
    lambda_n=Beta*(R+E)/M*Qpn,
    lambda_w=Beta*(R+E)/M*Qpw,
    lambda_u=Beta-lambda_r-lambda_e-lambda_n-lambda_w,
    lambda_u_w_ratio=lambda_u/lambda_w,
    lambda_u_n_ratio=lambda_u/lambda_n,
    lambda_u_wn_ratio=lambda_u/(lambda_w+lambda_n),
    loss=(R+E)*(1-Qus),
    rbc=E+R,
    varN=N/rbc,
    perRetic=R/rbc,
    lik=exp(loglik-max(loglik))
  ) |>
  ungroup() |>
  select(-loglik,-Beta,-dose) |>
  gather(variable,value,-rep,-time,-mouse,-box,-lik,-mouseid) |>
  filter(is.finite(value)) %>%
  group_by(box,time,variable) |>
  dplyr::reframe(
    value=pomp::wquant(x=value,probs=c(0.05,0.5,0.95),weights=lik),
    name=c("lo","med","hi")
  ) |>
  ungroup() |>
  pivot_wider() -> group_traj

group_traj$pABA <- factor(group_traj$box,levels=c("05","04","03","02","01"),
                          labels=c("Control","0%","0.05%","0.5%","5%"))

#Remove estimates for W for control mice
group_traj <- group_traj |> dplyr::slice(-which(group_traj$box=="05"&group_traj$variable=="W"))

###############################
#### Control verus pABA 0% ####
###############################
#Plot group-level trajectories for reticulocytes
group_traj |> filter(pABA%in%c("0%"),variable=="K") -> tmp
peak.day <- tmp$time[which(tmp$med==max(tmp$med))]

(RvsDay_UninfInf<-group_traj |>
  filter(time<=20,variable%in%c("R"),pABA%in%c("Control","0%")) |>
  mutate(variable=case_match(variable,
                             "R"~"Reticulocytes"),
         type=case_match(pABA,
                         "Control"~"Uninfected",
                         "0%"~"Infected")) |>
  ggplot()+
  geom_line(aes(x=time,y=med,color=type),linewidth=2)+
  geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=type),alpha=0.2)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+ 
  scale_colour_manual(values=cbPalette[c(2,1)],name=NULL)+
  scale_fill_manual(values=cbPalette[c(2,1)],name=NULL)+
  xlab("Day")+ylab("Reticulocytes (density per microlitre)")+
  theme_bw()+
  theme(strip.background=element_blank(),
        panel.grid = element_blank(),
        legend.position=c(0.2,0.8),
        axis.title=element_text(size=15),
        axis.text=element_text(size=13),
        legend.text=element_text(size=14)
  )
)
ggsave("ReticsvsDay_UninfectedInfected.png",plot=RvsDay_UninfInf,
       path="~/Documents/GitHub/bdd/nw11_hier/final_files/POMP_GroupLevel_figures/",
       width=15,height=12,units="cm",dpi=600)

(EvsDay_UninfInf<-group_traj |>
    filter(time<=20,variable%in%c("E"),pABA%in%c("Control","0%")) |>
    mutate(variable=case_match(variable,
                               "E"~"Erythrocytes"),
           type=case_match(pABA,
                           "Control"~"Uninfected",
                           "0%"~"Infected")) |>
    ggplot()+
    geom_line(aes(x=time,y=med,color=type),linewidth=2)+
    geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=type),alpha=0.2)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)))+ 
    scale_colour_manual(values=cbPalette[c(2,1)],name=NULL)+
    scale_fill_manual(values=cbPalette[c(2,1)],name=NULL)+
    xlab("Day")+ylab("Erythrocytes (density per microlitre)")+
    theme_bw()+
    theme(strip.background=element_blank(),
          panel.grid = element_blank(),
          legend.position=c(0.2,0.3),
          axis.title=element_text(size=15),
          axis.text=element_text(size=13),
          legend.text=element_text(size=14)
    )
)
ggsave("ErythsvsDay_UninfectedInfected.png",plot=EvsDay_UninfInf,
       path="~/Documents/GitHub/bdd/nw11_hier/final_files/POMP_GroupLevel_figures/",
       width=15,height=12,units="cm",dpi=600)

(KvsDay_UninfInf<-group_traj |>
    filter(time<=20,variable%in%c("K"),pABA%in%c("0%")) |>
    mutate(variable=case_match(variable,
                               "K"~"Parasites"),
           type=case_match(pABA,
                           "Control"~"Uninfected",
                           "0%"~"Infected")) |>
    ggplot()+
    geom_line(aes(x=time,y=med,color=type),linewidth=2)+
    geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=type),alpha=0.2)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)))+ 
    scale_colour_manual(values=cbPalette[c(2,1)],name=NULL)+
    scale_fill_manual(values=cbPalette[c(2,1)],name=NULL)+
    xlab("Day")+ylab("Parasites (density per microlitre)")+
    theme_bw()+
    theme(strip.background=element_blank(),
          panel.grid = element_blank(),
          legend.position="none",
          axis.title=element_text(size=15),
          axis.text=element_text(size=13),
          legend.text=element_text(size=14)
    )
)
ggsave("ParavsDay_UninfectedInfected.png",plot=KvsDay_UninfInf,
       path="~/Documents/GitHub/bdd/nw11_hier/final_files/POMP_GroupLevel_figures/",
       width=15,height=12,units="cm",dpi=600)

group_traj |> filter(pABA%in%c("0%"),variable=="K") -> tmp
peak.day <- tmp$time[which(tmp$med==max(tmp$med))]

Rtmp<-RvsDay_UninfInf+geom_vline(xintercept=peak.day,linetype="dashed",col="grey")+geom_vline(xintercept=10,linetype="dashed",col="grey")+ylab("Reticulocytes\n(density per microlitre)")+theme(legend.background = element_blank())
Etmp<-EvsDay_UninfInf+geom_vline(xintercept=peak.day,linetype="dashed",col="grey")+geom_vline(xintercept=10,linetype="dashed",col="grey")+theme(legend.position="none")+ylab("Erythrocytes\n(density per microlitre)")
Ktmp<-KvsDay_UninfInf+geom_vline(xintercept=peak.day,linetype="dashed",col="grey")+geom_vline(xintercept=10,linetype="dashed",col="grey")+ylab("Parasites\n(density per microlitre)")
REKvsDay_UninfInf<-grid.arrange(Rtmp,Etmp,Ktmp,nrow=3)
ggsave("REKvsDay_UninfectedInfected.png",plot=REKvsDay_UninfInf,
       path="~/Documents/GitHub/bdd/nw11_hier/final_files/POMP_GroupLevel_figures/",
       width=12,height=25,units="cm",dpi=600)


#Read in data of weighted individual trajectories
data_weighted <- read_csv("data_weighted.csv")
data_weighted$W[which(data_weighted$box=="05")]<-0
data_weighted$pABA <- factor(data_weighted$box,levels=c("05","04","03","02","01"),
                             labels=c("Control","0%","0.05%","0.5%","5%"))
data_weighted_path <- data_weighted |>
  filter(time<=20,pABA%in%c("Control","0%")) |>
  mutate(
         type=case_match(pABA,
                         "Control"~"Uninfected",
                         "0%"~"Infected"))

group_traj |>
  filter(time<=20,variable%in%c("R","E"),pABA%in%c("Control","0%")) |>
  mutate(variable=case_match(variable,
                             "R"~"Reticulocytes",
                             "E"~"Erythrocytes"),
         type=case_match(pABA,
                         "Control"~"Uninfected",
                         "0%"~"Infected")) |>
  select(-lo,-hi) |>
  pivot_wider(names_from="variable",values_from="med") -> group_traj_path



(ReticvsEryth_UninfInf<-group_traj_path |>
  ggplot()+
  geom_path(data=data_weighted_path,aes(x=E,y=R,group=interaction(rep,rep2),col=type),alpha=0.01)+
  geom_path(aes(x=Erythrocytes,y=Reticulocytes,col=type),linewidth=2)+
  geom_text(data=group_traj_path |> filter(time%in%c(0,6,10,20)),aes(x=Erythrocytes,y=Reticulocytes,label=time),col="black",size=5)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+ 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)))+ 
  scale_colour_manual(values=cbPalette[c(2,1)],name=NULL)+
  xlab("Erythrocytes (density per microlitre)")+ylab("Reticulocytes (density per microlitre)")+
  theme_bw()+
  theme(strip.background=element_blank(),
        panel.grid = element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=13),
        legend.text=element_text(size=14)
  )
)
ggsave("ReticvsEryth_UninfectedInfected_indiv.png",plot=ReticvsEryth_UninfInf,
       path="~/Documents/GitHub/bdd/nw11_hier/final_files/POMP_GroupLevel_figures/",
       width=15,height=12,units="cm",dpi=600)

#####################
#### Growth rate ####
#####################
#Plot group-level trajectories for reticulocytes
(RvsDay_pABA<-group_traj |>
    filter(time<=20,variable%in%c("R")) |>
    ggplot()+
    geom_line(aes(x=time,y=med,color=pABA),linewidth=2)+
    geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)))+ 
    scale_colour_manual(values=cbPalette)+
    scale_fill_manual(values=cbPalette)+
    xlab("Day")+ylab("Reticulocytes (density per microlitre)")+
    theme_bw()+
    theme(strip.background=element_blank(),
          panel.grid = element_blank(),
          legend.position=c(0.2,0.7),
          axis.title=element_text(size=15),
          axis.text=element_text(size=13),
          legend.text=element_text(size=14),
          legend.title=element_text(size=15)
    )
)
ggsave("ReticsvsDay_pABA.png",plot=RvsDay_pABA,
       path="~/Documents/GitHub/bdd/nw11_hier/final_files/POMP_GroupLevel_figures/",
       width=15,height=12,units="cm",dpi=600)

#Plot group-level trajectories for erythrocytes
(EvsDay_pABA<-group_traj |>
    filter(time<=20,variable%in%c("E")) |>
    ggplot()+
    geom_line(aes(x=time,y=med,color=pABA),linewidth=2)+
    geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)))+ 
    scale_colour_manual(values=cbPalette)+
    scale_fill_manual(values=cbPalette)+
    xlab("Day")+ylab("Erythrocytes (density per microlitre)")+
    theme_bw()+
    theme(strip.background=element_blank(),
          panel.grid = element_blank(),
          legend.position=c(0.2,0.3),
          axis.title=element_text(size=15),
          axis.text=element_text(size=13),
          legend.text=element_text(size=14),
          legend.title=element_text(size=15)
    )
)
ggsave("ErythsvsDay_pABA.png",plot=EvsDay_pABA,
       path="~/Documents/GitHub/bdd/nw11_hier/final_files/POMP_GroupLevel_figures/",
       width=15,height=12,units="cm",dpi=600)

#Plot group-level trajectories for parasites
(KvsDay_pABA<-group_traj |>
    filter(time<=20,variable%in%c("K"),pABA!="Control") |>
    ggplot()+
    geom_line(aes(x=time,y=med,color=pABA),linewidth=2)+
    geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                  labels=trans_format('log10',math_format(10^.x)))+ 
    scale_colour_manual(values=cbPalette[2:5])+
    scale_fill_manual(values=cbPalette[2:5])+
    xlab("Day")+ylab("Parasites (density per microlitre)")+
    theme_bw()+
    theme(strip.background=element_blank(),
          panel.grid = element_blank(),
          legend.position=c(0.2,0.3),
          axis.title=element_text(size=15),
          axis.text=element_text(size=13),
          legend.text=element_text(size=14),
          legend.title=element_text(size=15)
    )
)
ggsave("ParavsDay_pABA.png",plot=KvsDay_pABA,
       path="~/Documents/GitHub/bdd/nw11_hier/final_files/POMP_GroupLevel_figures/",
       width=15,height=12,units="cm",dpi=600)
