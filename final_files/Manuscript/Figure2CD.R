library(pomp)
library(panelPomp)
library(ggpubr)
library(scales)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

source("POMP_GroupLevel_DataPrep.R")

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

sm1_sub <- sm1 |>
  as_tibble() |>
  filter(mouseid!="01-02",mouseid!="02-03",time<=20) |> #remove underdosed mice
  pivot_wider(names_from=variable,values_from=value) |>
  separate_wider_delim(cols="mouseid",delim="-",names=c("box","mouse"),cols_remove=FALSE) |>
  filter(box!="05") |>
  group_by(mouseid) |>
  mutate(
    SM=exp(-M/(R+E)),
    SN=exp(-N/(R+E)),
    Qun=SM*(1-SN),
    lik=exp(loglik-max(loglik))
  ) |>
  ungroup()
sm1_sub$pABA <- factor(sm1_sub$box,levels=c("05","04","03","02","01"),
                       labels=c("Uninfected","Unsupplemented","Low","Medium","High"))

ann_text <- as.data.frame(
  matrix(
    c(
      c("04","03","02","01"),
      c("n=3","n=3","n=2","n=2")
    ),
    nrow=4,ncol=2,byrow=FALSE
  )
)
names(ann_text) <- c("box","lab")
ann_text$pABA <- factor(ann_text$box,levels=c("05","04","03","02","01"),
                       labels=c("Uninfected","Unsupplemented","Low","Medium","High"))



facet <- sm1_sub |>
  ggplot()+
  geom_line(aes(x=time,y=R,group=interaction(rep,mouse),col=pABA),alpha=0.025)+
  geom_text(data=ann_text,aes(x=17,y=10^6.9,label=lab),size=5)+
  facet_wrap(pABA~.)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(10^5,10^7))+
  scale_colour_manual(values=cbPalette[2:5])+
  theme_bw()+
  xlab("Day post-infection")+ylab("Reticulocyte supply (density per µL)")+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    plot.title=element_text(hjust=0.5,face="bold",size=15),
    panel.grid=element_blank(),
    legend.position="none",
    legend.title=element_text(size=15),
    legend.text=element_text(size=13),
    strip.background=element_blank(),
    strip.text=element_text(size=13)
  )

group <- group_traj |>
  filter(time<=20,pABA!="Uninfected",variable=="R") |>
  ggplot()+
  geom_line(aes(x=time,y=med,col=pABA),linewidth=2)+
  geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=pABA),alpha=0.2)+
  geom_text(aes(x=17,y=10^6.9,label="n=10"),size=5)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels=trans_format('log10',math_format(10^.x)),
                limits=c(10^5,10^7))+
  scale_colour_manual(values=cbPalette[2:5])+
  scale_fill_manual(values=cbPalette[2:5])+
  theme_bw()+
  xlab("Day post-infection")+ylab("Reticulocyte supply (density per µL)")+
  labs(colour="Parasite nutrient (pABA)",fill="Parasite nutrient (pABA)")+
  theme(
    axis.title=element_text(size=15),
    axis.text=element_text(size=13),
    plot.title=element_text(hjust=0.5,face="bold",size=15),
    panel.grid=element_blank(),
    legend.position=c(0.3,0.7),
    legend.title=element_text(size=10),
    legend.text=element_text(size=8),
    strip.background=element_blank(),
    strip.text=element_text(size=13),
    legend.background=element_blank()
  )

joint_plot <- ggarrange(facet, group, nrow = 1, labels = c("A","B"))

ggsave("Figure3.png",plot=joint_plot,
       width=25,height=12.5,units="cm")
