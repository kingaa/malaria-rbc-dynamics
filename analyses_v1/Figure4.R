#SET DIRECTORY TO SOURCE FILE LOCATION

#Load required packages
library(tidyverse)
library(pomp)
library(foreach)
library(iterators)
library(doFuture)
library(aakmisc)

#Set up for parallelisation 
plan(multisession) #for faster results, use multicore outside of RStudio

#Set seed (important because random sampling occurs below)
seed_choice <- 851657743
set.seed(seed_choice)

##Load in smooth distribution samples from PNAS work
sm1name <- "m5sm1_mod.rds"
sm1 <- readRDS(sm1name)

##Create sm1_mod tibble, with columns rep, mouse, mousid, box, time, E, R, lik and key
##From sm1_mod, we will be sampling individual distribution samples
sm1 |>
  as_tibble() |>
  filter(time<=20,mouseid!="01-02",mouseid!="02-03") |> #remove underdosed mice
  pivot_wider(names_from=variable,values_from=value) |>
  group_by(mouseid) |>
  mutate(
    RBC = R+E,
    lik=exp(loglik-max(loglik))
  ) |>
  ungroup() |>
  filter(box!="05") |> #remove control mice
  select(rep,mouse,mouseid,box,time,E,R,RBC,lik) |>
  unite("key",c(rep,box,mouse),sep="_",remove=FALSE) -> sm1_mod

##Create data frame with breakpoints for each of the four pABA boxes
box_list <- unique(sm1_mod$box)
breakpoint_grid <- as_tibble(
  matrix(
    data=rep(c(NA,8,9,10,11),each=4),
    ncol=4,byrow=TRUE,
    dimnames=list(NULL,box_list)
  )
)

#Loop over reps
rep_num <- 1000

bake(file="results_df_RBC_lag2.rds",{
  foreach (
    r=1:rep_num,
    .combine=rbind,
    .options.future = list(seed = seed_choice)
  ) %dofuture% {
    
    mouse_id_list <- unique(sm1_mod$mouseid)
    
    joint_mouse_df <- data.frame()
    for (id in mouse_id_list){
      
      key_opts <- sm1_mod |> filter(mouseid==id) |> select(key,lik) |> unique()
      
      key_choice <- sample(key_opts$key,size=1,prob=key_opts$lik)
      
      sm1_mod_mouse <- sm1_mod |> 
        filter(key==key_choice) |>
        mutate(lagRBC=lag(RBC,2)) |>
        na.omit()
      
      joint_mouse_df <- rbind(joint_mouse_df,sm1_mod_mouse)
      
    }
    
    joint_AIC_df <- data.frame()
    for (i in seq(1,nrow(breakpoint_grid),1)){
      
      bp <- breakpoint_grid[i,]
      
      if (is.na(bp[1])){
        
        joint_mouse_df -> df
        
        ##Models
        models <- list(
          m1=lm(R~poly(lagRBC,2,raw=F):(box-1)+(box-1),data=df),
          m2=lm(R~(poly(lagRBC,2,raw=F)+box-1),data=df),
          m3=lm(R~poly(lagRBC,2,raw=F),data=df),
          m4=lm(R~lagRBC:(box-1)+(box-1),data=df),
          m5=lm(R~(lagRBC+box-1),data=df),
          m6=lm(R~lagRBC,data=df)
        )
        
        models |>
          lapply(\(m) tibble(AIC=AIC(m),BIC=BIC(m),bp)) |>
          bind_rows(.id="model") -> model_df
        
      } else {
        
        joint_mouse_df |>
          mutate(
            phase=as.factor(if_else(time<=unlist(bp)[box],1,2))
          ) -> df
        
        ##Models
        models <- list(
          m1=lm(R~poly(lagRBC,2,raw=F):(box-1):(phase-1)+(box-1):(phase-1),data=df),
          m2=lm(R~(poly(lagRBC,2,raw=F)+box-1):(phase-1)+(phase-1),data=df),
          m3=lm(R~poly(lagRBC,2,raw=F):(phase-1)+(phase-1),data=df),
          m4=lm(R~lagRBC:(box-1):(phase-1)+(box-1):(phase-1),data=df),
          m5=lm(R~(lagRBC+box-1):(phase-1)+(phase-1),data=df),
          m6=lm(R~lagRBC:(phase-1)+(phase-1),data=df)
        )
        
        models |>
          lapply(\(m) tibble(AIC=AIC(m),BIC=BIC(m),bp)) |>
          bind_rows(.id="model") -> model_df
        
      }
      
      joint_AIC_df <- rbind(joint_AIC_df,model_df)
    }
    
    ##Extract breakpoints from best model (lowest AIC)
    joint_AIC_df |>
      filter(AIC==min(AIC)) -> best
    
    ##Obtain data frame with above breakpoints specified
    joint_mouse_df |>
      mutate(
        phase=as.factor(if_else(time<=as.integer(unlist(best)[box]),1,2))
      ) -> df
    
    
    ##Extract best model based on lowest AIC
    if (is.na(best$`01`)){
      models <- list(
        m1=lm(R~poly(lagRBC,2,raw=F):(box-1)+(box-1),data=df),
        m2=lm(R~(poly(lagRBC,2,raw=F)+box-1),data=df),
        m3=lm(R~poly(lagRBC,2,raw=F),data=df),
        m4=lm(R~lagRBC:(box-1)+(box-1),data=df),
        m5=lm(R~(lagRBC+box-1),data=df),
        m6=lm(R~lagRBC,data=df)
      )
    } else {
      models <- list(
        m1=lm(R~poly(lagRBC,2,raw=F):(box-1):(phase-1)+(box-1):(phase-1),data=df),
        m2=lm(R~(poly(lagRBC,2,raw=F)+box-1):(phase-1)+(phase-1),data=df),
        m3=lm(R~poly(lagRBC,2,raw=F):(phase-1)+(phase-1),data=df),
        m4=lm(R~lagRBC:(box-1):(phase-1)+(box-1):(phase-1),data=df),
        m5=lm(R~(lagRBC+box-1):(phase-1)+(phase-1),data=df),
        m6=lm(R~lagRBC:(phase-1)+(phase-1),data=df)
      )
    }
    
    chosen_model <- models[[best$model]]
    coefs <- chosen_model$coefficients
    
    pred <- model.matrix(chosen_model) %*% coefs |> as.data.frame() |> select(pred=V1)
    df <- cbind(df,pred)
    df$r <- r
    df$b <- if (is.na(best$`01`)){"None"} else {best$`01`}
    df$model <- best$model
    
    df
    
  } -> results_df
  
  results_df
  
}) -> results_df_RBC

results_bar <- results_df_RBC |> 
  select(r,b,model) |> 
  unique() |>
  unite("model_b",c(model,b),sep=", ",remove=FALSE)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

table(results_bar$b)/1000
table(results_bar$model_b) |> 
  as.data.frame() |>
  arrange(desc(Freq))

xaxis_labs <- c("Model B,\nbreakpoint 9",
                "Model B,\nbreakpoint 10",
                "Model E,\nbreakpoint 10",
                "Model A,\nbreakpoint 10",
                "Model C,\nbreakpoint 10",
                "Model C,\nbreakpoint 9",
                "Model F,\nbreakpoint 10",
                "Model D,\nbreakpoint 10",
                "Model E,\nbreakpoint 9",
                "Model A,\nbreakpoint 9",
                "Linear,\nuniphasic")

freq_df <- c(sort(table(results_bar$model_b)/1000,decreasing=TRUE),0) |> 
  as.data.frame() |> 
  rownames_to_column()
names(freq_df) <- c("model","freq")
freq_df$labs <- xaxis_labs

(bar_plot <- freq_df |>
    ggplot()+
    geom_bar(aes(x=reorder(labs, -freq),y=freq),stat="identity")+
    xlab("Model, breakpoint")+ylab("Frequency model selected")+
    scale_x_discrete(labels=xaxis_labs,breaks=xaxis_labs)+
    
    annotate("text",x=11,y=0,label="X",size=5)+
    
    annotate("text",x=3,y=0.4,label=expression(paste("Model B: ",
                                                      R[t]%~%RBC[t-2]^{2}:phase+RBC[t-2]:phase+pABA:phase+phase)),size=5,hjust=0,parse=T)+ 
    #annotate("text",x=2,y=0.35,label=expression(paste("Model C: ",
                                                      #R[t]%~%E[t-3]^{2}:phase+E[t-3]:phase+phase)),size=5,hjust=0,parse=T)+ 
    annotate("text",x=3,y=0.3,label=expression(paste("Model E: ",
                                                      R[t]%~%RBC[t-2]:phase+pABA:phase+phase)),size=5,hjust=0,parse=T)+
    annotate("text",x=3,y=0.2,label=expression(paste("Model A: ",
                                                      R[t]%~%RBC[t-2]^{2}:pABA:phase+RBC[t-2]:pABA:phase+pABA:phase)),size=5,hjust=0,parse=T)+ 
    theme_bw()+
    theme(
      axis.title=element_text(size=13),
      axis.text=element_text(size=11),
      axis.text.x=element_text(angle=45,vjust=0.7),
      axis.title.x=element_blank(),
      panel.grid=element_blank(),
      legend.position="none",
      legend.title=element_text(size=15),
      legend.text=element_text(size=13),
      strip.text=element_text(size=12),
      strip.background=element_blank()
      
    )
)


results_df_RBC$pABA <- factor(results_df_RBC$box,
                              levels=c("04","03","02","01"),
                              labels=c("Unsupplemented","Low","Medium","High"))


best_model_plot <- results_df_RBC |>
  unite("model_b",c(model,b),sep=", ",remove=FALSE) |>
  filter(model_b=="m2, 9") |>
  ggplot()+
  geom_line(aes(x=lagRBC,y=pred,group=interaction(r,phase,pABA),col=pABA),alpha=0.05)+
  scale_x_continuous(labels = aakmisc::scinot)+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4000000))+
  annotate("text",x=7000000,y=3200000,label="Phase 2",size=5)+
  annotate("text",x=6000000,y=1100000,label="Phase 1",size=5)+
  geom_segment(aes(x=6000000,xend=6000000,y=1000000,yend=800000))+
  geom_segment(aes(x=7000000,xend=6100000,y=3000000,yend=2600000))+
  xlab("")+ylab("Reticulocyte supply (t)")+
  scale_colour_manual(values=cbPalette[2:5])+
  geom_text(aes(x=7500000,y=3900000,label="Model B, breakpoint 9"))+
  theme_bw()+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  labs(colour="Parasite nutrient (pABA)")+
  theme(
    axis.title=element_blank(),
    axis.text=element_text(size=12),
    axis.text.x=element_text(colour="white"),
    panel.grid.minor=element_blank(),
    legend.position=c(0.225,0.25),
    legend.background=element_blank(),
    legend.title=element_text(size=10),
    legend.text=element_text(size=8),
    strip.text=element_text(size=12),
    strip.background=element_blank(),
    plot.title=element_text(size=14,hjust=0.5,face="bold"),
    panel.border = element_rect(linewidth=4)
    
  )

second_model_plot <- results_df_RBC |>
  unite("model_b",c(model,b),sep=", ",remove=FALSE) |>
  filter(model_b=="m2, 10") |>
  ggplot()+
  geom_line(aes(x=lagRBC,y=pred,group=interaction(r,phase,pABA),col=pABA),alpha=0.05)+
  xlab("")+ylab("")+
  scale_x_continuous(labels = aakmisc::scinot)+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4000000))+
  scale_colour_manual(values=cbPalette[2:5])+
  geom_text(aes(x=7500000,y=3900000,label="Model B, breakpoint 10"))+
  theme_bw()+
  theme(
    axis.title=element_blank(),
    axis.text=element_text(size=12,colour="white"),
    panel.grid.minor=element_blank(),
    legend.position="none",
    legend.title=element_text(size=10),
    legend.text=element_text(size=18),
    strip.text=element_text(size=12),
    strip.background=element_blank(),
    plot.title=element_text(size=14,hjust=0.5)
    
  )

third_model_plot <- results_df_RBC |>
  unite("model_b",c(model,b),sep=", ",remove=FALSE) |>
  filter(model_b=="m5, 10") |>
  ggplot()+
  geom_line(aes(x=lagRBC,y=pred,group=interaction(r,phase,box),col=pABA),alpha=0.15)+
  xlab("RBC density (t-1)")+ylab("Reticulocyte supply (t)")+
  scale_x_continuous(labels = aakmisc::scinot)+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4000000))+
  scale_colour_manual(values=cbPalette[c(2:5)])+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  geom_text(aes(x=7500000,y=3900000,label="Model E, breakpoint 10"))+
  labs(colour="Parasite nutrient (pABA)")+
  theme_bw()+
  theme(
    axis.title=element_blank(),
    axis.text=element_text(size=12),
    panel.grid.minor=element_blank(),
    legend.position="none",
    legend.background=element_blank(),
    legend.title=element_text(size=10),
    legend.text=element_text(size=8),
    strip.text=element_text(size=12),
    strip.background=element_blank(),
    plot.title=element_text(size=14,hjust=0.5)
    
  )

fourth_model_plot <- results_df_RBC |>
  unite("model_b",c(model,b),sep=", ",remove=FALSE) |>
  filter(model_b=="m1, 10") |>
  ggplot()+
  geom_line(aes(x=lagRBC,y=pred,group=interaction(r,phase,box),col=pABA),alpha=0.2)+
  xlab("Erythrocyte density (t-1)")+ylab("")+
  scale_x_continuous(labels = aakmisc::scinot)+
  scale_y_continuous(labels = aakmisc::scinot,limits=c(0,4000000))+
  scale_colour_manual(values=cbPalette[c(2:5)])+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  geom_text(aes(x=7500000,y=3900000,label="Model A, breakpoint 10"))+
  labs(colour="Parasite nutrient\n(pABA)")+
  theme_bw()+
  theme(
    axis.title=element_blank(),
    axis.text=element_text(size=12),
    axis.text.y=element_text(colour="white"),
    panel.grid.minor=element_blank(),
    legend.position="none",
    legend.background=element_blank(),
    legend.title=element_text(size=10),
    legend.text=element_text(size=8),
    strip.text=element_text(size=12),
    strip.background=element_blank(),
    plot.title=element_text(size=14,hjust=0.5)
    
  )

library("gridExtra")
library("cowplot")
library(ggpubr)
# Arrange plots using arrangeGrob
# returns a gtable (gt)

gt1 <- bar_plot
gt2 <- arrangeGrob(best_model_plot, second_model_plot, third_model_plot, fourth_model_plot,
                   ncol = 2, nrow = 2, 
                   layout_matrix = cbind(c(1,3), c(2,4)),
                   left="Reticulocyte supply (t)",
                   bottom="RBC density (t-2)")



gt <- arrangeGrob(gt1,gt2,
                  nrow=2,
                  heights=c(0.66,1))

# Add labels to the arranged plots
p <- as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C","D","E"), size = 15,
                  x = c(0, 0.05,0.55,0.05,0.55), y = c(1, 0.625, 0.625,0.33,0.33))

ggsave("Figure4.jpeg",width=25,height=30,units="cm")