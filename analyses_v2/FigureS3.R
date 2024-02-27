library(tidyverse)
library(gridExtra)

###################
#### 1-day lag ####
###################
df1 <- readRDS("results_df_RBC_lag1.rds")

results_bar <- df1 |> 
  select(r,b,model) |> 
  unique() |>
  unite("model_b",c(model,b),sep=", ",remove=FALSE)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

table(results_bar$b)/1000
table(results_bar$model_b) |> 
  as.data.frame() |>
  arrange(desc(Freq))

xaxis_labs <- c("Model B,\nbreakpoint 9",
                "Model C,\nbreakpoint 9",
                "Model E,\nbreakpoint 10",
                "Model E,\nbreakpoint 9",
                "Model A,\nbreakpoint 9",
                "Model F,\nbreakpoint 10",
                "Model B,\nbreakpoint 10",
                "Model D,\nbreakpoint 10",
                "Model D,\nbreakpoint 9",
                "Model F,\nbreakpoint 9",
                "Linear,\nuniphasic")

freq_df <- c(sort(table(results_bar$model_b)/1000,decreasing=TRUE),0) |> 
  as.data.frame() |> 
  rownames_to_column()
names(freq_df) <- c("model","freq")
freq_df$labs <- xaxis_labs

(bar_plot1 <- freq_df |>
    ggplot()+
    geom_bar(aes(x=reorder(labs, -freq),y=freq),stat="identity")+
    xlab("Model, breakpoint")+ylab("Frequency model selected")+
    scale_x_discrete(labels=xaxis_labs,breaks=xaxis_labs)+
    
    annotate("text",x=11,y=0,label="X",size=5)+
    
    annotate("text",x=2,y=0.8,label=expression(paste("Model B: ",
                                                     R[t]%~%RBC[t-1]^{2}:phase+RBC[t-1]:phase+pABA:phase+phase)),size=5,hjust=0,parse=T)+ 
    annotate("text",x=2,y=0.6,label=expression(paste("Model C: ",
                                                     R[t]%~%RBC[t-1]^{2}:phase+RBC[t-1]:phase+phase)),size=5,hjust=0,parse=T)+ 
    annotate("text",x=2,y=0.4,label=expression(paste("Model E: ",
                                                     R[t]%~%RBC[t-1]:phase+pABA:phase+phase)),size=5,hjust=0,parse=T)+
    #annotate("text",x=2,y=0.2,label=expression(paste("Model A: ",
    #R[t]%~%RBC[t-2]^{2}:pABA:phase+RBC[t-2]:pABA:phase+pABA:phase)),size=5,hjust=0,parse=T)+ 
    ggtitle("1-day reticulocyte response lag")+
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

###################
#### 2-day lag ####
###################
df2 <- readRDS("results_df_RBC_lag2.rds")

results_bar <- df2 |> 
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

(bar_plot2 <- freq_df |>
    ggplot()+
    geom_bar(aes(x=reorder(labs, -freq),y=freq),stat="identity")+
    xlab("Model, breakpoint")+ylab("Frequency model selected")+
    scale_x_discrete(labels=xaxis_labs,breaks=xaxis_labs)+
    
    annotate("text",x=11,y=0,label="X",size=5)+
    
    annotate("text",x=3,y=0.4,label=expression(paste("Model B: ",
                                                     R[t]%~%RBC[t-2]^{2}:phase+RBC[t-2]:phase+pABA:phase+phase)),size=5,hjust=0,parse=T)+ 
    annotate("text",x=3,y=0.3,label=expression(paste("Model E: ",
                                                     R[t]%~%RBC[t-2]:phase+pABA:phase+phase)),size=5,hjust=0,parse=T)+
    annotate("text",x=3,y=0.2,label=expression(paste("Model A: ",
                                                     R[t]%~%RBC[t-2]^{2}:pABA:phase+RBC[t-2]:pABA:phase+pABA:phase)),size=5,hjust=0,parse=T)+ 
    ggtitle("2-day reticulocyte response lag")+
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

###################
#### 3-day lag ####
###################
df3 <- readRDS("results_df_RBC_lag3.rds")

results_bar <- df3 |> 
  select(r,b,model) |> 
  unique() |>
  unite("model_b",c(model,b),sep=", ",remove=FALSE)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

table(results_bar$b)/1000
table(results_bar$model_b) |> 
  as.data.frame() |>
  arrange(desc(Freq))

xaxis_labs <- c("Model B,\nbreakpoint 10",
                "Model E,\nbreakpoint 10",
                "Model B,\nbreakpoint 9",
                "Model E,\nbreakpoint 9",
                "Model B,\nbreakpoint 11",
                "Model F,\nbreakpoint 9",
                "Model F,\nbreakpoint 10",
                "Model A,\nbreakpoint 11",
                "Model C,\nbreakpoint 9",
                "Model C,\nbreakpoint 10",
                "Model A,\nbreakpoint 10",
                "Model D,\nbreakpoint 9",
                "Model E,\nbreakpoint 11",
                "Linear,\nuniphasic")

freq_df <- c(sort(table(results_bar$model_b)/1000,decreasing=TRUE),0) |> 
  as.data.frame() |> 
  rownames_to_column()
names(freq_df) <- c("model","freq")
freq_df$labs <- xaxis_labs

(bar_plot3 <- freq_df |>
    ggplot()+
    geom_bar(aes(x=reorder(labs, -freq),y=freq),stat="identity")+
    xlab("Model, breakpoint")+ylab("Frequency model selected")+
    scale_x_discrete(labels=xaxis_labs,breaks=xaxis_labs)+
    #annotate("text",x=5,y=0.55,label="Model F: R~(lagged E):phase+phase",size=5,hjust=0,fontface=3)+
    
    annotate("text",x=14,y=0,label="X",size=5)+
    
    annotate("text",x=2,y=0.45,label=expression(paste("Model B: ",
                                                      R[t]%~%RBC[t-3]^{2}:phase+RBC[t-3]:phase+pABA:phase+phase)),size=5,hjust=0,parse=T)+ 
    annotate("text",x=2,y=0.35,label=expression(paste("Model E: ",
                                                      R[t]%~%RBC[t-3]:phase+pABA:phase+phase)),size=5,hjust=0,parse=T)+
    ggtitle("3-day reticulocyte response lag")+
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

bar_panel <- arrangeGrob(bar_plot1,bar_plot2,bar_plot3,nrow=3)

# Add labels to the arranged plots
p <- as_ggplot(bar_panel) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0, 0), y = c(1, 0.66, 0.33))
p
ggsave("FigureS3.jpeg",width=20,height=35,units="cm")
