library(tidyverse)

df50<-readRDS("jAICdf50.rds")
df100<-readRDS("jAICdf100.rds")
df500<-readRDS("jAICdf500.rds")
df1000<-readRDS("jAICdf1000.rds")
df2000<-readRDS("jAICdf2000.rds")

df50[which(df50$AIC==min(df50$AIC)),]
df100[which(df100$AIC==min(df100$AIC)),]
df500[which(df500$AIC==min(df500$AIC)),]
df1000[which(df1000$AIC==min(df1000$AIC)),]
df2000[which(df2000$AIC==min(df2000$AIC)),]

df50ord<-df50[order(df50$AIC,decreasing=FALSE),]
df100ord<-df100[order(df100$AIC,decreasing=FALSE),]
df500ord<-df500[order(df500$AIC,decreasing=FALSE),]
df1000ord<-df1000[order(df1000$AIC,decreasing=FALSE),]
df2000ord<-df2000[order(df2000$AIC,decreasing=FALSE),]

diff50<-df50ord$AIC[1]-df50ord$AIC[2]
diff100<-df100ord$AIC[1]-df100ord$AIC[2]
diff500<-df500ord$AIC[1]-df500ord$AIC[2]
diff1000<-df1000ord$AIC[1]-df1000ord$AIC[2]
diff2000<-df2000ord$AIC[1]-df2000ord$AIC[2]

diff_df<-as.data.frame(matrix(data=c(
                            50,100,500,1000,2000,
                            diff50,diff100,diff500,diff1000,diff2000
                            ),nrow=5,byrow=FALSE))

ggplot(data=diff_df)+
  geom_line(aes(x=V1,y=V2))+
  geom_point(aes(x=V1,y=V2))

model<-lm(V2~V1,data=diff_df)
summary(model)

new_data<-as.data.frame(matrix(data=c(1),nrow=1,byrow=FALSE))
predict(model,newdata=new_data)
