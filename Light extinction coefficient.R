# attach LI data sheet

attach(Final_Lights_interception_data)

PAR_LI<-read.csv(file ="/Users/pjayasin/Documents/Litriture/Supervisors documents/From Dave/Gatton Dairy Research data/Final Lights interception data.csv")
library(dplyr)















Brachiaria_LI<-filter(Final_Lights_interception_data, Pasture=="Brachiaria")
Brachiaria_LI

New_Brachiaria_LI<-Brachiaria_LI[-c(5,6,12), ]
New_Brachiaria_LI

library(tidyverse)
library(broom)
library(ggplot2)
library(ggrepel)

attach(LI_RR)


par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(LAI_Total,Light_fractionused,xlim = c(0,max(LAI_Total)+2), ylim=c(0,max(Light_fractionused)+0.5))

LI_nls_RR<-nls(LI_RR_fraction~ 1- exp(-k*LAI_Total), start = list(k=0.2))   
LI_nls_RR            
summary(LI_nls_RR)
LI_nls_coef<-coef(LI_nls_RR)
LI_nls_coef
function(x)  {(exp(-(LI_nls_coef[2])*LAI_Total))}

#######################

# Brachiaria ###########

Graph_LI_BM<-New_Brachiaria_LI$LightSoilcorrected
Graph_LAI_BM<-New_Brachiaria_LI$TotalLAI
Graph_BM<-data.frame(Graph_LI_BM,Graph_LAI_BM)
Graph_BM

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(TotalLAI,LightSoilcorrected,xlim = c(0,max(TotalLAI)+2), ylim=c(0,max(LightSoilcorrected)+0.5))

BrachiariaLI<-nls(Graph_LI_BM~ 1- exp(-k*Graph_LAI_BM), start = list(k=0.2))   
BrachiariaLI   
summary(BrachiariaLI)

BrachiariaLI_coef<-coef(BrachiariaLI)
BrachiariaLI_coef
function(x)  {(exp(-(BrachiariaLI_coef[2])*Graph_LAI_BM))}
library(ggplot2)
ggplot(Graph_BM, aes(Graph_LAI_BM, Graph_LI_BM )) + ylim(0,1) + xlim(0, 8)+
  geom_point(color="limegreen", size=3, alpha=0.8) + xlab("LAI") +ylab("Fraction of light intercepted (LI)")+
  geom_smooth(method = "nls", method.args = list(formula = y ~ 1- exp(-k*x),
                                                 start = list(k=0.94)),  
              data = Graph_BM,colour = "green",
              se = FALSE, fullrange = TRUE )

text( LI = 1-(e^(-k*LAI)))

summary(LI_nls_RR)
LI_nls_coef<-coef(LI_nls_RR)
LI_nls_coef
par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(LAI_Total,Light_fractionused,xlim = c(0,max(LAI_Total)+2), ylim=c(0,max(Light_fractionused)+ 0.5))
curve(1- exp(-(LI_nls_coef[2])*LAI_Total), col="red", lwd= 2)


###########
predict_LI_RR<-predict(LI_nls_RR)
predict_LI_RR
plot(predict_LI_RR, Light_used)
qplot(LAI_Total, Light_fractionused, data = augment(LI_nls_RR)) + geom_line(aes(y = .fitted), colour= "brown") + xlim(0,max(LAI_Total)+1)+ ylim(0,max(Light_fractionused)+0.5)

ggplot(LI_RR, aes(LAI_Total, Light_fractionused )) + ylim(0,1) + xlim(0, 8)+
  geom_point () + xlab("LAI") +ylab("Fraction of light intercepted (LI)")+
  geom_smooth(method = "nls", method.args = list(formula = y ~ 1- exp(-k*x),
                                               start = list(k=0.4)),  
              data = LI_RR,colour = "red",
              se = FALSE, fullrange = TRUE )
text( LI = 1-(e^(-k*LAI)))

# Gatton panic ##########

Gatton_LI<-filter(Final_Lights_interception_data, Pasture=="Gatton panic")
Gatton_LI
New_Gatton_LI<-Gatton_LI[-c(5,8), ]
New_Gatton_LI



Graph_LI_GP<-New_Gatton_LI$LightSoilcorrected
Graph_LAI_GP<-New_Gatton_LI$TotalLAI
Graph_GP<-data.frame(Graph_LI_GP,Graph_LAI_GP)
Graph_GP

GattonLI<-nls(Graph_LI_GP~ 1- exp(-k*Graph_LAI_GP), start = list(k=0.2))   
GattonLI   
summary(GattonLI)

GattonLI_coef<-coef(GattonLI)
GattonLI_coef
function(x)  {(exp(-(GattonLI_coef[2])*Graph_LAI_GP))}
library(ggplot2)
ggplot(Graph_GP, aes(Graph_LAI_GP, Graph_LI_GP )) + ylim(0,1) + xlim(0, 8)+
  geom_point(color="navyblue", size=3, alpha=0.8) + xlab("LAI") +ylab("Fraction of light intercepted (LI)")+
  geom_smooth(method = "nls", method.args = list(formula = y ~ 1- exp(-k*x),
                                                 start = list(k=0.2)),  
              data = Graph_GP,colour = "blue",
              se = FALSE, fullrange = TRUE )

# Rhodes grass

Rhodes_LI<-filter(Final_Lights_interception_data, Pasture=="Rhodes grass")
Rhodes_LI
New_Rhodes_LI<-Rhodes_LI[-c(5,6,7,14,15,16), ]
New_Rhodes_LI



Graph_LI_RR<-New_Rhodes_LI$LightSoilcorrected
Graph_LAI_RR<-New_Rhodes_LI$TotalLAI
Graph_RR<-data.frame(Graph_LI_RR,Graph_LAI_RR)
Graph_RR

RhodesLI<-nls(Graph_LI_RR~ 1- exp(-k*Graph_LAI_RR), start = list(k=0.2))   
RhodesLI   
summary(RhodesLI)

RhodesLI_coef<-coef(RhodesLI)
RhodesLI_coef
function(x)  {(exp(-(RhodesLI_coef[2])*Graph_LAI_RR))}
library(ggplot2)
ggplot(Graph_RR, aes(Graph_LAI_RR, Graph_LI_RR )) + ylim(0,1) + xlim(0, 8)+
  geom_point(color="violetred2", size=3, alpha=0.8) + xlab("LAI") +ylab("Fraction of light intercepted (LI)")+
  geom_smooth(method = "nls", method.args = list(formula = y ~ 1- exp(-k*x),
                                                 start = list(k=0.2)),  
              data = Graph_RR,colour = "red",
              se = FALSE, fullrange = TRUE )
totalLI<- data.frame(Graph_LI_BM,Graph_LI_GP,Graph_LI_RR)

#plot three curves#

library(plyr) #join unequal tables

totalcurvedata<-rbind.fill(New_Brachiaria_LI,New_Gatton_LI,New_Rhodes_LI)
totalcurvedata
totalLI_grapgh<-data.frame(Graph_LI_BM,Graph_LI_GP,Graph_LI_RR)
totalLAI_graph<- data.frame(Graph_LAI_BM,Graph_LAI_GP,Graph_LAI_RR)
totalgraph_data<-data.frame(totalLI_grapgh,totalLAI_graph)

ggplot(totalcurvedata, aes(TotalLAI, LightSoilcorrected, color=Pasture )) + ylim(0,1) + xlim(0, 8)+
  geom_point( cosize=3, alpha=0.8) + xlab("LAI") + scale_color_manual(values = c("limegreen", "navyblue", "violetred"))+ ylab("Fraction of light intercepted (LI)")+
  geom_smooth(method = "nls", method.args = list(formula = y ~ 1- exp(-k*x),
                                                 start = list(k=0.2)),  
              data = totalcurvedata,
              se = FALSE, fullrange = TRUE)


###############

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(LAI_GP,LI_usedfraction  ,xlim = c(0,max(LAI_GP)+2), ylim=c(0,max(LI_usedfraction)+0.5))

lines(LAI_GP,LI_nls_coef[1]*(1- exp(-(LI_nls_coef[2])*LAI_Total)), col="red", lwd= 2)

LI_nls_GP<- nls(LI_GP_transmited~ 1- exp(-k*LAI_GP), start = list(k=0.2))   
LI_nls_GP
ggplot(LI_GP, aes(LAI_GP, LI_usedfraction )) + ylim(0,1) + xlim(0, 8)+
  geom_point () + xlab("LAI") +ylab("Fraction of light intercepted (LI)")+
  geom_smooth(method = "nls", method.args = list(formula = y ~ 1- exp(-k*x),
                                                 start = list(k=0.4)),  
              data = LI_GP,colour = "blue",
              se = FALSE, fullrange = TRUE )

LI_nls_coef_GP<-coef(LI_nls_GP)
LI_nls_coef_GP

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(LAI_GP,Light_GP  ,xlim = c(0,max(LAI_GP)+2), ylim=c(0,max(Light_GP)+500))
lines(LAI_GP, LI_nls_coef_GP[1]*(1- exp(-(LI_nls_coef_GP[2])*LAI_GP)), col="green", lwd= 2)
attach(LI_Brachiaria)

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(LAI_BM,Light_BM  ,xlim = c(0,max(LAI_BM)+2), ylim=c(0,max(Light_BM)+500))

LI_nls_BM<- nls(light_BM_usedfraction~ 1- exp(-k*LAI_BM), start = list(k=0.5))   
LI_nls_BM
summary(LI_nls_BM)
Plot_LI<-ggplot(LI_Brachiaria, aes(LAI_BM, light_BM_usedfraction )) + ylim(0,1) + xlim(0, 10)+
  geom_point () + xlab("LAI") +ylab("Fraction of light intercepted (LI)") +
  geom_smooth(method = "nls", method.args = list(formula = y ~ 1- exp(-k*x),
                                                 start = list(k=0.4)),  
              data = LI_Brachiaria,colour = "green",
              se = FALSE, fullrange = TRUE )
Plot_LI
Plot_LI+ geom_text(x=5, y=0.5, label= "LI = 1 - e^(-k*LAI) 
                                        k=0.47") 
LI_predict= 1- exp(-0.68*LAItot)
LI_predict
LI_BM_predict<-data.frame(LI_predict)
LI_BM_predict
BM_LI_predict<-LI_BM_predict$LI_predict
BMpredictdata<-data.frame(LI_predict, LAItot)
dataLAI<-rbind.fill(Graph_LI_BM, BM_LI_predict)
 write.csv(BMpredictdata, "predictLIBMdata.csv")
write.csv(New_Brachiaria_LI, "actuaLIBMdata.csv")
###

graph_predict<- data.frame(LI_predict, LAItot)

Plot_LIPredict<-ggplot(graph_predict, aes(LAItot, LI_predict )) + ylim(0,1) + xlim(0, 10)+
  geom_point (color="navyblue", size=3, alpha=0.8) + xlab("LAI") +ylab("Fraction of light intercepted (LI)") +
  geom_smooth(method = "nls", method.args = list(formula = y ~ 1- exp(-k*x),
                                                 start = list(k=0.68)),  
              data = graph_predict,colour = "red",
              se = FALSE, fullrange = TRUE )
Plot_LIPredict
####
LI_nls_BM_I<- nls(LI_BM_intercepted~ a* exp(-k*LAI_BM), start = list(a= max(LI_BM_intercepted), k=0.5))   
LI_nls_BM_I

Plot_LI<-ggplot(LI_Brachiaria, aes(LAI_BM, LI_BM_intercepted )) + ylim(0,1) + xlim(0, 10)+
  geom_point () + xlab("LAI") +ylab("PAR transmittance (PARt/ PARi)") +
  geom_smooth(method = "nls", method.args = list(formula = y ~ a* exp(-k*x),
                                                 start = list(a= 0.87, k=0.4)),  
              data = LI_Brachiaria,colour = "green",
              se = FALSE, fullrange = TRUE )
#predict vs actual
attach(predictLIBMdata)

ggplot(predictLIBMdata, aes(LAItot, LI_predict, color=method )) + ylim(0,1) + xlim(0, 8)+
  geom_point( cosize=3, alpha=0.8) + xlab("LAI") + scale_color_manual(values = c("limegreen", "navyblue"))+ ylab("Fraction of light intercepted (LI)")+
  geom_smooth(method = "nls", method.args = list(formula = y ~ 1- exp(-k*x),
                                                 start = list(k=0.68)),  
              data = predictLIBMdata,
              se = FALSE, fullrange = FALSE)
cor(LI_predict, Graph_LI_BM)

predicted<-filter(predictLIBMdata, method=="predicted")
actual<- filter(predictLIBMdata, method=="actual")

compareplot<-rbind.fill (predicted,actual)
predictvsactual<-ggplot(compareplot, aes(x=actual, y=predicted )) 
  geom_point() + xlab("LAI") + scale_color_manual(values = c("limegreen", "navyblue"))+ ylab("Fraction of light intercepted (LI)")+
  
attach(LAIData)

predict(LI_nls_BM, newdata = data.frame(x=LAItot))


#############use csv data to calculation of K###

BM_k_data<-filter(PAR_LI, Pasture=="Brachiaria")
BM_k_data

BM_k_LAI<- BM_k_data$TotalLAI
BM_k_PAR<- BM_k_data$LI_PAR

BM_k<-nls(BM_k_PAR~ 1- exp(-k_BM*BM_k_LAI), start = list(k_BM=0.2))   
BM_k  

plot_BM_k<-ggplot(BM_k_data, aes(BM_k_LAI, BM_k_PAR)) + geom_point(color= "red")
plot_BM_k
hist(BM_k_LAI)
hist(BM_k_PAR)
