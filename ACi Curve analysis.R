#CO2 curve analysis
meanCo2data<-read.csv(file = "/Users/pjayasin/OneDrive - Massey University/Meta Analysis data/PhD_Data_Analysis/meanCo2data.csv")
meanCo2data

attach(meanCo2data)
library(dplyr)

BM_aci<-filter(meanCo2data, species=="Brachiaria")
BM_aci
BM_aci_omit<-BM_aci[-c(4), ]
BM_aci_omit

BM_Ci<-BM_aci_omit$avgCi
BM_Ci_pa<-BM_aci_omit$AvgCipa
BM_photo<-BM_aci_omit$avgphoto_aci
BM_CO2R<-BM_aci_omit$avgCO2R
BM_Tleaf<-BM_aci_omit$avgTleaf
BM_SE_Photo<-BM_aci_omit$stdPhoto_aci
BM_co2_func<-BM_aci_omit$CO2function

BM_curve_aci<-data.frame(BM_Ci,BM_Ci_pa,BM_photo,BM_CO2R,BM_Tleaf,BM_SE_Photo,BM_co2_func)
BM_curve_aci

#Plot co2 response with photo

curve.nlslrc_BMc = nls(BM_photo ~ (1/(2*theta))*(AQY*BM_CO2R+Am-sqrt((AQY*BM_CO2R+Am)^2- 4*AQY*theta*Am*BM_CO2R))-Rd, start=list(Am= (max(BM_photo)- min(BM_photo)),AQY=0.05, Rd=-min(BM_photo),theta=1))   

summary(curve.nlslrc_BMc)

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(BM_CO2R,BM_photo ,xaxt = "n", yaxt = "n",frame.plot = FALSE, xlab="", xlim = c(0,max(BM_CO2R)+2), ylab="", ylim=c(0,max(BM_photo)+5),cex.lab=3.2,cex.axis=2.5,cex=1.5)
arrows(BM_CO2R, BM_photo-BM_SE_Photo, BM_CO2R, BM_photo+BM_SE_Photo, length = 0.05, angle = 90, code = 3)

mtext(expression(CO [2]~ ppm),side=1,line=2.5,cex=1, cex.axis=1)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1)
curve((1/(2*summary(curve.nlslrc_BMc)$coef[4,1]))*(summary(curve.nlslrc_BMc)$coef[2,1]*x+summary(curve.nlslrc_BMc)$coef[1,1]-sqrt((summary(curve.nlslrc_BMc)$coef[2,1]*x+summary(curve.nlslrc_BMc)$coef[1,1])^2-4*summary(curve.nlslrc_BMc)$coef[2,1]*summary(curve.nlslrc_BMc)$coef[4,1]*summary(curve.nlslrc_BMc)$coef[1,1]*x))-summary(curve.nlslrc_BMc)$coef[3,1],lwd=2,col="orange",add=T)
abline(h=0)
abline(v=0)
axis(1, pos = 0,cex.axis=1)
axis(2, pos = 0, cex.axis=1)
box(col= "black")

#CO2 compensation point
x<-function(x) {(1/(2*summary(curve.nlslrc_BMc)$coef[4,1]))*(summary(curve.nlslrc_BMc)$coef[2,1]*x+summary(curve.nlslrc_BMc)$coef[1,1]-sqrt((summary(curve.nlslrc_BMc)$coef[2,1]*x+summary(curve.nlslrc_BMc)$coef[1,1])^2-4*summary(curve.nlslrc_BMc)$coef[2,1]*summary(curve.nlslrc_BMc)$coef[4,1]*summary(curve.nlslrc_BMc)$coef[1,1]*x))-summary(curve.nlslrc_BMc)$coef[3,1]}
uniroot(x,c(0,1))$root 
CCP=Rd/0.13874
CCP

#CO2 saturation point 
x<-function(x) {(1/(2*summary(curve.nlslrc_BMc)$coef[4,1]))*(summary(curve.nlslrc_BMc)$coef[2,1]*x+summary(curve.nlslrc_BMc)$coef[1,1]-sqrt((summary(curve.nlslrc_BMc)$coef[2,1]*x+summary(curve.nlslrc_BMc)$coef[1,1])^2-4*summary(curve.nlslrc_BMc)$coef[2,1]*summary(curve.nlslrc_BMc)$coef[4,1]*summary(curve.nlslrc_BMc)$coef[1,1]*x))-summary(curve.nlslrc_BMc)$coef[3,1]-(0.90*summary(curve.nlslrc_BMc)$coef[1,1])+0.90*(summary(curve.nlslrc_BMc)$coef[3,1])}

uniroot(x,c(0,1000))$root


#Aci curve fitting

R_BM=0.008314 #(kJ mol^-1 K^-1)
Kc_BM=exp(35.9774-80.99/(R*(BM_curve_aci$BM_Tleaf+273.15))) #Michaelis-Menten constant for Rubisco for O2 (Pa)
Ko_BM=exp(12.3772-23.72/(R*(BM_curve_aci$BM_Tleaf+273.15))) #Michaelis-Menten constant for Rubisco for CO2 (kPa)
GammaStar_BM=exp(11.187-24.46/(R*(BM_curve_aci$BM_Tleaf+273.15))) #Photorespiration compensation point (Pa)
O_BM=21

BM_aci_fit<-nls(BM_photo~ifelse(((Vcmax*(BM_Ci-GammaStar_BM))/(BM_Ci+(Kc_BM*(1+(O_BM/Ko_BM)))))<((J*(BM_Ci-GammaStar_BM))/((4*BM_Ci)+(8*GammaStar_BM))),((Vcmax*(BM_Ci-GammaStar_BM))/(BM_Ci+(Kc_BM*(1+(O_BM/Ko_BM))))),((J*(BM_Ci-GammaStar_BM))/((4*BM_Ci)+(8*GammaStar_BM))))-Rd,start=list(Vcmax=50,J=100,Rd=5),data=BM_curve_aci) #if error: reconsider starting values, bad dataset? (too few points or response curve not clear)

summary(BM_aci_fit)

Vcmax_BM<-summary(BM_aci_fit)$coef[1,1]
Vcmax_BM
J_BM<-summary(BM_aci_fit)$coef[2,1]
J_BM
Rd_BM<-summary(BM_aci_fit)$coef[3,1]
Rd_BM

# ---Graph raw data with modeled curve---
par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(BM_Ci,BM_photo,ylab="", xlab="",cex.lab=1.2,cex.axis=1.5,cex=2)
arrows(BM_Ci, BM_photo-BM_SE_Photo, BM_Ci, BM_photo+ BM_SE_Photo, length = 0.05, angle = 90, code = 3)
mtext(expression("Ci (ppm)"),side=1,line=3.3,cex=1.5)
mtext(expression("Net photosynthetic rate  (umol  "* CO[2]*   m^-2*   s^-1*")"),side=2,line=2.5,cex=1.5)

Ci_BM<-curve(ifelse(((Vcmax_BM*(x-mean(GammaStar_BM)))/(x+(mean(Kc_BM)*(1+(O_BM/mean(Ko_BM))))))<((J*(x-mean(GammaStar_BM)))/((4*x)+(8*mean(GammaStar_BM)))),((Vcmax_BM*(x-mean(GammaStar_BM)))/(x+(mean(Kc_BM)*(1+(O_BM/mean(Ko_BM)))))),((J_BM*(x-mean(GammaStar_BM)))/((4*x)+(8*mean(GammaStar_BM)))))-Rd_BM,col="green",add=T) #Reasonable fit? Could check goodness of fit, model assumptions
Ci_BM

#Rhode grass

RR_aci<-filter(meanCo2data, species=="Rhodes grass")
RR_aci

RR_aci_omit<-RR_aci[-c(3,4,8), ]
RR_aci_omit

RR_Ci<-RR_aci_omit$avgCi
RR_Ci_pa<-RR_aci_omit$AvgCipa
RR_photo<-RR_aci_omit$avgphoto_aci
RR_CO2R<-RR_aci_omit$avgCO2R
RR_Tleaf<-RR_aci_omit$avgTleaf
RR_SE_Photo<-RR_aci_omit$stdPhoto_aci
RR_CO2func<-RR_aci_omit$CO2function

RR_curve_aci<-data.frame(RR_Ci,RR_Ci_pa,RR_photo,RR_CO2R,RR_Tleaf,RR_SE_Photo,RR_CO2func)
RR_curve_aci
write.csv(RR_curve_aci, "RR_CO2.csv")
#Rhodes grass CO2 response curve 

curve.nlslrc_RRc = nls(RR_photo ~ (1/(2*theta))*(AQY*RR_CO2R+Am-sqrt((AQY*RR_CO2R+Am)^2- 4*AQY*theta*Am*RR_CO2R))-Rd, start=list(Am= (max(RR_photo)- min(RR_photo)),AQY=0.05, Rd=-min(RR_photo),theta=1))   

summary(curve.nlslrc_RRc)


par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(RR_CO2R,RR_photo ,xaxt = "n", yaxt = "n",frame.plot = FALSE, xlab="", xlim = c(0,max(RR_CO2R)+2), ylab="", ylim=c(0,max(RR_photo)+10),cex.lab=1.2,cex.axis=1,cex=1.5)
arrows(RR_CO2R, RR_photo-RR_SE_Photo, RR_CO2R, RR_photo+RR_SE_Photo, length = 0.05, angle = 90, code = 3)
mtext(expression(CO [2]~ ppm),side=1,line=2.5,cex=1)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1)
curve((1/(2*summary(curve.nlslrc_RRc)$coef[4,1]))*(summary(curve.nlslrc_RRc)$coef[2,1]*x+summary(curve.nlslrc_RRc)$coef[1,1]-sqrt((summary(curve.nlslrc_RRc)$coef[2,1]*x+summary(curve.nlslrc_RRc)$coef[1,1])^2-4*summary(curve.nlslrc_RRc)$coef[2,1]*summary(curve.nlslrc_RRc)$coef[4,1]*summary(curve.nlslrc_RRc)$coef[1,1]*x))-summary(curve.nlslrc_RRc)$coef[3,1],lwd=2,col="brown",add=T)
abline(h=0)
abline(v=0)
axis(1, pos = 0)
axis(2, pos = 0)
box(col="black")




#ACi curve fitting for Rhodes grass

R_RR=0.008314 #(kJ mol^-1 K^-1)
Kc_RR=exp(35.9774-80.99/(R*(RR_curve_aci$RR_Tleaf+273.15))) #Michaelis-Menten constant for Rubisco for O2 (Pa)
Ko_RR=exp(12.3772-23.72/(R*(RR_curve_aci$RR_Tleaf+273.15))) #Michaelis-Menten constant for Rubisco for CO2 (kPa)
GammaStar_RR=exp(11.187-24.46/(R*(RR_curve_aci$RR_Tleaf+273.15))) #Photorespiration compensation point (Pa)
O_RR=21

RR_aci_fit<-nls(RR_photo~ifelse(((Vcmax*(RR_Ci-GammaStar_RR))/(RR_Ci+(Kc_RR*(1+(O_RR/Ko_RR)))))<((J*(RR_Ci-GammaStar_RR))/((4*RR_Ci)+(8*GammaStar_RR))),((Vcmax*(RR_Ci-GammaStar_RR))/(RR_Ci+(Kc_RR*(1+(O_RR/Ko_RR))))),((J*(RR_Ci-GammaStar_RR))/((4*RR_Ci)+(8*GammaStar_RR))))-Rd,start=list(Vcmax=50,J=80,Rd=0.5),data=RR_curve_aci) #if error: reconsider starting values, bad dataset? (too few points or response curve not clear)

summary(RR_aci_fit)

Vcmax_RR<-summary(RR_aci_fit)$coef[1,1]
Vcmax_RR
J_RR<-summary(RR_aci_fit)$coef[2,1]
J_RR
Rd_RR<-summary(RR_aci_fit)$coef[3,1]
Rd_RR

# ---Graph raw data with modeled curve---
par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(RR_Ci,RR_photo,ylab="", xlab="",cex.lab=1.2,cex.axis=1.5,cex=2)
arrows(RR_Ci, RR_photo-RR_SE_Photo, RR_Ci, RR_photo+RR_SE_Photo, length = 0.05, angle = 90, code = 3)
mtext(expression("Ci (ppm)"),side=1,line=3.3,cex=1.5)
mtext(expression("Net photosynthetic rate  (umol  "* CO[2]*   m^-2*   s^-1*")"),side=2,line=2.5,cex=1.5)

Ci_RR<-curve(ifelse(((Vcmax_RR*(x-mean(GammaStar_RR)))/(x+(mean(Kc_RR)*(1+(O_RR/mean(Ko_RR))))))<((J_RR*(x-mean(GammaStar_RR)))/((4*x)+(8*mean(GammaStar_RR)))),((Vcmax_RR*(x-mean(GammaStar_RR)))/(x+(mean(Kc)*(1+(O_RR/mean(Ko_RR)))))),((J_RR*(x-mean(GammaStar_RR)))/((4*x)+(8*mean(GammaStar_RR)))))-Rd_RR,col="red",add=T) #Reasonable fit? Could check goodness of fit, model assumptions
Ci_RR
#Gaton panic 

GP_aci<-filter(meanCo2data, species=="Gatton pannic")
GP_aci


GP_Ci<-GP_aci$avgCi
GP_Ci_pa<-GP_aci$AvgCipa
GP_photo<-GP_aci$avgphoto_aci
GP_CO2R<-GP_aci$avgCO2R
GP_Tleaf<-GP_aci$avgTleaf
GP_SE_Photo<-GP_aci$stdPhoto_aci
GP_CO2func<-GP_aci$CO2function
GP_curve_aci<-data.frame(GP_Ci,GP_Ci_pa, GP_photo,GP_CO2R, GP_Tleaf, GP_SE_Photo,GP_CO2func)
GP_curve_aci
#CO2 response 
par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(GP_CO2R,GP_photo,xlab="PAR", ylab="Assimilation")
mtext(expression(CO [2]~ ppm),side=1,line=2.5,cex=1)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1)
abline(h=0)
abline(v=0)
axis(1, pos = 0)
axis(2, pos = 0)

curve.nlslrc_GPc = nls(GP_photo ~ (1/(2*theta))*(AQY*GP_CO2R+Am-sqrt((AQY*GP_CO2R+Am)^2- 4*AQY*theta*Am*GP_CO2R))-Rd, start=list(Am= (max(GP_photo)- min(GP_photo)),AQY=0.05, Rd=-min(GP_photo),theta=1))   

summary(curve.nlslrc_GPc)

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(GP_CO2R,GP_photo ,xaxt = "n", yaxt = "n",frame.plot = FALSE, xlab="", xlim = c(-50,max(GP_CO2R)+2), ylab="", ylim=c(-5,max(GP_photo)+10),cex.lab=1.2,cex.axis=1,cex=1.5)
arrows(GP_CO2R, GP_photo-GP_SE_Photo, GP_CO2R, GP_photo+GP_SE_Photo, length = 0.05, angle = 90, code = 3)

mtext(expression(CO [2]~ ppm),side=1,line=2.5,cex=1)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1)
abline(h=0)
abline(v=0)
axis(1, pos = 0)
axis(2, pos = 0)
box(col="black")

curve((1/(2*summary(curve.nlslrc_GPc)$coef[4,1]))*(summary(curve.nlslrc_GPc)$coef[2,1]*x+summary(curve.nlslrc_GPc)$coef[1,1]-sqrt((summary(curve.nlslrc_GPc)$coef[2,1]*x+summary(curve.nlslrc_GPc)$coef[1,1])^2-4*summary(curve.nlslrc_GPc)$coef[2,1]*summary(curve.nlslrc_GPc)$coef[4,1]*summary(curve.nlslrc_GPc)$coef[1,1]*x))-summary(curve.nlslrc_GPc)$coef[3,1],lwd=2,col="purple",add=T)

#CO2 compensation point
x<-function(x) {(1/(2*summary(curve.nlslrc_GPc)$coef[4,1]))*(summary(curve.nlslrc_GPc)$coef[2,1]*x+summary(curve.nlslrc_GPc)$coef[1,1]-sqrt((summary(curve.nlslrc_GPc)$coef[2,1]*x+summary(curve.nlslrc_GPc)$coef[1,1])^2-4*summary(curve.nlslrc_GPc)$coef[2,1]*summary(curve.nlslrc_GPc)$coef[4,1]*summary(curve.nlslrc_GPc)$coef[1,1]*x))-summary(curve.nlslrc_GPc)$coef[3,1]}
uniroot(x,c(0,50))$root 
CCP=Rd/0.13874
CCP

#CO2 saturation point 
x<-function(x) {(1/(2*summary(curve.nlslrc_GPc)$coef[4,1]))*(summary(curve.nlslrc_GPc)$coef[2,1]*x+summary(curve.nlslrc_GPc)$coef[1,1]-sqrt((summary(curve.nlslrc_GPc)$coef[2,1]*x+summary(curve.nlslrc_GPc)$coef[1,1])^2-4*summary(curve.nlslrc_GPc)$coef[2,1]*summary(curve.nlslrc_GPc)$coef[4,1]*summary(curve.nlslrc_GPc)$coef[1,1]*x))-summary(curve.nlslrc_GPc)$coef[3,1]-(0.90*summary(curve.nlslrc_GPc)$coef[1,1])+0.90*(summary(curve.nlslrc_GPc)$coef[3,1])}
uniroot(x,c(0,1000))$root


R_GP=0.008314 #(kJ mol^-1 K^-1)
Kc_GP=exp(35.9774-80.99/(R*(GP_curve_aci$GP_Tleaf+273.15))) #Michaelis-Menten constant for Rubisco for O2 (Pa)
Ko_GP=exp(12.3772-23.72/(R*(GP_curve_aci$GP_Tleaf+273.15))) #Michaelis-Menten constant for Rubisco for CO2 (kPa)
GammaStar_GP=exp(11.187-24.46/(R*(GP_curve_aci$GP_Tleaf+273.15))) #Photorespiration compensation point (Pa)
O_GP=21

GP_aci_fit<-nls(GP_photo~ifelse(((Vcmax*(GP_Ci-GammaStar_GP))/(GP_Ci+(Kc_GP*(1+(O_GP/Ko_GP)))))<((J*(GP_Ci-GammaStar_GP))/((4*GP_Ci)+(8*GammaStar_GP))),((Vcmax*(GP_Ci-GammaStar_GP))/(GP_Ci+(Kc_GP*(1+(O_GP/Ko_GP))))),((J*(GP_Ci-GammaStar))/((4*GP_Ci)+(8*GammaStar_GP))))-Rd,start=list(Vcmax=50,J=80,Rd=0.5),data=GP_curve_aci) #if error: reconsider starting values, bad dataset? (too few points or response curve not clear)

summary(GP_aci_fit)

Vcmax_GP<-summary(GP_aci_fit)$coef[1,1]
Vcmax_GP
J_GP<-summary(GP_aci_fit)$coef[2,1]
J_GP
Rd_GP<-summary(GP_aci_fit)$coef[3,1]
Rd_GP

# ---Graph raw data with modeled curve---
par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(GP_Ci,GP_photo,ylab="", xlab="",cex.lab=1.2,cex.axis=1.5,cex=2)
arrows(GP_Ci, GP_photo-GP_SE_Photo, GP_Ci, GP_photo+GP_SE_Photo, length = 0.05, angle = 90, code = 3)
mtext(expression("Ci (ppm)"),side=1,line=3.3,cex=1.5)
mtext(expression("Net photosynthetic rate  (umol  "* CO[2]*   m^-2*   s^-1*")"),side=2,line=2.5,cex=1.5)

Ci_GP<-curve(ifelse(((Vcmax_GP*(x-mean(GammaStar_GP)))/(x+(mean(Kc_GP)*(1+(O_GP/mean(Ko_GP))))))<((J_GP*(x-mean(GammaStar_GP)))/((4*x)+(8*mean(GammaStar_GP)))),((Vcmax_GP*(x-mean(GammaStar_GP)))/(x+(mean(Kc_GP)*(1+(O_GP/mean(Ko_GP)))))),((J_GP*(x-mean(GammaStar_GP)))/((4*x)+(8*mean(GammaStar_GP)))))-Rd_GP,col="blue",add=T) #Reasonable fit? Could check goodness of fit, model assumptions
Ci_GP
#all three A-Ci cures

Total_aci_photo<-c(RR_photo, BM_photo, GP_photo)
Total_aci_ci<-c(RR_Ci, BM_Ci, GP_Ci)

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(Total_aci_ci,Total_aci_photo,xaxt = "n", yaxt = "n",frame.plot = FALSE, xlab="", xlim = c(0,max(Total_aci_ci)+20), ylab="", ylim=c(0,max(Total_aci_photo)+5),cex.lab=1.2,cex.axis=1,cex=1.5)

arrows(RR_Ci, RR_photo-RR_SE_Photo, RR_Ci, RR_photo+RR_SE_Photo, length = 0.05, angle = 90, code = 3)
curve(ifelse(((Vcmax_RR*(x-mean(GammaStar_RR)))/(x+(mean(Kc_RR)*(1+(O_RR/mean(Ko_RR))))))<((J_RR*(x-mean(GammaStar_RR)))/((4*x)+(8*mean(GammaStar_RR)))),((Vcmax_RR*(x-mean(GammaStar_RR)))/(x+(mean(Kc_RR)*(1+(O_RR/mean(Ko_RR)))))),((J_RR*(x-mean(GammaStar_RR)))/((4*x)+(8*mean(GammaStar_RR)))))-Rd_RR,lwd=2, col="deeppink3",add=T) #Reasonable fit? Could check goodness of fit, model assumptions

arrows(BM_Ci, BM_photo-BM_SE_Photo, BM_Ci, BM_photo+ BM_SE_Photo, length = 0.05, angle = 90, code = 3)
curve(ifelse(((Vcmax_BM*(x-mean(GammaStar_BM)))/(x+(mean(Kc_BM)*(1+(O_BM/mean(Ko_BM))))))<((J_BM*(x-mean(GammaStar_BM)))/((4*x)+(8*mean(GammaStar_BM)))),((Vcmax_BM*(x-mean(GammaStar_BM)))/(x+(mean(Kc_BM)*(1+(O_BM/mean(Ko_BM)))))),((J_BM*(x-mean(GammaStar_BM)))/((4*x)+(8*mean(GammaStar_BM)))))-Rd_BM,lwd=2, col="darkolivegreen",add=T) #Reasonable fit? Could check goodness of fit, model assumptions

arrows(GP_Ci, GP_photo-GP_SE_Photo, GP_Ci, GP_photo+GP_SE_Photo, length = 0.05, angle = 90, code = 3)
curve(ifelse(((Vcmax_GP*(x-mean(GammaStar_GP)))/(x+(mean(Kc_GP)*(1+(O_GP/mean(Ko_GP))))))<((J_GP*(x-mean(GammaStar_GP)))/((4*x)+(8*mean(GammaStar_GP)))),((Vcmax_GP*(x-mean(GammaStar_GP)))/(x+(mean(Kc_GP)*(1+(O_GP/mean(Ko_GP)))))),((J_GP*(x-mean(GammaStar_GP)))/((4*x)+(8*mean(GammaStar_GP)))))-Rd_GP,lwd=2, col="cornflowerblue",add=T) #Reasonable fit? Could check goodness of fit, model assumptions

mtext(expression("Ci (ppm)"),side=1,line=3.3,cex=1.5)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1.5)
abline(h=0)
abline(v=0)
axis(1, pos = 0)
axis(2, pos = 0)

legend(15,40, legend = c("Brachiaria mulato II", "Gatton panic", "Reclaimer rhodes"), title= "Grass species",c("darkolivegreen", "cornflowerblue", "deeppink3"),
       cex=0.8,box.lty=0)

#CO2 response curves all in one graph


CO2_all_photo<- c(BM_photo,RR_photo, GP_photo)
CO2_all_photo
CO2_all_R<-c(BM_CO2R, RR_CO2R,GP_CO2R)

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(CO2_all_R,CO2_all_photo ,xaxt = "n", yaxt = "n",frame.plot = FALSE, xlab="", xlim = c(0,max(CO2_all_R)+2), ylab="", ylim=c(-10,max(CO2_all_photo)+10),cex.lab=1.2,cex.axis=1,cex=1.5)
arrows(BM_CO2R, BM_photo-BM_SE_Photo, BM_CO2R, BM_photo+BM_SE_Photo, length = 0.05, angle = 90, code = 3)
curve((1/(2*summary(curve.nlslrc_BMc)$coef[4,1]))*(summary(curve.nlslrc_BMc)$coef[2,1]*x+summary(curve.nlslrc_BMc)$coef[1,1]-sqrt((summary(curve.nlslrc_BMc)$coef[2,1]*x+summary(curve.nlslrc_BMc)$coef[1,1])^2-4*summary(curve.nlslrc_BMc)$coef[2,1]*summary(curve.nlslrc_BMc)$coef[4,1]*summary(curve.nlslrc_BMc)$coef[1,1]*x))-summary(curve.nlslrc_BMc)$coef[3,1],lwd=2,col="orange",add=T)


arrows(GP_CO2R, GP_photo-GP_SE_Photo, GP_CO2R, GP_photo+GP_SE_Photo, length = 0.05, angle = 90, code = 3)
curve((1/(2*summary(curve.nlslrc_GPc)$coef[4,1]))*(summary(curve.nlslrc_GPc)$coef[2,1]*x+summary(curve.nlslrc_GPc)$coef[1,1]-sqrt((summary(curve.nlslrc_GPc)$coef[2,1]*x+summary(curve.nlslrc_GPc)$coef[1,1])^2-4*summary(curve.nlslrc_GPc)$coef[2,1]*summary(curve.nlslrc_GPc)$coef[4,1]*summary(curve.nlslrc_GPc)$coef[1,1]*x))-summary(curve.nlslrc_GPc)$coef[3,1],lwd=2,col="purple",add=T)

arrows(RR_CO2R, RR_photo-RR_SE_Photo, RR_CO2R, RR_photo+RR_SE_Photo, length = 0.05, angle = 90, code = 3)
curve((1/(2*summary(curve.nlslrc_RRc)$coef[4,1]))*(summary(curve.nlslrc_RRc)$coef[2,1]*x+summary(curve.nlslrc_RRc)$coef[1,1]-sqrt((summary(curve.nlslrc_RRc)$coef[2,1]*x+summary(curve.nlslrc_RRc)$coef[1,1])^2-4*summary(curve.nlslrc_RRc)$coef[2,1]*summary(curve.nlslrc_RRc)$coef[4,1]*summary(curve.nlslrc_RRc)$coef[1,1]*x))-summary(curve.nlslrc_RRc)$coef[3,1],lwd=2,col="brown",add=T)


abline(h=0)
abline(v=0)
axis(1, pos = 0)
axis(2, pos = 0)
mtext(expression(CO [2]~ ppm),side=1,line=2.0,cex=1.5)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1.5)
legend(95,45, legend = c("BM", "GP", "RR"), title= "Grass species",c("orange", "purple", "brown"),
       col=c(curve.nlslrc_BMc, curve.nlslrc_RRc, curve.nlslrc_GPc), lty=1:2, cex=1,
       box.lty=0)

#CO2 function curve for DairyMod

#Brachiaria mulato II

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(BM_CO2R,BM_co2_func ,xaxt = "n", yaxt = "n",frame.plot = FALSE, xlab="", xlim = c(0,max(BM_CO2R)+2), ylab="", ylim=c(0,max(BM_co2_func)+0.5),cex.lab=1.2,cex.axis=1,cex=1.5)

mtext(expression(CO[2] ~ppm ),side=1,line=3.3,cex=1.5)
mtext(expression(CO[2] ~ "fc (C)"),side=2,line=2.5,cex=1.5)
abline(h=0)
abline(v=0)
axis(1, pos = 0)
axis(2, pos = 0)


curve.nlslrc_BMd = nls(BM_co2_func ~ (1/(2*theta))*(AQY*BM_CO2R+Am-sqrt((AQY*BM_CO2R+Am)^2- 4*AQY*theta*Am*BM_CO2R)), start=list(Am= (max(BM_co2_func)- min(BM_co2_func)),AQY=0.005,theta=1))   

summary(curve.nlslrc_BMd)


BM_co2_func_double= (1/(2 * 0.8319316)) * (0.0049088 * 800 + 1.1791921 - sqrt((0.0049088 * 
                                                                    800 + 1.1791921)^2 - 4 * 0.0049088 * 0.8319316 * 1.1791921  * 800))
BM_co2_func_double


par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(BM_CO2R,BM_co2_func ,xaxt = "n", yaxt = "n",frame.plot = FALSE, xlab="", xlim = c(0,max(BM_CO2R)+2), ylab="", ylim=c(0,max(BM_co2_func)+0.5),cex.lab=1.2,cex.axis=1,cex=1.5)
curve((1/(2*summary(curve.nlslrc_BMd)$coef[3,1]))*(summary(curve.nlslrc_BMd)$coef[2,1]*x+summary(curve.nlslrc_BMd)$coef[1,1]-sqrt((summary(curve.nlslrc_BMd)$coef[2,1]*x+summary(curve.nlslrc_BMd)$coef[1,1])^2-4*summary(curve.nlslrc_BMd)$coef[2,1]*summary(curve.nlslrc_BMd)$coef[3,1]*summary(curve.nlslrc_BMd)$coef[1,1]*x)),lwd=2,col="Orange",add=T)
mtext(expression(CO[2] ~ppm ),side=1,line=3.3,cex=1.5)
mtext(expression(CO[2] ~ "fc (C)"),side=2,line=2.5,cex=1.5)
abline(h=0)
abline(v=0)
axis(1, pos = 0)
axis(2, pos = 0)
box(col="black")

#Gatton panic

curve.nlslrc_GPd = nls(GP_CO2func ~ (1/(2*theta))*(AQY*GP_CO2R+Am-sqrt((AQY*GP_CO2R+Am)^2- 4*AQY*theta*Am*GP_CO2R)), start=list(Am= (max(GP_CO2func)- min(GP_CO2func)),AQY=0.005,theta=1))   

summary(curve.nlslrc_GPd)

GP_CO2func_double= (1/(2 * 0.7946003)) * (0.0061752 * 800 + 1.1904513 - sqrt((0.0061752 * 
                                                                   800 + 1.1904513)^2 - 4 * 0.0061752 * 0.7946003 * 1.1904513 * 800))
GP_CO2func_double

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(GP_CO2R,GP_CO2func ,xaxt = "n", yaxt = "n",frame.plot = FALSE, xlab="", xlim = c(0,max(GP_CO2R)+2), ylab="", ylim=c(0,max(GP_CO2func)+0.5),cex.lab=1.2,cex.axis=1,cex=1.5)
curve((1/(2*summary(curve.nlslrc_GPd)$coef[3,1]))*(summary(curve.nlslrc_GPd)$coef[2,1]*x+summary(curve.nlslrc_GPd)$coef[1,1]-sqrt((summary(curve.nlslrc_GPd)$coef[2,1]*x+summary(curve.nlslrc_GPd)$coef[1,1])^2-4*summary(curve.nlslrc_GPd)$coef[2,1]*summary(curve.nlslrc_GPd)$coef[3,1]*summary(curve.nlslrc_GPd)$coef[1,1]*x)),lwd=2,col="Red",add=T)
mtext(expression(CO[2] ~ppm ),side=1,line=3.3,cex=1.5)
mtext(expression(CO[2] ~ "fc (C)"),side=2,line=2.5,cex=1.5)
abline(h=0)
abline(v=0)
axis(1, pos = 0)
axis(2, pos = 0)
box(col="black")

#Rhodes grass

attach(RR_CO2)
curve.nlslrc_Rhodes = nls(Rhodes_CO2func ~ (1/(2*theta))*(AQY*Rhodes_CO2+Am-sqrt((AQY*Rhodes_CO2+Am)^2- 4*AQY*theta*Am*Rhodes_CO2)), start=list(Am= (max(Rhodes_CO2func)- min(Rhodes_CO2func)),AQY=0.005,theta=1))   

summary(curve.nlslrc_Rhodes)

RR_co2func_double = (1/(2 *  0.311350)) * (0.006596 * 800 +  1.403067 - sqrt(( 0.006596 * 
                                                                   800 +  1.403067)^2 - 4 *  0.006596 * 0.311350 * 1.403067 * 800))
RR_co2func_double

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(Rhodes_CO2,Rhodes_CO2func ,xaxt = "n", yaxt = "n",frame.plot = FALSE, xlab="", xlim = c(0,max(Rhodes_CO2)+2), ylab="", ylim=c(0,max(Rhodes_CO2func)+0.5),cex.lab=1.2,cex.axis=1,cex=1.5)
curve((1/(2*summary(curve.nlslrc_Rhodes)$coef[3,1]))*(summary(curve.nlslrc_Rhodes)$coef[2,1]*x+summary(curve.nlslrc_Rhodes)$coef[1,1]-sqrt((summary(curve.nlslrc_Rhodes)$coef[2,1]*x+summary(curve.nlslrc_Rhodes)$coef[1,1])^2-4*summary(curve.nlslrc_Rhodes)$coef[2,1]*summary(curve.nlslrc_Rhodes)$coef[3,1]*summary(curve.nlslrc_Rhodes)$coef[1,1]*x)),lwd=2,col="darkolivegreen1",add=T)
mtext(expression(CO[2] ~ppm ),side=1,line=3.3,cex=1.5)
mtext(expression(CO[2] ~ "fc (C)"),side=2,line=2.5,cex=1.5)
abline(h=0)
abline(v=0)
axis(1, pos = 0)
axis(2, pos = 0)
box(col="black")

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(RR_CO2R,RR_CO2func, xaxt = "n", yaxt = "n",frame.plot = FALSE, xlab="", xlim = c(0,max(RR_CO2R)+2), ylab="", ylim=c(0,max(RR_CO2func)+0.5),cex.lab=1.2,cex.axis=1,cex=1.5)
curve((1/(2*summary(curve.nlslrc_RRd)$coef[3,1]))*(summary(curve.nlslrc_RRd)$coef[2,1]*x+summary(curve.nlslrc_RRd)$coef[1,1]-sqrt((summary(curve.nlslrc_RRd)$coef[2,1]*x+summary(curve.nlslrc_RRd)$coef[1,1])^2-4*summary(curve.nlslrc_RRd)$coef[2,1]*summary(curve.nlslrc_RRd)$coef[3,1]*summary(curve.nlslrc_RRd)$coef[1,1]*x)),lwd=2,col="Blue",add=T)
mtext(expression(CO[2] ~ppm ),side=1,line=3.3,cex=1.5)
mtext(expression(CO[2] ~ "fc (C)"),side=2,line=2.5,cex=1.5)
abline(h=0)
abline(v=0)
axis(1, pos = 0)
axis(2, pos = 0)
box(col="black")





summary(RR_CO2func)


# three co2 functions in one graph

Total_co2func<-c(BM_co2_func,RR_co2func,GP_CO2func)
Total_co2R<-c(BM_CO2R,RR_CO2R,GP_CO2R)


par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(Total_co2R,Total_co2func ,xaxt = "n", yaxt = "n",frame.plot = FALSE, xlab="", xlim = c(0,max(Total_co2R)+2), ylab="", ylim=c(0,max(Total_co2func)+0.5),cex.lab=1.2,cex.axis=1,cex=1.5)
curve((1/(2*summary(curve.nlslrc_BMd)$coef[3,1]))*(summary(curve.nlslrc_BMd)$coef[2,1]*x+summary(curve.nlslrc_BMd)$coef[1,1]-sqrt((summary(curve.nlslrc_BMd)$coef[2,1]*x+summary(curve.nlslrc_BMd)$coef[1,1])^2-4*summary(curve.nlslrc_BMd)$coef[2,1]*summary(curve.nlslrc_BMd)$coef[3,1]*summary(curve.nlslrc_BMd)$coef[1,1]*x)),lwd=2,col="aquamarine3",add=T)
curve((1/(2*summary(curve.nlslrc_GPd)$coef[3,1]))*(summary(curve.nlslrc_GPd)$coef[2,1]*x+summary(curve.nlslrc_GPd)$coef[1,1]-sqrt((summary(curve.nlslrc_GPd)$coef[2,1]*x+summary(curve.nlslrc_GPd)$coef[1,1])^2-4*summary(curve.nlslrc_GPd)$coef[2,1]*summary(curve.nlslrc_GPd)$coef[3,1]*summary(curve.nlslrc_GPd)$coef[1,1]*x)),lwd=2,col="coral2",add=T)
curve((1/(2*summary(curve.nlslrc_RRd)$coef[3,1]))*(summary(curve.nlslrc_RRd)$coef[2,1]*x+summary(curve.nlslrc_RRd)$coef[1,1]-sqrt((summary(curve.nlslrc_RRd)$coef[2,1]*x+summary(curve.nlslrc_RRd)$coef[1,1])^2-4*summary(curve.nlslrc_RRd)$coef[2,1]*summary(curve.nlslrc_RRd)$coef[3,1]*summary(curve.nlslrc_RRd)$coef[1,1]*x)),lwd=2,col="darkolivegreen1",add=T)
abline(h=0)
abline(v=0)
axis(1, pos = 0)
axis(2, pos = 0)
mtext(expression(CO[2] ~ppm ),side=1,line=3.3,cex=1.5)
mtext(expression(CO[2] ~ "fc (C)"),side=2,line=2.5,cex=1.5)
legend(100,1.5, legend = c("Brachiaria mulato II", "Gatton panic", "Rhodes reclaimer"), title= "Grass species",c("aquamarine3", "coral2", "darkolivegreen1"),
       col=c(curve.nlslrc_BMd, curve.nlslrc_RRd, curve.nlslrc_GPd), lty=1:2, cex=1,
       box.lty=0)

#A-Ci curves

CO2_all_photo<- c(BM_photo,RR_photo, GP_photo)
CO2_all_photo
CO2_all_R<-c(BM_CO2R, RR_CO2R,GP_CO2R)

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(CO2_all_R,CO2_all_photo ,xaxt = "n", yaxt = "n",frame.plot = FALSE, xlab="", xlim = c(0,max(CO2_all_R)+2), ylab="", ylim=c(-10,max(CO2_all_photo)+10),cex.lab=1.2,cex.axis=1,cex=1.5)
arrows(BM_CO2R, BM_photo-BM_SE_Photo, BM_CO2R, BM_photo+BM_SE_Photo, length = 0.05, angle = 90, code = 3)

Ci_BM<-curve(ifelse(((Vcmax_BM*(x-mean(GammaStar_BM)))/(x+(mean(Kc_BM)*(1+(O_BM/mean(Ko_BM))))))<((J*(x-mean(GammaStar_BM)))/((4*x)+(8*mean(GammaStar_BM)))),((Vcmax_BM*(x-mean(GammaStar_BM)))/(x+(mean(Kc_BM)*(1+(O_BM/mean(Ko_BM)))))),((J_BM*(x-mean(GammaStar_BM)))/((4*x)+(8*mean(GammaStar_BM)))))-Rd_BM,col="green",add=T) #Reasonable fit? Could check goodness of fit, model assumptions
Ci_BM

Ci_RR<-curve(ifelse(((Vcmax_RR*(x-mean(GammaStar_RR)))/(x+(mean(Kc_RR)*(1+(O_RR/mean(Ko_RR))))))<((J_RR*(x-mean(GammaStar_RR)))/((4*x)+(8*mean(GammaStar_RR)))),((Vcmax_RR*(x-mean(GammaStar_RR)))/(x+(mean(Kc)*(1+(O_RR/mean(Ko_RR)))))),((J_RR*(x-mean(GammaStar_RR)))/((4*x)+(8*mean(GammaStar_RR)))))-Rd_RR,col="red",add=T) #Reasonable fit? Could check goodness of fit, model assumptions
Ci_RR

Ci_GP<-curve(ifelse(((Vcmax_GP*(x-mean(GammaStar_GP)))/(x+(mean(Kc_GP)*(1+(O_GP/mean(Ko_GP))))))<((J_GP*(x-mean(GammaStar_GP)))/((4*x)+(8*mean(GammaStar_GP)))),((Vcmax_GP*(x-mean(GammaStar_GP)))/(x+(mean(Kc_GP)*(1+(O_GP/mean(Ko_GP)))))),((J_GP*(x-mean(GammaStar_GP)))/((4*x)+(8*mean(GammaStar_GP)))))-Rd_GP,col="blue",add=T) #Reasonable fit? Could check goodness of fit, model assumptions
Ci_GP

abline(h=0)
abline(v=0)
axis(1, pos = 0)
axis(2, pos = 0)
mtext(expression(CO [2]~ ppm),side=1,line=2.0,cex=1.5)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1.5)
legend(95,45, legend = c("Brachiaria mulato II", "Gatton panic", "Rhodes reclaimer"), title= "Grass species",c("orange", "purple", "brown"),
       col=c(curve.nlslrc_BMc, curve.nlslrc_RRc, curve.nlslrc_GPc), lty=1:2, cex=1,
       box.lty=0)



par(mfrow=c(1,2))
totalPARi<-omit_lightmean$avgPARi
totalphoto<-omit_lightmean$avgphoto
par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(totalPARi,totalphoto, xaxt = "n", yaxt = "n",frame.plot = FALSE, xlab="", xlim = c(-500,max(totalPARi)+2), ylab="", ylim=c(-5,max(totalphoto)+5),cex.lab=1.2,cex.axis=1,cex=1.3)
arrows(PARlrcBM, photolrcBM-SEphoto, PARlrcBM, photolrcBM+SEphoto, length = 0.05, angle = 90, code = 3)
curve((1/(2*summary(curve.nlslrc_BM)$coef[4,1]))*(summary(curve.nlslrc_BM)$coef[2,1]*x+summary(curve.nlslrc_BM)$coef[1,1]-sqrt((summary(curve.nlslrc_BM)$coef[2,1]*x+summary(curve.nlslrc_BM)$coef[1,1])^2-4*summary(curve.nlslrc_BM)$coef[2,1]*summary(curve.nlslrc_BM)$coef[4,1]*summary(curve.nlslrc_BM)$coef[1,1]*x))-summary(curve.nlslrc_BM)$coef[3,1],lwd=2,col="Green",add=T)

arrows(PARlrcRR, photolrcRR-SEphotoRR, PARlrcRR, photolrcRR+SEphotoRR, length = 0.05, angle = 90, code = 3)
curve((1/(2*summary(curve.nlslrc_RR)$coef[4,1]))*(summary(curve.nlslrc_RR)$coef[2,1]*x+summary(curve.nlslrc_RR)$coef[1,1]-sqrt((summary(curve.nlslrc_RR)$coef[2,1]*x+summary(curve.nlslrc_RR)$coef[1,1])^2-4*summary(curve.nlslrc_RR)$coef[2,1]*summary(curve.nlslrc_RR)$coef[4,1]*summary(curve.nlslrc_RR)$coef[1,1]*x))-summary(curve.nlslrc_RR)$coef[3,1],lwd=2,col="Red",add=T)

arrows(PARlrcGP, photolrcGP-SEphotoGP, PARlrcGP, photolrcGP+SEphotoGP, length = 0.05, angle = 90, code = 3)
curve((1/(2*summary(curve.nlslrc_GP)$coef[4,1]))*(summary(curve.nlslrc_GP)$coef[2,1]*x+summary(curve.nlslrc_GP)$coef[1,1]-sqrt((summary(curve.nlslrc_GP)$coef[2,1]*x+summary(curve.nlslrc_GP)$coef[1,1])^2-4*summary(curve.nlslrc_GP)$coef[2,1]*summary(curve.nlslrc_GP)$coef[4,1]*summary(curve.nlslrc_GP)$coef[1,1]*x))-summary(curve.nlslrc_GP)$coef[3,1],lwd=2,col="blue",add=T)
abline(h=0)
abline(v=0)
axis(1, pos = 0)
axis(2, pos = 0)
mtext(expression("PPFD ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=1,cex=1.5)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=1,cex=1.5)
legend(100,31, legend = c("BM", "RR", "GP"), title= "Grass species",c("green", "red", "blue"),
       col=c(curve.nlslrc_BM, curve.nlslrc_RR, curve.nlslrc_GP), lty=1:2, cex=1.0,
       box.lty=0)

CO2_all_photo<- c(BM_photo,RR_photo, GP_photo)
CO2_all_photo
CO2_all_R<-c(BM_CO2R, RR_CO2R,GP_CO2R)

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(CO2_all_R,CO2_all_photo ,xaxt = "n", yaxt = "n",frame.plot = FALSE, xlab="", xlim = c(0,max(CO2_all_R)+2), ylab="", ylim=c(-10,max(CO2_all_photo)+10),cex.lab=1.2,cex.axis=1,cex=1.5)
arrows(BM_CO2R, BM_photo-BM_SE_Photo, BM_CO2R, BM_photo+BM_SE_Photo, length = 0.05, angle = 90, code = 3)
curve((1/(2*summary(curve.nlslrc_BMc)$coef[4,1]))*(summary(curve.nlslrc_BMc)$coef[2,1]*x+summary(curve.nlslrc_BMc)$coef[1,1]-sqrt((summary(curve.nlslrc_BMc)$coef[2,1]*x+summary(curve.nlslrc_BMc)$coef[1,1])^2-4*summary(curve.nlslrc_BMc)$coef[2,1]*summary(curve.nlslrc_BMc)$coef[4,1]*summary(curve.nlslrc_BMc)$coef[1,1]*x))-summary(curve.nlslrc_BMc)$coef[3,1],lwd=2,col="orange",add=T)


arrows(GP_CO2R, GP_photo-GP_SE_Photo, GP_CO2R, GP_photo+GP_SE_Photo, length = 0.05, angle = 90, code = 3)
curve((1/(2*summary(curve.nlslrc_GPc)$coef[4,1]))*(summary(curve.nlslrc_GPc)$coef[2,1]*x+summary(curve.nlslrc_GPc)$coef[1,1]-sqrt((summary(curve.nlslrc_GPc)$coef[2,1]*x+summary(curve.nlslrc_GPc)$coef[1,1])^2-4*summary(curve.nlslrc_GPc)$coef[2,1]*summary(curve.nlslrc_GPc)$coef[4,1]*summary(curve.nlslrc_GPc)$coef[1,1]*x))-summary(curve.nlslrc_GPc)$coef[3,1],lwd=2,col="purple",add=T)

arrows(RR_CO2R, RR_photo-RR_SE_Photo, RR_CO2R, RR_photo+RR_SE_Photo, length = 0.05, angle = 90, code = 3)
curve((1/(2*summary(curve.nlslrc_RRc)$coef[4,1]))*(summary(curve.nlslrc_RRc)$coef[2,1]*x+summary(curve.nlslrc_RRc)$coef[1,1]-sqrt((summary(curve.nlslrc_RRc)$coef[2,1]*x+summary(curve.nlslrc_RRc)$coef[1,1])^2-4*summary(curve.nlslrc_RRc)$coef[2,1]*summary(curve.nlslrc_RRc)$coef[4,1]*summary(curve.nlslrc_RRc)$coef[1,1]*x))-summary(curve.nlslrc_RRc)$coef[3,1],lwd=2,col="brown",add=T)


abline(h=0)
abline(v=0)
axis(1, pos = 0)
axis(2, pos = 0)
mtext(expression(CO [2]~ ppm),side=1,line=2.0,cex=1.5)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1.5)
legend(95,45, legend = c("BM", "GP", "RR"), title= "Grass species",c("orange", "purple", "brown"),
       col=c(curve.nlslrc_BMc, curve.nlslrc_RRc, curve.nlslrc_GPc), lty=1:2, cex=1,
       box.lty=0)

#Temperature response function#

attach(Tresponse)

Tfunc<- nls(GPTfunc~max(GPTfunc)*exp((Tempfunc-TO)/B), start = list(TO=30, B=5))

TfunctionG<- nls(GPTfunc~ (Tempfunc-min(Tempfunc))/(25-min(Tempfunc))^Q*((max(Tempfunc)-Tempfunc)/(max(Tempfunc)-25)), start = list(Q=1))
TfunctionG

TfunctionGPa<- nls(GPTfunc~ (Tempfunc-Tmin)/(25- Tmin)^Q*((Tmax-Tempfunc)/(Tmax-25)), start = list(Tmin=Tmax=GPTfunc=0, Q=1))
TfunctionGPa

topt<-(min(Tempfunc)+2.02*max(Tempfunc))/(1+2.02)
topt
GPTresp<-nls(GPTfunc~ ((Tempfunc- min(Tempfunc))/(25 -min(Tempfunc))^q*(((1+q)*Topt-min(Tempfunc)-q*25))/((1+q)*Topt-(min(Tempfunc)-q*25))), start = list(q=1, Topt=max(GPTfunc)))


