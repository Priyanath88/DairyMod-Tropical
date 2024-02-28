# November 2013, JM Heberling (jmheberling@gmail.com)
# Modeling CO2 response curves (A/Ci)
# Farquhar-von Caemmerer-Berry (FvCB) model (1980) Planta 149:78-90
# Simultaneous estimation detailed in Dubois et al. (2007) New Phyt 176:402-414
#---

 # Read in text file from Licor 6400
library(readxl)
aciGRS<- read_excel("CO2Data.xlsx")

acicurve<-filter(CO2Data, species=="Brachiaria")
acicurve
aci_RCi<-acicurve$Ci
aci_RCiPa<-acicurve$CiPa
aci_Rphoto<-acicurve$Photo
aci_RCO2R<-acicurve$CO2R
aci_rep<-acicurve$Repblock  

aci_Brachiaria<-data.frame(aci_RCi,aci_RCiPa,aci_Rphoto,aci_RCO2R,aci_rep)
acicurve=data.table(acicurve)

aci_brachiaria_output=acicurve[,. ("meanPhoto"=mean(Photo),
                "std"=std.error(Photo)),
                by=c("CO2 setting")]
aci_brachiaria_output
aci_brachiaria_outputCi_pa=acicurve[,. ("meanCipa"=mean(CiPa),
                                   "std"=std.error(CiPa)),
                               by=c("CO2 setting")]
aci_brachiaria_outputCi_pa
write.csv(aci_brachiaria_outputCi_pa,"Brachiaria_cipa.csv")
write.csv(aci_brachiaria_output, "acibrachiaria.csv")
aci_brachiaria_outputC02R=acicurve[,. ("meanCO2R"=mean(CO2R),
                                        "std"=std.error(CO2R)),
                                    by=c("CO2 setting")]
aci_brachiaria_outputC02R
aci_brachiaria_outputCi=acicurve[,. ("meanCi"=mean(Ci),
                                       "std"=std.error(Ci)),
                                   by=c("CO2 setting")]
aci_brachiaria_outputCi
aci_brachiaria_outputTleaf=acicurve[,. ("meanleafT"=mean(Tleaf),
                                     "std"=std.error(Tleaf)),
                                 by=c("CO2 setting")]
aci_brachiaria_outputTleaf
Aci_B<-data.frame(aci_brachiaria_output,aci_brachiaria_outputCi,aci_brachiaria_outputCi_pa,aci_brachiaria_outputC02R,aci_brachiaria_outputTleaf)
Aci_B
write.csv(Aci_B,"aci_b.csv")
read.csv("aci_b.csv")
Brachiaria_acicurve<- read_excel("Brachiaria_acicurve.xlsx")
View(Brachiaria_acicurve)
attach(Brachiaria_acicurve)

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(meanCipa,meanPhoto,ylab="", xlab="",cex.lab=1.2,cex.axis=1.5,cex=2)
mtext(expression("Intercellular "*CO[2]*" Pressure (Pa)"),side=1,line=3.3,cex=1.5)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1.5)


# ---Temperature adjusted coefficients:

# Constants published in Sharkey et al (2007) Plant Cell Env 30: 1035-1040 
# Measured using transgenic tobacco (ASSUMED to be similar across higher plants)
# Ci units in Pa; Sharkey et al (2007) recommend partial pressures
# **Be sure units are correct for your input data** (Ci is in Pa or ppm?)

R=0.008314 #(kJ mol^-1 K^-1)
Kc=exp(35.9774-80.99/(R*(Brachiaria_acicurve$meanleafT+273.15))) #Michaelis-Menten constant for Rubisco for O2 (Pa)
Ko=exp(12.3772-23.72/(R*(Brachiaria_acicurve$meanleafT+273.15))) #Michaelis-Menten constant for Rubisco for CO2 (kPa)
GammaStar=exp(11.187-24.46/(R*(Brachiaria_acicurve$meanleafT+273.15))) #Photorespiration compensation point (Pa)
O=21 #oxygen (O2) partial pressure (kPa)

# Alternative constants from Bernacchi et al. (2001) Plant Cell Env 24: 253-259 
# Ci units in ppm (but can be converted to Pa by atmospheric pressure)
#R=0.008314 #(kJ mol^-1 K^-1)
#aci$Kc=exp(38.05-79.43/(R*(aci$Tleaf+273.15))) #umol mol-1
#aci$Ko=exp(20.30-36.38/(R*(aci$Tleaf+273.15))) #mmol mol-1
#aci$GammaStar=exp(19.02-37.83/(R*(aci$Tleaf+273.15))) #umol mol-1
#O=210

# ---RuBisCO limited portion---
#(Vcmax*(Ci_Pa-GammaStar))/(Ci_Pa+(Kc*(1+(O/Ko)))))-Rd
	
# ---RUBP limited portion---
#((J*(Ci_Pa-GammaStar))/((4*Ci_Pa)+(8*GammaStar)))-Rd
	
# ---TPU.limited portion 
#3*TPU-Rd
# Few studies find triose phosphate limitation (under natural conditions), but could easily be added to model below to test (but with loss of statistical power); it is often not considered, but depends on the dataset
# Miao et al (2009) Plant, Cell, Env 32:109-122 recommend fitting when possible or removing TPU limited points from dataset before fitting


library(minpack.lm)
# Simultaneous estimation method described by Dubois et al. 2007 New Phyt 176:402-414
# Could change optimization algorithm (default here is Gauss-Newton)
# Could also do a "grid search" if estimates are sensitive to starting values

Brahiaria_aci.fit<-nls(meanPhoto~ifelse(((Vcmax*(meanCipa-GammaStar))/(meanCipa+(Kc*(1+(O/Ko)))))<((J*(meanCipa-GammaStar))/((4*meanCipa)+(8*GammaStar))),((Vcmax*(meanCipa-GammaStar))/(meanCipa+(Kc*(1+(O/Ko))))),((J*(meanCipa-GammaStar))/((4*meanCipa)+(8*GammaStar))))-Rd,start=list(Vcmax=50,J=100,Rd=5),data=Brachiaria_acicurve) #if error: reconsider starting values, bad dataset? (too few points or response curve not clear)

summary(Brahiaria_aci.fit)

Vcmax<-summary(Brahiaria_aci.fit)$coef[1,1]
Vcmax
J<-summary(Brahiaria_aci.fit)$coef[2,1]
Rd<-summary(Brahiaria_aci.fit)$coef[3,1]

# ---Graph raw data with modeled curve---

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(meanCipa,meanPhoto,ylab="", xlab="",cex.lab=1.2,cex.axis=1.5,cex=2)
mtext(expression("Intercellular "*CO[2]*" Pressure (Pa)"),side=1,line=3.3,cex=1.5)
mtext(expression("Net photosynthetic rate  (umol  "* CO[2]*   m^-2*   s^-1*")"),side=2,line=2.5,cex=1.5)

curve(ifelse(((Vcmax*(x-mean(GammaStar)))/(x+(mean(Kc)*(1+(O/mean(Ko))))))<((J*(x-mean(GammaStar)))/((4*x)+(8*mean(GammaStar)))),((Vcmax*(x-mean(GammaStar)))/(x+(mean(Kc)*(1+(O/mean(Ko)))))),((J*(x-mean(GammaStar)))/((4*x)+(8*mean(GammaStar)))))-Rd,col="green",add=T) #Reasonable fit? Could check goodness of fit, model assumptions

# ---Frequently used alternative: "disjunct segment estimation" (sensu Dubois et al 2007)---
# See Sharkey et al (2007) Plant Cell Env 30: 1035-1040
# This method entails selecting which points on the curve are limited by rubisco, RUBP, or TPU a priori (or subjectively). The data is subsetted by these cutoff points and each segment modeled separately.


# ---Considering mesophyll conductance (gm) ---

# The original Farquhar-von Caemmerer-Berry (1980) model was developed for Cc (chlorplastic CO2 concentration) not Ci (intercellular CO2 concentration).  Ci equals Cc, assuming mesophyll conductance (gm) is infinite, which might be an unfair assumption in many cases.  See reviews on estimating mesophyll conductance and the effect on A/Cc versus A/Ci modeling (eg.,  Niinemets et al, 2009 J Exp Bot 60:2271-2282; Flexas et al., 2013 J Exp Bot 64:3965-3981)

# Cc = Ci - A/gm

# ---Revised model with gm as a fitted parameter (modified "Ethier method") --- 
#  See Ethier & Livingston (2004) Plant, Cell, Env 27:137-153
#  Some authors highlight the need to fit gm (e.g., Niinemets et al, 2009 J Exp Bot 60:2271-2282) while others suggest gm estimation based on gas exchange data curve fitting alone is not preferred (Pons et al. 2009 J Exp Bot 60:2217-2234)
#  Techniques for measuring gm are still under development (see Flexas et al., 2013 J Exp Bot 64:3965-3981)
#  Ideally, gm would be independently measured and Cc would be known

aci.fit<-nlsLM(Photo~ifelse(((Vcmax*((Ci_Pa-(Photo*gmInv))-GammaStar))/((Ci_Pa-(Photo*gmInv))+(Kc*(1+(O/Ko)))))<((J*((Ci_Pa-(Photo*gmInv))-GammaStar))/((4*(Ci_Pa-(Photo*gmInv)))+(8*GammaStar))),((Vcmax*((Ci_Pa-(Photo*gmInv))-GammaStar))/((Ci_Pa-(Photo*gmInv))+(Kc*(1+(O/Ko))))),((J*((Ci_Pa-(Photo*gmInv))-GammaStar))/((4*(Ci_Pa-(Photo*gmInv)))+(8*GammaStar))))-Rd,start=list(Vcmax=50,J=100,Rd=0.5,gmInv=0),data=aci) #gm cannot be negative; gm=1/gmInv; fitting gm can be statistically difficult for many datasets
aci.fit
# ---Other notes: ---
# Miao et al (2009) Plant, Cell, Env 32:109-122 review A/Cc fitting methods and recommend fitting full model (including gm and TPU) combining grid search and two stage nonlinear least square regression 
# Dubois et al (2007) do not fit gm and warn it may not be reliable. However, if decided to model, they suggest fitting 1/gm if you do since gm can approach zero; They also suggest using a grid search to avoid local minima (sensitivity to initial starting values).  I have not found this to be an issue, but a grid search could be incorporated in this code.

# use Ci value instead Cipa##

Brahiaria_aci.fit1<-nls(meanPhoto~ifelse(((Vcmax*(meanCi-GammaStar))/(meanCi+(Kc*(1+(O/Ko)))))<((J*(meanCi-GammaStar))/((4*meanCi)+(8*GammaStar))),((Vcmax*(meanCi-GammaStar))/(meanCi+(Kc*(1+(O/Ko))))),((J*(meanCi-GammaStar))/((4*meanCi)+(8*GammaStar))))-Rd,start=list(Vcmax=50,J=100,Rd=5),data=Brachiaria_acicurve) #if error: reconsider starting values, bad dataset? (too few points or response curve not clear)

summary(Brahiaria_aci.fit1)

Vcmax<-summary(Brahiaria_aci.fit1)$coef[1,1]
Vcmax
J<-summary(Brahiaria_aci.fit1)$coef[2,1]
Rd<-summary(Brahiaria_aci.fit1)$coef[3,1]

# ---Graph raw data with modeled curve---

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(meanCi,meanPhoto,xlab="",xlim = c(0,max(meanCi)+5),ylab="", ylim=c(0,max(meanPhoto)+10), cex.lab=1.2,cex.axis=1.5,cex=2)
arrows(meanCi, meanPhoto-Photostd, meanCi, meanPhoto+Photostd, length = 0.05, angle = 90, code = 3)
mtext(expression("Ci (ppm)"),side=1,line=3.3,cex=1.5)
mtext(expression("Net photosynthetic rate  (umol  "* CO[2]*   m^-2*   s^-1*")"),side=2,line=2.5,cex=1.5)

curve(ifelse(((Vcmax*(x-mean(GammaStar)))/(x+(mean(Kc)*(1+(O/mean(Ko))))))<((J*(x-mean(GammaStar)))/((4*x)+(8*mean(GammaStar)))),((Vcmax*(x-mean(GammaStar)))/(x+(mean(Kc)*(1+(O/mean(Ko)))))),((J*(x-mean(GammaStar)))/((4*x)+(8*mean(GammaStar)))))-Rd,col="green",add=T) #Reasonable fit? Could check goodness of fit, model assumptions


aci.fit1<-nls(meanPhoto~ifelse(((Vcmax*((meanCi-(meanPhoto*gmInv))-GammaStar))/((meanCi-(meanPhoto*gmInv))+(Kc*(1+(O/Ko)))))<((J*((meanCi-(meanPhoto*gmInv))-GammaStar))/((4*(meanCi-(meanPhoto*gmInv)))+(8*GammaStar))),((Vcmax*((meanCi-(meanPhoto*gmInv))-GammaStar))/((meanCi-(meanPhoto*gmInv))+(Kc*(1+(O/Ko))))),((J*((meanCi-(meanPhoto*gmInv))-GammaStar))/((4*(meanCi-(Photo*gmInv)))+(8*GammaStar))))-Rd,start=list(Vcmax=50,J=100,Rd=0.5,gmInv=0),data=Brachiaria_acicurve) #gm cannot be negative; gm=1/gmInv; fitting gm can be statistically difficult for many datasets
aci.fit1
library(ggplot2)
ggplot(Brachiaria_acicurve, aes(x=meanCi, y=meanPhoto))+ geom_pointrange(aes(ymin=meanPhoto- Photostd, ymax=meanPhoto+ Photostd))+
  geom_smooth(se=FALSE, method = "nls", formula = y~x(log(x)))

attach(B_acibook)

alldata_ci_B<- B_acibook$Ci
alldata_cipa_B<-B_acibook$CiPa
alldata_CO2R<-B_acibook$CO2R
alldata_photo_B<-B_acibook$Photo
alldata_Tleaf_B<-B_acibook$Tleaf
Alldata_B_aci<-data.frame(alldata_ci_B,alldata_cipa_B,alldata_CO2R,alldata_photo_B,alldata_Tleaf_B)

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(alldata_ci_B,alldata_photo_B,ylab="", xlab="",cex.lab=1.2,cex.axis=1.5,cex=2)
mtext(expression("Intercellular "*CO[2]*" Pressure (Pa)"),side=1,line=3.3,cex=1.5)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1.5)

R=0.008314 #(kJ mol^-1 K^-1)
Kc=exp(35.9774-80.99/(R*(Alldata_B_aci$alldata_Tleaf_B+273.15))) #Michaelis-Menten constant for Rubisco for O2 (Pa)
Ko=exp(12.3772-23.72/(R*(Alldata_B_aci$alldata_Tleaf_B+273.15))) #Michaelis-Menten constant for Rubisco for CO2 (kPa)
GammaStar=exp(11.187-24.46/(R*(Alldata_B_aci$alldata_Tleaf_B+273.15))) #Photorespiration compensation point (Pa)
O=21

Brahiaria_aci.fit3<-nls(alldata_photo_B~ifelse(((Vcmax*(alldata_ci_B-GammaStar))/(alldata_ci_B+(Kc*(1+(O/Ko)))))<((J*(alldata_ci_B-GammaStar))/((4*alldata_ci_B)+(8*GammaStar))),((Vcmax*(alldata_ci_B-GammaStar))/(alldata_ci_B+(Kc*(1+(O/Ko))))),((J*(alldata_ci_B-GammaStar))/((4*alldata_ci_B)+(8*GammaStar))))-Rd,start=list(Vcmax=50,J=100,Rd=5),data=Alldata_B_aci) #if error: reconsider starting values, bad dataset? (too few points or response curve not clear)

summary(Brahiaria_aci.fit3)

Vcmax<-summary(Brahiaria_aci.fit3)$coef[1,1]
Vcmax
J<-summary(Brahiaria_aci.fit3)$coef[2,1]
J
Rd<-summary(Brahiaria_aci.fit3)$coef[3,1]
Rd

# ---Graph raw data with modeled curve---

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(alldata_ci_B,alldata_photo_B,ylab="", xlab="",cex.lab=1.2,cex.axis=1.5,cex=2)
mtext(expression("Ci (ppm)"),side=1,line=3.3,cex=1.5)
mtext(expression("Net photosynthetic rate  (umol  "* CO[2]*   m^-2*   s^-1*")"),side=2,line=2.5,cex=1.5)

curve(ifelse(((Vcmax*(x-mean(GammaStar)))/(x+(mean(Kc)*(1+(O/mean(Ko))))))<((J*(x-mean(GammaStar)))/((4*x)+(8*mean(GammaStar)))),((Vcmax*(x-mean(GammaStar)))/(x+(mean(Kc)*(1+(O/mean(Ko)))))),((J*(x-mean(GammaStar)))/((4*x)+(8*mean(GammaStar)))))-Rd,col="light green",add=T) #Reasonable fit? Could check goodness of fit, model assumptions

####Rhodes grass ACi curves analysis

attach(R_acibook)
Rhodegrassaci<-data.frame(R_acibook$Ci, R_acibook$Photo,R_acibook$CiPa, R_acibook$CO2R, R_acibook$Tleaf, R_acibook$`CO2 setting`)
Rhodegrassaci
RR_ci<-R_acibook$Ci
RR_cipa<-R_acibook$CiPa
RR_photo<-R_acibook$Photo
R
R_acibook=data.table(R_acibook)

aci_Rhodes_Ci=R_acibook[,. ("meanCi"=mean(Ci),
                            "std"=std.error(Ci)),
                        by=c("CO2 setting")]
aci_Rhodes_Ci



aci_Rhodes_photo=R_acibook[,. ("meanPhotoR"=mean(Photo),
                                   "std"=std.error(Photo)),
                            by=c("CO2 setting")]
aci_Rhodes_photo                

aci_Rhodes_ci_pa=R_acibook[,. ("meanCipa"=mean(CiPa),
                                        "std"=std.error(CiPa)),
                                    by=c("CO2 setting")]
aci_Rhodes_ci_pa

aci_Rhodes_CO2R=R_acibook[,. ("meanCO2R"=mean(CO2R),
                                       "std"=std.error(CO2R)),
                                   by=c("CO2 setting")]
aci_Rhodes_CO2R



aci_Rhodes_Tleaf=R_acibook[,. ("meanleafT"=mean(Tleaf),
                                        "std"=std.error(Tleaf)),
                                    by=c("CO2 setting")]
aci_Rhodes_Tleaf

Aci_RR<-data.frame(aci_Rhodes_Ci,aci_Rhodes_ci_pa,aci_Rhodes_photo,aci_Rhodes_Tleaf,aci_Rhodes_CO2R)
Aci_RR
write.csv(Aci_RR,"aci_RR1.csv")
read.csv("aci_RR.csv")
attach(aci_RR)
par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(mean_CiR,meanPhotoR,ylab="", xlab="",cex.lab=1.2,cex.axis=1.5,cex=2)
mtext(expression("Intercellular "*CO[2]*" Pressure (Pa)"),side=1,line=3.3,cex=1.5)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1.5)

R=0.008314 #(kJ mol^-1 K^-1)
Kc=exp(35.9774-80.99/(R*(aci_RR$TleafstdR+273.15))) #Michaelis-Menten constant for Rubisco for O2 (Pa)
Ko=exp(12.3772-23.72/(R*(aci_RR$TleafstdR+273.15))) #Michaelis-Menten constant for Rubisco for CO2 (kPa)
GammaStar=exp(11.187-24.46/(R*(aci_RR$TleafstdR+273.15))) #Photorespiration compensation point (Pa)
O=21

RR_aci.fit<-nls(meanPhotoR~ifelse(((Vcmax*(mean_CiR-GammaStar))/(mean_CiR+(Kc*(1+(O/Ko)))))<((J*(mean_CiR-GammaStar))/((4*mean_CiR)+(8*GammaStar))),((Vcmax*(mean_CiR-GammaStar))/(mean_CiR+(Kc*(1+(O/Ko))))),((J*(mean_CiR-GammaStar))/((4*mean_CiR)+(8*GammaStar))))-Rd,start=list(Vcmax=40,J=90,Rd=0.5),data=aci_RR) #if error: reconsider starting values, bad dataset? (too few points or response curve not clear)

summary(RR_aci.fit)

Vcmax<-summary(Brahiaria_aci.fit3)$coef[1,1]
Vcmax
J<-summary(Brahiaria_aci.fit3)$coef[2,1]
J
Rd<-summary(Brahiaria_aci.fit3)$coef[3,1]
Rd

# ---Graph raw data with modeled curve---

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(alldata_ci_B,alldata_photo_B,ylab="", xlab="",cex.lab=1.2,cex.axis=1.5,cex=2)
mtext(expression("Ci (ppm)"),side=1,line=3.3,cex=1.5)
mtext(expression("Net photosynthetic rate  (umol  "* CO[2]*   m^-2*   s^-1*")"),side=2,line=2.5,cex=1.5)

curve(ifelse(((Vcmax*(x-mean(GammaStar)))/(x+(mean(Kc)*(1+(O/mean(Ko))))))<((J*(x-mean(GammaStar)))/((4*x)+(8*mean(GammaStar)))),((Vcmax*(x-mean(GammaStar)))/(x+(mean(Kc)*(1+(O/mean(Ko)))))),((J*(x-mean(GammaStar)))/((4*x)+(8*mean(GammaStar)))))-Rd,col="light green",add=T) #Reasonable fit? Could check goodness of fit, model assumptions

