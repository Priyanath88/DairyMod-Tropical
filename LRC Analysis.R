attach(lightmean)
omit_lightmean<-lightmean[-c(15), ]
omit_lightmean
##Brahchiaria grass##

library(dplyr)
BB_light<-filter(lightmean, species=="Brachiaria")
BB_light
View(BB_light)

PARlrcBM<-BB_light$avgPARi #PAR (aka PPFD or Q)
photolrcBM<-BB_light$avgphoto #net photosynthetic rate (Anet)
SEphoto<-BB_light$stdPhoto

curvelrcBM<-data.frame(PARlrcBM, photolrcBM,SEphoto)
curvelrcBM # *inspect raw data and check notebook (data reasonable or need edited/discarded?)

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(PARlrcBM,photolrcBM,xlab="PAR", ylab="Assimilation")

mtext(expression("PPFD ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1.5)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1.5)

curve.nlslrc_BM = nls(photolrcBM ~ (1/(2*theta))*(AQY*PARlrcBM+Am-sqrt((AQY*PARlrcBM+Am)^2- 4*AQY*theta*Am*PARlrcBM))-Rd, start=list(Am= (max(photolrcBM)- min(photolrcBM)),AQY=0.05, Rd=-min(photolrcBM),theta=1))   

summary(curve.nlslrc_BM)

plot<- par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(PARlrcBM,photolrcBM ,xaxt = "n", yaxt = "n",frame.plot = FALSE, xlab="", xlim = c(-500,max(PARlrcBM)+2), ylab="", ylim=c(-5,max(photolrcBM)+10),cex.lab=1.5,cex.axis=3,cex=1.5)
arrows(PARlrcBM, photolrcBM-SEphoto, PARlrcBM, photolrcBM+SEphoto, length = 0.05, angle = 90, code = 3)

mtext(expression("PPFD ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=2.5,cex=1.5)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1.5)
curve((1/(2*summary(curve.nlslrc_BM)$coef[4,1]))*(summary(curve.nlslrc_BM)$coef[2,1]*x+summary(curve.nlslrc_BM)$coef[1,1]-sqrt((summary(curve.nlslrc_BM)$coef[2,1]*x+summary(curve.nlslrc_BM)$coef[1,1])^2-4*summary(curve.nlslrc_BM)$coef[2,1]*summary(curve.nlslrc_BM)$coef[4,1]*summary(curve.nlslrc_BM)$coef[1,1]*x))-summary(curve.nlslrc_BM)$coef[3,1],lwd=2,col="Orange",add=T)
abline(h=0)
abline(v=0)
axis(1, pos = 0,cex.axis=1)
axis(2, pos = 0, cex.axis=1)
plot
box(col= "black")

# ---Solve for light compensation point (LCPT), PPFD where Anet=0 ---
x<-function(x) {(1/(2*summary(curve.nlslrc_BM)$coef[4,1]))*(summary(curve.nlslrc_BM)$coef[2,1]*x+summary(curve.nlslrc_BM)$coef[1,1]-sqrt((summary(curve.nlslrc_BM)$coef[2,1]*x+summary(curve.nlslrc_BM)$coef[1,1])^2-4*summary(curve.nlslrc_BM)$coef[2,1]*summary(curve.nlslrc_BM)$coef[4,1]*summary(curve.nlslrc_BM)$coef[1,1]*x))-summary(curve.nlslrc_BM)$coef[3,1]}

uniroot(x,c(0,50))$root #LCPT


# ---Solve for light saturation point (LSP), PPFD where 75% of Amax is achieved (75% is arbitrary - cutoff could be changed)
x<-function(x) {(1/(2*summary(curve.nlslrc_BM)$coef[4,1]))*(summary(curve.nlslrc_BM)$coef[2,1]*x+summary(curve.nlslrc_BM)$coef[1,1]-sqrt((summary(curve.nlslrc_BM)$coef[2,1]*x+summary(curve.nlslrc_BM)$coef[1,1])^2-4*summary(curve.nlslrc_BM)$coef[2,1]*summary(curve.nlslrc_BM)$coef[4,1]*summary(curve.nlslrc_BM)$coef[1,1]*x))-summary(curve.nlslrc_BM)$coef[3,1]-(0.90*summary(curve.nlslrc_BM)$coef[1,1])+0.90*(summary(curve.nlslrc_BM)$coef[3,1])}

uniroot(x,c(0,1500))$root #LSP 


#Rhodes Grass###
RR_light<-filter(lightmean, species=="Rhodes grass")
RR_light
View(RR_light)

PARlrcRR<-RR_light$avgPARi #PAR (aka PPFD or Q)
photolrcRR<-RR_light$avgphoto #net photosynthetic rate (Anet)
SEphotoRR<-RR_light$stdPhoto

curvelrcRR<-data.frame(PARlrcRR,photolrcRR,SEphotoRR)
curvelrcRR

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(PARlrcRR,photolrcRR,xlab="PAR", ylab="Assimilation")

mtext(expression("PPFD ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1.5)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1.5)

curve.nlslrc_RR = nls(photolrcRR ~ (1/(2*theta))*(AQY*PARlrcRR+Am-sqrt((AQY*PARlrcRR+Am)^2- 4*AQY*theta*Am*PARlrcRR))-Rd, start=list(Am= (max(photolrcRR)- min(photolrcRR)),AQY=0.05, Rd=-min(photolrcRR),theta=1))   

summary(curve.nlslrc_RR)

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(PARlrcRR,photolrcRR ,xaxt = "n", yaxt = "n",frame.plot = FALSE, xlab="", xlim = c(-500,max(PARlrcBM)+2), ylab="", ylim=c(-5,max(photolrcBM)+10),cex.lab=1.2,cex.axis=1,cex=1.5)
arrows(PARlrcRR, photolrcRR-SEphotoRR, PARlrcRR, photolrcRR+SEphotoRR, length = 0.05, angle = 90, code = 3)

mtext(expression("PPFD ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=2.5,cex=1.5)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1.5)
curve((1/(2*summary(curve.nlslrc_RR)$coef[4,1]))*(summary(curve.nlslrc_RR)$coef[2,1]*x+summary(curve.nlslrc_RR)$coef[1,1]-sqrt((summary(curve.nlslrc_RR)$coef[2,1]*x+summary(curve.nlslrc_RR)$coef[1,1])^2-4*summary(curve.nlslrc_RR)$coef[2,1]*summary(curve.nlslrc_RR)$coef[4,1]*summary(curve.nlslrc_RR)$coef[1,1]*x))-summary(curve.nlslrc_RR)$coef[3,1],lwd=2,col="Blue",add=T)
abline(h=0)
abline(v=0)
axis(1, pos = 0)
axis(2, pos = 0)
box(col= "black")

# ---Solve for light compensation point (LCPT), PPFD where Anet=0 ---
x<-function(x) {(1/(2*summary(curve.nlslrc_RR)$coef[4,1]))*(summary(curve.nlslrc_RR)$coef[2,1]*x+summary(curve.nlslrc_RR)$coef[1,1]-sqrt((summary(curve.nlslrc_RR)$coef[2,1]*x+summary(curve.nlslrc_RR)$coef[1,1])^2-4*summary(curve.nlslrc_RR)$coef[2,1]*summary(curve.nlslrc_RR)$coef[4,1]*summary(curve.nlslrc_RR)$coef[1,1]*x))-summary(curve.nlslrc_RR)$coef[3,1]}

uniroot(x,c(0,50))$root #LCPT

# ---Solve for light saturation point (LSP), PPFD where 75% of Amax is achieved (75% is arbitrary - cutoff could be changed)
x<-function(x) {(1/(2*summary(curve.nlslrc_RR)$coef[4,1]))*(summary(curve.nlslrc_RR)$coef[2,1]*x+summary(curve.nlslrc_RR)$coef[1,1]-sqrt((summary(curve.nlslrc_RR)$coef[2,1]*x+summary(curve.nlslrc_RR)$coef[1,1])^2-4*summary(curve.nlslrc_RR)$coef[2,1]*summary(curve.nlslrc_RR)$coef[4,1]*summary(curve.nlslrc_RR)$coef[1,1]*x))-summary(curve.nlslrc_RR)$coef[3,1]-(0.90*summary(curve.nlslrc_RR)$coef[1,1])+0.90*(summary(curve.nlslrc_RR)$coef[3,1])}

uniroot(x,c(0,1500))$root #LSP 


#Gatton panic##

GP_light<-filter(lightmean, species=="Gatton pannic")
GP_light
View(GP_light)
GP_light_omit<-GP_light[-c(5), ]  #remove one outlier in raw#5
GP_light_omit
PARlrcGP<-GP_light_omit$avgPARi #PAR (aka PPFD or Q)
photolrcGP<-GP_light_omit$avgphoto #net photosynthetic rate (Anet)
SEphotoGP<-GP_light_omit$stdPhoto

curvelrcGP<-data.frame(PARlrcGP,photolrcGP,SEphotoGP)
curvelrcGP


par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(PARlrcGP,photolrcGP,xlab="PAR", ylab="Assimilation")

mtext(expression("PPFD ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1.5)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1.5)

curve.nlslrc_GP = nls(photolrcGP ~ (1/(2*theta))*(AQY*PARlrcGP+Am-sqrt((AQY*PARlrcGP+Am)^2- 4*AQY*theta*Am*PARlrcGP))-Rd, start=list(Am= (max(photolrcGP)- min(photolrcGP)),AQY=0.05, Rd=-min(photolrcGP),theta=1))   

summary(curve.nlslrc_GP)

par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(PARlrcGP,photolrcGP ,xaxt = "n", yaxt = "n",frame.plot = FALSE, xlab="", xlim = c(-500,max(PARlrcGP)+2), ylab="", ylim=c(-5,max(photolrcGP)+10),cex.lab=1.2,cex.axis=1,cex=1.5)
arrows(PARlrcGP, photolrcGP-SEphotoGP, PARlrcGP, photolrcGP+SEphotoGP, length = 0.05, angle = 90, code = 3)

mtext(expression("PPFD ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=2.5,cex=1.5)
mtext(expression(A[net]*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2.5,cex=1.5)
curve((1/(2*summary(curve.nlslrc_GP)$coef[4,1]))*(summary(curve.nlslrc_GP)$coef[2,1]*x+summary(curve.nlslrc_GP)$coef[1,1]-sqrt((summary(curve.nlslrc_GP)$coef[2,1]*x+summary(curve.nlslrc_GP)$coef[1,1])^2-4*summary(curve.nlslrc_GP)$coef[2,1]*summary(curve.nlslrc_GP)$coef[4,1]*summary(curve.nlslrc_GP)$coef[1,1]*x))-summary(curve.nlslrc_GP)$coef[3,1],lwd=2,col="Red",add=T)
abline(h=0)
abline(v=0)
axis(1, pos = 0)
axis(2, pos = 0)
box(col= "black")

# ---Solve for light compensation point (LCPT), PPFD where Anet=0 ---
x<-function(x) {(1/(2*summary(curve.nlslrc_GP)$coef[4,1]))*(summary(curve.nlslrc_GP)$coef[2,1]*x+summary(curve.nlslrc_GP)$coef[1,1]-sqrt((summary(curve.nlslrc_GP)$coef[2,1]*x+summary(curve.nlslrc_GP)$coef[1,1])^2-4*summary(curve.nlslrc_GP)$coef[2,1]*summary(curve.nlslrc_GP)$coef[4,1]*summary(curve.nlslrc_GP)$coef[1,1]*x))-summary(curve.nlslrc_GP)$coef[3,1]}

uniroot(x,c(0,50))$root #LCPT

# ---Solve for light saturation point (LSP), PPFD where 75% of Amax is achieved (75% is arbitrary - cutoff could be changed)
x<-function(x) {(1/(2*summary(curve.nlslrc_GP)$coef[4,1]))*(summary(curve.nlslrc_GP)$coef[2,1]*x+summary(curve.nlslrc_GP)$coef[1,1]-sqrt((summary(curve.nlslrc_GP)$coef[2,1]*x+summary(curve.nlslrc_GP)$coef[1,1])^2-4*summary(curve.nlslrc_GP)$coef[2,1]*summary(curve.nlslrc_GP)$coef[4,1]*summary(curve.nlslrc_GP)$coef[1,1]*x))-summary(curve.nlslrc_GP)$coef[3,1]-(0.90*summary(curve.nlslrc_GP)$coef[1,1])+0.90*(summary(curve.nlslrc_GP)$coef[3,1])}

uniroot(x,c(0,1600))$root #LSP 


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

