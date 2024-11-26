attach(CO2response) #attach the CO2 data file

CO2fun_R<-CO2response$CO2R
co2fun<-CO2response$CO2response
co2function<-data.frame(CO2fun_R,co2fun)
co2function
par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(CO2fun_R,co2fun,xlab="", ylab="")

mtext(expression(CO[2] ~ppm ),side=1,line=3.3,cex=1.5)
mtext(expression(CO[2] ~ "fc (C)"),side=2,line=2.5,cex=1.5)

curve.CO2fun_RR = nls(co2fun ~ (1/(2*theta))*(AQY*CO2fun_R+Am-sqrt((AQY*CO2fun_R+Am)^2- 4*AQY*theta*Am*CO2fun_R)), start = list(Am= (max(co2fun)- min(co2fun)),AQY=0.01, theta= (0<= theta <= 1)))   

summary(curve.CO2fun_RR)
