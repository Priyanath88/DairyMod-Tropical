# Analyzing the root measurement data for tropical pastures#

#Load required packages#

#read excel root data file#

library(readxl)

Root_data <- read_excel("~/Litriture/Supervisors documents/From Dave/Gatton Dairy Research data/Root data.xlsx")
view(Root_data)

Root_data<- read.csv(file = "/Users/pjayasin/Documents/Litriture/Supervisors documents/From Dave/Gatton Dairy Research data/Root data.csv")
Mean_Root<- read.csv(file = "/Users/pjayasin/Documents/Litriture/Supervisors documents/From Dave/Gatton Dairy Research data/Mean_root_data.csv")

attach(Root_data)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(broom)
library(ggrepel)
library(ggpubr)
library(plotly)
library(fastmap)
update.packages("fastmap")

# Plot root distribution of Rhodes grass #

RR_rootdata<-filter(Root_data, pasture=="Rhodes grass")
RR_rootdata

RR_mean_root<-filter(Mean_Root, Pasture=="RR")
RR_mean_root

fitRR_weight<-RR_rootdata$cumulativeweight
fitRR_soil<-RR_rootdata$Soildepth
fitRR_relative<-RR_rootdata$relativeweight
RRnlsfit<-data.frame(fitRR_soil,fitRR_relative)
RRnlsfit

RRfit_relative<-nls(fitRR_soil~ 1/(1+((fitRR_relative/d50)^R)), start = list(d50=40, R=2))   
RRfit_relative

attach(root_test)
fit_test<-nls(Sdepth~ 1/((1+((Rweight/d50)^R))), start = list(d50=10, R=1))   
fit_test

summary(RRfit)


library(plotly)
plot_root_RR<-ggplot(RR_mean_root, aes(x= RR_mean_root$SoilDepth,y= RR_mean_root$Mean_Rel_Root))+ 
   geom_point(color="navyblue", size=5) + 
  geom_errorbar(aes(ymin=RR_mean_root$Mean_Rel_Root- RR_mean_root$STd_Root, ymax=RR_mean_root$Mean_Rel_Root+ RR_mean_root$STd_Root), width=2,
                position=position_dodge(.9))+ xlab("Soil Depth (cm) ") +ylab("Relative Root Weight")+
  geom_smooth(method = "nls", 
              method.args = list(formula = y ~ 1/(1+((x/d50)^R)),
                                                 start = list(d50=4.421 , R= 1.1766    )),  
              data = RR_mean_root,colour = "#1955B7",
              se = FALSE )+theme(text = element_text(size = 20))  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_rect(colour = "black", fill = NA, size=1), axis.line = element_line(colour = "black"))

plot_root_RR

Flip_RR<-plot_root_RR+coord_flip()


library(investr)

invest(RRrootfit1, y0=c(0,0.5,1))


# Gatton panic

GP_rootdata<-filter(Root_data, pasture=="Gatton panic")
GP_rootdata
GP_mean_root<-filter(Mean_Root,Pasture=="GP")
GP_mean_root

fitGP_weight<-GP_rootdata$cumulativeweight
fitGP_soil<-GP_rootdata$soildepth
fitGP_relative<-GP_rootdata$relativeweight
GPnlsfit<-data.frame(fitGP_weight,fitGP_relative, fitGP_soil)
GPnlsfit



GPfit_relative<-nls(fitGP_relative~ 1/(1+((fitGP_soil/d50)^R)), start = list(d50=10, R=1))   
GPfit_relative


summary(GPfit)



plot_root_GP1<-ggplot(GP_mean_root, aes(x= GP_mean_root$SoilDepth,y= GP_mean_root$Mean_Rel_Root)) + 
  geom_point(color="red", size=5) + 
  geom_errorbar(aes(ymin=GP_mean_root$Mean_Rel_Root- GP_mean_root$STd_Root, ymax=GP_mean_root$Mean_Rel_Root+ GP_mean_root$STd_Root), width=2,
                position=position_dodge(.9))+ xlab("Soil Depth (cm) ") +ylab("Relative Root Weight")+
  geom_smooth(method = "nls", 
              method.args = list(formula = y ~ 1/(1+((x/d50)^R)),
                                 start = list(d50= 11.694 , R= 3.322   )),  
              data = GP_mean_root,colour = "#FC717F",
              se = FALSE )+ theme(text = element_text(size = 20))  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                        panel.background = element_rect(colour = "black", fill = NA, size=1), axis.line = element_line(colour = "black"))
plot_root_GP1
Flip_GP<-plot_root_GP1+coord_flip()

#Brachiaria

BM_rootdata<-filter(Root_data, pasture=="Brachiaria")
BM_rootdata
view(BM_rootdata)
BM_mean_root<-filter(Mean_Root, Pasture=="BM")
BM_mean_root

fitBM_soil<-BM_rootdata$soildepth
fitBM_relative<-BM_rootdata$relativeweight
BMnlsfit<-data.frame(fitBM_soil,fitBM_relative)
BMnlsfit



BMfit_relative<-nls(fitBM_relative~ 1/(1+((fitBM_soil/d50)^R)), start = list(d50=10, R=1))   
BMfit_relative

plot_root_BM1<-ggplot(BM_mean_root, aes(x= BM_mean_root$SoilDepth,y= BM_mean_root$Mean_Rel_Root)) + 
  geom_point(color="Orange", size=5) +geom_errorbar(aes(ymin=BM_mean_root$Mean_Rel_Root- BM_mean_root$STd_Root, ymax=BM_mean_root$Mean_Rel_Root+ BM_mean_root$STd_Root), width=2,
                                                     position=position_dodge(.9))+ xlab("Soil Depth (cm) ") +ylab("Relative Root Weight")+
  geom_smooth(method = "nls", 
              method.args = list(formula = y ~ 1/(1+((x/d50)^R)),
                                 start = list(d50= 10.818 , R= 2.739  )),  
              data = BM_mean_root,colour = "Orange",
              se = FALSE )+theme(text = element_text(size = 20))  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                         panel.background = element_rect(colour = "black", fill = NA, size=1), axis.line = element_line(colour = "black"))

plot_root_BM1
FLip_BM

All_root_plot<-ggarrange(Flip_RR,Flip_GP, FLip_BM, nrow = 1, ncol = 3)
All_root_plot
ggsave("All_root_plot.png",All_root_plot,height = 4, width = 8)
