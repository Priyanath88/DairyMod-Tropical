##attach the required packages//

library(nasapower) // reading climate data from the NASAPower data
library(apsimx)  //Reading data from the APSIM database
library(soilDB) // reading soil data from U.S. Department of Agriculture Natural Resources Conservation Service (USDA-NRCS) and National Cooperative Soil Survey (NCSS) databases//
library(ggplot2)
library(tidyverse)
library(chirps) // reading soil data from CHIRPS database

CLimate_Australia<- get_power(
  community = "ag",
  lonlat = c(151.81, -27.48),
  pars = c("RH2M", "T2M_MAX", "T2M_MIN", "PRECTOTCORR"),
  dates = c("1990-01-01","2020-12-31"),
  temporal_api = "daily"
  )
CLimate_Australia
  
APSIM_CLimate_Australia<- get_power_apsim_met(
    lonlat = c(152.33, -27.54),
    dates = c("1990-01-01","2020-12-31"),
  )
APSIM_CLimate_Australia
plot(APSIM_CLimate_Australia, met.var = "rain", years = 1990: 2020, cumulative= TRUE, climatology = TRUE )
frost_days<- plot(APSIM_CLimate_Australia, summary= TRUE, compute.frost = TRUE, met.var= "frost_days")

US_soil_data<- get_ssurgo_soil_profile(lonlat = c(-82.35,27.26 ))
US_soil_data
plot(US_soil_data, property= "water")


Duri_metdata<-

  params = list(
    lat="-31.16",
    lon="150.52",
    start="20000101", 
    finish="20160131",
    format="csv",
    comment="RXN",
    username="john.doe@xyz.com.au",
    password="silo"
  )
res <- httr::GET("https://www.longpaddock.qld.gov.au/cgi-bin/silo/DataDrillDataset.php", query=params)
silodata <- read_csv(httr::content(res, as="text")) 
head(silodata)

write.csv(silodata, "DuriMet.csv")


soil_Duri_Aus<- get_isric_soil_profile(lonlat = c(150.8666, -31.2666))
soil_Duri_Aus
view(soil_Duri_Aus)
plot(soil_Duri_Aus, property = "water")
plot(soil_Aus, prperty= "soilorganicmatter")
plot(soil_Aus, prperty= "soil.KS")


library(tidyverse)
library(ggplot2)
library(sf)
library(ggpmisc)
library(maps)
require(viridis)

Validation_map<-read.csv(file = "/Users/pjayasin/OneDrive - Massey University/Thesis chapeters/DairyMod modelling chapter/validation_data_map.csv")

world <- map_data("world")
world
ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region)
  )  

validation <- ggplot() +
  geom_map(data = world, map = world,
    aes(long, lat, map_id = region),
    color = "black", fill = "white", size = 0.1
  )+coord_cartesian(ylim = c(-60, +90))+ geom_point(
    data = Validation_map,
    aes(Long, Lat, color = Pasture, shape = Pasture),
    alpha = 0.5, size = 1.5
  )+ylab("Latitude")+ xlab("Longitude")+ scale_shape_manual(values = c(16, 17, 15))+ scale_colour_manual(values = c("Orange", "Red", "Blue"))+
  theme(legend.position = "bottom")+ theme(panel.background = element_blank())+
theme(panel.border = element_rect(color = "black",
                                  fill = NA, size = 1))+
  theme(legend.position = c(0.15, 0.25))+ theme(legend.background=element_blank())+                                             # Apply guides function
  guides(color = guide_legend(override.aes = list(size = 4)))+ theme(legend.key=element_rect(fill="white"))
validation

ggsave(filename= "validation_data_map.png",width = 7, height = 4,  validation)


##Climate data retrieve for Atenas, Costa Rica##

CLimate_Atenas_Costa_Rica<- get_power(
  community = "ag",
  lonlat = c(-84.24, 9.57),
  pars = c("ALLSKY_SFC_SW_DWN", "T2M_MAX", "T2M_MIN", "PRECTOTCORR"),
  dates = c("1990-01-01","2020-12-31"),
  temporal_api = "daily"
)
CLimate_Atenas_Costa_Rica
write.csv(CLimate_Atenas_Costa_Rica, "Clamate_Costa Rica.csv")


soil_Atenas_Costa_Rica<- get_isric_soil_profile(lonlat = c(-84.24, 9.57))
soil_Atenas_Costa_Rica
view(soil_Atenas_Costa_Rica)


##SOil profile for Moora, Australia
soil_Moora_AUS<- get_isric_soil_profile(lonlat = c(115.55,-30.45))
soil_Moora_AUS
view(soil_Moora_AUS)
plot(soil_Moora_AUS,property= "water")


#soil for Vietnam#

soil_Vietnam<- get_isric_soil_profile(lonlat = c(108.29,14.10))
soil_Vietnam
view(soil_Vietnam)
plot(soil_Vietnam, property= "water")


#Soil Ecuador ##

Soil_ECU<- get_isric_soil_profile(lonlat = c(-79.2608,-0.2222))
Soil_ECU
view(Soil_ECU)
write.csv(Soil_ECU, "soil_ECU.csv")
plot(Soil_ECU,property= "water")


lonlat_EC <- data.frame(lon = 34.93, lat = -9.43)
dates_EC <- c("2010-01-01", "2012-12-31")
data_EC <- get_chirps(lonlat_EC, dates_EC, server = "ClimateSERV")

data_EC
write.csv(data_EC, file = "ecuador2.csv")


#soil Tanzania#

SOil_TANZ<- get_isric_soil_profile(lonlat = c(34.93, -9.43))
SOil_TANZ
view(SOil_TANZ)

#Tanzania rainfall data#

lonlat <- data.frame(lon = 34.93, lat = -9.43)
dates <- c("2015-01-01", "2021-12-31")
data <- get_chirps(lonlat, dates, server = "ClimateSERV")

data
write.csv(data, file = "TRANZ data.csv")

#Mtwango rainfall#

lonlat_MT <- data.frame(lon = 35.58, lat = -8.58)
dates_MT <- c("2015-01-01", "2021-12-31")
data_MT <- get_chirps(lonlat_MT, dates_MT, server = "ClimateSERV")

data_MT
write.csv(data_MT, file = "TRANZ_MT data.csv")


SOil_Mtwango<- get_isric_soil_profile(lonlat = c(35.58, -8.58))
SOil_Mtwango
view(SOil_Mtwango)


# Igarape_Brazil rainfall data#

lonlat_IG <- data.frame(lon = -47.36, lat = -1.01)
dates_IG <- c("2010-01-01", "2021-03-31")
data_IG <- get_chirps(lonlat, dates, server = "ClimateSERV")

data_IG
write.csv(data_IG, file = "Igarape_data.csv")


#Weerawila_LKA

data_WR<- data.frame(lon= 81.23, lat= 6.28)
dates_wee_RF<- c("1982-01-01", "1986-12-31")
data_wee_RF<- get_chirps(data_WR, dates_wee_RF, server = "ClimateSERV")
data_wee_RF
write.csv(data_wee_RF, file = "werawilaweatherdata.csv")
