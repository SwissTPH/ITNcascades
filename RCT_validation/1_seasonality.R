#######################
# Generate seasonality profiles for all countries at admin 1 level
#
# monica.golumbeanu@unibas.ch
# 21.10.2021
# modified by Clara Champagne in december 2022
######################
rm(list=ls())
library(chirps)
library(sf)
library(raster)
library(dplyr)
library(reshape2)
library(ggplot2)

output_dir=""

#####################
# Extract CHIRPS data for Mosha trial (Misungwi)
#####################

shp_TZA=getData('GADM', country='TZA', level=2)
shp_Misungwi=subset(shp_TZA, NAME_2=="Misungwi")

# sample three points within the Misungwi council
set.seed(1234)
tp_point <- st_sample(st_as_sf(shp_Misungwi), 10)
tp_point <- st_as_sf(tp_point)
dat <- get_chirps(tp_point,
                  dates = c("2018-01-01","2022-12-31"), 
                  server = "ClimateSERV")

# Format CHIRPS data
get_seasonality_chirps=function(mydata){
  mydata_month=mydata %>% 
    mutate(month=format(date, "%y-%m"))%>%
    group_by(month) %>%
    summarise(rainfall=sum(chirps))%>%
    mutate(date=as.Date(paste0("20", month, "-01")))
  
  # Aggregate CHIRP per month over all considered years
  mydata_year=mydata_month %>% 
    mutate(month=as.numeric(as.character(format(date, "%m"))))%>%
    group_by(month) %>%
    summarise(rainfall=mean(rainfall))
  
  # change to same format as Worldclim
  seasonality = mydata_year  %>% 
    mutate(month=ifelse(month==12, 1, month+1) )%>%   # shift by 1 month for OM
    pivot_wider(names_from = month, values_from = rainfall, names_prefix = "m")
  
  # Normalize
  seasonality[, paste0("m", 1:12)] = seasonality[,  paste0("m", 1:12)]/rowSums(seasonality[,  paste0("m", 1:12)])
  
  return(seasonality)
}
# normalized anomaly: (value-mean)/sd

# CHIRPS 2018-2021 (for Mosha trial)
seasonality_chirps=get_seasonality_chirps(dat)

write.csv(seasonality_chirps, file.path(output_dir, "seasonality_Misungwi.csv"), row.names = FALSE)



##########################################
# format all datasets similarly for plotting

seasons_chirps=seasonality_chirps %>% pivot_longer(cols=starts_with("m"),names_to = "month",values_to="intensity",names_prefix = "m") %>%
  mutate(month=as.numeric(month), dataset="CHIRPS (2018-2022)") 


rbind(seasons_chirps)%>%
  ggplot()+
  geom_line(aes(x=month,y=intensity, color=dataset))+
  theme_minimal()+
  scale_color_manual(values=c("black","red", "orange", "gold", "darkred", "blue"), name="")
ggsave(file.path(output_dir, "misungwi_chirps_average1822.png"))


###############################################################################################################
###############################################################################################################
# MULEBA district for, Protopopoff 2018 trial
shp_Muleba=subset(shp_TZA, NAME_2=="Muleba")
set.seed(1234)
tp_point_muleba <- st_sample(st_as_sf(shp_Muleba), 10)
tp_point_muleba <- st_as_sf(tp_point_muleba)
dat_muleba <- get_chirps(tp_point_muleba,
                         dates = c("2014-01-01","2017-12-31"), 
                         server = "ClimateSERV")


# CHIRPS 2014-2016 (for Protopopoff trial)
seasonality_muleba=get_seasonality_chirps(dat_muleba)
write.csv(seasonality_muleba, file.path(output_dir, "seasonality_Muleba.csv"), row.names = FALSE)


##########################################
# format all datasets similarly for plotting

seasons_chirps_muleba=seasonality_muleba %>% pivot_longer(cols=starts_with("m"),names_to = "month",values_to="intensity",names_prefix = "m") %>%
  mutate(month=as.numeric(month), dataset="CHIRPS (2014-2017)") 

rbind( seasons_chirps_muleba)%>%
  ggplot()+
  geom_line(aes(x=month,y=intensity, color=dataset))+
  theme_minimal()+
  scale_color_manual(values=c( "black"), name="")
ggsave(file.path(output_dir, "muleba_chirps_average1417.png"))

