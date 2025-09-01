library(AnophelesModel)
library(dplyr)

list_country_hu=c("Cote D'Ivoire", "Haiti", "Kenya",
                  "Mozambique",  "Tanzania","Zambia")
list_country_mosq=c("Benin", "Burkina Faso","Cameroon", "Democratic Republic of Congo",  
                    "Erithrea", "Gabon", "Ghana", "Haiti", "India", "Kenya",
                    "Mozambique", "Nigeria", "PNG", "Sao Tome", "Senegal",
                    "Tanzania", "Thailand", "The Gambia", "Uganda","Zambia", "Zimbabwe")

list_country_africa=c("Benin", "Burkina Faso","Cameroon", "Democratic Republic of Congo", "Equatorial Guinea", 
                      "Erithrea", "Gabon", "Ghana","Kenya", "Ethiopia",
                      "Mozambique","Cote D'Ivoire", "Nigeria",  "Sao Tome", "Senegal",
                      "Tanzania","The Gambia", "Uganda","Zambia", "Zimbabwe")


rhythms_function=function(country_hu, country_mosq,my_activity_patterns=activity_patterns){
  # extracting all activity data from one specific country from the AnophelesModel database
  if(! country_hu %in% list_country_hu){
    select_idx_hu =   my_activity_patterns$country %in% list_country_africa
  } else {
    select_idx_hu =   my_activity_patterns$country == country_hu
  }
  
  if(! country_mosq %in% list_country_mosq){
    select_idx_mosq =   my_activity_patterns$country %in% list_country_africa
  } else {
    select_idx_mosq =   my_activity_patterns$country == country_mosq
  }
  
  my_rhythms_hu = my_activity_patterns[select_idx_hu,] %>%
    filter(sampling %in% c("BED", "IND"))
  my_rhythms_mosqs = my_activity_patterns[select_idx_mosq,] %>%
    filter(sampling %in% c("HBI", "HBO"))
  
  my_rhythms=rbind(my_rhythms_hu, my_rhythms_mosqs)%>%
    mutate(indicator=ifelse(sampling=="BED", "Humans in bed",
                            ifelse(sampling=="IND", "Humans indoors",
                                   ifelse(sampling=="HBO", "Mosquitoes biting outdoors", "Mosquitoes biting indoors"))))
  
  
  
  labels = substr(unique(my_rhythms$hour), 1, 5)
  my_rhythms$hour = factor(labels, levels = labels )
  my_rhythms_mean=my_rhythms%>%
    #mutate(species2=ifelse(species=="Homo sapiens", "Homo sapiens","Anopheles"))%>%
    group_by(species, hour,indicator, sampling)%>%
    summarise(value=mean(value, na.rm = T))%>%
    filter(!is.na(value))
  
  return(list("my_rhythms"=my_rhythms, "my_rhythms_mean"= my_rhythms_mean))
}


activity_function=function(my_rhythms_mean){
  humans_in_bed =my_rhythms_mean$value[my_rhythms_mean$sampling=="BED"]
  if(any(my_rhythms_mean$sampling=="IND")==FALSE){
    humans_indoors = humans_in_bed
  } else {
    humans_indoors = my_rhythms_mean$value[my_rhythms_mean$sampling=="IND"]
  }
  HBI =my_rhythms_mean$value[my_rhythms_mean$sampling=="HBI"]
  HBO = my_rhythms_mean$value[my_rhythms_mean$sampling=="HBO"]
  activity_p = def_activity_patterns(as.data.frame(cbind(HBI, HBO, humans_indoors,humans_in_bed )))
  #get_in_out_exp(activity_cycles = activity_p, vec_p = def_vector_params("Anopheles gambiae"))
  return(activity_p)
}


ent_params = def_vector_params(mosquito_species = ("Anopheles gambiae"))

## Mosha and Protopopoff
rhythms=rhythms_function(country_hu ="Tanzania",
                         country_mosq ="Tanzania")

exposures=get_in_out_exp(activity_cycles = activity_function(rhythms$my_rhythms_mean), vec_p = ent_params)



## Accrombessi
rhythms=rhythms_function(country_hu ="Benin",
                         country_mosq ="Benin")

exposures=get_in_out_exp(activity_cycles = activity_function(rhythms$my_rhythms_mean), vec_p = ent_params)



## Staedke
rhythms=rhythms_function(country_hu ="Uganda",
                         country_mosq ="Uganda")

exposures=get_in_out_exp(activity_cycles = activity_function(rhythms$my_rhythms_mean), vec_p = ent_params)
