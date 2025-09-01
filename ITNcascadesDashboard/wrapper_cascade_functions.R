
rhythms_function=function(country_hu, country_mosq,my_activity_patterns){
  # extracting all activity data from one specific country from the AnophelesModel database
  select_idx_hu =   my_activity_patterns$country == country_hu
  select_idx_mosq =   my_activity_patterns$country == country_mosq
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
    group_by(species2, hour,indicator, sampling)%>%
    summarise(value=mean(value, na.rm = T))%>%
    filter(!is.na(value))
  
  return(list("my_rhythms"=my_rhythms, "my_rhythms_mean"= my_rhythms_mean))
}

activity_function=function(my_rhythms_mean, species){
  humans_in_bed =my_rhythms_mean$value[my_rhythms_mean$sampling=="BED"]
  if(any(my_rhythms_mean$sampling=="IND")==FALSE){
    humans_indoors = humans_in_bed
  } else {
    humans_indoors = my_rhythms_mean$value[my_rhythms_mean$sampling=="IND"]
  }
  HBI =my_rhythms_mean$value[my_rhythms_mean$sampling=="HBI"]
  HBO = my_rhythms_mean$value[my_rhythms_mean$sampling=="HBO"]
  activity_p = def_activity_patterns(as.data.frame(cbind(HBI, HBO, humans_indoors,humans_in_bed )))
  inout=get_in_out_exp(activity_cycles = activity_p, vec_p = def_vector_params(species))
  return(inout$Exposure_Indoor_whileinbed)
}

run_anophelesModelforCascade = function(ent_params, host_params, activity_noRhythms,outputs_posteriormax, usage,exposure, insecticide_id, attrition, kappa_attrition){
  
  model_noRhythms = AnophelesModel::build_model_obj(vec_p=ent_params, hosts_p= host_params, activity=activity_noRhythms, total_pop=2000)
  
  # bednet efficacy
  impacts_efficacy=interventions_vectorial_capacity(outputs_posteriormax=outputs_posteriormax,
                                             insecticide_decay=F,
                                             model_p = model_noRhythms, cov=c(1), L=3*365, kappa=2)
  
  # bednet efficacy + attrition
  impacts_EfficacyUsage=interventions_vectorial_capacity(outputs_posteriormax=outputs_posteriormax,
                                                         cov=c(usage), L=3*365, kappa=2,
                                                 insecticide_decay=F, model_p = model_noRhythms)
  
  
  # bednet efficacy + attrition
  impacts_EfficacyAttritionUsage=interventions_vectorial_capacity(outputs_posteriormax=outputs_posteriormax,
                                                 L=attrition*365,
                                                 kappa=kappa_attrition,
                                                 insecticide_decay=F, model_p = model_noRhythms, cov=c(usage))
  
  # bednet efficacy + attrition + insecticide decay
  impacts_EfficacyAttritionDecay=interventions_vectorial_capacity(outputs_posteriormax=outputs_posteriormax,
                                                      L=attrition*365,
                                                      kappa=kappa_attrition,
                                                      insecticide_decay=T, model_p = model_noRhythms, cov=c(usage))
  
  
  # bednet efficacy + attrition + insecticide decay + rhythms
  impacts_EfficacyAttritionDecayRhythms=interventions_vectorial_capacity(outputs_posteriormax=outputs_posteriormax,
                                                             L=attrition*365,
                                                             kappa=kappa_attrition,
                                                             insecticide_decay=T, model_p = model_noRhythms, cov=c(usage*exposure))
  
  
  print(impacts_EfficacyAttritionDecayRhythms)
  
  VCreduc_all=rbind(data.frame(value=impacts_efficacy ,variable="Efficacy"),
                    data.frame(value=impacts_EfficacyUsage, variable="EfficacyUsage"),
                    data.frame(value=impacts_EfficacyAttritionUsage, variable="EfficacyAttritionUsage"),
                    data.frame(value=impacts_EfficacyAttritionDecay ,variable="EfficacyAttritionUsageDecay"),
                    data.frame(value=impacts_EfficacyAttritionDecayRhythms ,variable="EfficacyAttritionUsageDecayRhythms")
  )%>% mutate(value=100*value)
  
  print(VCreduc_all)
  return(VCreduc_all)
}


rhythms_function_plot=function(my_rhythms, my_rhythms_mean){
  
  plt=ggplot(my_rhythms, aes(x=hour, y=value)) +
    geom_line(size=1, aes( group=id, color=citation, group=sampling), alpha=0.7) + 
    geom_line(data=my_rhythms_mean,aes(x=hour, y=value,group=sampling, color="Average"), size=2, linetype="dashed") + 
    labs(color="",
         x="Time of the day (hh.mm)", y="Activity")+
    facet_wrap(.~indicator, scales = "free", ncol=2) +
    scale_color_manual(values=c("black", "darkred", "gold", "cyan4", "lightblue", "darkorange", "dodgerblue", "grey", "violet", "forestgreen"))+ 
    theme_bw()+
    theme(legend.position = "bottom", 
          strip.text.x = element_text(size = 12), axis.title=element_text(size=16), 
          axis.text=element_text(size=10),axis.text.x=element_text(angle=45), strip.background =element_rect(fill="white"), legend.text =element_text(size=12))
  
  
  return(plt)
}
