calculate_impact_EHT=function(results, results_washed,L_vector=3, kappa_vector=2, washedDecay=T, uncertainty=F, vectormodel, myname, coverage_vector, insecticide_id){
  
  npoints=3*365
  
  interventions_vec <- interventions_vectorial_capacity_wrapper(results=results,results_washed20 =results_washed,decay = "weibull",
                                                                        model_p=vectormodel,
                                                                        insecticides=insecticide_id,L = L_vector*365, kappa=kappa_vector,
                                                                        names=myname,  npoints=npoints,washedDecay=washedDecay, uncertainty=uncertainty)
  
  pyrethroid_params=c(0, 0.05, 0)
  names(pyrethroid_params)=c("Deterrency","PrePrandial", "PostPrandial" )
  pyrethroid_interv=get_interv(myproba=pyrethroid_params,L=L_vector*365, kappa=kappa_vector, npoints=3*365,
                               model_p=vectormodel, name="pyrethroid", intervention_type="LLINs", decay="weibull",halflife_insecticide=NULL)
  
  # calculate the impact
  impacts = calculate_combined_impact(intervention1=interventions_vec[[1]],
                                      intervention2=pyrethroid_interv[[1]],
                                      cov_intervention1 = coverage_vector,
                                      cov_intervention2 = coverage_vector,
                                      cov_combination = coverage_vector,
                                      N_vec =10000)

  # 
  return(impacts)
  
}


extract_VCreduction=function(myimpact){
  myeffect=myimpact$interventions_vec
  VCred=data.frame()
  for(i in 1:length(myeffect)){
    this.data=data.frame(insecticide=myeffect[[i]]$description,
                         VCred=myeffect[[i]]$effects$avg_impact,
                         coverages=myeffect[[i]]$coverages)
    VCred=rbind(VCred, this.data)
  }
  return(VCred)
}


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

rhythms_function_old=function(country_hu, country_mosq){
  # extracting all activity data from one specific country from the AnophelesModel database
  select_idx_hu =   activity_patterns$country == country_hu
  select_idx_mosq =   activity_patterns$country == country_mosq
  my_rhythms_hu = activity_patterns[select_idx_hu,] %>%
    filter(sampling %in% c("BED", "IND"))
  my_rhythms_mosqs = activity_patterns[select_idx_mosq,] %>%
    filter(sampling %in% c("HBI", "HBO"))
  
  my_rhythms=rbind(my_rhythms_hu, my_rhythms_mosqs)%>%
    mutate(indicator=ifelse(sampling=="BED", "Humans in bed",
                            ifelse(sampling=="IND", "Humans indoors",
                                   ifelse(sampling=="HBO", "Mosquitoes biting outdoors", "Mosquitoes biting indoors"))))
  
  
  
  labels = substr(unique(my_rhythms$hour), 1, 5)
  my_rhythms$hour = factor(labels, levels = labels )
  my_rhythms_mean=my_rhythms%>%
    mutate(species2=ifelse(species=="Homo sapiens", "Homo sapiens","Anopheles"))%>%
    group_by(species2, hour, sampling)%>%
    summarise(value=mean(value, na.rm = T))%>%
    filter(!is.na(value))
  
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

run_anophelesModelforCascade = function(ent_params, host_params, activity_p , activity_noRhythms,results, results_washed, usage, insecticide_id, attrition, kappa_attrition){
  model_p = AnophelesModel::build_model_obj(vec_p=ent_params, hosts_p= host_params, activity=activity_p, total_pop=2000)
  model_noRhythms = AnophelesModel::build_model_obj(vec_p=ent_params, hosts_p= host_params, activity=activity_noRhythms, total_pop=2000)
  
  # bednet efficacy
  impacts_efficacy=calculate_impact_EHT(results=results, results_washed=results_washed,
                                        washedDecay=F, vectormodel = model_noRhythms, myname="bednet",coverage_vector=c(1), insecticide_id=insecticide_id)
  
  # bednet efficacy + attrition
  impacts_EfficacyAttrition=calculate_impact_EHT(results=results, results_washed=results_washed,
                                                 L_vector=attrition,
                                                 kappa_vector=kappa_attrition,
                                                 washedDecay=F, vectormodel = model_noRhythms, myname="bednet",coverage_vector=c(1), insecticide_id=insecticide_id)
  
  # bednet efficacy + attrition
  impacts_EfficacyAttritionUsage=calculate_impact_EHT(results=results, results_washed=results_washed,
                                                 L_vector=attrition,
                                                 kappa_vector=kappa_attrition,
                                                 washedDecay=F, vectormodel = model_noRhythms, myname="bednet",coverage_vector=c(usage), insecticide_id=insecticide_id)
  
  # bednet efficacy + attrition + insecticide decay
  impacts_EfficacyAttritionDecay=calculate_impact_EHT(results=results, results_washed=results_washed,
                                                      L_vector=attrition,
                                                      kappa_vector=kappa_attrition,
                                                      washedDecay=T, vectormodel = model_noRhythms, myname="bednet",coverage_vector=c(usage), insecticide_id=insecticide_id)
  
  # bednet efficacy + attrition + insecticide decay + rhythms
  impacts_EfficacyAttritionDecayRhythms=calculate_impact_EHT(results=results, results_washed=results_washed,
                                                             L_vector=attrition,
                                                             kappa_vector=kappa_attrition,
                                                             washedDecay=T, vectormodel = model_p, myname="bednet",coverage_vector=c(usage), insecticide_id=insecticide_id)
  
  
  VCreduc_all=rbind(data.frame(value=impacts_efficacy ,variable="Efficacy"),
                    data.frame(value=impacts_EfficacyAttrition, variable="EfficacyAttrition"),
                    data.frame(value=impacts_EfficacyAttritionUsage, variable="EfficacyAttritionUsage"),
                    data.frame(value=impacts_EfficacyAttritionDecay ,variable="EfficacyAttritionUsageDecay"),
                    data.frame(value=impacts_EfficacyAttritionDecayRhythms ,variable="EfficacyAttritionUsageDecayRhythms")
  )%>% mutate(value=100*value)
  
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
          strip.text.x = element_text(size = 16), axis.title=element_text(size=16), 
          axis.text=element_text(size=10),axis.text.x=element_text(angle=45), strip.background =element_rect(fill="white"), legend.text =element_text(size=16))
  
  
  return(plt)
}
