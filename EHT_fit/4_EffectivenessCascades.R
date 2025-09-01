rm(list=ls())
library(dplyr)
library(ggplot2)
library(readxl)
library(lubridate)
library(stringr)
library(cowplot)
library(tidyr)
library(AnophelesModel)
library(gt)

mainDir="."
scriptDir=file.path(".","EHT_fit")

stanDir=file.path(mainDir, "results/stan_outputs")
plotDir=file.path(mainDir, "plots")
dataDir=file.path(scriptDir, "processed_data")
outputsDir=file.path(mainDir, "results/csv_outputs")
source(file.path(scriptDir, "functions_compute_vectorialCapacity.R"))
source(file.path(scriptDir, "functions_cascade_LLIN.R"))
source(file.path(scriptDir, "functions_fit_stan.R"))

rerun=FALSE

#####################################################################
# COMPUTE VECTORIAL CAPACITY REDUCTION
####################################################################


# define entomological parameters and activity rhythms
ent_params = def_vector_params(mosquito_species = "Anopheles gambiae")

# asuming everyone indoors and in bed all the time for cascade calculations
activity_noRhythms =def_activity_patterns()
activity_noRhythms$humans_in_bed=rep(1, length(activity_noRhythms$HBI))
activity_noRhythms$humans_indoors=rep(1, length(activity_noRhythms$HBI))
get_in_out_exp(activity_cycles = activity_noRhythms, vec_p = ent_params)

#Compile the model
model_noRhythms = build_model_obj(vec_p=ent_params, hosts_p= def_host_params(), activity=activity_noRhythms, total_pop=2000)


results_kibondo_w_72<- readRDS(file.path(stanDir,"/stan_kibondo_washed_72.rds"))
results_kibondo_unw_72<- readRDS(file.path(stanDir,"/stan_kibondo_unwashed_72.rds"))
results_kibondo_w_24<- readRDS(file.path(stanDir,"/stan_kibondo_washed_24.rds"))
results_kibondo_unw_24<- readRDS(file.path(stanDir,"/stan_kibondo_unwashed_24.rds"))


results_bit055_w<- readRDS(file.path(stanDir,"/stan_BIT055_washed_24.rds"))
results_bit055_unw<- readRDS(file.path(stanDir,"/stan_BIT055_unwashed_24.rds"))


results_bit103_unw_72<- readRDS(file.path(stanDir,"/stan_assenga_unwashed_72_all.rds"))
results_bit103_w_72<- readRDS(file.path(stanDir,"/stan_assenga_washed_72_all.rds"))


results_bit059_w<- readRDS(file.path(stanDir,"/stan_odufuwa_washed_24.rds"))
results_bit059_unw<- readRDS(file.path(stanDir,"/stan_odufuwa_unwashed_24.rds"))

results_bit080_w_72<- readRDS(file.path(stanDir,"/stan_BIT080_washed_72.rds"))
results_bit080_unw_72<- readRDS(file.path(stanDir,"/stan_BIT080_unwashed_72.rds"))
results_bit080_w_24<- readRDS(file.path(stanDir,"/stan_BIT080_washed_24.rds"))
results_bit080_unw_24<- readRDS(file.path(stanDir,"/stan_BIT080_unwashed_24.rds"))


results_martin_w72<- readRDS(file.path(stanDir,"/stan_martin_36m_72_gambiae.rds"))
results_martin_unw72<- readRDS(file.path(stanDir,"/stan_martin_unwashed_72_gambiae.rds"))
results_martin_w72_funestus<- readRDS(file.path(stanDir,"/stan_martin_36m_72_funestus.rds"))
results_martin_unw72_funestus<- readRDS(file.path(stanDir,"/stan_martin_unwashed_72_funestus.rds"))



# from fitting survey data on LLIN usage from Mosha et al.
L_vector_attrition=c("ig2"=2.4, "pbo"=1.7)
kappa_vector_attrition=c("ig2"=2.4,"pbo"=1.9)
coverage_vector_usage = c("ig2"=0.69,  "pbo"=0.75)

#####################################################################
# SETUP ANOPHELES MODEL
####################################################################

# define entomological parameters and activity rhythms
ent_params = def_vector_params(mosquito_species = "Anopheles gambiae")

# extracting all activity data from Tanzania from the AnophelesModel database
select_idx =   activity_patterns$country == "Tanzania"
TZ_rhythms = activity_patterns[select_idx,] %>%
  mutate(indicator=ifelse(sampling=="BED", "Humans in bed",
                          ifelse(sampling=="IND", "Humans indoors",
                                 ifelse(sampling=="HBO", "Mosquitoes biting outdoors", "Mosquitoes biting indoors"))))

TZ_rhythms$ref=ifelse(TZ_rhythms$id %in% c(3, 25,59, 60,61,190, 191, 192), "Geissbuhler et al. 2007",
                      ifelse(TZ_rhythms$id %in% c(8,9,69,70, 200, 201 ),  "Huho et al. 2013",
                             ifelse(TZ_rhythms$id %in% c(21,122, 123, 124, 125, 126, 127, 253, 254, 255, 256, 257, 258 ), "Russell et al. 2011",
                                    ifelse(TZ_rhythms$id %in% c(26,27,81, 82, 212, 213 ), "Killeen et al. 2006",
                                           ifelse(TZ_rhythms$id %in% c(52,53, 183, 184 ), "Maia et al. 2016",NA)))))

Vlabels = substr(unique(TZ_rhythms$hour), 1, 5)
TZ_rhythms$hour = factor(substr(TZ_rhythms$hour, 1, 5), levels =Vlabels )


TZ_rhythms_mean=TZ_rhythms%>%
  mutate(species2=ifelse(species=="Homo sapiens", "Homo sapiens","Anopheles"))%>%
  group_by(species2, hour, indicator, sampling)%>%
  summarise(value=mean(value, na.rm = T))%>% ungroup()

ggplot(TZ_rhythms, aes(x=hour, y=value)) +
  geom_line(size=1, aes( group=id, color=ref), alpha=0.7) +
  geom_line(data=TZ_rhythms_mean,aes(x=hour, y=value,group=sampling, color="Average"), size=2) +
  scale_color_manual(values=c("black", "darkred", "gold", "cyan4", "lightblue", "darkorange"))+
  labs(color="",
       x="Time of the day (hh.mm)", y="Activity")+
  facet_wrap(.~indicator, scales = "free", ncol=1) +
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 16), axis.title=element_text(size=16),
        axis.text=element_text(size=14), strip.background =element_rect(fill="white"), legend.text =element_text(size=16))
ggsave(file.path(plotDir, "activity_patterns_TZ.png"), width=10, height=15)

humans_in_bed =TZ_rhythms_mean$value[TZ_rhythms_mean$sampling=="BED"]
humans_indoors = TZ_rhythms_mean$value[TZ_rhythms_mean$sampling=="IND"]
HBI =TZ_rhythms_mean$value[TZ_rhythms_mean$sampling=="HBI"]
HBO = TZ_rhythms_mean$value[TZ_rhythms_mean$sampling=="HBO"]
custom_params = as.data.frame(cbind(HBI, HBO, humans_indoors,humans_in_bed ))
activity_p = def_activity_patterns(custom_params)
get_in_out_exp(activity_cycles = activity_p, vec_p = ent_params)

# asuming everyone indoors and in bed all the time for cascade calculations
activity_noRhythms =activity_p
activity_noRhythms$humans_in_bed=rep(1, length(activity_noRhythms$HBI))
activity_noRhythms$humans_indoors=rep(1, length(activity_noRhythms$HBI))
get_in_out_exp(activity_cycles = activity_noRhythms, vec_p = ent_params)

#Compile the mode;
model_noRhythms = build_model_obj(vec_p=ent_params, hosts_p= def_host_params(), activity=activity_noRhythms, total_pop=2000)

#####################################################################
# COMPUTE VECTORIAL CAPACITY REDUCTION
####################################################################
calculate_impact_4EHT=function(L_vector=c("ig2"=3, "pbo"=3),
                               kappa_vector=c("ig2"=2,  "pbo"=2), washedDecay=T, uncertainty=F, vectormodel, nsamples=100,
                               coverage_vector=c("ig2"=1,  "pbo"=1), inbed_exposure=1){

  npoints=3*365

  intervention_list_kibondo <- interventions_vectorial_capacity_wrapper(results=results_kibondo_unw,results_washed20 =results_kibondo_w,decay = "weibull",
                                                                        model_p=vectormodel,
                                                                        insecticides=c(2),L = L_vector[["ig2"]]*365, kappa=kappa_vector[["ig2"]],inbed_exposure=inbed_exposure,
                                                                        names=c("IG2 (Kibondo)"),  npoints=npoints,washedDecay=washedDecay, uncertainty=uncertainty,
                                                                        cov=coverage_vector[["ig2"]], nsamples=nsamples)


  intervention_list_bit055 <- interventions_vectorial_capacity_wrapper(results=results_bit055_unw,results_washed20=results_bit055_w,decay = "weibull",
                                                                       model_p=vectormodel,
                                                                       insecticides=c(2),L =  L_vector[["pbo"]]*365 , kappa=kappa_vector[["pbo"]],inbed_exposure=inbed_exposure,
                                                                       names=c("OlysetPlus (BIT055)"),  npoints=npoints, washedDecay=washedDecay, uncertainty=uncertainty,
                                                                       cov=coverage_vector[["pbo"]], nsamples=nsamples)

  intervention_list_bit103_PBO <- interventions_vectorial_capacity_wrapper(results=results_bit103_unw,results_washed20=results_bit103_w,decay = "weibull",
                                                                           model_p=vectormodel,
                                                                           insecticides=c(2),L =  L_vector[["pbo"]]*365 , kappa=kappa_vector[["pbo"]],inbed_exposure=inbed_exposure,
                                                                           names=c( "OlysetPlus (BIT103)"),  npoints=npoints,washedDecay=washedDecay, uncertainty=uncertainty,
                                                                           cov=coverage_vector[["pbo"]], nsamples=nsamples)

  intervention_list_bit103_IG2 <- interventions_vectorial_capacity_wrapper(results=results_bit103_unw,results_washed20=results_bit103_w,decay = "weibull",
                                                                           model_p=vectormodel,
                                                                           insecticides=c(3),L =  L_vector[["ig2"]]*365 , kappa=kappa_vector[["ig2"]],inbed_exposure=inbed_exposure,
                                                                           names=c("IG2 (BIT103)"),  npoints=npoints,washedDecay=washedDecay, uncertainty=uncertainty,
                                                                           cov=coverage_vector[["ig2"]], nsamples=nsamples)

  intervention_list_bit059 <- interventions_vectorial_capacity_wrapper(results=results_bit059_unw,results_washed20=results_bit059_w,decay = "weibull",
                                                                       model_p=vectormodel,
                                                                       insecticides=c(2),L =  L_vector[["pbo"]]*365 , kappa=kappa_vector[["pbo"]],inbed_exposure=inbed_exposure,
                                                                       names=c("OlysetPlus (BIT059)"),  npoints=npoints, washedDecay=washedDecay, uncertainty=uncertainty,
                                                                       cov=coverage_vector[["pbo"]], nsamples=nsamples)

  intervention_list_bit080 <- interventions_vectorial_capacity_wrapper(results=results_bit080_unw,results_washed20=results_bit080_w,decay = "weibull",
                                                                       model_p=vectormodel,
                                                                       insecticides=c(2),L =  L_vector[["ig2"]]*365 , kappa=kappa_vector[["ig2"]],inbed_exposure=inbed_exposure,
                                                                       names=c("IG2 (BIT080)"),  npoints=npoints, washedDecay=washedDecay, uncertainty=uncertainty,
                                                                       cov=coverage_vector[["ig2"]], nsamples=nsamples)
  
  intervention_list_martin_PBO <- interventions_vectorial_capacity_wrapper(results=results_martin_unw,results_washed20=results_martin_w,decay = "weibull",
                                                                           model_p=vectormodel,
                                                                           insecticides=c(2),L =  L_vector[["pbo"]]*365 , kappa=kappa_vector[["pbo"]],inbed_exposure=inbed_exposure,
                                                                           names=c( "OlysetPlus (Martin)"),  npoints=npoints,washedDecay=washedDecay, uncertainty=uncertainty,
                                                                           cov=coverage_vector[["pbo"]], nsamples=nsamples)
  
  intervention_list_martin_IG2 <- interventions_vectorial_capacity_wrapper(results=results_martin_unw,results_washed20=results_martin_w,decay = "weibull",
                                                                           model_p=vectormodel,
                                                                           insecticides=c(3),L =  L_vector[["ig2"]]*365 , kappa=kappa_vector[["ig2"]],inbed_exposure=inbed_exposure,
                                                                           names=c("IG2 (Martin)"),  npoints=npoints,washedDecay=washedDecay, uncertainty=uncertainty,
                                                                           cov=coverage_vector[["ig2"]], nsamples=nsamples)
  
  intervention_list_martinf_PBO <- interventions_vectorial_capacity_wrapper(results=results_martinf_unw,results_washed20=results_martinf_w,decay = "weibull",
                                                                           model_p=vectormodel,
                                                                           insecticides=c(2),L =  L_vector[["pbo"]]*365 , kappa=kappa_vector[["pbo"]],inbed_exposure=inbed_exposure,
                                                                           names=c( "OlysetPlus (Martin, f)"),  npoints=npoints,washedDecay=washedDecay, uncertainty=uncertainty,
                                                                           cov=coverage_vector[["pbo"]], nsamples=nsamples)
  
  intervention_list_martinf_IG2 <- interventions_vectorial_capacity_wrapper(results=results_martinf_unw,results_washed20=results_martinf_w,decay = "weibull",
                                                                           model_p=vectormodel,
                                                                           insecticides=c(3),L =  L_vector[["ig2"]]*365 , kappa=kappa_vector[["ig2"]],inbed_exposure=inbed_exposure,
                                                                           names=c("IG2 (Martin, f)"),  npoints=npoints,washedDecay=washedDecay, uncertainty=uncertainty,
                                                                           cov=coverage_vector[["ig2"]], nsamples=nsamples)
  
  impacts=rbind(data.frame(VCred=intervention_list_kibondo) %>% mutate(EHT="Kibondo", insecticide="IG2"),
                  data.frame(VCred=intervention_list_bit055) %>% mutate(EHT="BIT055", insecticide="OlysetPlus"),
                  data.frame(VCred=intervention_list_bit103_IG2) %>% mutate(EHT="BIT103", insecticide="IG2"),
                  data.frame(VCred=intervention_list_bit103_PBO) %>% mutate(EHT="BIT103", insecticide="OlysetPlus"),
                  data.frame(VCred=intervention_list_bit059) %>% mutate(EHT="BIT059", insecticide="OlysetPlus"),
                  data.frame(VCred=intervention_list_bit080) %>% mutate(EHT="BIT080", insecticide="IG2"),
                data.frame(VCred=intervention_list_martin_IG2) %>% mutate(EHT="Martin", insecticide="IG2"),
                data.frame(VCred=intervention_list_martin_PBO) %>% mutate(EHT="Martin", insecticide="OlysetPlus"),
                data.frame(VCred=intervention_list_martinf_IG2) %>% mutate(EHT="Martin, f", insecticide="IG2"),
                data.frame(VCred=intervention_list_martinf_PBO) %>% mutate(EHT="Martin, f", insecticide="OlysetPlus")
    )

  return(impacts)

}

#####################################################################
# CASCADE WITHOUT UNCERTAINTY
####################################################################

results_kibondo_w=results_kibondo_w_72
results_kibondo_unw=results_kibondo_unw_72
results_bit103_unw=results_bit103_unw_72
results_bit103_w=results_bit103_w_72
results_bit080_unw=results_bit080_unw_72
results_bit080_w=results_bit080_w_72
results_martin_unw=results_martin_unw72
results_martin_w=results_martin_w72
results_martinf_unw=results_martin_unw72_funestus
results_martinf_w=results_martin_w72_funestus

exposure=get_in_out_exp(activity_cycles = activity_p, vec_p = ent_params)$Exposure_Indoor_whileinbed


# bednet efficacy
impacts_efficacy=calculate_impact_4EHT(washedDecay=F, vectormodel = model_noRhythms)


# bednet efficacy + usage
impacts_EfficacyUsage=calculate_impact_4EHT(washedDecay=F, vectormodel = model_noRhythms,
                                                     coverage_vector = coverage_vector_usage)

# bednet efficacy + usage + attrition
impacts_EfficacyUsageAttrition=calculate_impact_4EHT(L_vector=L_vector_attrition,
                                                kappa_vector=kappa_vector_attrition,
                                                washedDecay=F, vectormodel = model_noRhythms,
                                                coverage_vector = coverage_vector_usage)



# bednet efficacy + attrition + insecticide decay
impacts_EfficacyAttritionDecayUsage=calculate_impact_4EHT(L_vector=L_vector_attrition,
                                                     kappa_vector=kappa_vector_attrition,
                                                     washedDecay=T, vectormodel = model_noRhythms,
                                                     coverage_vector = coverage_vector_usage)


impacts_EfficacyAttritionDecayUsageRhythms=calculate_impact_4EHT(L_vector=L_vector_attrition,
                                                             kappa_vector=kappa_vector_attrition,
                                                             washedDecay=T, vectormodel = model_noRhythms,
                                                             coverage_vector = coverage_vector_usage, inbed_exposure = exposure)


VCreduc_all=rbind(impacts_efficacy %>% mutate(scenario="Efficacy") ,
                  impacts_EfficacyUsage %>% mutate(scenario="EfficacyUsage"),
                  impacts_EfficacyUsageAttrition %>% mutate(scenario="impacts_EfficacyUsageAttrition") ,
                  impacts_EfficacyAttritionDecayUsage %>% mutate(scenario="EfficacyAttritionUsageDecay") ,
                  impacts_EfficacyAttritionDecayUsageRhythms  %>% mutate(scenario="EfficacyAttritionUsageDecayRhythms")
)

cascade_ig2_bit103=plot_bars_effectiveness_LLIN(df_VCred=VCreduc_all %>% filter(EHT=="BIT103", insecticide=="IG2") %>%
                                                  select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                colorfinal="orange")
cascade_pbo_bit103=plot_bars_effectiveness_LLIN(df_VCred=VCreduc_all %>% filter(EHT=="BIT103", insecticide=="OlysetPlus") %>%
                                                  select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                colorfinal="dodgerblue")
cascade_pbo_bit055=plot_bars_effectiveness_LLIN(df_VCred=VCreduc_all %>% filter(EHT=="BIT055", insecticide=="OlysetPlus") %>%
                                                  select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                colorfinal="dodgerblue")
cascade_ig2_kibondo=plot_bars_effectiveness_LLIN(df_VCred=VCreduc_all %>% filter(EHT=="Kibondo", insecticide=="IG2") %>%
                                                   select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                 colorfinal="orange")
cascade_pbo_bit059=plot_bars_effectiveness_LLIN(df_VCred=VCreduc_all %>% filter(EHT=="BIT059", insecticide=="OlysetPlus") %>%
                                                  select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                colorfinal="dodgerblue")
cascade_ig2_bit080=plot_bars_effectiveness_LLIN(df_VCred=VCreduc_all %>% filter(EHT=="BIT080", insecticide=="IG2") %>%
                                                  select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                colorfinal="orange")
cascade_ig2_martin=plot_bars_effectiveness_LLIN(df_VCred=VCreduc_all %>% filter(EHT=="Martin", insecticide=="IG2") %>%
                                                  select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                colorfinal="orange")
cascade_pbo_martin=plot_bars_effectiveness_LLIN(df_VCred=VCreduc_all %>% filter(EHT=="Martin", insecticide=="OlysetPlus") %>%
                                                  select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                colorfinal="dodgerblue")

cascade_ig2_martinf=plot_bars_effectiveness_LLIN(df_VCred=VCreduc_all %>% filter(EHT=="Martin, f", insecticide=="IG2") %>%
                                                  select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                colorfinal="orange")
cascade_pbo_martinf=plot_bars_effectiveness_LLIN(df_VCred=VCreduc_all %>% filter(EHT=="Martin, f", insecticide=="OlysetPlus") %>%
                                                  select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                colorfinal="dodgerblue")
IG2=plot_grid(ggdraw() +draw_label(  "Interceptor G2", x=0.1, y=0.9, fontface = "bold" )+theme_classic(),
              cascade_ig2_bit103,cascade_ig2_kibondo,cascade_ig2_bit080,cascade_ig2_martin,cascade_ig2_martinf,
              ncol = 1, byrow = F,
              labels=c("","BIT103", "Kibondo","BIT080", "Martin, gambiae", "Martin, funestus"),
              label_size = 12, rel_heights = c(0.15, 1,1,1,1,1),
              hjust = 0, label_x = 0.85)#+draw_plot_label("Interceptor G2", vjust = 0.5, hjust = -0.3)

OlysetPlus=plot_grid(ggdraw() +draw_label(  "Olyset Plus", x=0.05, y=0.9, fontface = "bold" )+theme_classic(),
                     cascade_pbo_bit103,cascade_pbo_bit055,cascade_pbo_bit059,cascade_pbo_martin,cascade_pbo_martinf,
                     ncol = 1, byrow = F,
                     labels=c( "","BIT103","BIT055",  "Odufuwa", "Martin, gambiae", "Martin, funestus"),
                     label_size = 12, rel_heights = c(0.15, 1,1,1,1,1),
                     hjust = 0, label_x = 0.85)#+draw_plot_label("Olyset Plus", vjust = 0.5, hjust = -0.3)

plot_grid(IG2, OlysetPlus,
          ncol = 2, byrow = F,
          labels=c("Interceptor G2", "Olyset Plus"),
          label_size = 12,
          hjust = 0, label_x = 0.4)
ggsave(file.path(plotDir, "cascades_LLIN_noUncertainty_perEHT.png"), width=12, height=16)


#####################################################################
# CASCADE WITH UNCERTAINTY
####################################################################

results_kibondo_w=results_kibondo_w_72
results_kibondo_unw=results_kibondo_unw_72
results_bit103_unw=results_bit103_unw_72
results_bit103_w=results_bit103_w_72
results_bit080_unw=results_bit080_unw_72
results_bit080_w=results_bit080_w_72
results_martin_unw=results_martin_unw72
results_martin_w=results_martin_w72
results_martinf_unw=results_martin_unw72_funestus
results_martinf_w=results_martin_w72_funestus

nsamples=1000

if(rerun){
  # bednet efficacy
  impacts_efficacy_uncertainty=calculate_impact_4EHT(washedDecay=F, vectormodel = model_noRhythms,uncertainty = T, nsamples = nsamples)
  
  
  # bednet efficacy + attrition
  impacts_EfficacyUsage_uncertainty=calculate_impact_4EHT(washedDecay=F,coverage_vector = coverage_vector_usage,
                                                          vectormodel = model_noRhythms,uncertainty = T, nsamples = nsamples)
  
  # bednet efficacy + attrition + usage
  impacts_EfficacyAttritionUsage_uncertainty=calculate_impact_4EHT(L_vector=L_vector_attrition,
                                                                   kappa_vector=kappa_vector_attrition,
                                                                   washedDecay=F, vectormodel = model_noRhythms,uncertainty = T,
                                                                   coverage_vector = coverage_vector_usage, nsamples = nsamples)
  # bednet efficacy + attrition + insecticide decay
  impacts_EfficacyAttritionDecay_uncertainty=calculate_impact_4EHT(L_vector=L_vector_attrition,
                                                                   kappa_vector=kappa_vector_attrition,
                                                                   washedDecay=T, vectormodel = model_noRhythms,uncertainty = T,
                                                                   coverage_vector = coverage_vector_usage, nsamples = nsamples)
  
  impacts_EfficacyAttritionDecayRhythms_uncertainty=calculate_impact_4EHT(L_vector=L_vector_attrition,
                                                                          kappa_vector=kappa_vector_attrition, inbed_exposure = exposure,
                                                                          washedDecay=T, vectormodel = model_noRhythms,uncertainty = T,
                                                                          coverage_vector = coverage_vector_usage, nsamples = nsamples)
  
  VCreduc_all_uncertainty=rbind(impacts_efficacy_uncertainty %>% mutate(scenario="Efficacy") ,
                                impacts_EfficacyUsage_uncertainty %>% mutate(scenario="EfficacyAAUsage"),
                                impacts_EfficacyAttritionUsage_uncertainty %>% mutate(scenario="EfficacyAttritionUsage") ,
                                impacts_EfficacyAttritionDecay_uncertainty %>% mutate(scenario="EfficacyAttritionUsageDecay") ,
                                impacts_EfficacyAttritionDecayRhythms_uncertainty  %>% mutate(scenario="EfficacyAttritionUsageDecayRhythms")
  )
  write.csv(VCreduc_all_uncertainty, file = file.path(plotDir, "VCreduc_all_uncertainty.csv"), row.names = F)
}

VCreduc_all_uncertainty=read.csv( file = file.path(plotDir, "VCreduc_all_uncertainty.csv"))

df_VCreduc_all_uncertainty=VCreduc_all_uncertainty %>% group_by(insecticide,  scenario)%>%
  summarise(q025=quantile(VCred, probs=0.025),q975=quantile(VCred, probs=0.975),
            min=min(VCred),max=max(VCred),
            med=quantile(VCred, probs=0.5), mean=mean(VCred) ) %>% ungroup()


cascade_ig2_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_all_uncertainty %>% filter( insecticide=="IG2") %>%
                                                                   select(-insecticide, -min, -max, -med)%>%
                                                                   rename(value=mean, variable=scenario) %>% mutate(value=100*value),
                                                                 colorfinal="orange")
cascade_pbo_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_all_uncertainty %>% filter( insecticide=="OlysetPlus") %>%
                                                                   select(-insecticide, -min, -max, -med)%>%
                                                                   rename(value=mean, variable=scenario) %>% mutate(value=100*value),
                                                                 colorfinal="dodgerblue")

plot_grid(cascade_ig2_uncertainty,
          cascade_pbo_uncertainty,
          ncol = 1,
          labels=c("Interceptor G2", "Olyset Plus"),
          label_size = 12,
          hjust = 0, label_x = 0.75)
ggsave(file.path(plotDir, "cascades_LLIN_uncertainty.png"), width=6, height=8)

plot_grid(cascade_ig2_uncertainty,
          cascade_pbo_uncertainty,
          ncol = 2,
          labels=c("Interceptor G2", "Olyset Plus"),
          label_size = 12,
          hjust = 0, label_x = 0.75)
ggsave(file.path(plotDir, "cascades_LLIN_uncertainty_horizontal.png"), width=12, height=4)



### Cascades per EHT with uncertainty
df_VCreduc_EHT_uncertainty=VCreduc_all_uncertainty %>% group_by(insecticide,  scenario, EHT)%>%
  summarise(q025=quantile(VCred, probs=0.025),q975=quantile(VCred, probs=0.975),
            min=min(VCred),max=max(VCred),
            med=quantile(VCred, probs=0.5), mean=mean(VCred) ) %>% ungroup()


cascade_ig2_bit103_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_EHT_uncertainty %>% filter(EHT=="BIT103", insecticide=="IG2") %>%
                                                                          select(-insecticide, -min, -max, -EHT, -med)%>%
                                                                          rename(value=mean, variable=scenario) %>% mutate(value=100*value),
                                                                        colorfinal="orange")
cascade_pbo_bit103_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_EHT_uncertainty %>% filter(EHT=="BIT103", insecticide=="OlysetPlus") %>%
                                                                          select(-insecticide, -min, -max, -EHT, -med)%>%
                                                                          rename(value=mean, variable=scenario) %>% mutate(value=100*value),
                                                                        colorfinal="dodgerblue")
cascade_pbo_bit055_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_EHT_uncertainty %>% filter(EHT=="BIT055", insecticide=="OlysetPlus") %>%
                                                                          select(-insecticide, -min, -max, -EHT, -med)%>%
                                                                          rename(value=mean, variable=scenario) %>% mutate(value=100*value),
                                                                        colorfinal="dodgerblue")
cascade_ig2_kibondo_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_EHT_uncertainty %>% filter(EHT=="Kibondo", insecticide=="IG2") %>%
                                                                           select(-insecticide, -min, -max, -EHT, -med)%>%
                                                                           rename(value=mean, variable=scenario)%>% mutate(value=100*value),
                                                                         colorfinal="orange")

cascade_pbo_bit059_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_EHT_uncertainty %>% filter(EHT=="BIT059", insecticide=="OlysetPlus") %>%
                                                                          select(-insecticide, -min, -max, -EHT, -med)%>%
                                                                          rename(value=mean, variable=scenario) %>% mutate(value=100*value),
                                                                        colorfinal="dodgerblue")
cascade_ig2_bit080_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_EHT_uncertainty %>% filter(EHT=="BIT080", insecticide=="IG2") %>%
                                                                          select(-insecticide, -min, -max, -EHT, -med)%>%
                                                                          rename(value=mean, variable=scenario)%>% mutate(value=100*value),
                                                                        colorfinal="orange")

cascade_pbo_martin_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_EHT_uncertainty %>% filter(EHT=="Martin", insecticide=="OlysetPlus") %>%
                                                                          select(-insecticide, -min, -max, -EHT, -med)%>%
                                                                          rename(value=mean, variable=scenario) %>% mutate(value=100*value),
                                                                        colorfinal="dodgerblue")
cascade_ig2_martin_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_EHT_uncertainty %>% filter(EHT=="Martin", insecticide=="IG2") %>%
                                                                          select(-insecticide, -min, -max, -EHT, -med)%>%
                                                                          rename(value=mean, variable=scenario)%>% mutate(value=100*value),
                                                                        colorfinal="orange")

cascade_pbo_martinf_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_EHT_uncertainty %>% filter(EHT=="Martin, f", insecticide=="OlysetPlus") %>%
                                                                          select(-insecticide, -min, -max, -EHT, -med)%>%
                                                                          rename(value=mean, variable=scenario) %>% mutate(value=100*value),
                                                                        colorfinal="dodgerblue")
cascade_ig2_martinf_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_EHT_uncertainty %>% filter(EHT=="Martin, f", insecticide=="IG2") %>%
                                                                          select(-insecticide, -min, -max, -EHT, -med)%>%
                                                                          rename(value=mean, variable=scenario)%>% mutate(value=100*value),
                                                                        colorfinal="orange")

plot_grid(cascade_ig2_bit103_uncertainty,cascade_ig2_kibondo_uncertainty, cascade_ig2_bit080_uncertainty, cascade_ig2_martin_uncertainty, cascade_ig2_martinf_uncertainty,
          cascade_pbo_bit103_uncertainty,cascade_pbo_bit055_uncertainty,cascade_pbo_bit059_uncertainty,cascade_pbo_martin_uncertainty,cascade_pbo_martinf_uncertainty,
          ncol = 2, byrow = F,
          labels=c("Interceptor G2 (Assenga)", "Olyset Plus (Assenga)",
                   "Interceptor G2 (Kibondo)",  "Olyset Plus (BIT055)",
                   "Interceptor G2 (BIT080)",  "Olyset Plus (Odufuwa)",
                   "Interceptor G2 (Martin, An. gambiae)",  "Olyset Plus (Martin, An. gambiae)",
                   "Interceptor G2 (Martin, An. funestus)",  "Olyset Plus (Martin, An. funestus)"),
          label_size = 12,
          hjust = 0, label_x = 0.45)
ggsave(file.path(plotDir, "cascades_LLIN_perEHT_uncertainty.png"), width=12, height=18)


###############################
# SENSITIVITY ON CASCADE ORDER


scenarios=list()
scenarios$nb1=c("attrition", "usage", "exposure", "decay")
scenarios$nb2=c("attrition", "usage", "exposure", "decay")
scenarios$nb3=c("attrition", "usage", "exposure", "decay")
scenarios$nb4=c("attrition", "usage", "exposure", "decay")


full_factorial=expand.grid(scenarios) %>%
  filter(nb2 !=nb1, nb3 !=nb1, nb4!=nb1,
         nb3 !=nb2, nb4 !=nb2,
         nb4 != nb3)


calculate_impact_scenario=function(myvec){
  my_washedDecay=("decay" %in% myvec)
  my_inbed=ifelse("exposure" %in% myvec,exposure, 1)
  
  if("usage" %in% myvec){
    my_coverage_vector=coverage_vector_usage
  } else {
    my_coverage_vector=c("ig2"=1,  "pbo"=1)
  }
  
  if("attrition" %in% myvec){
    my_L_vector=L_vector_attrition
    my_kappa_vector=kappa_vector_attrition
  } else {
    my_L_vector=c("ig2"=3, "pbo"=3)
    my_kappa_vector=c("ig2"=2,  "pbo"=2)
  }
  
  out=calculate_impact_4EHT(L_vector=my_L_vector,
                            kappa_vector=my_kappa_vector, washedDecay=my_washedDecay, uncertainty=F, vectormodel=model_noRhythms, nsamples=100,
                            coverage_vector=my_coverage_vector, inbed_exposure=my_inbed)
  return(out)
}



vec=full_factorial[1,]
calculate_impact_cascade=function(vec){
  # step0
  impacts_efficacy=calculate_impact_scenario(c(""))
  
  # step1
  impacts_1=calculate_impact_scenario(vec[1])
  
  # step2
  impacts_2=calculate_impact_scenario(vec[1:2])
  
  # step3
  impacts_3=calculate_impact_scenario(vec[1:3])
  
  # step4
  impacts_4=calculate_impact_scenario(vec[1:4])
  
  
  VCreduc_all=rbind(impacts_efficacy %>% mutate(scenario="0") ,
                    impacts_1 %>% mutate(scenario="1"),
                    impacts_2 %>% mutate(scenario="2") ,
                    impacts_3 %>% mutate(scenario="3") ,
                    impacts_4  %>% mutate(scenario="4")
  )
  return(VCreduc_all)
}




plot_cascades_order=function(vec, path, name){
  
  
  vec_rename=case_match(t(vec),
                        "attrition"~"Functional\nsurvival",
                        "usage"~"Usage at\ndistribution",
                        "decay"~"Insecticidal\ndurability",
                        "exposure"~"In-bed\nexposure", .default = NA)
  
  VCreduc_all=calculate_impact_cascade(t(vec))
  
  df_VCreduc_all=VCreduc_all %>% group_by(insecticide,  scenario)%>%
    summarise(VCred=mean(VCred) ) %>% ungroup()
  
  
  cascade_ig2=plot_bars_effectiveness_LLIN_sensitivity(df_VCred=df_VCreduc_all %>% filter( insecticide=="IG2") %>%
                                               select(scenario, VCred)%>% 
                                               rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                             colorfinal="orange", order_cascade = vec_rename)
  cascade_pbo=plot_bars_effectiveness_LLIN_sensitivity(df_VCred=df_VCreduc_all %>% filter( insecticide=="OlysetPlus")  %>%
                                               select(scenario, VCred)%>% 
                                               rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                             colorfinal="dodgerblue", order_cascade = vec_rename)
  
  plot_grid(cascade_ig2$p,
            cascade_pbo$p,
            ncol = 1,
            labels=c("Interceptor G2", "Olyset Plus"),
            label_size = 12,
            hjust = 0, label_x = 0.75)
  ggsave(file.path(path, paste0("cascades_LLIN_",name,".png")), width=6, height=8)
  
  db_mean=rbind(cascade_ig2$db %>%mutate(net="IG2"), cascade_pbo$db %>%mutate(net="PBO"))
  
  cascade_ig2_bit103=plot_bars_effectiveness_LLIN_sensitivity(df_VCred=VCreduc_all %>% filter(EHT=="BIT103", insecticide=="IG2") %>%
                                                      select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                    colorfinal="orange", order_cascade = vec_rename)
  cascade_pbo_bit103=plot_bars_effectiveness_LLIN_sensitivity(df_VCred=VCreduc_all %>% filter(EHT=="BIT103", insecticide=="OlysetPlus") %>%
                                                      select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                    colorfinal="dodgerblue", order_cascade = vec_rename)
  cascade_pbo_bit055=plot_bars_effectiveness_LLIN_sensitivity(df_VCred=VCreduc_all %>% filter(EHT=="BIT055", insecticide=="OlysetPlus") %>%
                                                      select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                    colorfinal="dodgerblue", order_cascade = vec_rename)
  cascade_ig2_kibondo=plot_bars_effectiveness_LLIN_sensitivity(df_VCred=VCreduc_all %>% filter(EHT=="Kibondo", insecticide=="IG2") %>%
                                                       select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                     colorfinal="orange", order_cascade = vec_rename)
  cascade_pbo_bit059=plot_bars_effectiveness_LLIN_sensitivity(df_VCred=VCreduc_all %>% filter(EHT=="BIT059", insecticide=="OlysetPlus") %>%
                                                      select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                    colorfinal="dodgerblue", order_cascade = vec_rename)
  cascade_ig2_bit080=plot_bars_effectiveness_LLIN_sensitivity(df_VCred=VCreduc_all %>% filter(EHT=="BIT080", insecticide=="IG2") %>%
                                                      select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                    colorfinal="orange", order_cascade = vec_rename)
  cascade_ig2_martin=plot_bars_effectiveness_LLIN_sensitivity(df_VCred=VCreduc_all %>% filter(EHT=="Martin", insecticide=="IG2") %>%
                                                      select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                    colorfinal="orange", order_cascade = vec_rename)
  cascade_pbo_martin=plot_bars_effectiveness_LLIN_sensitivity(df_VCred=VCreduc_all %>% filter(EHT=="Martin", insecticide=="OlysetPlus") %>%
                                                      select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                    colorfinal="dodgerblue", order_cascade = vec_rename)
  
  cascade_ig2_martinf=plot_bars_effectiveness_LLIN_sensitivity(df_VCred=VCreduc_all %>% filter(EHT=="Martin, f", insecticide=="IG2") %>%
                                                       select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                     colorfinal="orange", order_cascade = vec_rename)
  cascade_pbo_martinf=plot_bars_effectiveness_LLIN_sensitivity(df_VCred=VCreduc_all %>% filter(EHT=="Martin, f", insecticide=="OlysetPlus") %>%
                                                       select(scenario, VCred)%>% rename(value=VCred, variable=scenario) %>% mutate(value=100*value),
                                                     colorfinal="dodgerblue", order_cascade = vec_rename)
  db_eht=rbind(cascade_ig2_bit103$db %>%mutate(net="IG2", eht="Assenga et al."), 
               cascade_pbo_bit103$db %>%mutate(net="PBO", eht="Assenga et al."), 
               cascade_pbo_bit055$db %>%mutate(net="PBO", eht="BIT055"), 
               cascade_ig2_kibondo$db %>%mutate(net="IG2", eht="Kibondo et al."), 
               cascade_pbo_bit059$db %>%mutate(net="PBO", eht="BIT059"), 
               cascade_ig2_bit080$db %>%mutate(net="IG2", eht="Odufuwa et al."), 
               cascade_ig2_martin$db %>%mutate(net="IG2", eht="Martin et al. An. gambiae"), 
               cascade_pbo_martin$db %>%mutate(net="PBO", eht="Martin et al. An. gambiae"), 
               cascade_ig2_martinf$db %>%mutate(net="IG2", eht="Martin et al. An. funestus"), 
               cascade_pbo_martinf$db %>%mutate(net="PBO", eht="Martin et al. An. funestus")
  )
  
  return(list(db_mean, db_eht))
}

df_mean=data.frame()
df_eht=data.frame()
for(i in 1:nrow(full_factorial)){
  this.df=plot_cascades_order(full_factorial[i,],path=file.path(plotDir, "sensitivity_cascade"), name = i)
  df_mean=rbind(df_mean, this.df[[1]] %>% mutate(permut=i))
  df_eht=rbind(df_eht, this.df[[2]]%>% mutate(permut=i))
}

df_mean2=df_mean %>%
  filter(!variable %in% c("Perfect\nvector control", "Effectiveness", "Entomological\nefficacy")) %>% 
  arrange(net, permut) %>% group_by(net, permut)%>% mutate(rank=rank(-point_difference, ties.method = "first"))%>% 
  select(net, variable, point_difference, rank)

df_mean2%>%
  group_by(variable, net) %>%
  summarise(mean=mean(point_difference), min=min(point_difference), max=max(point_difference),
            mean_rank=mean(rank), min_rank=min(rank), max_rank=max(rank))%>% 
  mutate(diff=max-min,diff_rank=max_rank-min_rank)

df_mean2%>%
  group_by(variable, net) %>%
  summarise(min=min(point_difference), max=max(point_difference),
            min_rank=min(rank), max_rank=max(rank)) %>%
  mutate("Point difference"=paste0(min,"-", max), 
         "Rank"=paste0(min_rank,"-", max_rank))%>%
  select(variable, net, "Point difference", Rank)%>%
  rename(" "=variable)%>%
  group_by(net)%>% gt()%>%
  gt::gtsave(filename = file.path(plotDir, "sensitivity_cascade/cascade_sensitivity_average.docx"))


df_eht2=df_eht  %>% select(net, eht, variable, permut, point_difference)%>% filter(!variable %in% c("Perfect\nvector control", "Effectiveness", "Entomological\nefficacy")) %>%
  arrange(net, eht, permut) %>% group_by(net, eht,permut)%>% mutate(rank=rank(-point_difference, ties.method = "first"))

df_eht2$eht=factor(df_eht2$eht, levels = c("Assenga et al.","Kibondo et al.", "Odufuwa et al." ,"BIT055" ,"BIT059", "Martin et al. An. funestus", "Martin et al. An. gambiae"  ))

df_eht2%>%
  group_by(net, variable, eht) %>%
  summarise(min=min(point_difference), max=max(point_difference),
            min_rank=min(rank), max_rank=max(rank)) %>%
  mutate("Point difference"=paste0(min,"-", max), 
         "Rank"=paste0(min_rank,"-", max_rank))%>%
  arrange(eht, variable) %>%
  select(eht, variable, "Point difference", Rank)%>%
  pivot_wider(id_cols = c(eht, variable), names_from = net, values_from = c("Point difference", Rank))%>%
  select(eht, variable, "Point difference_IG2", Rank_IG2, "Point difference_PBO", Rank_PBO )%>%
  group_by( eht)%>% gt()%>%
  gt::gtsave(filename = file.path(plotDir, "sensitivity_cascade/cascade_sensitivity_perEHT.docx"))


df_eht2%>%
  filter(net=="IG2")%>%
  group_by(variable, eht) %>%
  summarise(min=min(point_difference), max=max(point_difference),
            min_rank=min(rank), max_rank=max(rank)) %>%
  mutate("Point difference"=paste0(min,"-", max), 
         "Rank"=paste0(min_rank,"-", max_rank))%>%
  arrange(eht, variable) %>%
  select(eht, variable, "Point difference", Rank)%>%
  rename("IG2"=variable)%>%
  group_by(eht)%>% gt()%>%
  gt::gtsave(filename = file.path(plotDir, "sensitivity_cascade/cascade_sensitivity_perEHT_IG2.docx"))

df_eht2%>%
  filter(net=="PBO")%>%
  group_by(variable, eht) %>%
  summarise(min=min(point_difference), max=max(point_difference),
            min_rank=min(rank), max_rank=max(rank)) %>%
  mutate("Point difference"=paste0(min,"-", max), 
         "Rank"=paste0(min_rank,"-", max_rank))%>%
  arrange(eht, variable) %>%
  select(eht, variable, "Point difference", Rank)%>%
  rename("Olyset Plus"=variable)%>%
  group_by(eht)%>% gt()%>%
  gt::gtsave(filename = file.path(plotDir, "sensitivity_cascade/cascade_sensitivity_perEHT_PBO.docx"))

