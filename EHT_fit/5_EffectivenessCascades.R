rm(list=ls())
library(dplyr)
library(ggplot2)
library(readxl)
library(lubridate)
library(stringr)
library(cowplot)
library(tidyr)
library(AnophelesModel)

mainDir="."
stanDir=file.path(mainDir, "stan_outputs")
plotDir=file.path(mainDir, "plots")
dataDir=file.path(mainDir, "data")
outputsDir=file.path(mainDir, "csv_outputs")
scriptDir=file.path(mainDir,"EHT_fit")

source(file.path(scriptDir, "compute_vectorialCapacity.R"))
source(file.path(scriptDir, "cascade_LLIN.R"))


results_kibondo_w_72<- readRDS(file.path(stanDir,"/stan_kibondo_washed_controlUnw_72.rds"))
results_kibondo_unw_72<- readRDS(file.path(stanDir,"/stan_kibondo_unwashed_72.rds"))
results_kibondo_w_24<- readRDS(file.path(stanDir,"/stan_kibondo_washed_controlUnw_24.rds"))
results_kibondo_unw_24<- readRDS(file.path(stanDir,"/stan_kibondo_unwashed_24.rds"))


results_bit055_w<- readRDS(file.path(stanDir,"/stan_BIT055_washed_controlUnw_24.rds"))
results_bit055_unw<- readRDS(file.path(stanDir,"/stan_BIT055_unwashed_24.rds"))


results_bit103_unw_72<- readRDS(file.path(stanDir,"/stan_BIT103_unwashed_72_ifakara.rds"))
results_bit103_w_72<- readRDS(file.path(stanDir,"/stan_BIT103_washed_controlUnw_72_ifakara.rds"))
results_bit103_unw_24<- readRDS(file.path(stanDir,"/stan_BIT103_unwashed_24_ifakara.rds"))
results_bit103_w_24<- readRDS(file.path(stanDir,"/stan_BIT103_washed_controlUnw_24_ifakara.rds"))

results_bit059_w<- readRDS(file.path(stanDir,"/stan_BIT059_washed_controlUnw_24.rds"))
results_bit059_unw<- readRDS(file.path(stanDir,"/stan_BIT059_unwashed_24.rds"))

results_bit080_w_72<- readRDS(file.path(stanDir,"/stan_BIT080_washed_controlUnw_72.rds"))
results_bit080_unw_72<- readRDS(file.path(stanDir,"/stan_BIT080_unwashed_72.rds"))

# from fitting survey data on LLIN usage from Mosha et al.
L_vector_attrition=c("kibondo"=2.4, "bit055"=1.7,"bit103_PBO"=1.7, "bit103_IG2"=2.4, "bit059"=1.7, "bit080"=2.4)
kappa_vector_attrition=c("kibondo"=2.4,"bit055"=1.9,"bit103_PBO"=1.9, "bit103_IG2"=2.4, "bit059"=1.9, "bit080"=2.4)
coverage_vector_usage = c("kibondo"=0.69,  "bit055"=0.75,"bit103_PBO"=0.75, "bit103_IG2"=0.69, "bit059"=0.75, "bit080"=0.69)

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
model_p = build_model_obj(vec_p=ent_params, hosts_p= def_host_params(), activity=activity_p, total_pop=2000)
model_noRhythms = build_model_obj(vec_p=ent_params, hosts_p= def_host_params(), activity=activity_noRhythms, total_pop=2000)

#####################################################################
# COMPUTE VECTORIAL CAPACITY REDUCTION
####################################################################

calculate_impact_4EHT=function(L_vector=c("kibondo"=3, "bit055"=3, "bit059"=3, "bit080"=3,"bit103_PBO"=3, "bit103_IG2"=3), 
                               kappa_vector=c("kibondo"=2,  "bit055"=2,"bit059"=2,"bit080"=2,"bit103_PBO"=2, "bit103_IG2"=2), washedDecay=T, uncertainty=F, vectormodel=model_p,nsamples=100,
                               coverage_vector=c("kibondo"=1,  "bit055"=1,"bit059"=1,"bit080"=1,"bit103_PBO"=1, "bit103_IG2"=1)){
 
  npoints=3*365
  pyrethroid_params=c(0, 0.05, 0)
  names(pyrethroid_params)=c("Deterrency","PrePrandial", "PostPrandial" )
  
  pyrethroid_interv_kibondo=get_interv(myproba=pyrethroid_params,L=L_vector[["kibondo"]]*365, kappa=kappa_vector[["kibondo"]], npoints=3*365,
                                       model_p=vectormodel, name="pyrethroid", intervention_type="LLINs", decay="weibull",halflife_insecticide=NULL)
  
  pyrethroid_interv_bit055=get_interv(myproba=pyrethroid_params,L=L_vector[["bit055"]]*365, kappa=kappa_vector[["bit055"]], npoints=3*365,
                                      model_p=vectormodel, name="pyrethroid", intervention_type="LLINs", decay="weibull",halflife_insecticide=NULL)
  
  pyrethroid_interv_bit103_PBO=get_interv(myproba=pyrethroid_params,L=L_vector[["bit103_PBO"]]*365, kappa=kappa_vector[["bit103_PBO"]], npoints=3*365,
                                          model_p=vectormodel, name="pyrethroid", intervention_type="LLINs", decay="weibull",halflife_insecticide=NULL)
  
  pyrethroid_interv_bit103_IG2=get_interv(myproba=pyrethroid_params,L=L_vector[["bit103_IG2"]]*365, kappa=kappa_vector[["bit103_IG2"]], npoints=3*365,
                                          model_p=vectormodel, name="pyrethroid", intervention_type="LLINs", decay="weibull",halflife_insecticide=NULL)
  
  pyrethroid_interv_bit059=get_interv(myproba=pyrethroid_params,L=L_vector[["bit059"]]*365, kappa=kappa_vector[["bit059"]], npoints=3*365,
                                      model_p=vectormodel, name="pyrethroid", intervention_type="LLINs", decay="weibull",halflife_insecticide=NULL)
  
  pyrethroid_interv_bit080=get_interv(myproba=pyrethroid_params,L=L_vector[["bit080"]]*365, kappa=kappa_vector[["bit080"]], npoints=3*365,
                                      model_p=vectormodel, name="pyrethroid", intervention_type="LLINs", decay="weibull",halflife_insecticide=NULL)
  
  intervention_list_kibondo <- interventions_vectorial_capacity_wrapper(results=results_kibondo_unw,results_washed20 =results_kibondo_w,decay = "weibull",
                                                                        model_p=vectormodel,
                                                                        insecticides=c(2),L = L_vector[["kibondo"]]*365, kappa=kappa_vector[["kibondo"]],
                                                                        names=c("IG2 (Kibondo)"),  npoints=npoints,washedDecay=washedDecay, uncertainty=uncertainty, 
                                                                        cov=coverage_vector[["kibondo"]], pyrethroid = pyrethroid_interv_kibondo, nsamples=nsamples)
  
  
  intervention_list_bit055 <- interventions_vectorial_capacity_wrapper(results=results_bit055_unw,results_washed20=results_bit055_w,decay = "weibull",
                                                                       model_p=vectormodel,
                                                                       insecticides=c(2),L =  L_vector[["bit055"]]*365 , kappa=kappa_vector[["bit055"]],
                                                                       names=c("OlysetPlus (BIT055)"),  npoints=npoints, washedDecay=washedDecay, uncertainty=uncertainty,
                                                                       cov=coverage_vector[["bit055"]], pyrethroid = pyrethroid_interv_bit055, nsamples=nsamples)
  
  intervention_list_bit103_PBO <- interventions_vectorial_capacity_wrapper(results=results_bit103_unw,results_washed20=results_bit103_w,decay = "weibull",
                                                                           model_p=vectormodel,
                                                                           insecticides=c(2),L =  L_vector[["bit103_PBO"]]*365 , kappa=kappa_vector[["bit103_PBO"]],
                                                                           names=c( "OlysetPlus (BIT103)"),  npoints=npoints,washedDecay=washedDecay, uncertainty=uncertainty,
                                                                           cov=coverage_vector[["bit103_PBO"]], pyrethroid = pyrethroid_interv_bit103_PBO, nsamples=nsamples)
  
  intervention_list_bit103_IG2 <- interventions_vectorial_capacity_wrapper(results=results_bit103_unw,results_washed20=results_bit103_w,decay = "weibull",
                                                                           model_p=vectormodel,
                                                                           insecticides=c(3),L =  L_vector[["bit103_IG2"]]*365 , kappa=kappa_vector[["bit103_IG2"]],
                                                                           names=c("IG2 (BIT103)"),  npoints=npoints,washedDecay=washedDecay, uncertainty=uncertainty,
                                                                           cov=coverage_vector[["bit103_IG2"]], pyrethroid = pyrethroid_interv_bit103_IG2, nsamples=nsamples)
  
  intervention_list_bit059 <- interventions_vectorial_capacity_wrapper(results=results_bit059_unw,results_washed20=results_bit059_w,decay = "weibull",
                                                                       model_p=vectormodel,
                                                                       insecticides=c(2),L =  L_vector[["bit059"]]*365 , kappa=kappa_vector[["bit059"]],
                                                                       names=c("OlysetPlus (BIT059)"),  npoints=npoints, washedDecay=washedDecay, uncertainty=uncertainty,
                                                                       cov=coverage_vector[["bit059"]], pyrethroid = pyrethroid_interv_bit059, nsamples=nsamples)
  
  intervention_list_bit080 <- interventions_vectorial_capacity_wrapper(results=results_bit080_unw,results_washed20=results_bit080_w,decay = "weibull",
                                                                       model_p=vectormodel,
                                                                       insecticides=c(2),L =  L_vector[["bit080"]]*365 , kappa=kappa_vector[["bit080"]],
                                                                       names=c("IG2 (BIT080)"),  npoints=npoints, washedDecay=washedDecay, uncertainty=uncertainty,
                                                                       cov=coverage_vector[["bit080"]], pyrethroid = pyrethroid_interv_bit055, nsamples=nsamples)
  
  
  
  if(uncertainty==T){
    impacts=rbind(data.frame(VCred=intervention_list_kibondo) %>% mutate(EHT="Kibondo", insecticide="IG2"),
                  data.frame(VCred=intervention_list_bit055) %>% mutate(EHT="BIT055", insecticide="PBO"),
                  data.frame(VCred=intervention_list_bit103_IG2) %>% mutate(EHT="BIT103", insecticide="IG2"),
                  data.frame(VCred=intervention_list_bit103_PBO) %>% mutate(EHT="BIT103", insecticide="PBO"),
                  data.frame(VCred=intervention_list_bit059) %>% mutate(EHT="BIT059", insecticide="PBO"),
                  data.frame(VCred=intervention_list_bit080) %>% mutate(EHT="BIT080", insecticide="IG2")
    )
    
    
  } else{
    
    # calculate the impact
    impacts_kibondo =  calculate_combined_impact(intervention1=intervention_list_kibondo[[1]],
                                                intervention2=pyrethroid_interv_kibondo[[1]],
                                                cov_intervention1 =  coverage_vector[["kibondo"]],
                                                cov_intervention2 = coverage_vector[["kibondo"]],
                                                cov_combination = coverage_vector[["kibondo"]],
                                                N_vec =10000)
    
    impacts_bit055 = calculate_combined_impact(intervention1=intervention_list_bit055[[1]],
                                               intervention2=pyrethroid_interv_bit055[[1]],
                                               cov_intervention1 = coverage_vector[["bit055"]],
                                               cov_intervention2 = coverage_vector[["bit055"]],
                                               cov_combination = coverage_vector[["bit055"]],
                                               N_vec =10000)
    
    impacts_bit103_pbo = calculate_combined_impact(intervention1=intervention_list_bit103_PBO[[1]],
                                                   intervention2=pyrethroid_interv_bit103_PBO[[1]],
                                                   cov_intervention1 = coverage_vector[["bit103_PBO"]],
                                                   cov_intervention2 = coverage_vector[["bit103_PBO"]],
                                                   cov_combination = coverage_vector[["bit103_PBO"]],
                                                   N_vec =10000)
    
    impacts_bit103_ig2 = calculate_combined_impact(intervention1=intervention_list_bit103_IG2[[1]],
                                                   intervention2=pyrethroid_interv_bit103_IG2[[1]],
                                                   cov_intervention1 = coverage_vector[["bit103_IG2"]],
                                                   cov_intervention2 = coverage_vector[["bit103_IG2"]],
                                                   cov_combination = coverage_vector[["bit103_IG2"]],
                                                   N_vec =10000)
    
    impacts_bit059 = calculate_combined_impact(intervention1=intervention_list_bit059[[1]],
                                               intervention2=pyrethroid_interv_bit059[[1]],
                                               cov_intervention1 = coverage_vector[["bit059"]],
                                               cov_intervention2 = coverage_vector[["bit059"]],
                                               cov_combination = coverage_vector[["bit059"]],
                                               N_vec =10000)
    
    impacts_bit080 = calculate_combined_impact(intervention1=intervention_list_bit080[[1]],
                                               intervention2=pyrethroid_interv_bit080[[1]],
                                               cov_intervention1 = coverage_vector[["bit080"]],
                                               cov_intervention2 = coverage_vector[["bit055"]],
                                               cov_combination = coverage_vector[["bit080"]],
                                               N_vec =10000)
    
    impacts=data.frame(EHT=c("Kibondo", "BIT055", "BIT103", "BIT103", "BIT059", "BIT080"),
                       insecticide=c("IG2","OlysetPlus", "OlysetPlus", "IG2" , "OlysetPlus", "IG2"),
                       VCred=c(impacts_kibondo,impacts_bit055, impacts_bit103_pbo, impacts_bit103_ig2,impacts_bit059,impacts_bit080))
  }

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


# bednet efficacy
impacts_efficacy=calculate_impact_4EHT(washedDecay=F, vectormodel = model_noRhythms)

# bednet efficacy + attrition
impacts_EfficacyAttrition=calculate_impact_4EHT(L_vector=L_vector_attrition, 
                                                kappa_vector=kappa_vector_attrition,
                                                washedDecay=F, vectormodel = model_noRhythms)

impacts_EfficacyAttritionUsage=calculate_impact_4EHT(L_vector=L_vector_attrition, 
                                                     kappa_vector=kappa_vector_attrition,
                                                     washedDecay=F, vectormodel = model_noRhythms,
                                                     coverage_vector = coverage_vector_usage)


# bednet efficacy + attrition + insecticide decay
impacts_EfficacyAttritionDecay=calculate_impact_4EHT(L_vector=L_vector_attrition, 
                                                     kappa_vector=kappa_vector_attrition,
                                                     washedDecay=T, vectormodel = model_noRhythms,
                                                     coverage_vector = coverage_vector_usage)

# bednet efficacy + attrition + insecticide decay + rhythms
impacts_EfficacyAttritionDecayRhythms=calculate_impact_4EHT(L_vector=L_vector_attrition, 
                                                            kappa_vector=kappa_vector_attrition,
                                                            washedDecay=T, vectormodel = model_p,
                                                            coverage_vector = coverage_vector_usage)

VCreduc_all=rbind(impacts_efficacy %>% mutate(scenario="Efficacy") , 
                  impacts_EfficacyAttrition %>% mutate(scenario="EfficacyAttrition"), 
                  impacts_EfficacyAttritionUsage %>% mutate(scenario="EfficacyAttritionUsage") , 
                  impacts_EfficacyAttritionDecay %>% mutate(scenario="EfficacyAttritionUsageDecay") ,
                  impacts_EfficacyAttritionDecayRhythms  %>% mutate(scenario="EfficacyAttritionUsageDecayRhythms")
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

plot_grid(cascade_ig2_bit103,cascade_ig2_kibondo,
          cascade_pbo_bit103,cascade_pbo_bit055,
          ncol = 2, byrow = F, 
          labels=c("BIT103", "BIT103","Kibondo",  "BIT055"),
          label_size = 12,
          hjust = 0, label_x = 0.85)
ggsave(file.path(plotDir, "cascades_LLIN_noUncertainty_endpoint.png"), width=12, height=8)

IG2=plot_grid(ggdraw() +draw_label(  "Interceptor G2", x=0.1, y=0.9, fontface = "bold" )+theme_classic(),
              cascade_ig2_bit103,cascade_ig2_kibondo,cascade_ig2_bit080,
          ncol = 1, byrow = F, 
          labels=c("","BIT103", "Kibondo","BIT080"),
          label_size = 12, rel_heights = c(0.15, 1,1,1),
          hjust = 0, label_x = 0.85)#+draw_plot_label("Interceptor G2", vjust = 0.5, hjust = -0.3)

OlysetPlus=plot_grid(ggdraw() +draw_label(  "Olyset Plus", x=0.05, y=0.9, fontface = "bold" )+theme_classic(),
                     cascade_pbo_bit103,cascade_pbo_bit055,cascade_pbo_bit059,
              ncol = 1, byrow = F, 
              labels=c( "","BIT103","BIT055",  "Odufuwa"),
              label_size = 12, rel_heights = c(0.15, 1,1,1),
              hjust = 0, label_x = 0.85)#+draw_plot_label("Olyset Plus", vjust = 0.5, hjust = -0.3)

plot_grid(IG2, OlysetPlus,
          ncol = 2, byrow = F, 
          #labels=c("Interceptor G2", "Olyset Plus"),
          label_size = 12,
          hjust = 0, label_x = 0.4)
ggsave(file.path(plotDir, "cascades_LLIN_noUncertainty6_endpoint.png"), width=12, height=12)


#####################################################################
# CASCADE WITH UNCERTAINTY
####################################################################

results_kibondo_w=results_kibondo_w_72
results_kibondo_unw=results_kibondo_unw_72
results_bit103_unw=results_bit103_unw_72
results_bit103_w=results_bit103_w_72
results_bit080_unw=results_bit080_unw_72
results_bit080_w=results_bit080_w_72



# bednet efficacy
impacts_efficacy_uncertainty=calculate_impact_4EHT(washedDecay=F, vectormodel = model_noRhythms,uncertainty = T, nsamples = 1000)

# bednet efficacy + attrition
impacts_EfficacyAttrition_uncertainty=calculate_impact_4EHT(L_vector=L_vector_attrition, 
                                                            kappa_vector=kappa_vector_attrition,
                                                            washedDecay=F, vectormodel = model_noRhythms,uncertainty = T, nsamples = 1000)

# bednet efficacy + attrition + usage
impacts_EfficacyAttritionUsage_uncertainty=calculate_impact_4EHT(L_vector=L_vector_attrition, 
                                                                 kappa_vector=kappa_vector_attrition,
                                                                 washedDecay=F, vectormodel = model_noRhythms,uncertainty = T,
                                                                 coverage_vector = coverage_vector_usage, nsamples = 1000)
# bednet efficacy + attrition + insecticide decay
impacts_EfficacyAttritionDecay_uncertainty=calculate_impact_4EHT(L_vector=L_vector_attrition, 
                                                                 kappa_vector=kappa_vector_attrition,
                                                                 washedDecay=T, vectormodel = model_noRhythms,uncertainty = T,
                                                                 coverage_vector = coverage_vector_usage, nsamples = 1000)

# bednet efficacy + attrition + insecticide decay + rhythms
impacts_EfficacyAttritionDecayRhythms_uncertainty=calculate_impact_4EHT(L_vector=L_vector_attrition, 
                                                                        kappa_vector=kappa_vector_attrition,
                                                                        washedDecay=T, vectormodel = model_p,uncertainty = T,
                                                                        coverage_vector = coverage_vector_usage, nsamples = 1000)


VCreduc_all_uncertainty=rbind(impacts_efficacy_uncertainty %>% mutate(scenario="Efficacy") , 
                              impacts_EfficacyAttrition_uncertainty %>% mutate(scenario="EfficacyAttrition"), 
                              impacts_EfficacyAttritionUsage_uncertainty %>% mutate(scenario="EfficacyAttritionUsage") , 
                              impacts_EfficacyAttritionDecay_uncertainty %>% mutate(scenario="EfficacyAttritionUsageDecay") ,
                              impacts_EfficacyAttritionDecayRhythms_uncertainty  %>% mutate(scenario="EfficacyAttritionUsageDecayRhythms")
)
write.csv(VCreduc_all_uncertainty, file = file.path(plotDir, "VCreduc_all_uncertainty.csv"), row.names = F)

df_VCreduc_all_uncertainty=VCreduc_all_uncertainty %>% group_by(insecticide,  scenario)%>%
  summarise(q025=quantile(VCred, probs=0.025),q975=quantile(VCred, probs=0.975),
            min=min(VCred),max=max(VCred),
            med=quantile(VCred, probs=0.5), mean=mean(VCred) ) %>% ungroup()


cascade_ig2_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_all_uncertainty %>% filter( insecticide=="IG2") %>% 
                                                                   select(-insecticide, -min, -max, -med)%>% 
                                                                   rename(value=mean, variable=scenario) %>% mutate(value=100*value), 
                                                                 colorfinal="orange")
cascade_pbo_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_all_uncertainty %>% filter( insecticide=="PBO") %>% 
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
cascade_pbo_bit103_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_EHT_uncertainty %>% filter(EHT=="BIT103", insecticide=="PBO") %>% 
                                                                          select(-insecticide, -min, -max, -EHT, -med)%>% 
                                                                          rename(value=mean, variable=scenario) %>% mutate(value=100*value), 
                                                                        colorfinal="dodgerblue")
cascade_pbo_bit055_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_EHT_uncertainty %>% filter(EHT=="BIT055", insecticide=="PBO") %>% 
                                                                          select(-insecticide, -min, -max, -EHT, -med)%>% 
                                                                          rename(value=mean, variable=scenario) %>% mutate(value=100*value), 
                                                                        colorfinal="dodgerblue")
cascade_ig2_kibondo_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_EHT_uncertainty %>% filter(EHT=="Kibondo", insecticide=="IG2") %>% 
                                                                           select(-insecticide, -min, -max, -EHT, -med)%>% 
                                                                           rename(value=mean, variable=scenario)%>% mutate(value=100*value), 
                                                                         colorfinal="orange")

cascade_pbo_bit059_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_EHT_uncertainty %>% filter(EHT=="BIT059", insecticide=="PBO") %>% 
                                                                          select(-insecticide, -min, -max, -EHT, -med)%>% 
                                                                          rename(value=mean, variable=scenario) %>% mutate(value=100*value), 
                                                                        colorfinal="dodgerblue")
cascade_ig2_bit080_uncertainty=plot_bars_effectiveness_LLIN_uncertainty(df_VCred=df_VCreduc_EHT_uncertainty %>% filter(EHT=="BIT080", insecticide=="IG2") %>% 
                                                                           select(-insecticide, -min, -max, -EHT, -med)%>% 
                                                                           rename(value=mean, variable=scenario)%>% mutate(value=100*value), 
                                                                         colorfinal="orange")

plot_grid(cascade_ig2_bit103_uncertainty,cascade_ig2_kibondo_uncertainty, cascade_ig2_bit080_uncertainty,
          cascade_pbo_bit103_uncertainty,cascade_pbo_bit055_uncertainty,cascade_pbo_bit059_uncertainty,
          ncol = 2, byrow = F, 
          labels=c("Interceptor G2 (BIT103)", "Olyset Plus (BIT103)",
                   "Interceptor G2 (Kibondo)",  "Olyset Plus (BIT055)",
                   "Interceptor G2 (BIT080)",  "Olyset Plus (Odufuwa)"),
          label_size = 12,
          hjust = 0, label_x = 0.65)
ggsave(file.path(plotDir, "cascades_LLIN_perEHT_uncertainty.png"), width=12, height=12)


plot_grid(cascade_ig2_bit103_uncertainty,
          cascade_pbo_bit103_uncertainty,
          ncol = 1, byrow = F, 
          labels=c("Interceptor G2 (BIT103)", "Olyset Plus (BIT103)"),
          label_size = 12,
          hjust = 0, label_x = 0.65)
ggsave(file.path(plotDir, "cascades_LLIN_EHT103_uncertainty.png"), width=6, height=8)


