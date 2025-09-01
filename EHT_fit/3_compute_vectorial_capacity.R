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


results_bit055_w<- readRDS(file.path(stanDir,"/stan_BIT055_washed_24.rds"))
results_bit055_unw<- readRDS(file.path(stanDir,"/stan_BIT055_unwashed_24.rds"))

results_bit103_unw_72<- readRDS(file.path(stanDir,"/stan_assenga_unwashed_72_all.rds"))
results_bit103_w_72<- readRDS(file.path(stanDir,"/stan_assenga_washed_72_all.rds"))

results_bit103_unw_72_CI<- readRDS(file.path(stanDir,"/stan_assengaCotedIvoire_unwashed_72.rds"))
results_bit103_w_72_CI<- readRDS(file.path(stanDir,"/stan_assengaCotedIvoire_washed_72.rds"))

results_bit059_w<- readRDS(file.path(stanDir,"/stan_odufuwa_washed_24.rds"))
results_bit059_unw<- readRDS(file.path(stanDir,"/stan_odufuwa_unwashed_24.rds"))

results_bit080_w_72<- readRDS(file.path(stanDir,"/stan_BIT080_washed_72.rds"))
results_bit080_unw_72<- readRDS(file.path(stanDir,"/stan_BIT080_unwashed_72.rds"))


results_martin_w72<- readRDS(file.path(stanDir,"/stan_martin_36m_72_gambiae.rds"))
results_martin_unw72<- readRDS(file.path(stanDir,"/stan_martin_unwashed_72_gambiae.rds"))
results_martin_w72_funestus<- readRDS(file.path(stanDir,"/stan_martin_36m_72_funestus.rds"))
results_martin_unw72_funestus<- readRDS(file.path(stanDir,"/stan_martin_unwashed_72_funestus.rds"))

results_nguessan_unw<- readRDS(file.path(stanDir,"/stan_nguessan_unwashed.rds"))
results_nguessan_w<- readRDS(file.path(stanDir,"/stan_nguessan_washed.rds"))

results_sovegnonIG2_unw<- readRDS(file.path(stanDir,"/stan_sovegnonIG2_unwashed.rds"))
results_sovegnonIG2_w<- readRDS(file.path(stanDir,"/stan_sovegnonIG2_aged.rds"))
results_sovegnonIG1_unw<- readRDS(file.path(stanDir,"/stan_sovegnonIG1_unwashed.rds"))
results_sovegnonIG1_w<- readRDS(file.path(stanDir,"/stan_sovegnonIG1_aged.rds"))

calculate_VC_uncertainty=function(results_EHT, nsamples, names, insecticide_id){
  
  intervention_list <- interventions_vectorial_capacity_wrapper(results=results_EHT,results_washed20 =NULL,decay = "weibull",
                                                                model_p=model_noRhythms,
                                                                insecticides=insecticide_id,L = 3*365, kappa=2,
                                                                names=names,  npoints=3*365,washedDecay=FALSE, uncertainty=TRUE, 
                                                                cov=1, nsamples=nsamples, inbed_exposure = 1)
  
  
  
  
  
}

nsamples=1000

if(rerun){
  # IG2 unwashed
  VCred_kibondo=calculate_VC_uncertainty(results_EHT=results_kibondo_unw_72, nsamples=nsamples, names="IG2", insecticide_id=2)
  VCred_bit103_ig2=calculate_VC_uncertainty(results_EHT=results_bit103_unw_72, nsamples=nsamples, names="IG2", insecticide_id=3)
  VCred_bit103_ig2_CI=calculate_VC_uncertainty(results_EHT=results_bit103_unw_72_CI, nsamples=nsamples, names="IG2", insecticide_id=3)
  VCred_bit080=calculate_VC_uncertainty(results_EHT=results_bit080_unw_72, nsamples=nsamples, names="IG2", insecticide_id=2)
  VCred_martin_ig2=calculate_VC_uncertainty(results_EHT=results_martin_unw72, nsamples=nsamples, names="IG2", insecticide_id=3)
  VCred_martin_ig2_funestus=calculate_VC_uncertainty(results_EHT=results_martin_unw72_funestus, nsamples=nsamples, names="IG2", insecticide_id=3)
  VCred_sovegnon=calculate_VC_uncertainty(results_EHT=results_sovegnonIG2_unw, nsamples=nsamples, names="IG2", insecticide_id=1)
  VCred_nguessan=calculate_VC_uncertainty(results_EHT=results_nguessan_unw, nsamples=nsamples, names="IG2", insecticide_id=2)

  impacts_ig2=rbind(data.frame(VCred=VCred_kibondo) %>% mutate(EHT_short="Kibondo"),
                data.frame(VCred=VCred_bit103_ig2) %>% mutate(EHT_short="Assenga"),
                data.frame(VCred=VCred_bit103_ig2_CI) %>% mutate(EHT_short="Assenga, CI"),
                data.frame(VCred=VCred_bit080) %>% mutate(EHT_short="BIT080"),
                data.frame(VCred=VCred_martin_ig2) %>% mutate(EHT_short="Martin"),
                data.frame(VCred=VCred_martin_ig2_funestus) %>% mutate(EHT_short="Martin, f"),
                data.frame(VCred=VCred_sovegnon) %>% mutate(EHT_short="Sovegnon"),
                data.frame(VCred=VCred_nguessan) %>% mutate(EHT_short="Nguessan")
  )%>% mutate(washed_status="Unwashed", insecticide_name="Interceptor G2")
  
  
  # PBO unwashed
  VCred_bit055=calculate_VC_uncertainty(results_EHT=results_bit055_unw, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit103_pbo=calculate_VC_uncertainty(results_EHT=results_bit103_unw_72, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit103_pbo_CI=calculate_VC_uncertainty(results_EHT=results_bit103_unw_72_CI, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit059=calculate_VC_uncertainty(results_EHT=results_bit059_unw, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_martin_pbo=calculate_VC_uncertainty(results_EHT=results_martin_unw72, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_martin_pbo_funestus=calculate_VC_uncertainty(results_EHT=results_martin_unw72_funestus, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  

  impacts_pbo=rbind(data.frame(VCred=VCred_bit055) %>% mutate(EHT_short="BIT055"),
                data.frame(VCred=VCred_bit103_pbo) %>% mutate(EHT_short="Assenga"),
                data.frame(VCred=VCred_bit103_pbo_CI) %>% mutate(EHT_short="Assenga, CI"),
                data.frame(VCred=VCred_bit059) %>% mutate(EHT_short="Odufuwa"),
                data.frame(VCred=VCred_martin_pbo) %>% mutate(EHT_short="Martin"),
                data.frame(VCred=VCred_martin_pbo_funestus) %>% mutate(EHT_short="Martin, f")
  )%>% mutate(washed_status="Unwashed", insecticide_name="Olyset Plus")
  
  
  # IG2 washed
  VCred_kibondo_w=calculate_VC_uncertainty(results_EHT=results_kibondo_w_72, nsamples=nsamples, names="IG2", insecticide_id=2)
  VCred_bit103_ig2_w=calculate_VC_uncertainty(results_EHT=results_bit103_w_72, nsamples=nsamples, names="IG2", insecticide_id=3)
  VCred_bit103_ig2_CI_w=calculate_VC_uncertainty(results_EHT=results_bit103_w_72_CI, nsamples=nsamples, names="IG2", insecticide_id=3)
  VCred_bit080_w=calculate_VC_uncertainty(results_EHT=results_bit080_w_72, nsamples=nsamples, names="IG2", insecticide_id=2)
  VCred_martin_ig2_w=calculate_VC_uncertainty(results_EHT=results_martin_w72, nsamples=nsamples, names="IG2", insecticide_id=3)
  VCred_martin_ig2_w_funestus=calculate_VC_uncertainty(results_EHT=results_martin_w72_funestus, nsamples=nsamples, names="IG2", insecticide_id=2)
  VCred_sovegnon_w=calculate_VC_uncertainty(results_EHT=results_sovegnonIG2_w, nsamples=nsamples, names="IG2", insecticide_id=1)
  VCred_nguessan_w=calculate_VC_uncertainty(results_EHT=results_nguessan_w, nsamples=nsamples, names="IG2", insecticide_id=2)

  
  impacts_ig2_w=rbind(data.frame(VCred=VCred_kibondo_w) %>% mutate(EHT_short="Kibondo"),
                  data.frame(VCred=VCred_bit103_ig2_w) %>% mutate(EHT_short="Assenga"),
                  data.frame(VCred=VCred_bit103_ig2_CI_w) %>% mutate(EHT_short="Assenga, CI"),
                  data.frame(VCred=VCred_bit080_w) %>% mutate(EHT_short="BIT080"),
                  data.frame(VCred=VCred_martin_ig2_w) %>% mutate(EHT_short="Martin"),
                  data.frame(VCred=VCred_martin_ig2_w_funestus) %>% mutate(EHT_short="Martin, f"),
                  data.frame(VCred=VCred_sovegnon_w) %>% mutate(EHT_short="Sovegnon"),
                  data.frame(VCred=VCred_nguessan_w) %>% mutate(EHT_short="Nguessan")
  )%>% mutate(washed_status="Washed", insecticide_name="Interceptor G2")
  
  # PBO washed
  VCred_bit103_pbo_w=calculate_VC_uncertainty(results_EHT=results_bit103_w_72, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit103_pbo_CI_w=calculate_VC_uncertainty(results_EHT=results_bit103_w_72_CI, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit059_w=calculate_VC_uncertainty(results_EHT=results_bit059_w, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit055_w=calculate_VC_uncertainty(results_EHT=results_bit055_w, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_martin_pbo_w=calculate_VC_uncertainty(results_EHT=results_martin_w72, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_martin_pbo_w_funestus=calculate_VC_uncertainty(results_EHT=results_martin_w72_funestus, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  
  impacts_pbo_w=rbind(data.frame(VCred=VCred_bit103_pbo_w) %>% mutate(EHT_short="Assenga"),
                      data.frame(VCred=VCred_bit103_pbo_CI_w) %>% mutate(EHT_short="Assenga, CI"),
                      data.frame(VCred=VCred_bit059_w) %>% mutate(EHT_short="Odufuwa"),
                      data.frame(VCred=VCred_bit055_w) %>% mutate(EHT_short="BIT055"),
                      data.frame(VCred=VCred_martin_pbo_w) %>% mutate(EHT_short="Martin"),
                      data.frame(VCred=VCred_martin_pbo_w_funestus) %>% mutate(EHT_short="Martin, f")
  )%>% mutate(washed_status="Washed", insecticide_name="Olyset Plus")
  
  # Pyrethroid unwashed
  VCred_kibondo_pyr=calculate_VC_uncertainty(results_EHT=results_kibondo_unw_72, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit080_pyr=calculate_VC_uncertainty(results_EHT=results_bit080_unw_72, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit103_pyr=calculate_VC_uncertainty(results_EHT=results_bit103_unw_72, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit103_pyr_CI=calculate_VC_uncertainty(results_EHT=results_bit103_unw_72_CI, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit059_pyr=calculate_VC_uncertainty(results_EHT=results_bit059_unw, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit055_pyr=calculate_VC_uncertainty(results_EHT=results_bit055_unw, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_martin_pyr=calculate_VC_uncertainty(results_EHT=results_martin_unw72, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_martin_pyr_funestus=calculate_VC_uncertainty(results_EHT=results_martin_unw72_funestus, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_sovegnon_pyr=calculate_VC_uncertainty(results_EHT=results_sovegnonIG1_unw, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_nguessan_pyr=calculate_VC_uncertainty(results_EHT=results_nguessan_unw, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  
  impacts_pyr=rbind(data.frame(VCred=VCred_kibondo_pyr) %>% mutate(EHT_short="Kibondo"),
                    data.frame(VCred=VCred_bit055_pyr) %>% mutate(EHT_short="BIT055"),
                    data.frame(VCred=VCred_bit103_pyr) %>% mutate(EHT_short="Assenga"),
                    data.frame(VCred=VCred_bit103_pyr_CI) %>% mutate(EHT_short="Assenga, CI"),
                    data.frame(VCred=VCred_bit059_pyr) %>% mutate(EHT_short="Odufuwa"),
                    data.frame(VCred=VCred_bit080_pyr) %>% mutate(EHT_short="BIT080"),
                    data.frame(VCred=VCred_martin_pyr) %>% mutate(EHT_short="Martin"),
                    data.frame(VCred=VCred_martin_pyr_funestus) %>% mutate(EHT_short="Martin, f"),
                    data.frame(VCred=VCred_sovegnon_pyr) %>% mutate(EHT_short="Sovegnon"),
                    data.frame(VCred=VCred_nguessan_pyr) %>% mutate(EHT_short="Nguessan")
  )%>% mutate(washed_status="Unwashed", insecticide_name="Pyrethroid")
  
  
  
  # Pyrethroid washed
  VCred_kibondo_pyr_w=calculate_VC_uncertainty(results_EHT=results_kibondo_w_72, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit080_pyr_w=calculate_VC_uncertainty(results_EHT=results_bit080_w_72, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit103_pyr_w=calculate_VC_uncertainty(results_EHT=results_bit103_w_72, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit103_pyr_CI_w=calculate_VC_uncertainty(results_EHT=results_bit103_w_72_CI, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit059_pyr_w=calculate_VC_uncertainty(results_EHT=results_bit059_w, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit055_pyr_w=calculate_VC_uncertainty(results_EHT=results_bit055_w, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_martin_pyr_w=calculate_VC_uncertainty(results_EHT=results_martin_w72, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_martin_pyr_w_funestus=calculate_VC_uncertainty(results_EHT=results_martin_w72_funestus, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_sovegnon_pyr_w=calculate_VC_uncertainty(results_EHT=results_sovegnonIG1_w, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_nguessan_pyr_w=calculate_VC_uncertainty(results_EHT=results_nguessan_w, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  

  
  impacts_pyr_w=rbind(data.frame(VCred=VCred_kibondo_pyr_w) %>% mutate(EHT_short="Kibondo"),
                      data.frame(VCred=VCred_bit055_pyr_w) %>% mutate(EHT_short="BIT055"),
                      data.frame(VCred=VCred_bit103_pyr_w) %>% mutate(EHT_short="Assenga"),
                      data.frame(VCred=VCred_bit103_pyr_CI_w) %>% mutate(EHT_short="Assenga, CI"),
                      data.frame(VCred=VCred_bit059_pyr_w) %>% mutate(EHT_short="Odufuwa"),
                      data.frame(VCred=VCred_bit080_pyr_w) %>% mutate(EHT_short="BIT080"),
                      data.frame(VCred=VCred_martin_pyr_w) %>% mutate(EHT_short="Martin"),
                      data.frame(VCred=VCred_martin_pyr_w_funestus) %>% mutate(EHT_short="Martin, f"),
                      data.frame(VCred=VCred_sovegnon_pyr_w) %>% mutate(EHT_short="Sovegnon"),
                      data.frame(VCred=VCred_nguessan_pyr_w) %>% mutate(EHT_short="Nguessan")
  )%>% mutate(washed_status="Washed", insecticide_name="Pyrethroid")
  
  
  write.csv(rbind(impacts_ig2, impacts_pbo,impacts_ig2_w, impacts_pbo_w, impacts_pyr,impacts_pyr_w), file = file.path(outputsDir, "impacts_EHT_control.csv"))
  
  impacts_ig2 %>% group_by(EHT_short, insecticide_name,washed_status)%>% summarise(mean=mean(VCred), min=min(VCred), max=max(VCred))
  impacts_pbo %>% group_by(EHT_short, insecticide_name,washed_status)%>% summarise(mean=mean(VCred), min=min(VCred), max=max(VCred))
  impacts_pyr %>% group_by(EHT_short, insecticide_name,washed_status)%>% summarise(mean=mean(VCred), min=min(VCred), max=max(VCred))
  
}

#######################
# PLOT SUMMARY
######################
estimates_kibondo_72=read.csv(file.path(outputsDir, "estimates_kibondo_72.csv"))
estimates_kibondo_24=read.csv(file.path(outputsDir, "estimates_kibondo_24.csv"))

estimates_bit059=read.csv(file.path(outputsDir, "estimates_odufuwa_24.csv"))

estimates_bit055=read.csv(file.path(outputsDir, "estimates_bit055_24.csv"))

estimates_bit080_72=read.csv(file.path(outputsDir, "estimates_bit080_72.csv"))
estimates_bit080_24=read.csv(file.path(outputsDir, "estimates_bit080_24.csv"))

estimates_bit103_72=read.csv(file.path(outputsDir, "estimates_assenga_72_all.csv"))

estimates_bit103_72_CI=read.csv(file.path(outputsDir, "estimates_assengaCotedIvoire_72.csv"))


estimates_martin_72=read.csv(file.path(outputsDir, "estimates_martin_72_gambiae.csv"))
estimates_martin_72_funestus=read.csv(file.path(outputsDir, "estimates_martin_72_funestus.csv"))

estimates_sovegnon_72=read.csv(file.path(outputsDir, "estimates_sovegnon_72.csv"))
estimates_nguessan_72=read.csv(file.path(outputsDir, "estimates_nguessan_72.csv"))


impacts_72=read.csv(file.path(outputsDir, "impacts_EHT_control.csv"))

impacts_72_summary=impacts_72 %>% group_by(EHT_short, insecticide_name,washed_status)%>% 
  summarise(mean=mean(VCred), X2.5.=quantile(VCred, probs = 0.025), X97.5.=quantile(VCred, probs = 0.975))%>%
  mutate(param ="VCred", EHT=ifelse(EHT_short=="Martin","Martin et al. 2024",
                                           ifelse(EHT_short=="Nguessan","Nguessan et al. 2016", 
                                                  ifelse(EHT_short=="Sovegnon","Sovegnon et al. 2024",
                                                         ifelse(EHT_short=="Kibondo","Kibondo et al. 2022",
                                                                ifelse(EHT_short=="Odufuwa","Odufuwa et al. 2024", 
                                                                       ifelse(EHT_short=="Assenga","Assenga et al. 2025", 
                                                                              ifelse(EHT_short=="Martin, f","Martin et al. 2024, f", 
                                                                                     ifelse(EHT_short=="Assenga, CI","Assenga et al. 2025, CI", 
                                                                                            EHT_short ) )
                                                 )))))), X=NA)


#############################
# endpoint: longest available per trial

# formating the results
all_estimates=rbind(estimates_kibondo_72 %>%mutate(EHT="Kibondo et al. 2022", EHT_short="Kibondo"),
                    estimates_bit055 %>%mutate(EHT="BIT055", EHT_short=EHT),
                    estimates_bit059 %>%mutate(EHT="Odufuwa et al. 2024", EHT_short="Odufuwa"),
                    estimates_bit080_72 %>%mutate(EHT="BIT080", EHT_short=EHT),
                    estimates_bit103_72 %>%mutate(EHT="Assenga et al. 2025", EHT_short="Assenga"),
                    estimates_bit103_72_CI %>%mutate(EHT="Assenga et al. 2025, CI", EHT_short="Assenga, CI"),
                    estimates_martin_72 %>%mutate(EHT="Martin et al. 2024", EHT_short="Martin"),
                    estimates_martin_72_funestus %>%mutate(EHT="Martin et al. 2024, f", EHT_short="Martin, f"),
                    estimates_sovegnon_72 %>%mutate(EHT="Sovegnon et al. 2024", EHT_short="Sovegnon"),
                    estimates_nguessan_72 %>%mutate(EHT="Nguessan et al. 2016", EHT_short="Nguessan"),
                    impacts_72_summary 
)%>%
  mutate(insecticide_name=ifelse(insecticide_name %in% c("IG2", "Interceptor_G2","InterceptorG2","Interceptor®G2", "Interceptor\xaeG2", "Interceptor G2", "interceptorG2", "IG2.Aged"), "Interceptor G2",
                                 ifelse(insecticide_name %in% c("Olyset_Plus", "OlysetPlus", "Olyset Plus", "olysetplus"), "Olyset Plus", "Pyrethroid-only")),
         param1=gsub("Initial", "", gsub("Efficacy", "", gsub("Rate" ,"",  param))))%>%
  #mutate(EHT=gsub(" et al.", "\net al.", EHT), EHT_full=EHT)%>% 
  filter(param1 !="KillingDuringHostSeeking")%>%
  mutate(EHT=ifelse((EHT_short=="Sovegnon" & insecticide_name=="Pyrethroid-only"), "Sovegnon et al. 2024, 1",ifelse(EHT_short=="Sovegnon" & insecticide_name=="Interceptor G2", "Sovegnon et al. 2024, 2", EHT)))

all_estimates$param2=factor(all_estimates$param1, levels=c("Repellency", "Preprandialkilling", "Postprandialkilling", "VCred"),
                            labels=c(expression(paste("Reduction in host availability (", pi,")")), expression(paste("Pre-prandial killing effect (", phi,")")),expression(paste("Post-prandial killing effect (", xi,")")), expression(paste("Entomological efficacy"))))
all_estimates$param3=factor(all_estimates$param1, levels=c("Repellency", "Preprandialkilling", "Postprandialkilling", "VCred"),
                            labels=c("Reduction in host availability", "Pre-prandial killing effect","Post-prandial killing effect","Entomological efficacy"))
all_estimates$EHT=factor(all_estimates$EHT, levels=c("Assenga et al. 2025", "Kibondo et al. 2022", "BIT080", "BIT055", "Odufuwa et al. 2024", "Martin et al. 2024", "Martin et al. 2024, f", "Assenga et al. 2025, CI", "Nguessan et al. 2016","Sovegnon et al. 2024, 1","Sovegnon et al. 2024, 2"))
all_estimates$EHT_id= as.numeric(all_estimates$EHT)
all_estimates$country=ifelse(all_estimates$EHT_short %in% c("Sovegnon", "Nguessan", "Assenga, CI"), "Benin", "Tanzania")

all_estimates_simple=all_estimates %>%
  mutate(washed_status=ifelse(washed_status=="Washed20", "Washed", washed_status),
         label_height=ifelse(param2=="Increase in host seeking mortality",-0.07, -0.03))

all_estimates_simple_mean=all_estimates_simple%>% group_by(insecticide_name, param2, washed_status)%>% summarise(meanmean=mean(mean))%>%
  mutate( washed_status=ifelse(washed_status=="Unwashed", "Unwashed", "Washed 20x/Aged"))


# plot all estimates
all_estimates_simple %>%
  mutate( washed_status=ifelse(washed_status=="Unwashed", "Unwashed", "Washed 20x/Aged"))%>%
  ggplot()+
  geom_col(aes(x=insecticide_name, y=mean, fill=insecticide_name, group=interaction( washed_status, EHT), alpha=washed_status), stat="identity", position=position_dodge())+
  geom_errorbar(aes(x=insecticide_name, ymin=X2.5., ymax=X97.5., color=insecticide_name,  group=interaction( washed_status, EHT)),position=position_dodge(0.9), width=0.4)+
  geom_segment(data=all_estimates_simple_mean%>% filter(insecticide_name =="Interceptor G2"),
               aes(y = meanmean, yend=meanmean, x= 0.55, xend = 1.45, color=insecticide_name, linetype = washed_status), linewidth = 0.9) +
  geom_segment(data=all_estimates_simple_mean%>% filter(insecticide_name =="Olyset Plus"),
               aes(y = meanmean, yend=meanmean, x= 1.55, xend = 2.45, color=insecticide_name, linetype = washed_status), linewidth = 0.9) +
  geom_segment(data=all_estimates_simple_mean%>% filter(insecticide_name =="Pyrethroid-only"),
               aes(y = meanmean, yend=meanmean, x= 2.55, xend = 3.45, color=insecticide_name, linetype = washed_status), linewidth = 0.9) +
  facet_wrap( . ~param2, scales = "free_x", labeller = label_parsed)+
  theme_bw()+#ylim(0,1)+
  labs(x="", y="", color="", fill="", alpha="", linetype="")+
  scale_alpha_manual(values=c(  0.7, 0.4),na.translate=FALSE)+
  scale_color_manual(values=c("darkorange","dodgerblue","darkgrey" ))+
  scale_fill_manual(values=c("darkorange","dodgerblue","darkgrey"))+
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))+
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 14),
        axis.text=element_text(size=14), strip.background =element_rect(fill="white"), legend.text =element_text(size=14) )+ 
  geom_text(aes(label=EHT_id, x=insecticide_name, y=label_height,group=EHT),
            position = position_dodge(width=0.9),size=3.5)+#ylim(-0.05,1)+ 
  guides(fill = "none",color = "none")
ggsave(file.path(plotDir, "plot_EHTfit_pyr_all_new.png"), width=12, height=10)



final_table=all_estimates_simple %>%
  mutate(mean=round(mean, digits = 2),X2.5.=round(X2.5., digits = 2),X97.5.=round(X97.5., digits = 2))%>%
  select(insecticide_name, EHT, EHT_id, param3, washed_status, mean, X2.5., X97.5.)%>%
  rename(netType=insecticide_name, parameter=param3, q025=X2.5., q975=X97.5.)%>%
  tidyr::pivot_wider(id_cols =c( netType, EHT, parameter), names_from = washed_status, values_from = c(mean, q025, q975) )%>%
  mutate(halflife_insecticide=3/(1-mean_Washed/mean_Unwashed)
  )

write.csv(final_table, file.path(scriptDir, "fitted_parameters_pyr.csv"), row.names = F)

final_table %>%
  select(-halflife_insecticide)%>%
  pivot_longer(cols = c("mean_Washed" , "mean_Unwashed","q025_Washed"   ,       "q025_Unwashed"       
                        ,"q975_Washed"    ,      "q975_Unwashed"))%>%
  separate(name, into=c("param", "wash_status"))%>%
  ungroup()%>%
  pivot_wider(id_cols = c(wash_status, parameter, netType,EHT ), names_from = param, values_from = value)%>%
  mutate(value=paste0(mean, " (", q025, "-", q975, ")"))%>%
  select(-mean, -q025, -q975)%>%
  pivot_wider(id_cols = c(parameter, netType,EHT ), names_from = wash_status, values_from = value)%>%
  arrange(parameter, netType, EHT)%>%
  select(parameter, netType, EHT, Unwashed, Washed)%>%
  rename("Net type"=netType, "Washed 20x"=Washed)%>%
  gt()%>%
  gt::gtsave(filename = file.path(plotDir, "outputs_ento_efficacy_all.docx"))



final_table %>% filter(parameter=="Entomological efficacy")%>%
  mutate(reduction_wash=100*(mean_Unwashed-mean_Washed)/mean_Unwashed)%>%
  View()

final_table %>% filter(parameter=="Entomological efficacy", netType=="Pyrethroid-only")%>%
  mutate(reduction_wash=100*(mean_Unwashed-mean_Washed)/mean_Unwashed)%>%
  View()


impacts_72_summary%>% 
  pivot_wider(id_cols = c(EHT_short, insecticide_name), names_from = washed_status, values_from = mean)%>%
  mutate(red=(Unwashed-Washed)/Unwashed)%>%View()


###################################################
# COMPUTE DIFFERENCE BETWEEN PRODUCTS WITHIN EACH EHT

impacts_pivot=impacts_72 %>%
  filter(washed_status=="Unwashed")%>%
  mutate(insecticide_name=gsub(" ", "", insecticide_name), 
         Y=X%%1000, 
         VCred=VCred*100)%>%
  pivot_wider(id_cols = c(EHT_short,Y, washed_status), names_from = insecticide_name, values_from = VCred)%>%
  mutate(diff_IG2_pyr=InterceptorG2-Pyrethroid, 
         diff_OP_pyr=OlysetPlus-Pyrethroid,
         diff_IG2_OP=InterceptorG2-OlysetPlus)


impacts_pivot%>%
  group_by(EHT_short)%>%
  summarise(
    X2.5._IG2_pyr=quantile(diff_IG2_pyr, probs = 0.025, na.rm=T), 
    X97.5._IG2_pyr=quantile(diff_IG2_pyr, probs = 0.975, na.rm=T),
    X2.5._OP_pyr=quantile(diff_OP_pyr, probs = 0.025, na.rm=T), 
    X97.5._OP_pyr=quantile(diff_OP_pyr, probs = 0.975, na.rm=T),
    X2.5._IG2_OP=quantile(diff_IG2_OP, probs = 0.025, na.rm=T), 
    X97.5._IG2_OP=quantile(diff_IG2_OP, probs = 0.975, na.rm=T))


