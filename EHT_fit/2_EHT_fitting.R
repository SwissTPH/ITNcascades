rm(list=ls())
library(dplyr)
library(stringr)
library(tidyr)
library(rstan)

mainDir="."
scriptDir=file.path(".","EHT_fit")

stanDir=file.path(mainDir, "results/stan_outputs")
plotDir=file.path(mainDir, "plots")
dataDir=file.path(scriptDir, "processed_data")
outputsDir=file.path(mainDir, "results/csv_outputs")
source(file.path(scriptDir, "functions_fit_stan.R"))

niter=6000
nwarmup=3000
nchains=3
rerun=FALSE
######################################################
# KIBONDO 2022
######################################################

data_kibondo_72<-read.csv(file.path(dataDir,"kibondo_72.csv"))
data_kibondo_24<-read.csv(file.path(dataDir,"kibondo_24.csv"))

if(rerun){
  
  # unwashed
  fit_EHT_multinomial_stan(data=data_kibondo_72 %>% filter(wash_status %in% c("Unwashed", "Control") ),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_kibondo_unwashed_72.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  fit_EHT_multinomial_stan(data=data_kibondo_24 %>% filter(wash_status %in% c("Unwashed", "Control")),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_kibondo_unwashed_24.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  # washed
  fit_EHT_multinomial_stan(data=data_kibondo_72 %>% filter(wash_status %in% c("20x Washed", "Control")),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_kibondo_washed_72.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_kibondo_24 %>% filter(wash_status  %in% c("20x Washed", "Control")),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_kibondo_washed_24.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
}

results_kibondo_w72<- readRDS(file.path(stanDir,"/stan_kibondo_washed_72.rds"))
results_kibondo_unw72<- readRDS(file.path(stanDir,"/stan_kibondo_unwashed_72.rds"))
results_kibondo_w24<- readRDS(file.path(stanDir,"/stan_kibondo_washed_24.rds"))
results_kibondo_unw24<- readRDS(file.path(stanDir,"/stan_kibondo_unwashed_24.rds"))

estimates_kibondo_unwashed72=get_stan_summary_output(results_kibondo_unw72, decay="none",
                                                       save=FALSE,
                                                       path=NULL, data=data_kibondo_72) %>% 
  mutate(washed_status="Unwashed")


estimates_kibondo_washed72=get_stan_summary_output(results_kibondo_w72, decay="none",
                                                     save=FALSE,
                                                     path=NULL, data=data_kibondo_72) %>% 
  mutate(washed_status="Washed")

write.csv(rbind(estimates_kibondo_washed72, 
                estimates_kibondo_unwashed72), file.path(outputsDir, "estimates_kibondo_72.csv"))


estimates_kibondo_unwashed24=get_stan_summary_output(results_kibondo_unw24, decay="none",
                                                         save=FALSE,
                                                         path=NULL, data=data_kibondo_24) %>% 
  mutate(washed_status="Unwashed")


estimates_kibondo_washed24=get_stan_summary_output(results_kibondo_w24, decay="none",
                                                       save=FALSE,
                                                       path=NULL, data=data_kibondo_24) %>% 
  mutate(washed_status="Washed")

write.csv(rbind(estimates_kibondo_washed24, 
                estimates_kibondo_unwashed24), file.path(outputsDir, "estimates_kibondo_24.csv"))


######################################################
# BIT055
######################################################

data_bit055_24<-read.csv(file.path(dataDir,"bit055_24.csv"))

if(rerun){
  
  fit_EHT_multinomial_stan(data=data_bit055_24 %>% filter(washed==0 | washed=="Control"),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT055_unwashed_24.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
 fit_EHT_multinomial_stan(data=data_bit055_24 %>% filter(washed!=0),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT055_washed_24.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
}


results_bit055_w<- readRDS(file.path(stanDir,"/stan_BIT055_washed_24.rds"))
results_bit055_unw<- readRDS(file.path(stanDir,"/stan_BIT055_unwashed_24.rds"))

estimates_bit055_washed=get_stan_summary_output(results_bit055_w, decay="none",
                                                    save=FALSE,
                                                    path=NULL, data=data_bit055_24) %>% 
  mutate(washed_status="Washed")
estimates_bit055_unwashed=get_stan_summary_output(results_bit055_unw, decay="none",
                                                      save=FALSE,
                                                      path=NULL, data=data_bit055_24) %>% 
  mutate(washed_status="Unwashed")


write.csv(rbind(estimates_bit055_washed, 
                estimates_bit055_unwashed), file.path(outputsDir, "estimates_bit055_24.csv"))



######################################################
# BIT103
######################################################

data_bit103_72<-read.csv(file.path(dataDir,"assenga_all_72.csv"))
data_bit103_24<-read.csv(file.path(dataDir,"assenga_all_24.csv"))


if(rerun){

  fit_EHT_multinomial_stan(data=data_bit103_72 %>% dplyr::filter(wash_status==0 | wash_status=="Control"),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_assenga_unwashed_72_all.rds"),
                           stanpath = scriptDir, stanmodel  = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_bit103_72 %>% filter(wash_status==20| wash_status=="Control"),
                           iter=10000,
                           warmup=7000,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_assenga_washed_72_all.rds"),
                           stanpath = scriptDir, stanmodel  = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_bit103_24 %>% dplyr::filter(wash_status==0 | wash_status=="Control"),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_assenga_unwashed_24_all.rds"),
                           stanpath = scriptDir, stanmodel  = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_bit103_24 %>% filter(wash_status==20| wash_status=="Control"),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_assenga_washed_24_all.rds"),
                           stanpath = scriptDir, stanmodel  = "EHT_fitting_model.stan")
  
}

results_bit103_unw_72_all<- readRDS(file.path(stanDir,"/stan_assenga_unwashed_72_all.rds"))
results_bit103_w_72_all<- readRDS(file.path(stanDir,"/stan_assenga_washed_72_all.rds"))

results_bit103_unw_24_all<- readRDS(file.path(stanDir,"/stan_assenga_unwashed_24_all.rds"))
results_bit103_w_24_all<- readRDS(file.path(stanDir,"/stan_assenga_washed_24_all.rds"))


estimates_bit103_washed72_all=get_stan_summary_output(results_bit103_w_72_all, decay="none",
                                                      save=FALSE,
                                                      path=NULL, data=data_bit103_72) %>% 
  mutate(washed_status="Washed")
estimates_bit103_unwashed72_all=get_stan_summary_output(results_bit103_unw_72_all, decay="none",
                                                        save=FALSE,
                                                        path=NULL, data=data_bit103_72) %>% 
  mutate(washed_status="Unwashed")

write.csv(rbind(estimates_bit103_washed72_all, 
                estimates_bit103_unwashed72_all), file.path(outputsDir, "estimates_assenga_72_all.csv"))



estimates_bit103_washed24_all=get_stan_summary_output(results_bit103_w_24_all, decay="none",
                                                      save=FALSE,
                                                      path=NULL, data=data_bit103_24) %>% 
  mutate(washed_status="Washed")
estimates_bit103_unwashed24_all=get_stan_summary_output(results_bit103_unw_24_all, decay="none",
                                                        save=FALSE,
                                                        path=NULL, data=data_bit103_24) %>% 
  mutate(washed_status="Unwashed")

write.csv(rbind(estimates_bit103_washed24_all, 
                estimates_bit103_unwashed24_all), file.path(outputsDir, "estimates_assenga_24_all.csv"))



######################################################
# BIT103 Cote d'Ivoire
######################################################

data_bit103_CI_72<-read.csv(file.path(dataDir,"assenga_cotedivoire_72.csv"))
data_bit103_CI_24<-read.csv(file.path(dataDir,"assenga_cotedivoire_24.csv"))


if(rerun){
  
  fit_EHT_multinomial_stan(data=data_bit103_CI_72 %>% dplyr::filter(wash_status=="unwashed"),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_assengaCotedIvoire_unwashed_72.rds"),
                           stanpath = scriptDir, stanmodel  = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_bit103_CI_72 %>% filter(wash_status=="20x"),
                           iter=10000,
                           warmup=7000,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_assengaCotedIvoire_washed_72.rds"),
                           stanpath = scriptDir, stanmodel  = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_bit103_CI_24 %>% dplyr::filter(wash_status=="unwashed"),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_assengaCotedIvoire_unwashed_24.rds"),
                           stanpath = scriptDir, stanmodel  = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_bit103_CI_24 %>% filter(wash_status=="20x"),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_assengaCotedIvoire_washed_24.rds"),
                           stanpath = scriptDir, stanmodel  = "EHT_fitting_model.stan")
  
}

results_bit103_unw_72_CI<- readRDS(file.path(stanDir,"/stan_assengaCotedIvoire_unwashed_72.rds"))
results_bit103_w_72_CI<- readRDS(file.path(stanDir,"/stan_assengaCotedIvoire_washed_72.rds"))


estimates_bit103_washed72_CI=get_stan_summary_output(results_bit103_w_72_CI, decay="none",
                                                          save=FALSE,
                                                          path=NULL, data=data_bit103_CI_72) %>% 
  mutate(washed_status="Washed")
estimates_bit103_unwashed72_CI=get_stan_summary_output(results_bit103_unw_72_CI, decay="none",
                                                            save=FALSE,
                                                            path=NULL, data=data_bit103_CI_72) %>% 
  mutate(washed_status="Unwashed")

write.csv(rbind(estimates_bit103_washed72_CI, 
                estimates_bit103_unwashed72_CI), file.path(outputsDir, "estimates_assengaCotedIvoire_72.csv"))

## 24h
results_bit103_unw_24_CI<- readRDS(file.path(stanDir,"/stan_assengaCotedIvoire_unwashed_24.rds"))
results_bit103_w_24_CI<- readRDS(file.path(stanDir,"/stan_assengaCotedIvoire_washed_24.rds"))


estimates_bit103_washed24_CI=get_stan_summary_output(results_bit103_w_24_CI, decay="none",
                                                     save=FALSE,
                                                     path=NULL, data=data_bit103_CI_24) %>% 
  mutate(washed_status="Washed")
estimates_bit103_unwashed24_CI=get_stan_summary_output(results_bit103_unw_24_CI, decay="none",
                                                       save=FALSE,
                                                       path=NULL, data=data_bit103_CI_24) %>% 
  mutate(washed_status="Unwashed")

write.csv(rbind(estimates_bit103_washed24_CI, 
                estimates_bit103_unwashed24_CI), file.path(outputsDir, "estimates_assengaCotedIvoire_24.csv"))


######################################################
# BIT059, Odufuwa et al. 2024
######################################################

data_bit059<-read.csv(file.path(dataDir,"Odufuwa_24.csv"))

if(rerun){
  
  fit_EHT_multinomial_stan(data=data_bit059 %>% filter(wash_status==0 | wash_status=="Control"),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_odufuwa_unwashed_24.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_bit059 %>% filter(wash_status!=0),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_odufuwa_washed_24.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
}


results_bit059_w<- readRDS(file.path(stanDir,"/stan_odufuwa_washed_24.rds"))
results_bit059_unw<- readRDS(file.path(stanDir,"/stan_odufuwa_unwashed_24.rds"))

estimates_bit059_washed=get_stan_summary_output(results_bit059_w, decay="none",
                                                    save=FALSE,
                                                    path=NULL, data=data_bit059) %>% 
  mutate(washed_status="Washed")
estimates_bit059_unwashed=get_stan_summary_output(results_bit059_unw, decay="none",
                                                      save=FALSE,
                                                      path=NULL, data=data_bit059) %>% 
  mutate(washed_status="Unwashed")


write.csv(rbind(estimates_bit059_washed, 
                estimates_bit059_unwashed), file.path(outputsDir, "estimates_odufuwa_24.csv"))


######################################################
# BIT080
######################################################

data_bit080_72<-read.csv(file.path(dataDir,"bit080_72.csv"))
data_bit080_24<-read.csv(file.path(dataDir,"bit080_24.csv"))

if(rerun){
  
  fit_EHT_multinomial_stan(data=data_bit080_72 %>% filter(wash_status==0 | wash_status=="Control"),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT080_unwashed_72.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_bit080_72 %>% filter(wash_status!=0),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT080_washed_72.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_bit080_24 %>% filter(wash_status==0| wash_status=="Control"),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT080_unwashed_24.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_bit080_24 %>% filter(wash_status!=0),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT080_washed_24.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
}


results_bit080_w_72<- readRDS(file.path(stanDir,"/stan_BIT080_washed_72.rds"))
results_bit080_unw_72<- readRDS(file.path(stanDir,"/stan_BIT080_unwashed_72.rds"))

estimates_bit080_washed_72=get_stan_summary_output(results_bit080_w_72, decay="none",
                                                    save=FALSE,
                                                    path=NULL, data=data_bit080_72) %>% 
  mutate(washed_status="Washed")
estimates_bit080_unwashed_72=get_stan_summary_output(results_bit080_unw_72, decay="none",
                                                      save=FALSE,
                                                      path=NULL, data=data_bit080_72) %>% 
  mutate(washed_status="Unwashed")

results_bit080_w_24<- readRDS(file.path(stanDir,"/stan_BIT080_washed_24.rds"))
results_bit080_unw_24<- readRDS(file.path(stanDir,"/stan_BIT080_unwashed_24.rds"))

estimates_bit080_washed_24=get_stan_summary_output(results_bit080_w_24, decay="none",
                                                       save=FALSE,
                                                       path=NULL, data=data_bit080_24) %>% 
  mutate(washed_status="Washed")
estimates_bit080_unwashed_24=get_stan_summary_output(results_bit080_unw_24, decay="none",
                                                         save=FALSE,
                                                         path=NULL, data=data_bit080_24) %>% 
  mutate(washed_status="Unwashed")


write.csv(rbind(estimates_bit080_washed_72, 
                estimates_bit080_unwashed_72), file.path(outputsDir, "estimates_bit080_72.csv"))


write.csv(rbind(estimates_bit080_washed_24, 
                estimates_bit080_unwashed_24), file.path(outputsDir, "estimates_bit080_24.csv"))


######################################################
# Martin 2024
######################################################

data_martin_72<-read.csv(file.path(dataDir,"martin_gambiae_72.csv"))
data_martin_24<-read.csv(file.path(dataDir,"martin_gambiae_24.csv"))

data_martin_72_funestus<-read.csv(file.path(dataDir,"martin_funestus_72.csv"))
data_martin_24_funestus<-read.csv(file.path(dataDir,"martin_funestus_24.csv"))

# separate data per year, as robustness check
data_martin2022_72<-read.csv(file.path(dataDir,"martin2022_gambiae_72.csv"))
data_martin2020_72<-read.csv(file.path(dataDir,"martin2020_gambiae_72.csv"))

data_martin2022_72_funestus<-read.csv(file.path(dataDir,"martin2022_funestus_72.csv"))
data_martin2020_72_funestus<-read.csv(file.path(dataDir,"martin2020_funestus_72.csv"))




if(rerun){
  
  # unwashed
  fit_EHT_multinomial_stan(data=data_martin_72 %>% filter(wash_status %in% c("0", "Control") ),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_martin_unwashed_72_gambiae.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  # 36
  fit_EHT_multinomial_stan(data=data_martin_72 %>% filter(wash_status %in% c("36", "Control")),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_martin_36m_72_gambiae.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  
  # unwashed
  fit_EHT_multinomial_stan(data=data_martin_72_funestus %>% filter(wash_status %in% c("0", "Control") ),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_martin_unwashed_72_funestus.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  # 36
  fit_EHT_multinomial_stan(data=data_martin_72_funestus %>% filter(wash_status %in% c("36", "Control")),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_martin_36m_72_funestus.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  ###### 24h
  # unwashed
  fit_EHT_multinomial_stan(data=data_martin_24 %>% filter(wash_status %in% c("0", "Control") ),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_martin_unwashed_24_gambiae.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  # 36
  fit_EHT_multinomial_stan(data=data_martin_24 %>% filter(wash_status %in% c("36", "Control")),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_martin_36m_24_gambiae.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  
  # unwashed
  fit_EHT_multinomial_stan(data=data_martin_24_funestus %>% filter(wash_status %in% c("0", "Control") ),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_martin_unwashed_24_funestus.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  # 36
  fit_EHT_multinomial_stan(data=data_martin_24_funestus %>% filter(wash_status %in% c("36", "Control")),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_martin_36m_24_funestus.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  ###
  # separate the data per year
  
  # unwashed
  fit_EHT_multinomial_stan(data=data_martin2020_72 ,
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_martin2020_unwashed_72_gambiae.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  # 36
  fit_EHT_multinomial_stan(data=data_martin2022_72 ,
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_martin2022_36m_72_gambiae.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  
  # unwashed
  fit_EHT_multinomial_stan(data=data_martin2020_72_funestus ,
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_martin2020_unwashed_72_funestus.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  # 36
  fit_EHT_multinomial_stan(data=data_martin2022_72_funestus ,
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_martin2022_36m_72_funestus.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  
}

results_martin_w72<- readRDS(file.path(stanDir,"/stan_martin_36m_72_gambiae.rds"))
results_martin_unw72<- readRDS(file.path(stanDir,"/stan_martin_unwashed_72_gambiae.rds"))
results_martin_w24<- readRDS(file.path(stanDir,"/stan_martin_36m_24_gambiae.rds"))
results_martin_unw24<- readRDS(file.path(stanDir,"/stan_martin_unwashed_24_gambiae.rds"))


estimates_martin_unwashed72=get_stan_summary_output(results_martin_unw72, decay="none",
                                                        save=FALSE,
                                                        path=NULL, data=data_martin_72) %>% 
  mutate(washed_status="Unwashed")


estimates_martin_washed72=get_stan_summary_output(results_martin_w72, decay="none",
                                                      save=FALSE,
                                                      path=NULL, data=data_martin_72) %>% 
  mutate(washed_status="Washed")

write.csv(rbind(estimates_martin_washed72, 
                estimates_martin_unwashed72), file.path(outputsDir, "estimates_martin_72_gambiae.csv"))


results_martin_w72_funestus<- readRDS(file.path(stanDir,"/stan_martin_36m_72_funestus.rds"))
results_martin_unw72_funestus<- readRDS(file.path(stanDir,"/stan_martin_unwashed_72_funestus.rds"))

estimates_martin_unwashed72_funestus=get_stan_summary_output(results_martin_unw72_funestus, decay="none",
                                                        save=FALSE,
                                                        path=NULL, data=data_martin_72_funestus) %>% 
  mutate(washed_status="Unwashed")


estimates_martin_washed72_funestus=get_stan_summary_output(results_martin_w72_funestus, decay="none",
                                                      save=FALSE,
                                                      path=NULL, data=data_martin_72_funestus) %>% 
  mutate(washed_status="Washed")

write.csv(rbind(estimates_martin_washed72_funestus, 
                estimates_martin_unwashed72_funestus), file.path(outputsDir, "estimates_martin_72_funestus.csv"))


# 24h
estimates_martin_unwashed24=get_stan_summary_output(results_martin_unw24, decay="none",
                                                    save=FALSE,
                                                    path=NULL, data=data_martin_24) %>% 
  mutate(washed_status="Unwashed")


estimates_martin_washed24=get_stan_summary_output(results_martin_w24, decay="none",
                                                  save=FALSE,
                                                  path=NULL, data=data_martin_24) %>% 
  mutate(washed_status="Washed")

write.csv(rbind(estimates_martin_washed24, 
                estimates_martin_unwashed24), file.path(outputsDir, "estimates_martin_24_gambiae.csv"))


results_martin_w24_funestus<- readRDS(file.path(stanDir,"/stan_martin_36m_24_funestus.rds"))
results_martin_unw24_funestus<- readRDS(file.path(stanDir,"/stan_martin_unwashed_24_funestus.rds"))

estimates_martin_unwashed24_funestus=get_stan_summary_output(results_martin_unw24_funestus, decay="none",
                                                             save=FALSE,
                                                             path=NULL, data=data_martin_24_funestus) %>% 
  mutate(washed_status="Unwashed")


estimates_martin_washed24_funestus=get_stan_summary_output(results_martin_w24_funestus, decay="none",
                                                           save=FALSE,
                                                           path=NULL, data=data_martin_24_funestus) %>% 
  mutate(washed_status="Washed")

write.csv(rbind(estimates_martin_washed24_funestus, 
                estimates_martin_unwashed24_funestus), file.path(outputsDir, "estimates_martin_24_funestus.csv"))


######################################################
# Sovegnon 2024
######################################################

data_sovegnon_IG2<-read.csv(file.path(dataDir,"data_sovegnon_IG2_SRS_72.csv"))
data_sovegnon_IG1<-read.csv(file.path(dataDir,"data_sovegnon_IG1_LRS_72.csv"))

if(rerun){
  
  # unwashed
  fit_EHT_multinomial_stan(data=data_sovegnon_IG2%>% filter(wash_status %in% c("Unwashed", "Control")) ,
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_sovegnonIG2_unwashed.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_sovegnon_IG1%>% filter(wash_status %in% c("Unwashed", "Control")) ,
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_sovegnonIG1_unwashed.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  # washed
  
  fit_EHT_multinomial_stan(data=data_sovegnon_IG2%>% filter(wash_status %in% c("Aged", "Control")) ,
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_sovegnonIG2_aged.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_sovegnon_IG1%>% filter(wash_status %in% c("Aged", "Control")) ,
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_sovegnonIG1_aged.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
}

results_sovegnonIG2_unw<- readRDS(file.path(stanDir,"/stan_sovegnonIG2_unwashed.rds"))
results_sovegnonIG2_w<- readRDS(file.path(stanDir,"/stan_sovegnonIG2_aged.rds"))

results_sovegnonIG1_unw<- readRDS(file.path(stanDir,"/stan_sovegnonIG1_unwashed.rds"))
results_sovegnonIG1_w<- readRDS(file.path(stanDir,"/stan_sovegnonIG1_aged.rds"))


estimates_sovegnonIG2_unwashed72=get_stan_summary_output(results_sovegnonIG2_unw, decay="none",
                                                          save=FALSE,
                                                          path=NULL, data=data_sovegnon_IG2%>% filter(wash_status %in% c("Unwashed", "Control")) 
)%>% 
  mutate(washed_status="Unwashed")


estimates_sovegnonIG1_unwashed72=get_stan_summary_output(results_sovegnonIG1_unw, decay="none",
                                                         save=FALSE,
                                                         path=NULL, data=data_sovegnon_IG1%>% filter(wash_status %in% c("Unwashed", "Control")) 
)%>% 
  mutate(washed_status="Unwashed")

estimates_sovegnonIG2_washed72=get_stan_summary_output(results_sovegnonIG2_w, decay="none",
                                                        save=FALSE,
                                                        path=NULL, data=data_sovegnon_IG2 %>% filter(wash_status %in% c("Control", "Aged"))
) %>% 
  mutate(washed_status="Washed")


estimates_sovegnonIG1_washed72=get_stan_summary_output(results_sovegnonIG1_w, decay="none",
                                                       save=FALSE,
                                                       path=NULL, data=data_sovegnon_IG1 %>% filter(wash_status %in% c("Control", "Aged"))
) %>% 
  mutate(washed_status="Washed")

write.csv(rbind(estimates_sovegnonIG2_unwashed72, estimates_sovegnonIG1_unwashed72, 
                estimates_sovegnonIG2_washed72, estimates_sovegnonIG1_washed72), file.path(outputsDir, "estimates_sovegnon_72.csv"))


######################################################
# Nguessan 2016
######################################################

data_nguessan<-read.csv(file.path(dataDir,"nguessan_72.csv"))%>%
  mutate(dead=FD+UD)

if(rerun){
  
  # unwashed
  fit_EHT_multinomial_stan(data=data_nguessan%>% filter(wash_status %in% c("unwashed", "Control")) ,
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_nguessan_unwashed.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  # unwashed
  fit_EHT_multinomial_stan(data=data_nguessan%>% filter(wash_status %in% c("Control", "20x")) ,
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_nguessan_washed.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  
}

results_nguessan_unw<- readRDS(file.path(stanDir,"/stan_nguessan_unwashed.rds"))
results_nguessan_w<- readRDS(file.path(stanDir,"/stan_nguessan_washed.rds"))


estimates_nguessan_unwashed72=get_stan_summary_output(results_nguessan_unw, decay="none",
                                                          save=FALSE,
                                                          path=NULL, data=data_nguessan) %>% 
  mutate(washed_status="Unwashed")


estimates_nguessan_washed72=get_stan_summary_output(results_nguessan_w, decay="none",
                                                        save=FALSE,
                                                        path=NULL, data=data_nguessan) %>% 
  mutate(washed_status="Washed")

write.csv(rbind(estimates_nguessan_washed72, 
                estimates_nguessan_unwashed72), file.path(outputsDir, "estimates_nguessan_72.csv"))


##########################
# EXTRACT POSTERIOR MAXIMUM


map_kibondo_unwashed72=extract_stan_posteriormax(res=results_kibondo_unw72, data=data_kibondo_72) %>% 
  mutate(washed_status="Unwashed", EHT="Kibondo", reference="Kibondo et al.")
map_kibondo_washed72=extract_stan_posteriormax(res=results_kibondo_w72, data=data_kibondo_72) %>% 
  mutate(washed_status="Washed", EHT="Kibondo", reference="Kibondo et al.")
#

map_bit055_washed=extract_stan_posteriormax(results_bit055_w,data=data_bit055_24) %>% 
  mutate(washed_status="Washed", EHT="BIT055", reference="BIT055")
map_bit055_unwashed=extract_stan_posteriormax(results_bit055_unw,data=data_bit055_24) %>% 
  mutate(washed_status="Unwashed", EHT="BIT055", reference="BIT055")



map_bit103_washed72_all=extract_stan_posteriormax(results_bit103_w_72_all, data=data_bit103_ifakara_72) %>% 
  mutate(washed_status="Washed", EHT="Assenga", reference="Assenga et al., Tanzania")
map_bit103_unwashed72_all=extract_stan_posteriormax(results_bit103_unw_72_all,  data=data_bit103_ifakara_72) %>% 
  mutate(washed_status="Unwashed", EHT="Assenga", reference="Assenga et al., Tanzania")

map_bit103_washed72_CI=extract_stan_posteriormax(results_bit103_w_72_CI, data=data_bit103_CI_72) %>% 
  mutate(washed_status="Washed", EHT="Assenga, CI", reference="Assenga et al., Côte d'Ivoire")
map_bit103_unwashed72_CI=extract_stan_posteriormax(results_bit103_unw_72_CI, data=data_bit103_CI_72) %>% 
  mutate(washed_status="Unwashed", EHT="Assenga, CI", reference="Assenga et al., Côte d'Ivoire")



map_bit059_washed=extract_stan_posteriormax(results_bit059_w, data=data_bit059) %>% 
  mutate(washed_status="Washed", EHT="Odufuwa", reference="Odufuwa et al.")
map_bit059_unwashed=extract_stan_posteriormax(results_bit059_unw,data=data_bit059) %>% 
  mutate(washed_status="Unwashed", EHT="Odufuwa", reference="Odufuwa et al.")



map_bit080_washed_72=extract_stan_posteriormax(results_bit080_w_72,  data=data_bit080_72) %>% 
  mutate(washed_status="Washed", EHT="BIT080", reference="BIT080")
map_bit080_unwashed_72=extract_stan_posteriormax(results_bit080_unw_72, data=data_bit080_72) %>% 
  mutate(washed_status="Unwashed", EHT="BIT080", reference="BIT080")



map_martin_unwashed72=extract_stan_posteriormax(results_martin_unw72,  data=data_martin_72) %>% 
  mutate(washed_status="Unwashed", EHT="Martin", reference="Martin et al., An. gambiae")
map_martin_washed72=extract_stan_posteriormax(results_martin_w72,data=data_martin_72) %>% 
  mutate(washed_status="Washed", EHT="Martin", reference="Martin et al., An. gambiae")

map_martin_unwashed72_funestus=extract_stan_posteriormax(results_martin_unw72_funestus, data=data_martin_72_funestus) %>% 
  mutate(washed_status="Unwashed", EHT="Martin, f", reference="Martin et al., An. funestus")
map_martin_washed72_funestus=extract_stan_posteriormax(results_martin_w72_funestus,  data=data_martin_72_funestus) %>% 
  mutate(washed_status="Washed", EHT="Martin, f", reference="Martin et al., An. funestus")


map_sovegnonIG2_unwashed72=extract_stan_posteriormax(results_sovegnonIG2_unw,  data=data_sovegnon_IG2%>% filter(wash_status %in% c("Unwashed", "Control")) )%>% 
  mutate(washed_status="Unwashed", EHT="Sovegnon", reference="Sovegnon et al.")
map_sovegnonIG1_unwashed72=extract_stan_posteriormax(results_sovegnonIG1_unw,  data=data_sovegnon_IG1%>% filter(wash_status %in% c("Unwashed", "Control")) )%>% 
  mutate(washed_status="Unwashed", EHT="Sovegnon", reference="Sovegnon et al.", treatment=ifelse(treatment==-1, -2, treatment))

map_sovegnonIG2_washed72=extract_stan_posteriormax(results_sovegnonIG2_w,  data=data_sovegnon_IG2%>% filter(wash_status %in% c("Aged", "Control")) )%>% 
  mutate(washed_status="Washed", EHT="Sovegnon", reference="Sovegnon et al.")
map_sovegnonIG1_washed72=extract_stan_posteriormax(results_sovegnonIG1_w,  data=data_sovegnon_IG1%>% filter(wash_status %in% c("Aged", "Control")) )%>% 
  mutate(washed_status="Washed", EHT="Sovegnon", reference="Sovegnon et al.", treatment=ifelse(treatment==-1, -2, treatment))

map_nguessan_unwashed72=extract_stan_posteriormax(results_nguessan_unw,  data=data_nguessan) %>% 
  mutate(washed_status="Unwashed", EHT="Nguessan", reference="Nguessan et al.")
map_nguessan_washed72=extract_stan_posteriormax(results_nguessan_w,  data=data_nguessan) %>% 
  mutate(washed_status="Washed", EHT="Nguessan", reference="Nguessan et al.")



# formating the results
all_posteriormax=rbind(map_kibondo_unwashed72, map_kibondo_washed72,
                    map_bit055_unwashed, map_bit055_washed,
                    map_bit103_unwashed72_all, map_bit103_washed72_all,
                    map_bit103_unwashed72_CI, map_bit103_washed72_CI,
                    map_bit059_unwashed, map_bit059_washed,
                    map_bit080_unwashed_72, map_bit080_washed_72,
                    map_martin_unwashed72, map_martin_washed72 ,
                    map_martin_unwashed72_funestus, map_martin_washed72_funestus ,
                    map_sovegnonIG2_unwashed72, map_sovegnonIG1_unwashed72,map_sovegnonIG2_washed72,map_sovegnonIG1_washed72 ,
                    map_nguessan_unwashed72, map_nguessan_washed72
)%>%
  mutate(insecticide_name=ifelse(insecticide_name %in% c("IG2", "Interceptor_G2","InterceptorG2","Interceptor®G2", "Interceptor\xaeG2", "Interceptor G2", "interceptorG2", "IG2.Aged"), "Interceptor G2",
                                 ifelse(insecticide_name %in% c("Olyset_Plus", "OlysetPlus", "Olyset Plus", "olysetplus"), "Olyset Plus", 
                                        ifelse(is.na(insecticide_name) , "control", "Pyrethroid-only"))),
         param1=gsub("Initial", "", gsub("Efficacy", "", gsub("Rate" ,"",  param))),
         aged_time=ifelse(EHT=="Sovegnon", 2, 3)) %>%
  mutate(value=round(value, digits = 2))%>%
  select(insecticide_name, EHT, reference, param1, washed_status, value, aged_time, treatment)%>%
  rename(netType=insecticide_name, parameter=param1)%>%
  tidyr::pivot_wider(id_cols =c( netType, EHT, reference, parameter, aged_time, treatment), names_from = washed_status, values_from = c(value) )%>%
  mutate(halflife_insecticide=aged_time/(1-Washed/Unwashed),
         halflife_insecticide_old=3/(1-Washed/Unwashed)
  )
write.csv(all_posteriormax, file.path(scriptDir, "fitted_parameters_posteriormax.csv"), row.names = F)
