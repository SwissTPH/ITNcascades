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
scriptDir=file.path(mainDir,"EHT_fit")

stanDir=file.path(mainDir, "stan_outputs")
plotDir=file.path(mainDir, "plots")
dataDir=file.path(mainDir, "data")
outputsDir=file.path(mainDir, "csv_outputs")
source(file.path(scriptDir, "fit_stan_function_multinomial.R"))
source(file.path(scriptDir, "compute_vectorialCapacity.R"))

niter=6000
nwarmup=3000
nchains=3
rerun=FALSE
######################################################
# KIBONDO 2022
######################################################

data_kibondo_72<-read.csv(file.path(outputsDir,"kibondo_reformat_72.csv"))
data_kibondo_24<-read.csv(file.path(outputsDir,"kibondo_reformat_24.csv"))

data_kibondo2_72=data_kibondo_72 %>%
  filter(! (treatment==0 & wash_status!="Unwashed"))%>%
  mutate(wash_status=ifelse(treatment==0, "Control", wash_status))

data_kibondo2_24=data_kibondo_24 %>%
  filter(! (treatment==0 & wash_status!="Unwashed"))%>%
  mutate(wash_status=ifelse(treatment==0, "Control", wash_status))

if(rerun){
  
  # unwashed
  fit_EHT_multinomial_stan(data=data_kibondo_72 %>% filter(wash_status =="Unwashed"),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_kibondo_unwashed_72.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  fit_EHT_multinomial_stan(data=data_kibondo_24 %>% filter(wash_status =="Unwashed"),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_kibondo_unwashed_24.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  # washed
  fit_EHT_multinomial_stan(data=data_kibondo2_72 %>% filter(wash_status !="Unwashed"),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_kibondo_washed_controlUnw_72.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_kibondo2_24 %>% filter(wash_status !="Unwashed"),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_kibondo_washed_controlUnw_24.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
}

results_kibondo_w72<- readRDS(file.path(stanDir,"/stan_kibondo_washed_controlUnw_72.rds"))
results_kibondo_unw72<- readRDS(file.path(stanDir,"/stan_kibondo_unwashed_72.rds"))
results_kibondo_w24<- readRDS(file.path(stanDir,"/stan_kibondo_washed_controlUnw_24.rds"))
results_kibondo_unw24<- readRDS(file.path(stanDir,"/stan_kibondo_unwashed_24.rds"))

estimates_kibondo_unwashed72=get_IRSmodel_summary_output(results_kibondo_unw72, decay="none",
                                                       save=FALSE,
                                                       path=NULL, data=data_kibondo_72) %>% 
  mutate(washed_status="Unwashed")


estimates_kibondo_washed72=get_IRSmodel_summary_output(results_kibondo_w72, decay="none",
                                                     save=FALSE,
                                                     path=NULL, data=data_kibondo_72) %>% 
  mutate(washed_status="Washed")

write.csv(rbind(estimates_kibondo_washed72, 
                estimates_kibondo_unwashed72), file.path(outputsDir, "estimates_kibondo_washed_72.csv"))


estimates_kibondo_unwashed24=get_IRSmodel_summary_output(results_kibondo_unw24, decay="none",
                                                         save=FALSE,
                                                         path=NULL, data=data_kibondo_24) %>% 
  mutate(washed_status="Unwashed")


estimates_kibondo_washed24=get_IRSmodel_summary_output(results_kibondo_w24, decay="none",
                                                       save=FALSE,
                                                       path=NULL, data=data_kibondo_24) %>% 
  mutate(washed_status="Washed")

write.csv(rbind(estimates_kibondo_washed24, 
                estimates_kibondo_unwashed24), file.path(outputsDir, "estimates_kibondo_washed_24.csv"))


######################################################
# BIT055
######################################################

data_bit055_24<-read.csv(file.path(outputsDir,"bit055_reformat_24.csv"))

data_bit055_24_2=data_bit055_24 %>%
  filter(! (treatment==0 & washed!=0))%>%
  mutate(washed=ifelse(treatment==0, 1, washed))

if(rerun){
  
  fit_EHT_multinomial_stan(data=data_bit055_24 %>% filter(washed==0),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT055_unwashed_24.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
 fit_EHT_multinomial_stan(data=data_bit055_24_2 %>% filter(washed!=0),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT055_washed_controlUnw_24.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
}


results_bit055_w<- readRDS(file.path(stanDir,"/stan_BIT055_washed_controlUnw_24.rds"))
results_bit055_unw<- readRDS(file.path(stanDir,"/stan_BIT055_unwashed_24.rds"))

estimates_bit055_washed=get_IRSmodel_summary_output(results_bit055_w, decay="none",
                                                    save=FALSE,
                                                    path=NULL, data=data_bit055_24) %>% 
  mutate(washed_status="Washed")
estimates_bit055_unwashed=get_IRSmodel_summary_output(results_bit055_unw, decay="none",
                                                      save=FALSE,
                                                      path=NULL, data=data_bit055_24) %>% 
  mutate(washed_status="Unwashed")


write.csv(rbind(estimates_bit055_washed, 
                estimates_bit055_unwashed), file.path(outputsDir, "estimates_bit055_washed_24.csv"))



######################################################
# BIT103
######################################################

data_bit103_ifakara_72<-read.csv(file.path(outputsDir,"bit103_reformat_new_ifakara_72.csv"))
data_bit103_ifakara_24<-read.csv(file.path(outputsDir,"bit103_reformat_new_ifakara_24.csv"))

data_bit103_ifakara_72_2=data_bit103_ifakara_72 %>%
  filter(! (treatment==0 & wash_status!="0"))%>%
  mutate(wash_status=ifelse(treatment==0, 20, wash_status))

data_bit103_ifakara_24_2=data_bit103_ifakara_24 %>%
  filter(! (treatment==0 & wash_status!="0"))%>%
  mutate(wash_status=ifelse(treatment==0, 20, wash_status))



if(rerun){

  fit_EHT_multinomial_stan(data=data_bit103_ifakara_72 %>% filter(wash_status==0),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT103_unwashed_72_ifakara.rds"),
                           stanpath = scriptDir, stanmodel  = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_bit103_ifakara_24 %>% filter(wash_status==0),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT103_unwashed_24_ifakara.rds"),
                           stanpath = scriptDir, stanmodel  = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_bit103_ifakara_72_2 %>% filter(wash_status==20),
                           iter=10000,
                           warmup=7000,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT103_washed_controlUnw_72_ifakara.rds"),
                           stanpath = scriptDir, stanmodel  = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_bit103_ifakara_24_2 %>% filter(wash_status==20),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT103_washed_controlUnw_24_ifakara.rds"),
                           stanpath = scriptDir, stanmodel  = "EHT_fitting_model.stan")
  
}

results_bit103_unw_72<- readRDS(file.path(stanDir,"/stan_BIT103_unwashed_72_ifakara.rds"))
results_bit103_w_72<- readRDS(file.path(stanDir,"/stan_BIT103_washed_controlUnw_72_ifakara.rds"))
results_bit103_unw_24<- readRDS(file.path(stanDir,"/stan_BIT103_unwashed_24_ifakara.rds"))
results_bit103_w_24<- readRDS(file.path(stanDir,"/stan_BIT103_washed_controlUnw_24_ifakara.rds"))

estimates_bit103_washed24=get_IRSmodel_summary_output(results_bit103_w_24, decay="none",
                                                    save=FALSE,
                                                    path=NULL, data=data_bit103_ifakara_24) %>% 
  mutate(washed_status="Washed")
estimates_bit103_unwashed24=get_IRSmodel_summary_output(results_bit103_unw_24, decay="none",
                                                      save=FALSE,
                                                      path=NULL, data=data_bit103_ifakara_24) %>% 
  mutate(washed_status="Unwashed")


write.csv(rbind(estimates_bit103_washed24, 
                estimates_bit103_unwashed24), file.path(outputsDir, "estimates_bit103_washed24.csv"))

estimates_bit103_washed72=get_IRSmodel_summary_output(results_bit103_w_72, decay="none",
                                                    save=FALSE,
                                                    path=NULL, data=data_bit103_ifakara_72_2) %>% 
  mutate(washed_status="Washed")
estimates_bit103_unwashed72=get_IRSmodel_summary_output(results_bit103_unw_72, decay="none",
                                                      save=FALSE,
                                                      path=NULL, data=data_bit103_ifakara_72_2) %>% 
  mutate(washed_status="Unwashed")

write.csv(rbind(estimates_bit103_washed72, 
                estimates_bit103_unwashed72), file.path(outputsDir, "estimates_bit103_washed72.csv"))


######################################################
# BIT059, Odufuwa et al. 2024
######################################################

data_bit059<-read.csv(file.path(outputsDir,"bit059_reformat_24.csv"))

data_bit059_2=data_bit059 %>%
  filter(! (treatment==0 & wash_status!=0))%>%
  mutate(wash_status=ifelse(treatment==0, 1, wash_status))

if(rerun){
  
  fit_EHT_multinomial_stan(data=data_bit059 %>% filter(wash_status==0),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT059_unwashed_24.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_bit059_2 %>% filter(wash_status!=0),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT059_washed_controlUnw_24.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
}


results_bit059_w<- readRDS(file.path(stanDir,"/stan_BIT059_washed_controlUnw_24.rds"))
results_bit059_unw<- readRDS(file.path(stanDir,"/stan_BIT059_unwashed_24.rds"))

estimates_bit059_washed=get_IRSmodel_summary_output(results_bit059_w, decay="none",
                                                    save=FALSE,
                                                    path=NULL, data=data_bit059) %>% 
  mutate(washed_status="Washed")
estimates_bit059_unwashed=get_IRSmodel_summary_output(results_bit059_unw, decay="none",
                                                      save=FALSE,
                                                      path=NULL, data=data_bit059) %>% 
  mutate(washed_status="Unwashed")


write.csv(rbind(estimates_bit059_washed, 
                estimates_bit059_unwashed), file.path(outputsDir, "estimates_bit059_washed_24.csv"))


######################################################
# BIT080
######################################################

data_bit080_72<-read.csv(file.path(outputsDir,"bit080_reformat_72.csv"))
data_bit080_24<-read.csv(file.path(outputsDir,"bit080_reformat_24.csv"))

data_bit080_72_2=data_bit080_72 %>%
  filter(! (treatment==0 & wash_status!=0))%>%
  mutate(wash_status=ifelse(treatment==0, 1, wash_status))

data_bit080_24_2=data_bit080_24 %>%
  filter(! (treatment==0 & wash_status!=0))%>%
  mutate(wash_status=ifelse(treatment==0, 1, wash_status))

if(rerun){
  
  fit_EHT_multinomial_stan(data=data_bit080_72 %>% filter(wash_status==0),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT080_unwashed_72.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_bit080_72_2 %>% filter(wash_status!=0),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT080_washed_controlUnw_72.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_bit080_24 %>% filter(wash_status==0),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT080_unwashed_24.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
  
  fit_EHT_multinomial_stan(data=data_bit080_24_2 %>% filter(wash_status!=0),
                           iter=niter,
                           warmup=nwarmup,
                           chains=nchains,
                           path=file.path(stanDir,"/stan_BIT080_washed_controlUnw_24.rds"),
                           stanpath = scriptDir, stanmodel = "EHT_fitting_model.stan")
}


results_bit080_w_72<- readRDS(file.path(stanDir,"/stan_BIT080_washed_controlUnw_72.rds"))
results_bit080_unw_72<- readRDS(file.path(stanDir,"/stan_BIT080_unwashed_72.rds"))

estimates_bit080_washed_72=get_IRSmodel_summary_output(results_bit080_w_72, decay="none",
                                                    save=FALSE,
                                                    path=NULL, data=data_bit080_72) %>% 
  mutate(washed_status="Washed")
estimates_bit080_unwashed_72=get_IRSmodel_summary_output(results_bit080_unw_72, decay="none",
                                                      save=FALSE,
                                                      path=NULL, data=data_bit080_72) %>% 
  mutate(washed_status="Unwashed")

results_bit080_w_24<- readRDS(file.path(stanDir,"/stan_BIT080_washed_controlUnw_24.rds"))
results_bit080_unw_24<- readRDS(file.path(stanDir,"/stan_BIT080_unwashed_24.rds"))

estimates_bit080_washed_24=get_IRSmodel_summary_output(results_bit080_w_24, decay="none",
                                                       save=FALSE,
                                                       path=NULL, data=data_bit080_24) %>% 
  mutate(washed_status="Washed")
estimates_bit080_unwashed_24=get_IRSmodel_summary_output(results_bit080_unw_24, decay="none",
                                                         save=FALSE,
                                                         path=NULL, data=data_bit080_24) %>% 
  mutate(washed_status="Unwashed")


write.csv(rbind(estimates_bit080_washed_72, 
                estimates_bit080_unwashed_72), file.path(outputsDir, "estimates_bit080_washed_72.csv"))


write.csv(rbind(estimates_bit080_washed_24, 
                estimates_bit080_unwashed_24), file.path(outputsDir, "estimates_bit080_washed_24.csv"))




#######################
# PLOT SUMMARY
######################
estimates_kibondo_72=read.csv(file.path(outputsDir, "estimates_kibondo_washed_72.csv"))
estimates_kibondo_24=read.csv(file.path(outputsDir, "estimates_kibondo_washed_24.csv"))

estimates_bit059=read.csv(file.path(outputsDir, "estimates_bit059_washed_24.csv"))

estimates_bit055=read.csv(file.path(outputsDir, "estimates_bit055_washed_24.csv"))

estimates_bit080_72=read.csv(file.path(outputsDir, "estimates_bit080_washed_72.csv"))
estimates_bit080_24=read.csv(file.path(outputsDir, "estimates_bit080_washed_24.csv"))

estimates_bit103_72=read.csv(file.path(outputsDir, "estimates_bit103_washed72.csv"))
estimates_bit103_24=read.csv(file.path(outputsDir, "estimates_bit103_washed24.csv"))

#############################
# endpoint: longest available per trial

# formating the results
all_estimates=rbind(estimates_kibondo_72 %>%mutate(EHT="Kibondo et al. 2022", EHT_short="Kibondo", insecticide_name="InterceptorG2"),
                    estimates_bit055 %>%mutate(EHT="BIT055", EHT_short=EHT),
                    estimates_bit059 %>%mutate(EHT="BIT059", EHT_short=EHT),
                    estimates_bit080_72 %>%mutate(EHT="BIT080", EHT_short=EHT),
                    estimates_bit103_72 %>%mutate(EHT="BIT103", EHT_short=EHT)
)%>%
  mutate(insecticide_name=ifelse(insecticide_name %in% c("IG2", "Interceptor_G2","InterceptorG2","Interceptor®G2", "Interceptor\xaeG2"), "Interceptor G2",
                                 ifelse(insecticide_name %in% c("Olyset_Plus", "OlysetPlus"), "Olyset Plus", insecticide_name)),
         param2=gsub("Initial", "", gsub("Efficacy", "", gsub("Rate" ,"",  param))))%>%
  mutate(EHT=gsub(" et al.", "\net al.", EHT))

all_estimates$param2=factor(all_estimates$param2, levels=c("Repellency",  "KillingDuringHostSeeking","Preprandialkilling", "Postprandialkilling"),
                            labels=c("Reduction in host availability",  "Increase in host seeking mortality",  "Pre-prandial killing effect","Post-prandial killing effect"))

all_estimates_simple=all_estimates %>%
  mutate(washed_status=ifelse(washed_status=="Washed20", "Washed", washed_status),
         label_height=ifelse(param2=="Increase in host seeking mortality",-0.07, -0.03))

whitepointforscale=all_estimates_simple%>% filter(EHT_short=="Kibondo", washed_status=="Washed")%>%
  mutate(mean=NA, X2.5.=NA, X97.5.=NA,
         test=1)

# plot all estimates
all_estimates_simple %>%
  mutate( washed_status=ifelse(washed_status=="Unwashed", "Unwashed", "Washed 20x"))%>%
  ggplot()+
  geom_col(aes(x=insecticide_name, y=mean, fill=insecticide_name, group=interaction( washed_status, EHT), alpha=washed_status), stat="identity", position=position_dodge(preserve = "single"))+
  geom_errorbar(aes(x=insecticide_name, ymin=X2.5., ymax=X97.5., color=insecticide_name,  group=interaction( washed_status, EHT)),position=position_dodge(0.9, preserve = "single"), width=0.4)+
  geom_point(data=whitepointforscale, aes(x=insecticide_name, y=test), width=0.4, color="white")+
  #facet_grid( . ~param2, scales = "free")+
  facet_wrap( . ~param2, scales = "free")+
  theme_bw()+#ylim(0,1)+
  labs(x="", y="", color="", fill="", alpha="")+
  scale_alpha_manual(values=c(  0.8, 0.5),na.translate=FALSE)+
  scale_color_manual(values=c("darkorange","dodgerblue","darkblue", "red", "orange"  ))+
  scale_fill_manual(values=c("darkorange","dodgerblue","darkblue", "red", "orange" ))+
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 14),
        axis.text=element_text(size=14), strip.background =element_rect(fill="white"), legend.text =element_text(size=14) )+ 
  geom_text(aes(label=EHT_short, x=insecticide_name, y=label_height,group=EHT),
            position = position_dodge(width=0.9),size=3.5)+#ylim(-0.05,1)+ 
  guides(fill = "none",color = "none")
ggsave(file.path(plotDir, "plot_EHTfit_all.png"), width=9, height=12)




final_table=all_estimates_simple %>%
  select(insecticide_name, EHT_short, param2, washed_status, mean, X2.5., X97.5.)%>%
  rename(netType=insecticide_name, parameter=param2, q025=X2.5., q975=X97.5., EHT=EHT_short)%>%
  tidyr::pivot_wider(id_cols =c( netType, EHT, parameter), names_from = washed_status, values_from = c(mean, q025, q975) )%>%
  mutate(halflife_insecticide=3/(1-mean_Washed/mean_Unwashed)
  )

write.csv(final_table, file.path(scriptDir, "../FOR_OM_USERS/fitted_parameters_all_new.csv"), row.names = F)

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



######################
# 24H EVERYWHERE

# formating the results
all_estimates_24=rbind(estimates_kibondo_24 %>%mutate(EHT="Kibondo et al. 2022", EHT_short="Kibondo", insecticide_name="InterceptorG2"),
                    estimates_bit055 %>%mutate(EHT="BIT055", EHT_short=EHT),
                    estimates_bit059 %>%mutate(EHT="BIT059", EHT_short=EHT),
                    estimates_bit080_24 %>%mutate(EHT="BIT080", EHT_short=EHT),
                    estimates_bit103_24 %>%mutate(EHT="BIT103", EHT_short=EHT)
)%>%
  mutate(insecticide_name=ifelse(insecticide_name %in% c("IG2", "Interceptor_G2", "InterceptorG2","Interceptor®G2", "Interceptor\xaeG2"), "Interceptor G2",
                                 ifelse(insecticide_name %in% c("Olyset_Plus", "OlysetPlus"), "Olyset Plus", insecticide_name)),
         param2=gsub("Initial", "", gsub("Efficacy", "", gsub("Rate" ,"",  param))))%>%
  mutate(EHT=gsub(" et al.", "\net al.", EHT))

all_estimates_24$param2=factor(all_estimates_24$param2, levels=c("Repellency",  "KillingDuringHostSeeking","Preprandialkilling", "Postprandialkilling"),
                            labels=c("Reduction in host availability",  "Increase in host seeking mortality",  "Pre-prandial killing effect","Post-prandial killing effect"))

all_estimates_24_simple=all_estimates_24 %>%
  filter(insecticide_name %in% c("Interceptor G2", "Olyset Plus"))%>%
  mutate(washed_status=ifelse(washed_status=="Washed20", "Washed", washed_status),
         label_height=ifelse(param2=="Increase in host seeking mortality",-0.07, -0.03))

whitepointforscale=all_estimates_24_simple%>% filter(EHT_short=="Kibondo", washed_status=="Washed")%>%
  mutate(mean=NA, X2.5.=NA, X97.5.=NA,
         test=c(1,1,1,NA))

# plot all estimates
all_estimates_24_simple %>%
  mutate( washed_status=ifelse(washed_status=="Unwashed", "Unwashed", "Washed 20x"))%>%
  ggplot()+
  geom_col(aes(x=insecticide_name, y=mean, fill=insecticide_name, group=interaction( washed_status, EHT), alpha=washed_status), stat="identity", position=position_dodge(preserve = "single"))+
  geom_errorbar(aes(x=insecticide_name, ymin=X2.5., ymax=X97.5., color=insecticide_name,  group=interaction( washed_status, EHT)),position=position_dodge(0.9, preserve = "single"), width=0.4)+
  geom_point(data=whitepointforscale, aes(x=insecticide_name, y=test), width=0.4, color="white")+
  #facet_grid( . ~param2, scales = "free")+
  facet_wrap( . ~param2, scales = "free")+
  theme_bw()+#ylim(0,1)+
  labs(x="", y="", color="", fill="", alpha="")+
  scale_alpha_manual(values=c(  0.8, 0.5),na.translate=FALSE)+
  scale_color_manual(values=c("darkorange","dodgerblue","darkblue", "red", "orange"  ))+
  scale_fill_manual(values=c("darkorange","dodgerblue","darkblue", "red", "orange" ))+
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 14),
        axis.text=element_text(size=14), strip.background =element_rect(fill="white"), legend.text =element_text(size=14) )+ 
  geom_text(aes(label=EHT_short, x=insecticide_name, y=label_height,group=EHT),
            position = position_dodge(width=0.9),size=3.5)+#ylim(-0.05,1)+ 
  guides(fill = "none",color = "none")
ggsave(file.path(plotDir, "plot_EHTfit_24_all.png"), width=9, height=12)




final_table_24=all_estimates_24_simple %>%
  select(insecticide_name, EHT_short, param2, washed_status, mean, X2.5., X97.5.)%>%
  rename(netType=insecticide_name, parameter=param2, q025=X2.5., q975=X97.5., EHT=EHT_short)%>%
  tidyr::pivot_wider(id_cols =c( netType, EHT, parameter), names_from = washed_status, values_from = c(mean, q025, q975) )%>%
  mutate(halflife_insecticide=3/(1-mean_Washed/mean_Unwashed)
  )

write.csv(final_table_24, file.path(scriptDir, "../FOR_OM_USERS/fitted_parameters_24_all_new.csv"), row.names = F)

final_table_24 %>%
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
  gt::gtsave(filename = file.path(plotDir, "outputs_ento_efficacy_24_all.docx"))





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

#Compile the mode;
model_noRhythms = build_model_obj(vec_p=ent_params, hosts_p= def_host_params(), activity=activity_noRhythms, total_pop=2000)


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
results_bit080_w_24<- readRDS(file.path(stanDir,"/stan_BIT080_washed_controlUnw_24.rds"))
results_bit080_unw_24<- readRDS(file.path(stanDir,"/stan_BIT080_unwashed_24.rds"))


calculate_VC_uncertainty=function(results_EHT, nsamples, names, insecticide_id){
  
  npoints=3*365
  L=3
  kappa=2
  pyrethroid_params=c(0, 0.05, 0)
  names(pyrethroid_params)=c("Deterrency","PrePrandial", "PostPrandial" )
  
  pyrethroid_interv=get_interv(myproba=pyrethroid_params,L=L*365, kappa=kappa, npoints=npoints,
                               model_p=model_noRhythms, name="pyrethroid", intervention_type="LLINs", decay="weibull",halflife_insecticide=NULL)
  
  
  intervention_list_kibondo <- interventions_vectorial_capacity_wrapper(results=results_EHT,results_washed20 =NULL,decay = "weibull",
                                                                        model_p=model_noRhythms,
                                                                        insecticides=insecticide_id,L = L*365, kappa=kappa,
                                                                        names=names,  npoints=npoints,washedDecay=FALSE, uncertainty=TRUE, 
                                                                        cov=1, pyrethroid = pyrethroid_interv, nsamples=nsamples)
  
  
  
  
  
}

nsamples=1000

if(rerun){
  # bednet efficacy
  VCred_kibondo=calculate_VC_uncertainty(results_EHT=results_kibondo_unw_72, nsamples=nsamples, names="IG2", insecticide_id=2)
  VCred_bit103_ig2=calculate_VC_uncertainty(results_EHT=results_bit103_unw_72, nsamples=nsamples, names="IG2", insecticide_id=3)
  VCred_bit080=calculate_VC_uncertainty(results_EHT=results_bit080_unw_72, nsamples=nsamples, names="IG2", insecticide_id=2)
  
  VCred_bit103_pbo=calculate_VC_uncertainty(results_EHT=results_bit103_unw_72, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit059=calculate_VC_uncertainty(results_EHT=results_bit059_unw, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit055=calculate_VC_uncertainty(results_EHT=results_bit055_unw, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  
  
  VCred_kibondo_w=calculate_VC_uncertainty(results_EHT=results_kibondo_w_72, nsamples=nsamples, names="IG2", insecticide_id=2)
  VCred_bit103_ig2_w=calculate_VC_uncertainty(results_EHT=results_bit103_w_72, nsamples=nsamples, names="IG2", insecticide_id=3)
  VCred_bit080_w=calculate_VC_uncertainty(results_EHT=results_bit080_w_72, nsamples=nsamples, names="IG2", insecticide_id=2)
  
  VCred_bit103_pbo_w=calculate_VC_uncertainty(results_EHT=results_bit103_w_72, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit059_w=calculate_VC_uncertainty(results_EHT=results_bit059_w, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit055_w=calculate_VC_uncertainty(results_EHT=results_bit055_w, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  
  
  impacts=rbind(data.frame(VCred=VCred_kibondo) %>% mutate(EHT_short="Kibondo", insecticide_name="Interceptor G2"),
                data.frame(VCred=VCred_bit055) %>% mutate(EHT_short="BIT055", insecticide_name="Olyset Plus"),
                data.frame(VCred=VCred_bit103_ig2) %>% mutate(EHT_short="BIT103", insecticide_name="Interceptor G2"),
                data.frame(VCred=VCred_bit103_pbo) %>% mutate(EHT_short="BIT103", insecticide_name="Olyset Plus"),
                data.frame(VCred=VCred_bit059) %>% mutate(EHT_short="BIT059", insecticide_name="Olyset Plus"),
                data.frame(VCred=VCred_bit080) %>% mutate(EHT_short="BIT080", insecticide_name="Interceptor G2")
  )%>% mutate(washed_status="Unwashed")
  
  
  impacts_w=rbind(data.frame(VCred=VCred_kibondo_w) %>% mutate(EHT_short="Kibondo", insecticide_name="Interceptor G2"),
                  data.frame(VCred=VCred_bit055_w) %>% mutate(EHT_short="BIT055", insecticide_name="Olyset Plus"),
                  data.frame(VCred=VCred_bit103_ig2_w) %>% mutate(EHT_short="BIT103", insecticide_name="Interceptor G2"),
                  data.frame(VCred=VCred_bit103_pbo_w) %>% mutate(EHT_short="BIT103", insecticide_name="Olyset Plus"),
                  data.frame(VCred=VCred_bit059_w) %>% mutate(EHT_short="BIT059", insecticide_name="Olyset Plus"),
                  data.frame(VCred=VCred_bit080_w) %>% mutate(EHT_short="BIT080", insecticide_name="Interceptor G2")
  )%>% mutate(washed_status="Washed")
  
  
  write.csv(rbind(impacts,impacts_w), file = file.path(outputsDir, "impacts_EHT.csv"))
  
  impacts %>% group_by(EHT_short, insecticide_name,washed_status)%>% summarise(mean=mean(VCred), min=min(VCred), max=max(VCred))
  
}


if(rerun){
  # bednet efficacy
  VCred24_kibondo=calculate_VC_uncertainty(results_EHT=results_kibondo_unw_24, nsamples=nsamples, names="IG2", insecticide_id=2)
  VCred24_bit103_ig2=calculate_VC_uncertainty(results_EHT=results_bit103_unw_24, nsamples=nsamples, names="IG2", insecticide_id=3)
  VCred24_bit080=calculate_VC_uncertainty(results_EHT=results_bit080_unw_24, nsamples=nsamples, names="IG2", insecticide_id=2)
  
  VCred24_bit103_pbo=calculate_VC_uncertainty(results_EHT=results_bit103_unw_24, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit059=calculate_VC_uncertainty(results_EHT=results_bit059_unw, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit055=calculate_VC_uncertainty(results_EHT=results_bit055_unw, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  
  
  VCred24_kibondo_w=calculate_VC_uncertainty(results_EHT=results_kibondo_w_24, nsamples=nsamples, names="IG2", insecticide_id=2)
  VCred24_bit103_ig2_w=calculate_VC_uncertainty(results_EHT=results_bit103_w_24, nsamples=nsamples, names="IG2", insecticide_id=3)
  VCred24_bit080_w=calculate_VC_uncertainty(results_EHT=results_bit080_w_24, nsamples=nsamples, names="IG2", insecticide_id=2)
  
  VCred24_bit103_pbo_w=calculate_VC_uncertainty(results_EHT=results_bit103_w_24, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit059_w=calculate_VC_uncertainty(results_EHT=results_bit059_w, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit055_w=calculate_VC_uncertainty(results_EHT=results_bit055_w, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  
  
  impacts24=rbind(data.frame(VCred=VCred24_kibondo) %>% mutate(EHT_short="Kibondo", insecticide_name="Interceptor G2"),
                data.frame(VCred=VCred_bit055) %>% mutate(EHT_short="BIT055", insecticide_name="Olyset Plus"),
                data.frame(VCred=VCred24_bit103_ig2) %>% mutate(EHT_short="BIT103", insecticide_name="Interceptor G2"),
                data.frame(VCred=VCred24_bit103_pbo) %>% mutate(EHT_short="BIT103", insecticide_name="Olyset Plus"),
                data.frame(VCred=VCred_bit059) %>% mutate(EHT_short="BIT059", insecticide_name="Olyset Plus"),
                data.frame(VCred=VCred24_bit080) %>% mutate(EHT_short="BIT080", insecticide_name="Interceptor G2")
  )%>% mutate(washed_status="Unwashed")
  
  
  impacts_w24=rbind(data.frame(VCred=VCred24_kibondo_w) %>% mutate(EHT_short="Kibondo", insecticide_name="Interceptor G2"),
                  data.frame(VCred=VCred_bit055_w) %>% mutate(EHT_short="BIT055", insecticide_name="Olyset Plus"),
                  data.frame(VCred=VCred24_bit103_ig2_w) %>% mutate(EHT_short="BIT103", insecticide_name="Interceptor G2"),
                  data.frame(VCred=VCred24_bit103_pbo_w) %>% mutate(EHT_short="BIT103", insecticide_name="Olyset Plus"),
                  data.frame(VCred=VCred_bit059_w) %>% mutate(EHT_short="BIT059", insecticide_name="Olyset Plus"),
                  data.frame(VCred=VCred24_bit080_w) %>% mutate(EHT_short="BIT080", insecticide_name="Interceptor G2")
  )%>% mutate(washed_status="Washed")
  
  
  write.csv(rbind(impacts24,impacts_w24), file = file.path(outputsDir, "impacts_EHT_24.csv"))
  
  impacts24 %>% group_by(EHT_short, insecticide_name,washed_status)%>% summarise(mean=mean(VCred), min=min(VCred), max=max(VCred))
  
}



#######################
# PLOT SUMMARY
######################
estimates_kibondo_72=read.csv(file.path(outputsDir, "estimates_kibondo_washed_72.csv"))
estimates_kibondo_24=read.csv(file.path(outputsDir, "estimates_kibondo_washed_24.csv"))

estimates_bit059=read.csv(file.path(outputsDir, "estimates_bit059_washed_24.csv"))

estimates_bit055=read.csv(file.path(outputsDir, "estimates_bit055_washed_24.csv"))

estimates_bit080_72=read.csv(file.path(outputsDir, "estimates_bit080_washed_72.csv"))
estimates_bit080_24=read.csv(file.path(outputsDir, "estimates_bit080_washed_24.csv"))

estimates_bit103_72=read.csv(file.path(outputsDir, "estimates_bit103_washed72.csv"))
estimates_bit103_24=read.csv(file.path(outputsDir, "estimates_bit103_washed24.csv"))


impacts_72=read.csv(file.path(outputsDir, "impacts_EHT.csv"))%>%
  mutate(EHT_short=ifelse(EHT_short=="BIT059","Odufuwa", EHT_short ))

impacts_72_summary=impacts_72 %>% group_by(EHT_short, insecticide_name,washed_status)%>% 
  summarise(mean=mean(VCred), X2.5.=quantile(VCred, probs = 0.025), X97.5.=quantile(VCred, probs = 0.975))%>%
  mutate(param ="VCred", EHT=ifelse(EHT_short=="Kibondo","Kibondo et al. 2022",
                                    ifelse(EHT_short=="Odufuwa","Odufuwa et al. 2024", EHT_short ) ), X=NA)


impacts_24=read.csv(file.path(outputsDir, "impacts_EHT_24.csv"))%>%
  mutate(EHT_short=ifelse(EHT_short=="BIT059","Odufuwa", EHT_short ))

impacts_24_summary=impacts_24 %>% group_by(EHT_short, insecticide_name,washed_status)%>% 
  summarise(mean=mean(VCred), X2.5.=quantile(VCred, probs = 0.025), X97.5.=quantile(VCred, probs = 0.975))%>%
  mutate(param ="VCred", EHT=ifelse(EHT_short=="Kibondo","Kibondo et al. 2022", 
                                    ifelse(EHT_short=="Odufuwa","Odufuwa et al. 2024", EHT_short ) ), X=NA)

#############################
# endpoint: longest available per trial

# formating the results
all_estimates=rbind(estimates_kibondo_72 %>%mutate(EHT="Kibondo et al. 2022", EHT_short="Kibondo", insecticide_name="InterceptorG2"),
                    estimates_bit055 %>%mutate(EHT="BIT055", EHT_short=EHT),
                    estimates_bit059 %>%mutate(EHT="Odufuwa et al. 2024", EHT_short="Odufuwa"),
                    estimates_bit080_72 %>%mutate(EHT="BIT080", EHT_short=EHT),
                    estimates_bit103_72 %>%mutate(EHT="BIT103", EHT_short=EHT),
                    impacts_72_summary
)%>%
  mutate(insecticide_name=ifelse(insecticide_name %in% c("IG2", "Interceptor_G2","InterceptorG2","Interceptor®G2", "Interceptor\xaeG2"), "Interceptor G2",
                                 ifelse(insecticide_name %in% c("Olyset_Plus", "OlysetPlus"), "Olyset Plus", insecticide_name)),
         param2=gsub("Initial", "", gsub("Efficacy", "", gsub("Rate" ,"",  param))))%>%
  mutate(EHT=gsub(" et al.", "\net al.", EHT))%>% 
  filter(param2 !="KillingDuringHostSeeking")

all_estimates$param2=factor(all_estimates$param2, levels=c("Repellency", "Preprandialkilling", "Postprandialkilling", "VCred"),
                            labels=c("Reduction in host availability",  "Pre-prandial killing effect","Post-prandial killing effect", "Reduction in vectorial capacity"))
all_estimates$EHT=factor(all_estimates$EHT, levels=c("BIT103", "Kibondo\net al. 2022", "BIT080", "BIT055", "Odufuwa\net al. 2024"))

all_estimates_simple=all_estimates %>%
  mutate(washed_status=ifelse(washed_status=="Washed20", "Washed", washed_status),
         label_height=ifelse(param2=="Increase in host seeking mortality",-0.07, -0.03))

whitepointforscale=all_estimates_simple%>% filter(EHT_short=="Kibondo", washed_status=="Washed")%>%
  mutate(mean=NA, X2.5.=NA, X97.5.=NA,
         test=1)

# plot all estimates
all_estimates_simple %>%
  mutate( washed_status=ifelse(washed_status=="Unwashed", "Unwashed", "Washed 20x"))%>%
  ggplot()+
  geom_col(aes(x=insecticide_name, y=mean, fill=insecticide_name, group=interaction( washed_status, EHT), alpha=washed_status), stat="identity", position=position_dodge(preserve = "single"))+
  geom_errorbar(aes(x=insecticide_name, ymin=X2.5., ymax=X97.5., color=insecticide_name,  group=interaction( washed_status, EHT)),position=position_dodge(0.9, preserve = "single"), width=0.4)+
  geom_point(data=whitepointforscale, aes(x=insecticide_name, y=test), width=0.4, color="white")+
  #facet_grid( . ~param2, scales = "free")+
  facet_wrap( . ~param2, scales = "free")+
  theme_bw()+#ylim(0,1)+
  labs(x="", y="", color="", fill="", alpha="")+
  scale_alpha_manual(values=c(  0.8, 0.5),na.translate=FALSE)+
  scale_color_manual(values=c("darkorange","dodgerblue","darkblue", "red", "orange"  ))+
  scale_fill_manual(values=c("darkorange","dodgerblue","darkblue", "red", "orange" ))+
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 14),
        axis.text=element_text(size=14), strip.background =element_rect(fill="white"), legend.text =element_text(size=14) )+ 
  geom_text(aes(label=EHT_short, x=insecticide_name, y=label_height,group=EHT),
            position = position_dodge(width=0.9),size=3.5)+#ylim(-0.05,1)+ 
  guides(fill = "none",color = "none")
ggsave(file.path(plotDir, "plot_EHTfit_all_VCred.png"), width=10, height=12)




final_table=all_estimates_simple %>%
  mutate(mean=round(mean, digits = 2),X2.5.=round(X2.5., digits = 2),X97.5.=round(X97.5., digits = 2))%>%
  select(insecticide_name, EHT_short, param2, washed_status, mean, X2.5., X97.5.)%>%
  rename(netType=insecticide_name, parameter=param2, q025=X2.5., q975=X97.5., EHT=EHT_short)%>%
  tidyr::pivot_wider(id_cols =c( netType, EHT, parameter), names_from = washed_status, values_from = c(mean, q025, q975) )%>%
  mutate(halflife_insecticide=3/(1-mean_Washed/mean_Unwashed)
  )

write.csv(final_table, file.path(scriptDir, "../FOR_OM_USERS/fitted_parameters_all_new_VCred.csv"), row.names = F)

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
  gt::gtsave(filename = file.path(plotDir, "outputs_ento_efficacy_all_VCred.docx"))


#########################################################################################
#############################
# endpoint: 24h everywhere

# formating the results
all_estimates=rbind(estimates_kibondo_24 %>%mutate(EHT="Kibondo et al. 2022", EHT_short="Kibondo", insecticide_name="InterceptorG2"),
                    estimates_bit055 %>%mutate(EHT="BIT055", EHT_short=EHT),
                    estimates_bit059 %>%mutate(EHT="Odufuwa et al. 2024", EHT_short=EHT),
                    estimates_bit080_24 %>%mutate(EHT="BIT080", EHT_short=EHT),
                    estimates_bit103_24 %>%mutate(EHT="BIT103", EHT_short=EHT),
                    impacts_24_summary
)%>%
  mutate(insecticide_name=ifelse(insecticide_name %in% c("IG2", "Interceptor_G2","InterceptorG2","Interceptor®G2", "Interceptor\xaeG2"), "Interceptor G2",
                                 ifelse(insecticide_name %in% c("Olyset_Plus", "OlysetPlus"), "Olyset Plus", insecticide_name)),
         param2=gsub("Initial", "", gsub("Efficacy", "", gsub("Rate" ,"",  param))))%>%
  mutate(EHT=gsub(" et al.", "\net al.", EHT))%>% 
  filter(param2 !="KillingDuringHostSeeking")

all_estimates$param2=factor(all_estimates$param2, levels=c("Repellency", "Preprandialkilling", "Postprandialkilling", "VCred"),
                            labels=c("Reduction in host availability",  "Pre-prandial killing effect","Post-prandial killing effect", "Reduction in vectorial capacity"))
all_estimates$EHT=factor(all_estimates$EHT, levels=c("BIT103", "Kibondo\net al. 2022", "BIT080", "BIT055", "Odufuwa\net al. 2024"))

all_estimates_simple=all_estimates %>%
  mutate(washed_status=ifelse(washed_status=="Washed20", "Washed", washed_status),
         label_height=ifelse(param2=="Increase in host seeking mortality",-0.07, -0.03))

whitepointforscale=all_estimates_simple%>% filter(EHT_short=="Kibondo", washed_status=="Washed")%>%
  mutate(mean=NA, X2.5.=NA, X97.5.=NA,
         test=1)

# plot all estimates
all_estimates_simple %>%
  mutate( washed_status=ifelse(washed_status=="Unwashed", "Unwashed", "Washed 20x"))%>%
  ggplot()+
  geom_col(aes(x=insecticide_name, y=mean, fill=insecticide_name, group=interaction( washed_status, EHT), alpha=washed_status), stat="identity", position=position_dodge(preserve = "single"))+
  geom_errorbar(aes(x=insecticide_name, ymin=X2.5., ymax=X97.5., color=insecticide_name,  group=interaction( washed_status, EHT)),position=position_dodge(0.9, preserve = "single"), width=0.4)+
  geom_point(data=whitepointforscale, aes(x=insecticide_name, y=test), width=0.4, color="white")+
  #facet_grid( . ~param2, scales = "free")+
  facet_wrap( . ~param2, scales = "free")+
  theme_bw()+#ylim(0,1)+
  labs(x="", y="", color="", fill="", alpha="")+
  scale_alpha_manual(values=c(  0.8, 0.5),na.translate=FALSE)+
  scale_color_manual(values=c("darkorange","dodgerblue","darkblue", "red", "orange"  ))+
  scale_fill_manual(values=c("darkorange","dodgerblue","darkblue", "red", "orange" ))+
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 14),
        axis.text=element_text(size=14), strip.background =element_rect(fill="white"), legend.text =element_text(size=14) )+ 
  geom_text(aes(label=EHT_short, x=insecticide_name, y=label_height,group=EHT),
            position = position_dodge(width=0.9),size=3.5)+#ylim(-0.05,1)+ 
  guides(fill = "none",color = "none")
ggsave(file.path(plotDir, "plot_EHTfit_all_VCred_24.png"), width=10, height=12)



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
  gt::gtsave(filename = file.path(plotDir, "outputs_ento_efficacy_all_VCred24.docx"))


