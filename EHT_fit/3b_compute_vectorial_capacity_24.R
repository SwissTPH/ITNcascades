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


results_kibondo_w_24<- readRDS(file.path(stanDir,"/stan_kibondo_washed_24.rds"))
results_kibondo_unw_24<- readRDS(file.path(stanDir,"/stan_kibondo_unwashed_24.rds"))


results_bit055_w<- readRDS(file.path(stanDir,"/stan_BIT055_washed_24.rds"))
results_bit055_unw<- readRDS(file.path(stanDir,"/stan_BIT055_unwashed_24.rds"))


results_bit103_unw_24<- readRDS(file.path(stanDir,"/stan_assenga_unwashed_24_all.rds"))
results_bit103_w_24<- readRDS(file.path(stanDir,"/stan_assenga_washed_24_all.rds"))

results_bit103_unw_24_CI<- readRDS(file.path(stanDir,"/stan_assengaCotedIvoire_unwashed_24.rds"))
results_bit103_w_24_CI<- readRDS(file.path(stanDir,"/stan_assengaCotedIvoire_washed_24.rds"))

results_bit059_w<- readRDS(file.path(stanDir,"/stan_odufuwa_washed_24.rds"))
results_bit059_unw<- readRDS(file.path(stanDir,"/stan_odufuwa_unwashed_24.rds"))

results_bit080_w_24<- readRDS(file.path(stanDir,"/stan_BIT080_washed_24.rds"))
results_bit080_unw_24<- readRDS(file.path(stanDir,"/stan_BIT080_unwashed_24.rds"))


results_martin_w24<- readRDS(file.path(stanDir,"/stan_martin_36m_24_gambiae.rds"))
results_martin_unw24<- readRDS(file.path(stanDir,"/stan_martin_unwashed_24_gambiae.rds"))
results_martin_w24_funestus<- readRDS(file.path(stanDir,"/stan_martin_36m_24_funestus.rds"))
results_martin_unw24_funestus<- readRDS(file.path(stanDir,"/stan_martin_unwashed_24_funestus.rds"))


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
  VCred_kibondo=calculate_VC_uncertainty(results_EHT=results_kibondo_unw_24, nsamples=nsamples, names="IG2", insecticide_id=2)
  VCred_bit103_ig2=calculate_VC_uncertainty(results_EHT=results_bit103_unw_24, nsamples=nsamples, names="IG2", insecticide_id=3)
  VCred_bit103_ig2_CI=calculate_VC_uncertainty(results_EHT=results_bit103_unw_24_CI, nsamples=nsamples, names="IG2", insecticide_id=3)
  VCred_bit080=calculate_VC_uncertainty(results_EHT=results_bit080_unw_24, nsamples=nsamples, names="IG2", insecticide_id=2)
  VCred_martin_ig2=calculate_VC_uncertainty(results_EHT=results_martin_unw24, nsamples=nsamples, names="IG2", insecticide_id=3)
  VCred_martin_ig2_funestus=calculate_VC_uncertainty(results_EHT=results_martin_unw24_funestus, nsamples=nsamples, names="IG2", insecticide_id=3)

  impacts_ig2=rbind(data.frame(VCred=VCred_kibondo) %>% mutate(EHT_short="Kibondo"),
                    data.frame(VCred=VCred_bit103_ig2) %>% mutate(EHT_short="Assenga"),
                    data.frame(VCred=VCred_bit103_ig2_CI) %>% mutate(EHT_short="Assenga, CI"),
                    data.frame(VCred=VCred_bit080) %>% mutate(EHT_short="BIT080"),
                    data.frame(VCred=VCred_martin_ig2) %>% mutate(EHT_short="Martin"),
                    data.frame(VCred=VCred_martin_ig2_funestus) %>% mutate(EHT_short="Martin, f")
  )%>% mutate(washed_status="Unwashed", insecticide_name="Interceptor G2")
  
  
  # PBO unwashed
  VCred_bit055=calculate_VC_uncertainty(results_EHT=results_bit055_unw, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit103_pbo=calculate_VC_uncertainty(results_EHT=results_bit103_unw_24, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit103_pbo_CI=calculate_VC_uncertainty(results_EHT=results_bit103_unw_24_CI, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit059=calculate_VC_uncertainty(results_EHT=results_bit059_unw, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_martin_pbo=calculate_VC_uncertainty(results_EHT=results_martin_unw24, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_martin_pbo_funestus=calculate_VC_uncertainty(results_EHT=results_martin_unw24_funestus, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  
  
  impacts_pbo=rbind(data.frame(VCred=VCred_bit055) %>% mutate(EHT_short="BIT055"),
                    data.frame(VCred=VCred_bit103_pbo) %>% mutate(EHT_short="Assenga"),
                    data.frame(VCred=VCred_bit103_pbo_CI) %>% mutate(EHT_short="Assenga, CI"),
                    data.frame(VCred=VCred_bit059) %>% mutate(EHT_short="Odufuwa"),
                    data.frame(VCred=VCred_martin_pbo) %>% mutate(EHT_short="Martin"),
                    data.frame(VCred=VCred_martin_pbo_funestus) %>% mutate(EHT_short="Martin, f")
  )%>% mutate(washed_status="Unwashed", insecticide_name="Olyset Plus")
  
  
  # IG2 washed
  VCred_kibondo_w=calculate_VC_uncertainty(results_EHT=results_kibondo_w_24, nsamples=nsamples, names="IG2", insecticide_id=2)
  VCred_bit103_ig2_w=calculate_VC_uncertainty(results_EHT=results_bit103_w_24, nsamples=nsamples, names="IG2", insecticide_id=3)
  VCred_bit103_ig2_CI_w=calculate_VC_uncertainty(results_EHT=results_bit103_w_24_CI, nsamples=nsamples, names="IG2", insecticide_id=3)
  VCred_bit080_w=calculate_VC_uncertainty(results_EHT=results_bit080_w_24, nsamples=nsamples, names="IG2", insecticide_id=2)
  VCred_martin_ig2_w=calculate_VC_uncertainty(results_EHT=results_martin_w24, nsamples=nsamples, names="IG2", insecticide_id=3)
  VCred_martin_ig2_w_funestus=calculate_VC_uncertainty(results_EHT=results_martin_w24_funestus, nsamples=nsamples, names="IG2", insecticide_id=2)

  
  impacts_ig2_w=rbind(data.frame(VCred=VCred_kibondo_w) %>% mutate(EHT_short="Kibondo"),
                      data.frame(VCred=VCred_bit103_ig2_w) %>% mutate(EHT_short="Assenga"),
                      data.frame(VCred=VCred_bit103_ig2_CI_w) %>% mutate(EHT_short="Assenga, CI"),
                      data.frame(VCred=VCred_bit080_w) %>% mutate(EHT_short="BIT080"),
                      data.frame(VCred=VCred_martin_ig2_w) %>% mutate(EHT_short="Martin"),
                      data.frame(VCred=VCred_martin_ig2_w_funestus) %>% mutate(EHT_short="Martin, f")
  )%>% mutate(washed_status="Washed", insecticide_name="Interceptor G2")
  
  # PBO washed
  VCred_bit103_pbo_w=calculate_VC_uncertainty(results_EHT=results_bit103_w_24, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit103_pbo_CI_w=calculate_VC_uncertainty(results_EHT=results_bit103_w_24_CI, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit059_w=calculate_VC_uncertainty(results_EHT=results_bit059_w, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_bit055_w=calculate_VC_uncertainty(results_EHT=results_bit055_w, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_martin_pbo_w=calculate_VC_uncertainty(results_EHT=results_martin_w24, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  VCred_martin_pbo_w_funestus=calculate_VC_uncertainty(results_EHT=results_martin_w24_funestus, nsamples=nsamples, names="Olyset Plus", insecticide_id=2)
  
  impacts_pbo_w=rbind(data.frame(VCred=VCred_bit103_pbo_w) %>% mutate(EHT_short="Assenga"),
                      data.frame(VCred=VCred_bit103_pbo_CI_w) %>% mutate(EHT_short="Assenga, CI"),
                      data.frame(VCred=VCred_bit059_w) %>% mutate(EHT_short="Odufuwa"),
                      data.frame(VCred=VCred_bit055_w) %>% mutate(EHT_short="BIT055"),
                      data.frame(VCred=VCred_martin_pbo_w) %>% mutate(EHT_short="Martin"),
                      data.frame(VCred=VCred_martin_pbo_w_funestus) %>% mutate(EHT_short="Martin, f")
  )%>% mutate(washed_status="Washed", insecticide_name="Olyset Plus")
  
  # Pyrethroid unwashed
  VCred_kibondo_pyr=calculate_VC_uncertainty(results_EHT=results_kibondo_unw_24, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit080_pyr=calculate_VC_uncertainty(results_EHT=results_bit080_unw_24, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit103_pyr=calculate_VC_uncertainty(results_EHT=results_bit103_unw_24, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit103_pyr_CI=calculate_VC_uncertainty(results_EHT=results_bit103_unw_24_CI, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit059_pyr=calculate_VC_uncertainty(results_EHT=results_bit059_unw, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit055_pyr=calculate_VC_uncertainty(results_EHT=results_bit055_unw, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_martin_pyr=calculate_VC_uncertainty(results_EHT=results_martin_unw24, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_martin_pyr_funestus=calculate_VC_uncertainty(results_EHT=results_martin_unw24_funestus, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)

  impacts_pyr=rbind(data.frame(VCred=VCred_kibondo_pyr) %>% mutate(EHT_short="Kibondo"),
                    data.frame(VCred=VCred_bit055_pyr) %>% mutate(EHT_short="BIT055"),
                    data.frame(VCred=VCred_bit103_pyr) %>% mutate(EHT_short="Assenga"),
                    data.frame(VCred=VCred_bit103_pyr_CI) %>% mutate(EHT_short="Assenga, CI"),
                    data.frame(VCred=VCred_bit059_pyr) %>% mutate(EHT_short="Odufuwa"),
                    data.frame(VCred=VCred_bit080_pyr) %>% mutate(EHT_short="BIT080"),
                    data.frame(VCred=VCred_martin_pyr) %>% mutate(EHT_short="Martin"),
                    data.frame(VCred=VCred_martin_pyr_funestus) %>% mutate(EHT_short="Martin, f")
  )%>% mutate(washed_status="Unwashed", insecticide_name="Pyrethroid")
  
  
  
  # Pyrethroid washed
  VCred_kibondo_pyr_w=calculate_VC_uncertainty(results_EHT=results_kibondo_w_24, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit080_pyr_w=calculate_VC_uncertainty(results_EHT=results_bit080_w_24, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit103_pyr_w=calculate_VC_uncertainty(results_EHT=results_bit103_w_24, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit103_pyr_CI_w=calculate_VC_uncertainty(results_EHT=results_bit103_w_24_CI, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit059_pyr_w=calculate_VC_uncertainty(results_EHT=results_bit059_w, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_bit055_pyr_w=calculate_VC_uncertainty(results_EHT=results_bit055_w, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_martin_pyr_w=calculate_VC_uncertainty(results_EHT=results_martin_w24, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)
  VCred_martin_pyr_w_funestus=calculate_VC_uncertainty(results_EHT=results_martin_w24_funestus, nsamples=nsamples, names="Pyrethroid", insecticide_id=1)

  
  impacts_pyr_w=rbind(data.frame(VCred=VCred_kibondo_pyr_w) %>% mutate(EHT_short="Kibondo"),
                      data.frame(VCred=VCred_bit055_pyr_w) %>% mutate(EHT_short="BIT055"),
                      data.frame(VCred=VCred_bit103_pyr_w) %>% mutate(EHT_short="Assenga"),
                      data.frame(VCred=VCred_bit103_pyr_CI_w) %>% mutate(EHT_short="Assenga, CI"),
                      data.frame(VCred=VCred_bit059_pyr_w) %>% mutate(EHT_short="Odufuwa"),
                      data.frame(VCred=VCred_bit080_pyr_w) %>% mutate(EHT_short="BIT080"),
                      data.frame(VCred=VCred_martin_pyr_w) %>% mutate(EHT_short="Martin"),
                      data.frame(VCred=VCred_martin_pyr_w_funestus) %>% mutate(EHT_short="Martin, f")
  )%>% mutate(washed_status="Washed", insecticide_name="Pyrethroid")
  
  
  write.csv(rbind(impacts_ig2, impacts_pbo,impacts_ig2_w, impacts_pbo_w, impacts_pyr,impacts_pyr_w), file = file.path(outputsDir, "impacts24_EHT_control.csv"))
  
  impacts_ig2 %>% group_by(EHT_short, insecticide_name,washed_status)%>% summarise(mean=mean(VCred), min=min(VCred), max=max(VCred))
  impacts_pbo %>% group_by(EHT_short, insecticide_name,washed_status)%>% summarise(mean=mean(VCred), min=min(VCred), max=max(VCred))
  impacts_pyr %>% group_by(EHT_short, insecticide_name,washed_status)%>% summarise(mean=mean(VCred), min=min(VCred), max=max(VCred))
  
}


#######################
# PLOT SUMMARY
######################
seed=1
list_eht=read.csv(file.path(scriptDir, "EHTlist.csv"))

all_estimates0=data.frame()

nsamples=1000

for (i in list_eht$id){
  
  eht=list_eht[list_eht$id==i,]
  if(eht$stan_name_unwashed24 !="NULL"){
    results_unw<- readRDS(file.path(stanDir,eht$stan_name_unwashed24))
    results_w<- readRDS(file.path(stanDir,eht$stan_name_washed24))
    data=read.csv(file.path(dataDir, eht$data_name24)) %>% select(insecticide_name, treatment)%>% unique()
    
    estimates_unwashed=extract_stan_uncertainty(res=results_unw,nsamples=nsamples, seed=seed)%>% 
      mutate(washed_status="Unwashed") %>% rename(treatment=insecticide) %>%
      left_join(data)
    
    estimates_washed=extract_stan_uncertainty(res=results_w,nsamples=nsamples, seed=seed)%>% 
      mutate(washed_status="Washed")%>% rename(treatment=insecticide) %>%
      left_join(data)
    
    estimates=rbind(estimates_unwashed, estimates_washed)%>%mutate(EHT=eht$EHT, EHT_short=eht$EHT_short, id=eht$id)
    
    all_estimates0=rbind(all_estimates0, estimates)
  }
  
}

impacts_24=read.csv(file.path(outputsDir, "impacts24_EHT_control.csv"))%>%
  rename(nsample=X, mean=VCred)%>% mutate(param="VCred")%>%
  mutate(EHT_short0=EHT_short)%>%
  left_join(list_eht %>% select(EHT_short, EHT, id))%>% select(-EHT_short0)



#############################
# endpoint: 24h everywhere

# formating the results
all_estimates=rbind(all_estimates0 %>% select( -treatment),impacts_24 )%>%
  mutate(insecticide_name=ifelse(insecticide_name %in% c("IG2", "Interceptor_G2","InterceptorG2","Interceptor®G2", "Interceptor\xaeG2", "Interceptor G2", "interceptorG2", "IG2.Aged"), "Interceptor G2",
                                 ifelse(insecticide_name %in% c("Olyset_Plus", "OlysetPlus", "Olyset Plus", "olysetplus"), "Olyset Plus", 
                                        ifelse(insecticide_name %in% c("IG1", "IG1.Aged", "ILN", "interceptor", "Interceptor®", "Magnet", "MagNet_ITN", "MiraNet", "Olyset", "PermaNet2.0","Pyrethroid"), "Pyrethroid-only", "Control"))),
         param1=gsub("Initial", "", gsub("Efficacy", "", gsub("Rate" ,"",  param))))%>%
  filter(!param1 %in% c("KillingDuringHostSeeking", "alpha_0", "mu_0"))%>%
  filter(insecticide_name !="Control")

all_estimates$param2=factor(all_estimates$param1, levels=c("Repellent", "Preprandialkilling", "Postprandialkilling", "VCred"),
                            labels=c(expression(paste("Reduction in host availability (", pi,")")), expression(paste("Pre-prandial killing effect (", phi,")")),expression(paste("Post-prandial killing effect (", xi,")")), expression(paste("Vectorial capacity reduction"))))
all_estimates$param3=factor(all_estimates$param1, levels=c("Repellent", "Preprandialkilling", "Postprandialkilling", "VCred"),
                            labels=c("Reduction in host availability", "Pre-prandial killing effect","Post-prandial killing effect","Vectorial capacity reduction"))
all_estimates$EHT_id= factor(all_estimates$id,levels=c("1", "2", "3", "4","5", "6", "7", "8"))
all_estimates$country=ifelse(all_estimates$EHT_short %in% c("Sovegnon1","Sovegnon2", "Nguessan", "Assenga, CI"), "Benin", "Tanzania")

all_estimates_simple=all_estimates %>%
  mutate(washed_status=ifelse(washed_status=="Washed20", "Washed", washed_status))%>%
  rename(value=mean)

all_estimates_simple_mean=all_estimates_simple%>% group_by(insecticide_name, param2, washed_status)%>% summarise(meanmean=mean(value))%>%
  mutate( washed_status=ifelse(washed_status=="Unwashed", "Unwashed", "Washed 20x/Aged"))


# plot all estimates
all_estimates_simple %>%
  mutate( washed_status=ifelse(washed_status=="Unwashed", "Unwashed", "Washed 20x/Aged"))%>%
  mutate( washed_status=factor(washed_status, levels=c("Unwashed", "Washed 20x/Aged")))%>%
  mutate( washed_status2=factor(washed_status, levels=c("Unwashed", "Washed 20x/Aged"), labels=c("Unw.", "20x/Aged")))%>%
  ggplot()+
  geom_boxplot(aes(x=insecticide_name, y=value, fill=insecticide_name, color=insecticide_name,  group=interaction( EHT_id, washed_status, insecticide_name), alpha=washed_status), position=position_dodge(0.9), width=0.5)+
  geom_segment(data=all_estimates_simple_mean%>% filter(insecticide_name =="Interceptor G2",washed_status=="Unwashed"),
               aes(y = meanmean, yend=meanmean, x= 0.55, xend = .95, linetype = washed_status), linewidth = 0.9) +
  geom_segment(data=all_estimates_simple_mean%>% filter(insecticide_name =="Olyset Plus",washed_status=="Unwashed"),
               aes(y = meanmean, yend=meanmean, x= 1.55, xend = 1.95, linetype = washed_status), linewidth = 0.9) +
  geom_segment(data=all_estimates_simple_mean%>% filter(insecticide_name =="Pyrethroid-only",washed_status=="Unwashed"),
               aes(y = meanmean, yend=meanmean, x= 2.55, xend = 2.95, linetype = washed_status), linewidth = 0.9) +
  geom_segment(data=all_estimates_simple_mean%>% filter(insecticide_name =="Interceptor G2",washed_status!="Unwashed"),
               aes(y = meanmean, yend=meanmean, x= 1, xend = 1.45, linetype = washed_status), linewidth = 0.9) +
  geom_segment(data=all_estimates_simple_mean%>% filter(insecticide_name =="Olyset Plus",washed_status!="Unwashed"),
               aes(y = meanmean, yend=meanmean, x= 2, xend = 2.45, linetype = washed_status), linewidth = 0.9) +
  geom_segment(data=all_estimates_simple_mean%>% filter(insecticide_name =="Pyrethroid-only",washed_status!="Unwashed"),
               aes(y = meanmean, yend=meanmean, x= 3, xend = 3.45, linetype = washed_status), linewidth = 0.9) +
  facet_wrap( . ~param2, scales = "free_x", labeller = label_parsed)+
  theme_bw()+#ylim(0,1)+
  labs(x="", y="", color="", fill="", alpha="", linetype="")+
  scale_alpha_manual(values=c(  0.7, 0.2),na.translate=FALSE)+
  scale_color_manual(values=c("darkorange","dodgerblue","darkgrey" ))+
  scale_fill_manual(values=c("darkorange","dodgerblue","darkgrey"))+
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))+
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 14),
        axis.text=element_text(size=14), strip.background =element_rect(fill="white"), legend.text =element_text(size=14) )+ 
  geom_text(aes(label=EHT_id, x=insecticide_name, y=-0.03,group=interaction(EHT_id,washed_status)),
            position = position_dodge(width=0.9),size=2.3)+#ylim(-0.05,1)+ 
  geom_text(aes(label=washed_status2, x=insecticide_name, y=-0.1,group=washed_status),
            position = position_dodge(width=0.9),size=3.5)+#ylim(-0.05,1)+ 
  guides(fill = "none",color = "none", alpha="none")
ggsave(file.path(plotDir, "Figure1_24.png"), width=13, height=10)
