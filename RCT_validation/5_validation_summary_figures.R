rm(list=ls())
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(cowplot)


inputDir="RCT_validation/csv_inputs"
outputDir=""

# load results protopopoff 2018
results_protopopoff_detail=read.csv(file.path(inputDir, "results_protopopoff_detail.csv"))%>% 
  filter(year<2018, year>2013,
         age=="0.5-14", 
         is.na(PR_middle) ==F) %>% 
  mutate(RCT="Protopopoff et al.")
fitdat_protopopoff=read.csv(file.path(inputDir, "fitdat_protopopoff.csv")) %>% 
  mutate(RCT="Protopopoff et al.") 
validat_protopopoff=read.csv(file.path(inputDir, "validat_protopopoff.csv")) %>% 
  mutate(RCT="Protopopoff et al.") 

# load results mosha 2022
results_mosha_detail=read.csv(file.path(inputDir, "results_mosha_detail.csv")) %>% 
  filter(year>2017, year>2013) %>% 
  mutate(RCT="Mosha et al.") 
fitdat_mosha=read.csv(file.path(inputDir, "fitdat_mosha.csv")) %>% 
  mutate(RCT="Mosha et al.")
validat_mosha=read.csv(file.path(inputDir, "validat_mosha.csv")) %>% 
  mutate(RCT="Mosha et al.")


# load results Staedke 2020
results_staedke_detail=read.csv(file.path(inputDir, "results_staedke_detail.csv"))%>% 
  mutate(RCT="Staedke et al.", date=as.Date(paste0(year, "-",month, "-01")))%>%
  filter(date<"2019-02-15")%>% select(-date)
fitdat_staedke=read.csv(file.path(inputDir, "fitdat_staedke.csv")) %>% 
  mutate(RCT="Staedke et al.")%>%mutate(pos_children=NA)
validat_staedke=read.csv(file.path(inputDir, "validat_staedke.csv")) %>% 
  mutate(RCT="Staedke et al.")


# load results mosha 2022
results_accrombessi_detail=read.csv(file.path(inputDir, "results_accrombessi_detail.csv")) %>% 
  mutate(RCT="Accrombessi et al.", date=as.Date(paste0(year, "-",month, "-01")))%>%
  filter(date>"2019-09-15",date<"2023-01-01")%>% select(-date)
fitdat_accrombessi=read.csv(file.path(inputDir, "fitdat_accrombessi.csv")) %>% 
  mutate(RCT="Accrombessi et al.")
validat_accrombessi=read.csv(file.path(inputDir, "validat_accrombessi.csv")) %>% 
  mutate(RCT="Accrombessi et al.")


# combine datasets
agePR="0.5-14"

results_detail=rbind(results_protopopoff_detail , 
                     results_mosha_detail,
                     results_staedke_detail,
                    results_accrombessi_detail
                     )
validat=rbind(validat_protopopoff, validat_mosha,
              validat_staedke, 
              validat_accrombessi
              )
fitdat=rbind(fitdat_protopopoff %>% select(-trial, -arm)%>% unique(), 
             fitdat_mosha,
             fitdat_staedke, 
             fitdat_accrombessi
             )


min_dates=fitdat %>% mutate(date=as.Date(paste0(year,"-", month, "-01")))%>%
  group_by(setting, age, RCT)%>%
  summarise(min_date=min(date))

fitdat=fitdat%>% left_join(min_dates)%>% 
  mutate(date=as.Date(paste0(year,"-", month, "-01")),
         is_calib=(date==min_date))

fitdat_pooledAll= fitdat %>%
  filter(is_calib)%>%
  mutate(setting0=setting,
         setting=ifelse(str_detect( setting0, pattern = ("IG2")), "IG2",
                        ifelse(str_detect( setting, pattern = ("PBO" )), "OlysetPlus",
                               ifelse(str_detect( setting, pattern = ("OlysetPlus" )), "OlysetPlus", "Pyrethroid"))),
         status="Calibration points")%>%
  unique()%>% select(-date,-min_date, -is_calib)

validat_pooledAll= rbind(validat, fitdat %>%
                           filter(!is_calib)%>% select(-date,-min_date, -is_calib)) %>%
  mutate(setting0=setting,
         setting=ifelse(str_detect( setting0, pattern = ("IG2")), "IG2",
                        ifelse(str_detect( setting, pattern = ("PBO" )), "OlysetPlus",
                               ifelse(str_detect( setting, pattern = ("OlysetPlus" )), "OlysetPlus", "Pyrethroid"))),
         status="Validation points")%>%
  unique()
# 
fitdat_pooledAll$setting=factor(fitdat_pooledAll$setting, levels=c("Pyrethroid", "OlysetPlus", "IG2"))
validat_pooledAll$setting=factor(validat_pooledAll$setting, levels=c("Pyrethroid", "OlysetPlus", "IG2"))


####################################
# pool all data
results_pooled=results_detail%>% 
  mutate(setting=ifelse(setting=="pyrethroid",  "pyrethroid_.",setting))%>%
  separate(setting, into = c("setting","EHT"))%>%
  group_by(setting, RCT, year, month, age) %>% 
  summarise(PR_mini=min(PR_mini),PR_maxi=max(PR_maxi),PR_middle=median(PR_middle) )%>%    
  mutate(distribution_time=ifelse(RCT=="Mosha et al.",  2019+27/365,  
                                  ifelse(RCT=="Protopopoff et al.", 2015+27/365, 
                                         ifelse(RCT=="Staedke et al.", 2017+(7-1)/12,
                                         2020+3/12
                                         ))))


# Match the pyrethroid net type
results_pooled_perPyrethroid=results_detail%>% 
  mutate(setting=ifelse(setting=="pyrethroid",  "pyrethroid_.",setting))%>%
  separate(setting, into = c("setting","EHT"))%>%
  group_by(setting, RCT, year, month, age, pyrethroid) %>% 
  summarise(PR_mini=min(PR_mini),PR_maxi=max(PR_maxi),PR_middle=median(PR_middle) )%>%    
  mutate(distribution_time=ifelse(RCT=="Mosha et al.",  2019+27/365,  
                                  ifelse(RCT=="Protopopoff et al.", 2015+27/365, 
                                         ifelse(RCT=="Staedke et al.", 2017+(7-1)/12,
                                                2020+3/12
                                         ))))


results_pooled$setting=factor(results_pooled$setting,
                                          levels=c("pyrethroid", "PBO", "IG2"),
                                          labels=c("Pyrethroid-only", "Olyset Plus", "Interceptor G2"))

results_pooled_perPyrethroid$setting=factor(results_pooled_perPyrethroid$setting,
                              levels=c("pyrethroid", "PBO", "IG2"),
                              labels=c("Pyrethroid-only", "Olyset Plus", "Interceptor G2"))

fitdat_pooledAll$setting=factor(fitdat_pooledAll$setting,
                                          levels=c("Pyrethroid", "OlysetPlus", "IG2"),
                                          labels=c("Pyrethroid-only", "Olyset Plus", "Interceptor G2"))

validat_pooledAll$setting=factor(validat_pooledAll$setting,
                                          levels=c("Pyrethroid", "OlysetPlus", "IG2"),
                                          labels=c("Pyrethroid-only", "Olyset Plus", "Interceptor G2"))


plot_results_function=function(myRCT, agegroup, colors=c( "darkgrey", "dodgerblue","orange"), title, legend=F, results=results_pooled, fitdat=fitdat_pooledAll, validat=validat_pooledAll, breaks=NULL){
  p=ggplot(results %>% filter(RCT==myRCT), aes(x=year+(month-1)/12))+
    geom_ribbon(aes(ymin=PR_mini,ymax=PR_maxi, fill=setting),alpha=0.4)+
    geom_line(aes(y=PR_middle, color=setting),size=1)+
    geom_vline(aes(xintercept=distribution_time,linetype="Net distribution"),size=0.8)+
    geom_pointrange(data=fitdat %>% filter(RCT==myRCT),aes(y=PR_obs,ymin=PR_obs_min,ymax=PR_obs_max,shape=status),
                    size=0.7,stroke=1,color="black")+
    geom_pointrange(data=validat %>% filter(RCT==myRCT),aes(y=PR_obs,ymin=PR_obs_min,ymax=PR_obs_max,shape=status),
                    size=0.7,color="black")+
    scale_y_continuous(labels=scales::percent, limits = c(0,0.85))+
    scale_shape_manual(values=c( 1, 19))+
    scale_fill_manual(values=colors)+
    scale_color_manual(values=colors )+
    labs(x="",y=paste0("PfPR \n",agegroup),fill="Trial arm",linetype="", color="Trial arm", shape="")+
    facet_grid(RCT~setting, scale="free_x")+
    theme_bw()+
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_rect(colour = "black"), 
      legend.margin = margin(t = 0, b = 1, r=3, l=1),
      axis.text=element_text(size=12),
      axis.title=element_text(size=14), legend.text = element_text(size=15),
      legend.title = element_text(size=15),
      plot.title = element_text(size = 18),panel.spacing = unit(2, "lines")
    )+ggtitle(title)
  
  if(! legend){
    p=p+
      theme(legend.position = "none")
  }
  
  if(! is.null(breaks)){
    p=p+
      scale_x_continuous(breaks = breaks)
  }
  
  
  return(p)
}


# plots
plot_mosha_full_match=plot_results_function(myRCT="Mosha et al.", agegroup="(6 months to 14 years old)",title="A. Mosha et al. 2022, 2024",legend=T)

plot_mosha_match=plot_mosha_full_match+
  theme(legend.position = "none")

plot_protopopoff_match=plot_results_function(myRCT="Protopopoff et al.", agegroup="(6 months to 14 years old)",title="B. Protopopoff et al. 2018, 2023")
plot_staedke_match=plot_results_function(myRCT="Staedke et al.", agegroup="(2 to 10 years old)",title="D. Staedke et al. 2020", breaks=c(2017, 2018, 2019))
plot_accrombessi_match=plot_results_function(myRCT="Accrombessi et al.", agegroup="(all ages)",title="C. Accrombessi et al. 2023, 2024", colors=c( "darkgrey","orange"))


grobs <- ggplotGrob(plot_mosha_full_match)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]



plot_line2=plot_grid(plot_protopopoff_match, legend, ncol=2, rel_widths = c(2/3, 1/3))
plot_line3=plot_grid(plot_accrombessi_match, plot_staedke_match, ncol=2, rel_widths = c(1, 1))
plot_all_match=plot_grid(plot_mosha_match, plot_line2, plot_line3, ncol=1)
plot_all_match

ggsave(plot_all_match,filename=file.path(outputDir,"plot_all_match.png")
       ,width = 12,height=12)




############################
# Diagnostics 



# ON PREVALENCE
validat_pooledAll %>%
  ungroup()%>%unique()%>%
  left_join( results_pooled) %>%
  mutate(maxmin=pmax(PR_mini,PR_obs_min),
         minmax=pmin(PR_maxi,PR_obs_max),
         overlap=(maxmin<= minmax))%>%
  View()

regdat_match=validat_pooledAll %>%
  left_join( results_pooled) 

# check whether intervals overlap
interval_overlap=regdat_match %>% 
  mutate(overlap1=(PR_obs_min<=PR_maxi & PR_obs_min>=PR_mini), 
         overlap2=(PR_obs_max<=PR_maxi & PR_obs_max>=PR_mini),
         overlap=(overlap1|overlap2))%>%
  select(RCT, setting, year, month, overlap) 

interval_overlap %>% group_by(setting)%>% summarise(true=sum(as.numeric(!overlap)), tot=n())

# regression model
model_PR <- lm(PR_obs ~ PR_middle+0,
            data=regdat_match)
summary(model_PR)
Rsquared <- round(summary(model_PR)$adj.r.squared,2)
Rsquared


slope_PR <- round(summary(model_PR)$coefficients[1],2)


scatter_PR=validat_pooledAll %>%
  left_join( results_pooled) %>%
  ggplot()+
  geom_point(aes(x=PR_obs*100, y=PR_middle*100, color=setting, shape=RCT), size=2)+
  geom_linerange(aes(y=PR_middle*100, xmin=PR_obs_min*100, xmax=PR_obs_max*100, color=setting))+
  geom_linerange(aes(ymin=PR_mini*100, x=PR_obs*100, ymax=PR_maxi*100, color=setting))+
  scale_color_manual(values=c( "darkgrey", "dodgerblue","orange"))+
  geom_abline(slope=1, intercept = 0, linetype="dashed")+
  #geom_abline(slope=slope_PR, intercept = 0)+
  theme_minimal()+xlim(0,100)+ylim(0,100)+labs(color="", y='Modelled prevalence', x="Observed prevalence", shape="")+
  theme(legend.position = "right", legend.box="vertical")

ggsave(scatter_PR,filename=file.path(outputDir,"scatter_PR_match.png")
       ,width = 6,height=6)

##################
# on effect sizes

observation_dates_validation=validat_pooledAll %>% select(RCT, year, month)%>% unique()

results_filter=results_pooled_perPyrethroid %>% right_join(observation_dates_validation)

eff_size_model_match=results_filter %>% filter(setting =="Pyrethroid-only") %>% rename("PR_mini_control"="PR_mini","PR_maxi_control"="PR_maxi","PR_middle_control"="PR_middle") %>% 
  ungroup()%>% dplyr::select(-setting)%>%
  left_join(results_filter %>% filter(setting !="Pyrethroid-only"))%>%
  mutate(effect_size=100*(PR_middle_control-PR_middle )/PR_middle_control)


eff_size_model_match_pool=eff_size_model_match%>% group_by(RCT, setting, year, month)%>%
  summarise(effect_size=median(effect_size, na.rm = T) )

eff_size_data=validat_pooledAll %>% filter(setting =="Pyrethroid-only") %>% rename("PR_obs_min_control"="PR_obs_min","PR_obs_max_control"="PR_obs_max","PR_obs_control"="PR_obs") %>% ungroup()%>% dplyr::select(age, year, month, RCT, PR_obs_min_control, PR_obs_max_control, PR_obs_control)%>%
  unique()%>%
  left_join(validat_pooledAll %>% filter(setting !="Pyrethroid-only") %>% select(-setting0)%>% unique())%>%
  mutate(effect_size_obs=100*(PR_obs_control-PR_obs)/PR_obs_control)%>%
  filter(status=="Validation points")



reg_effsize_match=eff_size_data %>% left_join(eff_size_model_match_pool) %>% select(RCT, year, month, setting,effect_size, effect_size_obs )

model_effsize <- lm(effect_size_obs ~ effect_size +0,
            data=reg_effsize_match)
summary(model_effsize)
Rsquared_effsize <- round(summary(model_effsize)$adj.r.squared,2)
Rsquared_effsize

slope_effsize<- round(summary(model_effsize)$coefficients[1],2)

scatter_effsize=eff_size_data %>% left_join(eff_size_model_match_pool) %>% 
  ggplot()+
  geom_point(aes(y=effect_size, x=effect_size_obs, color=setting, shape=RCT), size=2)+
  scale_color_manual(values=c( "dodgerblue","orange"))+xlim(-10,100)+ ylim(0,100)+
  geom_abline(slope=1, intercept = 0, linetype="dashed")+
  #geom_abline(slope=slope_effsize, intercept = 0)+
  theme_minimal()+xlim(0,100)+ylim(0,100)+
  labs(color="", y='Modelled effect size', x="Observed effect size",shape="")+
  theme(legend.position = "right")
ggsave(filename=file.path(outputDir,"scatter_effectSize_match.png")
       ,width = 6,height=6)

##################
# on effect sizes

observation_dates_all=rbind(fitdat_pooledAll %>% mutate(status="baseline"),validat_pooledAll %>%mutate(status="post")) %>% select(RCT, setting, year, month, status)%>% unique()

results_filter_all=results_pooled_perPyrethroid %>% right_join(observation_dates_all)

before_after_match=results_filter_all %>% filter(status =="baseline") %>% rename("PR_mini_baseline"="PR_mini","PR_maxi_baseline"="PR_maxi","PR_middle_baseline"="PR_middle") %>% 
  ungroup()%>% dplyr::select(-month, -year, -status)%>%
  left_join(results_filter_all %>% filter(status=="post"))%>%
  mutate(effect_size=100*(PR_middle_baseline-PR_middle )/PR_middle_baseline)

before_after_match_pool=before_after_match%>% group_by(RCT, setting, year, month)%>%
  summarise(effect_size=median(effect_size, na.rm = T) )

before_after_data=fitdat_pooledAll %>% rename("PR_obs_min_baseline"="PR_obs_min","PR_obs_max_baseline"="PR_obs_max","PR_obs_baseline"="PR_obs") %>%
  ungroup()%>% dplyr::select(age, setting, RCT, PR_obs_min_baseline, PR_obs_max_baseline, PR_obs_baseline)%>%
  unique()%>%
  left_join(validat_pooledAll  %>% select(-setting0)%>% unique())%>%
  mutate(effect_size_obs=100*(PR_obs_baseline-PR_obs)/PR_obs_baseline)



reg_before_after_match=before_after_data %>% left_join(before_after_match_pool) %>% select(RCT, year, month, setting,effect_size, effect_size_obs )

model_before_after <- lm(effect_size_obs ~ effect_size +0,
                    data=reg_before_after_match)
summary(model_before_after)
Rsquared_before_after <- round(summary(model_before_after)$adj.r.squared,2)
Rsquared_before_after

slope_before_after<- round(summary(model_before_after)$coefficients[1],2)


scatter_before_after=before_after_data %>% left_join(before_after_match_pool) %>% 
  ggplot()+
  geom_point(aes(y=effect_size, x=effect_size_obs, color=setting, shape=RCT), size=2)+
  scale_color_manual(values=c("darkgrey", "dodgerblue","orange"))+#xlim(-10,100)+ ylim(0,100)+
  geom_abline(slope=1, intercept = 0, linetype="dashed")+
  #geom_abline(slope=slope_before_after, intercept = 0)+
  theme_minimal()+#xlim(0,100)+ylim(0,100)+
  labs(color="", y='Modelled reduction to baseline', x="Observed reduction to baseline",shape="")+
  theme(legend.position = "right")
ggsave(filename=file.path(outputDir,"scatter_before_after_match.png")
       ,width = 6,height=6)




grobs_scatter <- ggplotGrob(scatter_PR)$grobs
legend_scatter <- grobs_scatter[[which(sapply(grobs_scatter, function(x) x$name) == "guide-box")]]


plot_all_scatter=plot_grid(scatter_PR+theme(legend.position = "none"), 
                            scatter_effsize+theme(legend.position = "none"), 
                            scatter_before_after+theme(legend.position = "none"),
                            ncol=1, rel_widths = c(1,1,1), scale = 0.95, labels=c("A", "B", "C"))
plot_all_scatter
plot_all_scatter=plot_grid(plot_all_scatter,
                           legend_scatter, ncol=2, scale=c(1,0.8), rel_widths  = c(1,0.3))
plot_all_scatter
ggsave(plot_all_scatter,filename=file.path(outputDir,"plot_scatter.png") ,width = 6,height=8)


