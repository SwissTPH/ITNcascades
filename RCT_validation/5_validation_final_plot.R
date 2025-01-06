rm(list=ls())
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(cowplot)


mainDir="."
scriptDir=file.path(mainDir, "RCTvalidation/")
source(file.path(scriptDir, "helpfunctions_rct.R"))

# load results protopopoff 2018
results_protopopoff=read.csv(file.path(mainDir, "results_protopopoff.csv"))%>% 
  filter(year<2018, year>2013) %>% 
  mutate(RCT="Protopopoff et al.")
results_protopopoff_detail=read.csv(file.path(mainDir, "results_protopopoff_detail.csv"))%>% 
  filter(year<2018, year>2013) %>% 
  mutate(RCT="Protopopoff et al.")
fitdat_protopopoff=read.csv(file.path(mainDir, "fitdat_protopopoff.csv")) %>% 
  mutate(RCT="Protopopoff et al.") 
validat_protopopoff=read.csv(file.path(mainDir, "validat_protopopoff.csv")) %>% 
  mutate(RCT="Protopopoff et al.") 

fitdat_protopopoff_plot=fitdat_protopopoff%>% full_join(validat_protopopoff %>% filter(setting=="Pyrethroid" & PR_obs<0.8))
validat_protopopoff_plot=validat_protopopoff%>% filter(setting!="Pyrethroid" | PR_obs>0.8)

# load results mosha 2022
results_mosha=read.csv(file.path(mainDir, "results_mosha.csv")) %>% 
  filter(year>2017, year>2013) %>% 
  mutate(RCT="Mosha et al.")
results_mosha_detail=read.csv(file.path(mainDir, "results_mosha_detail.csv")) %>% 
  filter(year>2017, year>2013) %>% 
  mutate(RCT="Mosha et al.")
fitdat_mosha=read.csv(file.path(mainDir, "fitdat_mosha.csv")) %>% 
  mutate(RCT="Mosha et al.")
validat_mosha=read.csv(file.path(mainDir, "validat_mosha.csv")) %>% 
  mutate(RCT="Mosha et al.")

fitdat_mosha_plot=fitdat_mosha%>% full_join(validat_mosha %>% filter(setting=="pyrethroid" & PR_obs<0.8))
validat_mosha_plot=validat_mosha%>% filter(setting!="pyrethroid" | PR_obs>0.8)


# combine datasets
agePR="0.5-14"

results=rbind(results_protopopoff %>% filter(age==agePR), results_mosha%>% filter(age==agePR))
results_detail=rbind(results_protopopoff_detail %>% filter(age==agePR), results_mosha_detail%>% filter(age==agePR))
validat=rbind(validat_protopopoff_plot, validat_mosha_plot)
fitdat=rbind(fitdat_protopopoff_plot %>% select(-trial), fitdat_mosha_plot)


fitdat_pooledAll= fitdat %>%
  mutate(setting0=setting,
         setting=ifelse(str_detect( setting0, pattern = ("IG2")), "IG2",
                        ifelse(str_detect( setting, pattern = ("PBO" )), "Olyset Plus",
                               ifelse(str_detect( setting, pattern = ("OlysetPlus" )), "OlysetPlus", "Pyrethroid"))),
         status="Calibration points")%>%
  unique()

validat_pooledAll= validat %>%
  mutate(setting0=setting,
         setting=ifelse(str_detect( setting0, pattern = ("IG2")), "IG2",
                        ifelse(str_detect( setting, pattern = ("PBO" )), "OlysetPlus",
                               ifelse(str_detect( setting, pattern = ("OlysetPlus" )), "OlysetPlus", "Pyrethroid"))),
         status="Validation points")%>%
  unique()
# 
fitdat_pooledAll$setting=factor(fitdat_pooledAll$setting, levels=c("Pyrethroid", "OlysetPlus", "IG2"))
validat_pooledAll$setting=factor(validat_pooledAll$setting, levels=c("Pyrethroid", "OlysetPlus", "IG2"))


###############
# Match the pyrethroid net type

results_PBO_match=rbind(
  results_mosha_match=results_mosha_detail %>%
    filter(setting =="PBO_BIT103", age==agePR),
  results_protopopoff_detail %>%
    filter(setting =="OlysetPlus_BIT055", age==agePR))%>%
  mutate(setting="OlysetPlus")
row.names(results_PBO_match)=NULL

results%>% filter(setting=="OlysetPlus")





results_match_pooledModels= rbind(results %>% filter(setting !="OlysetPlus"),
                                  results_PBO_match)%>%    
  mutate(distribution_time=ifelse(RCT=="Mosha et al.",  2019+27/365,  2015+27/365))


results_match_pooledModels$setting=factor(results_match_pooledModels$setting,
                                          levels=c("Pyrethroid", "OlysetPlus", "IG2"),
                                          labels=c("Pyrethroid-only", "Olyset Plus", "Interceptor G2"))

fitdat_pooledAll$setting=factor(fitdat_pooledAll$setting,
                                          levels=c("Pyrethroid", "OlysetPlus", "IG2"),
                                          labels=c("Pyrethroid-only", "Olyset Plus", "Interceptor G2"))

validat_pooledAll$setting=factor(validat_pooledAll$setting,
                                          levels=c("Pyrethroid", "OlysetPlus", "IG2"),
                                          labels=c("Pyrethroid-only", "Olyset Plus", "Interceptor G2"))

# plots
plot_mosha_full_match=ggplot(results_match_pooledModels %>% filter(RCT=="Mosha et al."), aes(x=year+(month-1)/12))+
  geom_ribbon(aes(ymin=PR_mini,ymax=PR_maxi, fill=setting),alpha=0.4)+
  geom_line(aes(y=PR_middle, color=setting),size=1)+
  geom_vline(aes(xintercept=distribution_time,linetype="Net distribution"),size=0.8)+
  geom_pointrange(data=fitdat_pooledAll %>% filter(RCT=="Mosha et al."),aes(y=PR_obs,ymin=PR_obs_min,ymax=PR_obs_max,shape=status),
                  size=0.7,stroke=1,color="black")+
  geom_pointrange(data=validat_pooledAll %>% filter(RCT=="Mosha et al."),aes(y=PR_obs,ymin=PR_obs_min,ymax=PR_obs_max,shape=status),
                  size=0.7,color="black")+
  scale_y_continuous(labels=scales::percent, limits = c(0,0.85))+
  scale_shape_manual(values=c( 1, 19))+
  scale_fill_manual(values=c( "darkgrey", "dodgerblue","orange"))+
  scale_color_manual(values=c("darkgrey", "dodgerblue","orange") )+
  #scale_x_continuous(labels=function(x){return(paste0("Jan\n",x))})+
  labs(x="",y="PfPR \n(6 months to 14 years old)",fill="Trial arm",linetype="", color="Trial arm", shape="")+
  facet_grid(RCT~setting, scale="free_x")+
  theme_minimal()+#guides(fill=guide_legend(ncol=3)) + 
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
  )+ggtitle("A. Mosha et al. 2022, 2024")

plot_mosha_match=plot_mosha_full_match+
  theme(legend.position = "none")



plot_protopopoff_match=ggplot(results_match_pooledModels %>% filter(RCT=="Protopopoff et al."), aes(x=year+(month-1)/12))+
  geom_ribbon(aes(ymin=PR_mini,ymax=PR_maxi, fill=setting),alpha=0.4)+
  geom_line(aes(y=PR_middle, color=setting),size=1)+
  geom_vline(aes(xintercept=distribution_time,linetype="Net distribution"),size=0.8)+
  geom_pointrange(data=fitdat_pooledAll %>% filter(RCT=="Protopopoff et al."),aes(y=PR_obs,ymin=PR_obs_min,ymax=PR_obs_max,shape=status),
                  size=0.7,stroke=1,color="black")+
  geom_pointrange(data=validat_pooledAll %>% filter(RCT=="Protopopoff et al."),aes(y=PR_obs,ymin=PR_obs_min,ymax=PR_obs_max,shape=status),
                  size=0.7,color="black")+
  scale_y_continuous(labels=scales::percent, limits = c(0,0.85))+
  scale_shape_manual(values=c( 1, 19))+
  scale_fill_manual(values=c( "darkgrey", "dodgerblue","orange"))+
  scale_color_manual(values=c("darkgrey", "dodgerblue","orange") )+
  # scale_x_continuous(labels=function(x){return(paste0("Jan\n",x))})+
  labs(y="PfPR \n(6 months to 14 years old)",fill="",linetype="", color="", shape="", x="")+
  facet_grid(RCT~setting, scale="free_x")+
  theme_minimal()+
  theme(legend.position = "none") + 
  theme(  strip.background = element_blank(),    strip.text = element_blank(),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14), legend.text = element_text(size=15),
          legend.title = element_text(size=15),
          plot.title = element_text(size = 18),panel.spacing = unit(2, "lines")
  )+ggtitle("B. Protopopoff et al. 2018, 2023")


grobs <- ggplotGrob(plot_mosha_full_match)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

plot0_match=plot_grid(plot_protopopoff_match, legend, ncol=2, rel_widths = c(2/3, 1/3))
plot_all_match=plot_grid(plot_mosha_match, plot0_match, ncol=1)
plot_all_match

ggsave(filename=file.path(mainDir,"plot_all_match.png")
       ,width = 12,height=12)



############################
# Diagnostics / Goodness of fit

# on prevalence
validat_pooledAll %>%
  ungroup()%>%unique()%>%
  left_join( results_match_pooledModels) %>%
  mutate(maxmin=pmax(PR_mini,PR_obs_min),
         minmax=pmin(PR_maxi,PR_obs_max),
         overlap=(maxmin<= minmax))%>%
  View()

validat_pooledAll %>%
  left_join( results_match_pooledModels) %>%
  ggplot()+
  geom_point(aes(x=PR_middle, y=PR_obs, color=setting, shape=RCT), size=2)+
  geom_linerange(aes(x=PR_middle, ymin=PR_obs_min, ymax=PR_obs_max, color=setting))+
  geom_linerange(aes(xmin=PR_mini, y=PR_obs, xmax=PR_maxi, color=setting))+
  scale_color_manual(values=c( "darkgrey", "dodgerblue","orange"))+
  geom_abline(slope=1, intercept = 0)+
  theme_minimal()+xlim(0,1)+ylim(0,1)+labs(color="", x='Modelled prevalence', y="Observed prevalence", shape="")+
  theme(legend.position = "bottom")
ggsave(filename=file.path(mainDir,"scatter_PR_match.png")
       ,width = 6,height=6)

regdat_match=validat_pooledAll %>%
  left_join( results_match_pooledModels) 
model <- lm(PR_obs ~ PR_middle,
            data=regdat_match)
summary(model)
Rsquared <- round(summary(model)$r.squared,2)

cor.test(regdat_match$PR_obs, regdat_match$PR_middle)

##################
# on effect sizes
eff_size_model_match=results_match_pooledModels %>% filter(setting =="Pyrethroid-only") %>% rename("PR_mini_control"="PR_mini","PR_maxi_control"="PR_maxi","PR_middle_control"="PR_middle") %>% ungroup()%>% dplyr::select(-setting)%>%
  left_join(results_match_pooledModels %>% filter(setting !="Pyrethroid-only"))%>%
  mutate(effect_size=100*(PR_middle_control-PR_middle )/PR_middle_control,
         effect_size_min=pmin(100*(PR_mini_control-PR_mini)/PR_mini_control,100*(PR_maxi_control-PR_maxi)/ PR_maxi_control),
         effect_size_max=pmax(100*(PR_mini_control-PR_mini)/PR_mini_control,100*(PR_maxi_control-PR_maxi)/ PR_maxi_control))


eff_size_data=rbind(fitdat_pooledAll,validat_pooledAll) %>% filter(setting =="Pyrethroid-only") %>% rename("PR_obs_min_control"="PR_obs_min","PR_obs_max_control"="PR_obs_max","PR_obs_control"="PR_obs") %>% ungroup()%>% dplyr::select(age, year, month, RCT, PR_obs_min_control, PR_obs_max_control, PR_obs_control)%>%
  unique()%>%
  left_join(rbind(fitdat_pooledAll,validat_pooledAll) %>% filter(setting !="Pyrethroid-only") %>% select(-setting0)%>% unique())%>%
  mutate(effect_size_obs=100*(PR_obs_control-PR_obs)/PR_obs_control,
         effect_size_obs_min=pmin(100*(PR_obs_min_control-PR_obs_min)/PR_obs_min_control,100*(PR_obs_max_control-PR_obs_max)/ PR_obs_max_control),
         effect_size_obs_max=pmax(100*(PR_obs_min_control-PR_obs_min)/PR_obs_min_control,100*(PR_obs_max_control-PR_obs_max)/ PR_obs_max_control))%>%
  filter(status=="Validation points")


eff_size_data %>% left_join(eff_size_model_match) %>%
  ggplot()+
  geom_point(aes(x=effect_size, y=effect_size_obs, color=setting, shape=RCT), size=2)+
  geom_linerange(aes(x=effect_size, ymin=effect_size_obs_min, ymax=effect_size_obs_max, color=setting))+
  geom_linerange(aes(xmin=effect_size_min, y=effect_size_obs, xmax=effect_size_max, color=setting))+
  scale_color_manual(values=c( "dodgerblue","orange"))+
  geom_abline(slope=1, intercept = 0)+
  theme_minimal()+xlim(0,100)+ylim(0,100)+labs(color="", x='Modelled effect size', y="Observed effect size",shape="")+
  theme(legend.position = "bottom")
ggsave(filename=file.path(mainDir,"scatter_effectSize_match.png")
       ,width = 6,height=6)



reg_effsize_match=eff_size_data %>% left_join(eff_size_model_match)
cor.test(reg_effsize_match$effect_size, reg_effsize_match$effect_size_obs)
model <- lm(effect_size_obs ~ effect_size ,
            data=reg_effsize_match)
summary(model)
Rsquared <- round(summary(model)$r.squared,2)
