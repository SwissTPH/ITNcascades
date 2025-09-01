rm(list=ls())
library(dplyr)
library(ggplot2)
mainDir="."
outputsDir=file.path(mainDir, "csv_outputs")
plotDir=file.path(mainDir, "plots")

all_ITNcov=data.frame()



optimise_decay_param=function(param, df_ITN_cov){
  #if(param[3]<1){
  out=data.frame(
    time=df_ITN_cov%>% pull(time),
    decay_obs=df_ITN_cov  %>% pull(ITNcov)
  )%>%
    mutate(decay=param[3]*exp( -(time/(param[1]))^param[2] * log(2) ),
           diff=(decay-decay_obs)^2)%>%
    summarise(diff=sum(diff))
  # } else {
  #   out=Inf
  # }
  return(out)
}


calculate_param=function(df_ITN_cov){
  opti=optim(c(0.448, 1.11, 0.7), optimise_decay_param, df_ITN_cov=df_ITN_cov , method = "L-BFGS-B", lower=c(0, 0, 0), upper = c(3, Inf, 1) )$par
  names(opti)=c("L", "k", "a")
  opti$setting=unique(df_ITN_cov$setting)
  return(opti)
}


#######################################################
#' Table S1 of Protopopoff et al. 2018 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(18)30427-6/fulltext
#' Table S2 of Protopopoff et al. 2024 https://malariajournal.biomedcentral.com/articles/10.1186/s12936-023-04727-8#Sec12
#' # proportion of residents declaring using a study net the previous night
ITN_covtrial=data.frame(arm=c(rep("Pyrethroid",6),rep("PBO",6))
                        ,time=rep(c(0.3, 0.75, 1.3, 1.75,2.3, 2.75),2)
                        ,ITNcov=c(1915/2567,2012/2510,1818/2669,1610/2595, 1386/2614, 784/2289,
                                  1997/2559,1820/2335,1571/2687,1445/2647, 1166/2459, 790/2513),
                        RCT="2018Protopopoff")

all_ITNcov=rbind(all_ITNcov, ITN_covtrial)

########################################################
## Use values not disagggregated by arm, but have similar trends (Supplemental file 2)
## Staedke 2020 page 1298

ITN_covtrial_1 = data.frame(time = c(.5,1,1.5)
                            ,ITNcov = c(.85,.79,.73),
                            RCT="2020Staedke", 
                            arm="Pyrethroid")
ITN_covtrial=rbind(ITN_covtrial_1, ITN_covtrial_1%>% mutate(arm="PBO"))


all_ITNcov=rbind(all_ITNcov, ITN_covtrial)



######################################################################################
#' Mosha 2022 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
#' Abtract, Findings

#' Mosha 2022 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
#' Appendix page 4, explicit dates in main paper page 1230
#' "Proportion of participants reporting using a net the night before"

ITN_covtrial_1=data.frame(arm=rep(c("Pyrethroid","IG2","PBO"),each=4),
                          time=rep(c(0.25,1,1.5,2),3),
                          ITNcov=c(868/1130,2848/4627,2606/4991,2488/5029,
                                   752/1099,3157/4833,2676/5194,2585/5576,
                                   771/1046,2506/4247,1885/4633,1534/5186)
)

# mosha 2024
#' Appendix S2 figure. Raw data not available. Only among selected children
ITN_covtrial_2=data.frame(arm=rep(c("Pyrethroid","IG2","PBO"),each=2),
                          time=rep(c(2.5,3),3),
                          ITNcov=c(0.33,0.25,
                                   0.29,0.21,
                                   0.16,0.09)
)

ITN_covtrial=rbind(ITN_covtrial_1, ITN_covtrial_2)%>% mutate(RCT="2022Mosha")
all_ITNcov=rbind(all_ITNcov, ITN_covtrial)

######################################################################################
#' Accrombessi 2023
#' Accrombessi 2024 Appendix 2.1
#' "study net coverage"

ITN_covtrial=data.frame(arm=rep(c("Pyrethroid","IG2"),each=6),
                        time=rep(c(0.5,0.75, 1.5, 2, 2.5, 3),2),
                        ITNcov=c(2127/3461,1867/2422,2023/3182,1430/2307, 2031/3329, 1400/2430,
                                 2287/3392,1935/2345,1951/2993,1440/2102, 1837/3204, 1278/2450),
                        RCT="2023Accrombessi"
) %>% filter(time!=0.5)
all_ITNcov=rbind(all_ITNcov, ITN_covtrial)

#########################################################################
#########################################################################
all_ITNcov$setting=paste0(all_ITNcov$RCT,"_" ,all_ITNcov$arm)

all_param=data.frame()
for(i in unique(all_ITNcov$setting)){
  print(i)
  all_param=rbind(all_param, calculate_param( df_ITN_cov =all_ITNcov %>% filter(setting==i)))
}

all_param=all_param%>%
  mutate(setting2=setting)%>%
  separate(setting2, into=c("RCT", "arm"))
write.csv(all_param, file.path(plotDir, "ITN_usage_fit.csv"), row.names = F)

##############
# Visualisation
curve_plot=expand.grid(list(time=seq(0, 3, 0.1),setting= unique(all_ITNcov$setting))) %>%
  left_join(all_param)%>%
  mutate(decay_calc=a*exp( -(time/(L))^k * log(2) )
  )%>%
  separate(setting, into=c("RCT", "arm"))

all_ITNcov%>%
  ggplot()+
  geom_point(aes(x=time, y=ITNcov))+
  geom_line(data=curve_plot, aes(x=time, y=decay_calc), col="dodgerblue")+
  facet_wrap(RCT~arm)+labs(x="Year after distribution", y="LLIN usage")+
  theme_bw()+ylim(0,1)+
  scale_color_manual(values=c("black", "grey"), name="")+
  geom_text(data=all_param, aes(label=paste0("Half-life: ",round(L,1), " years"), group=interaction(RCT, arm), x=1.5, y=0.9),size=4.5)+
  geom_text(data=all_param, aes(label=paste0("Initial usage: ",100*round(a,2), " %"), group=interaction(RCT, arm), x=0.6, y=0.5),size=4.5)+
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 14), axis.title=element_text(size=14), 
        axis.text=element_text(size=14), strip.background =element_rect(fill="white"), legend.text =element_text(size=14))
ggsave(file.path(plotDir, "all_trials_LLINusage.png"), width=12, height = 12)

