rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)


mainDir="."
scriptDir=file.path(mainDir, "RCTvalidation")

source(file.path(scriptDir, "calculate_decay_ITNs.R"))
mainDir="/Users/chamcl/switchdrive/AIM//5. Collaborations/External/Pie & Sarah/EHT_Fitting/"
plotDir=file.path(mainDir, "plots")



#' Mosha 2022 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
#' Abtract, Findings

#' Mosha 2022 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)02499-5/fulltext#sec1
#' Appendix page 4, explicit dates in main paper page 1230
#' "Proportion of participants reporting using a net the night before"

ITN_cov=data.frame(setting=rep(c("pyrethroid","pyriproxyfen","chlorfenapyr","PBO"),each=4),
                   time=rep(c(0.25,1,1.5,2),4),
                   ITNcov=c(868/1130,2848/4627,2606/4991,2488/5029,
                            764/1103,2646/4358,2124/4645,2087/5455,
                            752/1099,3157/4833,2676/5194,2585/5576,
                            771/1046,2506/4247,1885/4633,1534/5186)
)


#' "Proportion of selected children reporting using a net the night before"
ITN_cov_c=data.frame(setting=rep(c("pyrethroid","pyriproxyfen","chlorfenapyr","PBO"),each=3),
                     time=rep(c(1,1.5,2),4),
                     ITNcov=c(714/1200,671/1277,587/1246,
                              659/1149,566/1213,461/1306,
                              787/1230, 671/1281,615/1316,
                              659/1138,483/1205,370/1309)
)
# 

# mosha 2024
#' Appendix S2 figure. Raw data not available. Only among selected children
ITN_cov_after=data.frame(setting=rep(c("pyrethroid","pyriproxyfen","chlorfenapyr","PBO"),each=2),
                         time=rep(c(2.5,3),4),
                         ITNcov=c(0.33,0.25,
                                  0.25,0.14,
                                  0.29,0.21,
                                  0.16,0.09)
)

ITN_cov_before=data.frame(setting=rep(c("pyrethroid","pyriproxyfen","chlorfenapyr","PBO"),each=1),
                          time=rep(c(0.25),4),
                          ITNcov=c(0.76,
                                   0.68,
                                   0.67,
                                   0.73)
)


ITN_cov=rbind(ITN_cov, ITN_cov_after)
ITN_cov_children=rbind(ITN_cov_c, ITN_cov_after, ITN_cov_before)



all_param=data.frame()
for(i in unique(ITN_cov$setting)){
  all_param=rbind(all_param, calculate_param(i))
}

curve_plot=expand.grid(list(time=seq(0, 3, 0.1),setting= unique(ITN_cov$setting))) %>%
  left_join(all_param)%>%
  mutate(decay_calc=a*exp( -(time/(L))^k * log(2) )
  )

ITN_cov%>%
  filter(setting !="pyriproxyfen")%>%
  ggplot()+
  geom_point(aes(x=time, y=ITNcov, col="all age groups"))+
  geom_point(data=ITN_cov_children%>%
               filter(setting !="pyriproxyfen"),aes(x=time, y=ITNcov, col="selected children"))+
  geom_line(data=curve_plot%>%
              filter(setting !="pyriproxyfen"), aes(x=time, y=decay_calc), col="dodgerblue")+
  facet_wrap(.~setting)+labs(x="Year after distribution", y="LLIN usage")+
  theme_bw()+ylim(0,1)+
  scale_color_manual(values=c("black", "grey"), name="")+
  geom_text(data=all_param%>%
              filter(setting !="pyriproxyfen"), aes(label=paste0("Half-life: ",round(L,1), " years"), group=setting, x=1.5, y=0.9),size=4.5)+
  geom_text(data=all_param%>%
              filter(setting !="pyriproxyfen"), aes(label=paste0("Initial usage: ",100*round(a,2), " %"), group=setting, x=0.6, y=a+0.05),size=4.5)+
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 14), axis.title=element_text(size=14), 
        axis.text=element_text(size=14), strip.background =element_rect(fill="white"), legend.text =element_text(size=14))
ggsave(file.path(plotDir, "Mosha_llinUsage.png"), width=12, height = 5)


