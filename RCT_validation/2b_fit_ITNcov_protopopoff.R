rm(list=ls())
library(dplyr)
library(ggplot2)

mainDir="."
scriptDir=file.path(mainDir, "RCTvalidation")

source(file.path(scriptDir, "calculate_decay_ITNs.R"))
mainDir="/Users/chamcl/switchdrive/AIM//5. Collaborations/External/Pie & Sarah/EHT_Fitting/"
plotDir=file.path(mainDir, "plots")



#' Table S1 of Protopopoff et al. 2018 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(18)30427-6/fulltext
#' Table S2 of Protopopoff et al. 2024 https://malariajournal.biomedcentral.com/articles/10.1186/s12936-023-04727-8#Sec12
#' # proportion of residents declaring using a study net the previous night
ITN_cov=data.frame(setting=c(rep("Pyrethroid",6),rep("PBO",6),rep("Pyrethroid+IRS",6),rep("PBO+IRS",6))
                   ,time=rep(c(0.3, 0.75, 1.3, 1.75,2.3, 2.75),4)
                   ,ITNcov_0=c(1915/2567,2012/2510,1818/2669,1610/2595, 1386/2614, 784/2289,
                               1997/2559,1820/2335,1571/2687,1445/2647, 1166/2459, 790/2513,
                               1893/2493,1761/2333,1629/2652,1368/2570, 1134/2584,713/2375,
                               2002/2533,1941/2522,1580/2546,1482/2739, 1124/2575,706/2397),
                   prop_studyNets=c(1219/1329, 1144/1199, 1034/1126, 844/1025,1,1,
                                    1241/1337, 1043/1110, 981/1081, 785/1021,1,1,
                                    1175/1241, 1042/1122, 1018/1131, 789/1004,1,1,
                                    1214/1291, 1132/1216, 913/1013, 771/1122, 1, 1
                   ))%>%
  mutate(#ITNcov=ITNcov_0*prop_studyNets
    ITNcov=ITNcov_0
  )
#ITN_cov=ITN_cov%>%filter(time<2)

# ITN_cov=data.frame(setting=c(rep("Pyrethroid",6),rep("PBO",6))
#                    ,time=rep(c(0.3, 0.75, 1.3, 1.75,2.3, 2.75),4)
#                    ,ITNcov=c(0.6, 0.71, 0.44, 0.39, 0.4, 0.20,
#                              0.72, 0.70, 0.4, 0.43, 0.38, 0.23))
# 

all_param=data.frame()
for(i in unique(ITN_cov$setting)){
  all_param=rbind(all_param, calculate_param(i))
}


curve_plot=expand.grid(list(time=seq(0, 3, 0.1),setting= unique(ITN_cov$setting))) %>%
  left_join(all_param)%>%
  mutate(decay_calc=a*exp( -(time/(L))^k * log(2) )
  )


ggplot(ITN_cov %>% filter(setting %in% c("PBO", "Pyrethroid")))+
  geom_point(aes(x=time, y=ITNcov))+
  geom_line(data=curve_plot %>% filter(setting %in% c("PBO", "Pyrethroid")), aes(x=time, y=decay_calc), col="dodgerblue")+
  facet_wrap(.~setting)+labs(x="Year after distribution", y="LLIN usage")+
  theme_bw()+ylim(0,1)+
  geom_text(data=all_param %>% filter(setting %in% c("PBO", "Pyrethroid")), aes(label=paste0("Half-life: ",round(L,1), " years"), group=setting, x=1.7, y=0.95),size=4.5)+
  geom_text(data=all_param %>% filter(setting %in% c("PBO", "Pyrethroid")), aes(label=paste0("Initial usage: ",100*round(a,2), " %"), group=setting, x=0.6, y=a+0.08),size=4.5)+
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 14), axis.title=element_text(size=14), 
        axis.text=element_text(size=14), strip.background =element_rect(fill="white"), legend.text =element_text(size=14))
ggsave(file.path(plotDir, "Protopopoff_llinUsage.png"), width=8, height = 5)

